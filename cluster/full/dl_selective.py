"""Download only the CEOS members we keep, instead of the whole ALOS-2 zip.

Each ALOS-2 ScanSAR L1.1 granule is ~56.6 GiB, but the pipeline keeps only
VOL/LED/TRL and IMG-HH frames 2-4 (~16.7 GB): the rest is other polarisations
and frames. ASF serves the granule over a signed URL that honours HTTP range
requests, so the zip's central directory can be read first and only the needed
members fetched. That turns 56.6 GiB per date into roughly what we keep.

  dl_selective.py probe <url>
  dl_selective.py get <url> <outdir> [threads]

Falls back cleanly: any parse or range failure exits non-zero so the caller can
run the whole-zip downloader instead.
"""
import io
import os
import zlib
import re
import sys
import zipfile
import concurrent.futures as cf

import requests

KEEP = re.compile(r"^(VOL-|LED-|TRL-|IMG-HH-.*-F[234]$)")
CHUNK = 1 << 22


def signed_url(url):
    import asf_search as asf
    s = asf.ASFSession()
    r = s.get(url, stream=True, allow_redirects=True, timeout=120)
    r.raise_for_status()
    u, total = r.url, int(r.headers["Content-Length"])
    r.close()
    if not total:
        raise SystemExit("dl_selective: server gave no Content-Length")
    return u, total


class RangeFile(io.RawIOBase):
    """Seekable read-only file over HTTP ranges, enough for zipfile."""

    def __init__(self, url, size, session=None):
        self.url, self.size, self.pos = url, size, 0
        self.s = session or requests.Session()

    def seekable(self):
        return True

    def readable(self):
        return True

    def seek(self, off, whence=0):
        self.pos = off if whence == 0 else (self.pos + off if whence == 1 else self.size + off)
        return self.pos

    def tell(self):
        return self.pos

    def read(self, n=-1):
        if n < 0:
            n = self.size - self.pos
        n = min(n, self.size - self.pos)
        if n <= 0:
            return b""
        a, b = self.pos, self.pos + n - 1
        r = self.s.get(self.url, headers={"Range": f"bytes={a}-{b}"}, timeout=300)
        r.raise_for_status()
        d = r.content
        self.pos += len(d)
        return d


def members(url):
    u, total = signed_url(url)
    zf = zipfile.ZipFile(RangeFile(u, total))
    return u, total, zf, [i for i in zf.infolist() if KEEP.match(os.path.basename(i.filename))]


def probe(url):
    u, total, zf, keep = members(url)
    all_n = len(zf.infolist())
    comp = sum(i.compress_size for i in keep)
    full = sum(i.file_size for i in keep)
    print(f"granule {total / 2**30:.1f} GiB, {all_n} members")
    for i in keep:
        print(f"  keep {os.path.basename(i.filename):55s} "
              f"{i.compress_size / 2**20:9.1f} MiB -> {i.file_size / 2**20:9.1f} MiB "
              f"method {i.compress_type}")
    print(f"selective download {comp / 2**30:.2f} GiB of {total / 2**30:.1f} GiB "
          f"({100 * comp / total:.1f} percent), expands to {full / 2**30:.2f} GiB")


def _data_offset(rf, info):
    """Start of a member's data: the local header is 30 bytes plus the name and
    extra fields, whose lengths live in the local header, not the central one."""
    rf.seek(info.header_offset)
    h = rf.read(30)
    if h[:4] != b"PK\x03\x04":
        raise RuntimeError(f"{info.filename}: no local header at {info.header_offset}")
    n = int.from_bytes(h[26:28], "little")
    e = int.from_bytes(h[28:30], "little")
    return info.header_offset + 30 + n + e


def get(url, outdir, threads=12):
    u, total, zf, keep = members(url)
    if not keep:
        raise SystemExit("dl_selective: no matching members in the granule")
    for i in keep:
        if i.compress_type != zipfile.ZIP_STORED:
            raise SystemExit(f"dl_selective: {i.filename} is compressed "
                             f"(method {i.compress_type}); use the whole-zip downloader")
    os.makedirs(outdir, exist_ok=True)
    probe_rf = RangeFile(u, total)
    plan = []
    for i in keep:
        start = _data_offset(probe_rf, i)
        plan.append((os.path.basename(i.filename), start, i.file_size, i.CRC))
    comp = sum(p[2] for p in plan)
    print(f"fetching {len(plan)} members, {comp / 2**30:.2f} GiB of {total / 2**30:.1f} GiB "
          f"({100 * comp / total:.1f} percent)", flush=True)

    # split every member into contiguous pieces and pull them with one shared
    # pool, so a single 5 GB image still uses the full stream budget
    PIECE = 1 << 29        # 512 MiB
    jobs = []
    for name, start, size, _ in plan:
        dst = os.path.join(outdir, name)
        with open(dst + ".part", "wb") as f:
            f.truncate(size)
        for off in range(0, size, PIECE):
            jobs.append((dst + ".part", start + off, off, min(PIECE, size - off)))

    sess = requests.Session()

    def pull(j):
        path, src_off, dst_off, n = j
        for attempt in (1, 2, 3):
            try:
                r = sess.get(u, headers={"Range": f"bytes={src_off}-{src_off + n - 1}"},
                             stream=True, timeout=600)
                r.raise_for_status()
                got = 0
                with open(path, "r+b") as f:
                    f.seek(dst_off)
                    for blk in r.iter_content(CHUNK):
                        f.write(blk)
                        got += len(blk)
                if got != n:
                    raise RuntimeError(f"short read {got} != {n}")
                return n
            except Exception as e:
                if attempt == 3:
                    raise
                print(f"  retry {os.path.basename(path)} +{dst_off}: {e}", flush=True)
        return 0

    with cf.ThreadPoolExecutor(threads) as ex:
        got = sum(ex.map(pull, jobs))
    if got != comp:
        raise SystemExit(f"dl_selective: got {got} bytes, expected {comp}")
    # every member carries a CRC-32 in the zip directory, so the fetched bytes
    # are checked against the archive itself rather than against a local copy
    # the pipeline may already have consumed
    for name, _, size, crc in plan:
        dst = os.path.join(outdir, name)
        part = dst + ".part"
        if os.path.getsize(part) != size:
            raise SystemExit(f"dl_selective: {name} is {os.path.getsize(part)} B, expected {size}")
        h = 0
        with open(part, "rb") as f:
            while True:
                b = f.read(1 << 24)
                if not b:
                    break
                h = zlib.crc32(b, h)
        if h != crc:
            os.remove(part)
            raise SystemExit(f"dl_selective: {name} CRC {h:08x} != {crc:08x} from the zip directory")
        os.replace(part, dst)
        print(f"  {name} {size / 2**20:.1f} MiB crc {h:08x} ok", flush=True)
    print("SELECTIVE_OK", flush=True)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    if sys.argv[1] == "probe":
        probe(sys.argv[2])
    elif sys.argv[1] == "get":
        get(sys.argv[2], sys.argv[3], int(sys.argv[4]) if len(sys.argv) > 4 else 6)
    else:
        raise SystemExit(__doc__)
