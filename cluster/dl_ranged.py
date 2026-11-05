import sys, os, math, requests, concurrent.futures as cf
import asf_search as asf

url, out, nthreads = sys.argv[1], sys.argv[2], int(sys.argv[3] if len(sys.argv) > 3 else 8)
s = asf.ASFSession()
r = s.get(url, stream=True, allow_redirects=True, timeout=120)
r.raise_for_status()
signed = r.url
total = int(r.headers["Content-Length"])
r.close()
print(f"size {total/2**30:.1f} GiB, {nthreads} range streams", flush=True)
chunk = math.ceil(total / nthreads)
with open(out, "wb") as f:
    f.truncate(total)

def pull(i):
    a, b = i * chunk, min((i + 1) * chunk, total) - 1
    if a > b: return 0
    want = b - a + 1
    err = None
    for attempt in range(2):
        try:
            rr = requests.get(signed, headers={"Range": f"bytes={a}-{b}"}, stream=True, timeout=300)
            rr.raise_for_status()
            if rr.status_code != 206:
                # A 200 means the server ignored the Range header; writing the
                # full body at offset a would overwrite every other chunk.
                raise IOError(f"range {a}-{b}: got HTTP {rr.status_code}, need 206")
            n = 0
            with open(out, "r+b") as f:
                f.seek(a)
                for blk in rr.iter_content(1 << 22):
                    f.write(blk); n += len(blk)
            if n != want:
                raise IOError(f"range {a}-{b}: short body {n} != {want}")
            return n
        except Exception as e:
            err = e
            if attempt == 0:
                print(f"range {a}-{b} attempt 1 failed ({e}), retrying", flush=True)
    raise err

with cf.ThreadPoolExecutor(nthreads) as ex:
    got = sum(ex.map(pull, range(nthreads)))
assert got == total, f"short read {got} != {total}"
print("RANGED_OK", flush=True)
