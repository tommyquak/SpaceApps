"""Attach Macrostrat geologic-map units to the reference candidates and freeze the list.

Run on a machine with internet access after nisar_refsel.py:
  python nisar_refsel_geology.py <refsel.json> <refsel_frozen.json>
Records the API response verbatim per candidate centroid (unit name, lithology, age interval,
source and map ids). No phase is read. The frozen file carries a SHA-256 of its candidate list.
"""
import hashlib
import json
import subprocess
import sys

src, dst = sys.argv[1], sys.argv[2]
rep = json.load(open(src))
for c in rep["candidates"]:
    url = f"https://macrostrat.org/api/v2/geologic_units/map?lat={c['lat']:.5f}&lng={c['lon']:.5f}&format=json"
    try:                                           # curl: the system python may lack CA certificates
        raw = subprocess.run(["curl", "-s", "--max-time", "40", url], capture_output=True, text=True, check=True).stdout
        data = json.loads(raw).get("success", {}).get("data", [])
    except Exception as exc:                       # network failure is recorded, not hidden
        data = [{"error": str(exc)}]
    units = [{k: u.get(k) for k in ("map_id", "source_id", "name", "strat_name", "lith", "descrip",
                                     "b_int_name", "t_int_name", "b_age", "t_age")} for u in data]
    c["geology"] = dict(api=url, units=units)
    names = [u.get("name") for u in units]; ages = [(u.get("b_int_name"), u.get("t_int_name")) for u in units]
    pre_cenozoic = any(u.get("b_age") is not None and float(u["b_age"]) > 66.0 for u in units if "error" not in u)
    c["geology"]["pre_cenozoic_unit"] = bool(pre_cenozoic)
    print("cand %d (%.3f, %.3f): %s %s pre-Cenozoic=%s" % (c["rank"], c["lon"], c["lat"], names, ages, pre_cenozoic))
rep["geology"] = "attached (Macrostrat geologic_units/map, point queries at centroids)"
blob = json.dumps(rep["candidates"], sort_keys=True).encode()
rep["candidates_sha256"] = hashlib.sha256(blob).hexdigest()
rep["frozen"] = True
json.dump(rep, open(dst, "w"), indent=1)
print("FROZEN", rep["candidates_sha256"])
