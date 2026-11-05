import httpx, pathlib
CYCLES = {"D": 2005, "E": 2007, "F": 2009, "G": 2011, "H": 2013, "I": 2015, "J": 2017}
FILES = ["BMX", "SMQ", "ALB_CR", "HDL", "TCHOL", "BPQ", "BPX"]
OUT = pathlib.Path("data/raw")
for suf, yr in CYCLES.items():
    for f in FILES:
        name = "BPXO" if (f == "BPX" and suf == "J") else f   # 2017-18 uses a different BP file
        dest = OUT / f"{name}_{suf}.xpt"
        if dest.exists() and dest.stat().st_size > 20000:
            continue
        url = f"https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{yr}/DataFiles/{name}_{suf}.xpt"
        r = httpx.get(url, follow_redirects=True, timeout=120)
        ok = r.status_code == 200 and r.content[:6] == b"HEADER"
        if ok:
            dest.write_bytes(r.content)
        print(("OK   " if ok else "FAIL ") + dest.name)
