import httpx, pathlib
CYCLES = {"D": 2005, "E": 2007, "F": 2009, "G": 2011, "H": 2013, "I": 2015, "J": 2017}
FILES = ["DEMO", "GHB", "BIOPRO", "RXQ_RX", "DIQ", "MCQ"]

OUT = pathlib.Path("data/raw"); OUT.mkdir(parents=True, exist_ok=True)
for suf, yr in CYCLES.items():
    for f in FILES:
        dest = OUT / f"{f}_{suf}.xpt"
        if dest.exists() and dest.stat().st_size > 20000:
            continue
        url = f"https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{yr}/DataFiles/{f}_{suf}.xpt"
        r = httpx.get(url, follow_redirects=True, timeout=120)
        ok = r.status_code == 200 and r.content[:6] == b"HEADER"
        if ok:
            dest.write_bytes(r.content)
        print(("OK   " if ok else "FAIL ") + dest.name)
