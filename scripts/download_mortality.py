import httpx, pathlib
BASE = "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality/"
CYCLES = ["2005_2006", "2007_2008", "2009_2010", "2011_2012",
          "2013_2014", "2015_2016", "2017_2018"]
OUT = pathlib.Path("data/raw")
for c in CYCLES:
    name = f"NHANES_{c}_MORT_2019_PUBLIC.dat"
    r = httpx.get(BASE + name, follow_redirects=True, timeout=120)
    head = r.content[:40]
    ok = r.status_code == 200 and not head.lstrip().startswith(b"<")
    if ok:
        (OUT / name).write_bytes(r.content)
    print(("OK   " if ok else "FAIL ") + name, len(r.content))
