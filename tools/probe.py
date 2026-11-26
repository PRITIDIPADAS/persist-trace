import time, traceback, urllib.request, shutil, gzip, os
os.makedirs("out", exist_ok=True)
url = "https://zenodo.org/records/8192729/files/CollecTRI_regulons.csv?download=1"
for i in range(8):
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "PersistTrace"})
        with urllib.request.urlopen(req, timeout=120) as r, gzip.open("out/collectri_regulons.csv.gz", "wb") as f:
            shutil.copyfileobj(r, f)
        print("ok", os.path.getsize("out/collectri_regulons.csv.gz")); break
    except Exception as e:
        print("try", i, e); time.sleep(30)
with gzip.open("out/collectri_regulons.csv.gz", "rt") as f:
    for _ in range(3): print(f.readline().strip())
