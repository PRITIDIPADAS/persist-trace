import os, urllib.request, shutil
os.makedirs("out",exist_ok=True)
u="https://commons.wikimedia.org/wiki/Special:FilePath/Mucinous_lung_adenocarcinoma_and_airway_-_alt_--_high_mag.jpg?width=1800"
req=urllib.request.Request(u,headers={"User-Agent":"PersistTrace/1.0 (science fair; github.com/PRITIDIPADAS/persist-trace)"})
with urllib.request.urlopen(req,timeout=300) as r, open("out/histology.jpg","wb") as f: shutil.copyfileobj(r,f)
print(os.path.getsize("out/histology.jpg"))
