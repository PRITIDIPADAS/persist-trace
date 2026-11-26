import os, re, urllib.request, json, gzip, shutil
UA={"User-Agent":"PersistTrace probe"}
def get(url,dest):
    with urllib.request.urlopen(urllib.request.Request(url,headers=UA),timeout=600) as r, open(dest,"wb") as f: shutil.copyfileobj(r,f)
    return dest
os.makedirs("out",exist_ok=True); o={}
for g in ["GSM8472663","GSM8472664"]:
    u="https://ftp.ncbi.nlm.nih.gov/geo/samples/%snnn/%s/suppl/"%(g[:-3],g)
    try:
        get(u,"idx.html"); names=re.findall(r'href="([^"?/][^"]*)"',open("idx.html").read()); o[g]=names
        for n in names:
            if g=="GSM8472664" and n.endswith(".gz"):
                p=get(u+n,"out/"+n); o[n]=[l[:200] for l in gzip.open(p,"rt",errors="replace").read(3000).splitlines()[:10]]; o[n+"_bytes"]=os.path.getsize(p)
    except Exception as e: o[g]=str(e)
for pdb in ["4ZAU","6JX0"]:
    get("https://files.rcsb.org/download/%s.pdb"%pdb,"out/%s.pdb"%pdb)
json.dump(o,open("out/probe2.json","w"),indent=1); print(json.dumps(o,indent=1)[:4000])
