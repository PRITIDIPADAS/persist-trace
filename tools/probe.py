import gzip, io, os, re, subprocess, tarfile, urllib.request, json
UA={"User-Agent":"PersistTrace probe"}
def get(url, dest):
    req=urllib.request.Request(url,headers=UA)
    with urllib.request.urlopen(req,timeout=600) as r, open(dest,"wb") as f:
        while True:
            b=r.read(1<<20)
            if not b: break
            f.write(b)
    return dest
def head_gz(p,n=8):
    with gzip.open(p,"rt",errors="replace") as f:
        return [next(f,"").rstrip("\n")[:300] for _ in range(n)]
out={}
os.makedirs("dl",exist_ok=True)
base="https://ftp.ncbi.nlm.nih.gov/geo/series/"
def series(acc): return base+acc[:-3]+"nnn/"+acc+"/"
for acc in ["GSE275199","GSE243562","GSE335846","GSE335847","GSE335848","GSE193259"]:
    try:
        idx=get(series(acc)+"suppl/","dl/%s_idx.html"%acc)
        out[acc+"_suppl"]=re.findall(r'href="([^"?/][^"]*)"',open(idx).read())
    except Exception as e: out[acc+"_suppl"]=str(e)
    try:
        m=get(series(acc)+"matrix/%s_series_matrix.txt.gz"%acc,"dl/%s_m.gz"%acc)
        lines=[l for l in gzip.open(m,"rt",errors="replace") if l.startswith(("!Sample_title","!Sample_geo","!Sample_characteristics","!Sample_source","!Series_sample_id","!Sample_description"))]
        out[acc+"_matrix"]=[l[:600] for l in lines]
    except Exception as e: out[acc+"_matrix"]=str(e)
s=series("GSE275199")+"suppl/"
for f in ["GSE275199_Samples_hashtags_barcodesequences.txt.gz","GSE275199_features.tsv.gz","GSE275199_barcodes.tsv.gz"]:
    p=get(s+f,"dl/"+f); out[f]=head_gz(p,12)
    if "features" in f:
        rows=[l.rstrip("\n").split("\t") for l in gzip.open(p,"rt")]
        out["features_types"]={t:sum(1 for r in rows if r[-1]==t) for t in set(r[-1] for r in rows)}
        out["features_nonGene"]=[r for r in rows if r[-1]!="Gene Expression"][:20]
p=get(s+"GSE275199_matrix.mtx.gz","dl/m.mtx.gz"); out["mtx_head"]=head_gz(p,4); out["mtx_bytes"]=os.path.getsize(p)
p=get(series("GSE243562")+"suppl/GSE243562_RAW.tar","dl/raw.tar")
with tarfile.open(p) as t: out["GSE243562_tar"]=[(m.name,m.size) for m in t.getmembers()]
with tarfile.open(p) as t:
    for m in t.getmembers():
        if "features" in m.name or "genes" in m.name:
            out["GSE243562_features_head"]=gzip.decompress(t.extractfile(m).read()).decode()[:1200]; break
p=get(series("GSE335846")+"suppl/GSE335846_RNA-seq_for_GEO.txt.gz","dl/c.txt.gz"); out["GSE335846_head"]=head_gz(p,6)
for pdb in ["4ZAU","6JX0","6JXT","7UKV"]:
    try:
        p=get("https://files.rcsb.org/header/%s.pdb"%pdb,"dl/%s.txt"%pdb); out["pdb_"+pdb]=[l[:90] for l in open(p) if l.startswith(("HEADER","TITLE","COMPND   2","HETNAM"))][:8]
    except Exception as e: out["pdb_"+pdb]=str(e)
json.dump(out,open("probe_out.json","w"),indent=1)
print(json.dumps(out,indent=1)[:3000])
