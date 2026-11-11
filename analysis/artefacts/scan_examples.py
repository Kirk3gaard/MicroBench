import json, gzip, random, datetime, numpy as np, pysam, pandas as pd
random.seed(5)
S=json.load(open("scan_stats.json")); starts=np.array(S["scan_start_h"])*3600
t0=pd.read_csv("full.reads.tsv.gz",sep="\t",usecols=["st"]).st; t0=t0[t0>0].min()
pool=[]
for line in gzip.open("full.artefacts.jsonl.gz","rt"):
    r=json.loads(line)
    if "low_complexity" not in r["cls"] or not (1500<=r["len"]<=30000): continue
    a=[x for x in r["lowc"] if x[1]-x[0]>=100 and x[0]<0.05*r["len"]]
    if not a: continue
    st=datetime.datetime.fromisoformat(r["st"]).timestamp()-t0
    d=st-starts[np.argmin(np.abs(starts-st))]
    if 0<=d<=330: r["since_scan_s"]=round(float(d),1); pool.append(r)
pick=random.sample(pool,min(12,len(pool)))
want={r["id"]:r for r in pick}; seqs={}
for a in pysam.AlignmentFile("sup520.bam",check_sq=False,threads=4):
    if a.query_name in want:
        seqs[a.query_name]=(a.query_sequence,list(a.query_qualities))
        if len(seqs)==len(want): break
exec(open("build_full_page.py").read().split("def thin")[0].split("P, BAM, STATS, OUT")[0])
src=open("build_full_page.py").read()
exec("def thin"+src.split("def thin")[1].split("examples = {")[0])
ex=[]
for r in pick:
    e=example(r); e["since_scan_s"]=r["since_scan_s"]; ex.append(e)
json.dump({"n_pool":len(pool),"examples":ex},open("scan_examples.json","w"),separators=(",",":"))
print(len(pool),len(ex))
