import json, gzip, numpy as np, datetime, collections
S=json.load(open("scan_stats.json")); starts=np.array(S["scan_start_h"])*3600; dur=np.array(S["scan_dur_s"])
import pandas as pd
t0=pd.read_csv("full.reads.tsv.gz",sep="\t",usecols=["st"]).st
t0=t0[t0>0].min()
rel=[]; relend=[]; pos=[]; n=0; tot=0
for line in gzip.open("full.artefacts.jsonl.gz","rt"):
    r=json.loads(line)
    if "low_complexity" not in r["cls"]: continue
    st=datetime.datetime.fromisoformat(r["st"]).timestamp()-t0
    for a,b in r["lowc"]:
        if b-a<100: continue
        t=st+r["du"]*a/r["len"]
        i=np.argmin(np.abs(starts-t)); rel.append(t-starts[i]); pos.append(a/r["len"])
        tot+=1
rel=np.array(rel); pos=np.array(pos)
print("lowc runs>=100bp",tot)
h=np.histogram(rel,bins=[-3600,-600,-60,-10,0,10,30,60,120,210,300,600,3600])
for e,c in zip(h[1],h[0]): print(f"{e:>6}s {c}")
near=(rel>=-5)&(rel<=215)
print("frac within scan window",near.mean().round(3),"window share of run", (210*27/(41.57*3600)))
print("position in read: near-scan median",np.median(pos[near]).round(3),"other",np.median(pos[~near]).round(3))
json.dump({"hist_edges":h[1].tolist(),"hist":h[0].tolist(),"frac_in_window":float(near.mean()),"window_share":210*27/(41.57*3600),
 "pos_near":float(np.median(pos[near])),"pos_other":float(np.median(pos[~near])),"n_runs":int(tot)},open("junk_time.json","w"))
