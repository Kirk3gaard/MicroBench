import sys
P = sys.argv[1] if len(sys.argv) > 1 else "full"
BAM = sys.argv[2] if len(sys.argv) > 2 else BAM
TAG = sys.argv[3] if len(sys.argv) > 3 else ""
import pandas as pd, numpy as np, json
R=pd.read_csv(P + ".reads.tsv.gz",sep="\t",usecols=["id","len","qs","ch","st","du","split","aln","del","mm","ins","n_del10","cls"])
t0=R.st[R.st>0].min(); R["t"]=R.st-t0; R["end"]=R.t+R.du
for c in ["chimera_interspecies","chimera_intraspecies","foldback","internal_adapter","low_complexity","skip_prone","unmapped","end_adapter"]:
    R[c]=R.cls.str.contains(c,regex=False)
# mass-end seconds -> cluster into events (gap > 10 min)
ends=np.bincount(np.clip(R.end.astype(int),0,None))
act=np.zeros(len(ends)+2); np.add.at(act,R.t.clip(lower=0).astype(int),1); np.add.at(act,(R.end.astype(int)+1),-1); act=np.cumsum(act)[:len(ends)]
rate=ends.sum()/act.sum(); exp=act*rate
from scipy.stats import poisson
p=poisson.sf(ends-1,np.maximum(exp,1e-9)); mass=np.where((p<0.01/len(ends))&(ends>3*exp)&(ends>=5))[0]
clusters=[]; 
for s in mass:
    if not clusters or s-clusters[-1][-1]>600: clusters.append([s])
    else: clusters[-1].append(s)
starts=np.array([c[0] for c in clusters]); stops=np.array([c[-1] for c in clusters])
print("scan events",len(clusters),"start h",np.round(starts/3600,3).tolist()); print("durations s",(stops-starts).tolist())
# time since end of last scan for each read start
idx=np.searchsorted(stops,R.t.values)-1
since=np.where(idx>=0,R.t.values-stops[np.clip(idx,0,None)],np.nan)
R["since_min"]=since/60
L=R[(R["len"]>=1000)]
L=L.assign(b=pd.cut(L.since_min,[-1,1,2,5,10,20,40,60,90,200]))
g=L.groupby("b",observed=True).agg(n=("id","size"),lowc=("low_complexity","mean"),ia=("internal_adapter","mean"),inter=("chimera_interspecies","mean"),fold=("foldback","mean"),skip=("skip_prone","mean"),qs=("qs","median"),medlen=("len","median"))
L2=L[L.aln>0]; g["delkb"]=L2.groupby("b",observed=True).apply(lambda x:(1000*x["del"]/x.aln).median())
print(g.round(4).to_string())
json.dump({"scan_start_h":(starts/3600).round(4).tolist(),"scan_dur_s":(stops-starts).tolist(),"since":g.reset_index().astype({"b":str}).round(5).to_dict("records")},open(TAG + "scan_stats.json","w"))
