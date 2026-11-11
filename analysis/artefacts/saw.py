import sys
P = sys.argv[1] if len(sys.argv) > 1 else "full"
BAM = sys.argv[2] if len(sys.argv) > 2 else BAM
TAG = sys.argv[3] if len(sys.argv) > 3 else ""
import pandas as pd, numpy as np, json
R=pd.read_csv(P + ".reads.tsv.gz",sep="\t",usecols=["len","qs","ch","st","du","aln","del","mm","ins","cls"])
t0=R.st[R.st>0].min(); R["t"]=R.st-t0
L=R[(R["len"]>=1000)&(R.aln>0)].copy()
L["delkb"]=1000*L["del"]/L.aln; L["mmkb"]=1000*L.mm/L.aln; L["inskb"]=1000*L.ins/L.aln
L["lowc"]=L.cls.str.contains("low_complexity",regex=False); L["ia"]=L.cls.str.contains("internal_adapter",regex=False)
L["fold"]=L.cls.str.contains("foldback",regex=False)
L["bin"]=(L.t//300).astype(int)
g=L.groupby("bin").agg(n=("qs","size"),ch=("ch","nunique"),dels=("delkb","median"),mm=("mmkb","median"),ins=("inskb","median"),qs=("qs","median"),lowc=("lowc","mean"),ia=("ia","mean"),fold=("fold","mean")).reset_index()
g["h"]=g.bin/12
g=g[g.n>=200]
json.dump(g.drop(columns="bin").round(4).to_dict("records"),open(TAG + "saw.json","w"))
S=json.load(open(TAG + "scan_stats.json")); st=np.array(S["scan_start_h"])
# within-cycle: per-read del vs minutes since scan, per cycle, normalised
print(g[["h","n","ch","dels","mm","lowc"]].iloc[:40].round(3).to_string())
