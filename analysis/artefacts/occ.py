import pandas as pd, numpy as np, json
R=pd.read_csv("full.reads.tsv.gz",sep="\t",usecols=["id","len","qs","ch","st","du","split","aln","cls"])
t0=R.st[R.st>0].min(); R["t"]=R.st-t0
for c in ["chimera_interspecies","chimera_intraspecies","foldback","internal_adapter","low_complexity"]:
    R[c]=R.cls.str.contains(c,regex=False)
R["bin"]=(R.t//900).astype(int)
L=R[(R["len"]>=1000)&(R.aln>0)]
g=L.groupby("bin").agg(n=("id","size"),ch=("ch","nunique"),inter=("chimera_interspecies","mean"),intra=("chimera_intraspecies","mean"),fold=("foldback","mean"),ia=("internal_adapter","mean"),lowc=("low_complexity","mean"),split=("split","mean"),medlen=("len","median"))
g=g[g.n>500]
print(g.iloc[::8].round(4).to_string())
print("corr with active channels:", g[["inter","intra","fold","ia","lowc","split"]].corrwith(g.ch).round(3).to_dict())
C=L.groupby("ch").agg(n=("id","size"),inter=("chimera_interspecies","mean"),fold=("foldback","mean"),lowc=("low_complexity","mean"),split=("split","mean"))
C=C[C.n>=300]
print("channel corr with reads/channel:", C[["inter","fold","lowc","split"]].corrwith(C.n,method="spearman").round(3).to_dict())
q=pd.qcut(C.n,5); print(C.groupby(q)[["inter","fold","lowc","split"]].mean().round(4).to_string())
