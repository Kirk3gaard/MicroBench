"""Separate chimera junctions at recurrent genome positions from real chimeras.

A ligation chimera joins two random molecules, so its junction sits at a random position. When
many reads break at the same position, the cause is the reference, not the library: a structural
difference between the sequenced strain and the reference, a repeat (rRNA operon, IS element)
that the aligner placed in another copy, or the end of a linear contig/plasmid in the reference.

For every contig the expected number of junctions per 200-bp bin is computed from the contig's
junction density; bins with far more junctions than chance (Poisson, Bonferroni over all bins,
and at least 5 reads) are recurrent sites. Reads whose every chimera junction lies at a recurrent
site lose their chimera class and get "recurrent_junction" instead.

usage: refine_chimeras.py <prefix> ref.fa      (rewrites <prefix>.reads.tsv.gz and
       <prefix>.artefacts.jsonl.gz; originals kept as *.raw.*)
"""
import sys, os, gzip, json, shutil
import numpy as np, pandas as pd, mappy
from scipy.stats import poisson

P, REF = sys.argv[1], sys.argv[2]
BIN = 200
raw_r, raw_a = P + ".reads.raw.tsv.gz", P + ".artefacts.raw.jsonl.gz"
if not os.path.exists(raw_r): shutil.copy(P + ".reads.tsv.gz", raw_r)
if not os.path.exists(raw_a): shutil.copy(P + ".artefacts.jsonl.gz", raw_a)
clen = {n: len(s) for n, s, _ in mappy.fastx_read(REF)}

E = pd.read_csv(P + ".events.tsv.gz", sep="\t", usecols=["id", "type", "ctg", "rpos"])
J = E[E.type.isin(["junction_chimera", "junction_foldback"])].copy()
J["bin"] = (J.rpos // BIN).astype(int)
site = J.groupby(["ctg", "bin"]).id.nunique().rename("n").reset_index()
dens = J.groupby("ctg").id.nunique()
nbins = sum(max(1, L // BIN) for L in clen.values())
lam = site.ctg.map(lambda c: dens.get(c, 0) * BIN / max(BIN, clen.get(c, BIN)))
site["p"] = poisson.sf(site.n - 1, lam.clip(lower=1e-9))
site["recurrent"] = (site.p < 0.01 / nbins) & (site.n >= 5)
rec = set(map(tuple, site.loc[site.recurrent, ["ctg", "bin"]].values))
J["rec"] = [(c, b) in rec for c, b in zip(J.ctg, J.bin)]
chim_j = J[J.type == "junction_chimera"]
all_rec = chim_j.groupby("id").rec.all()
move = set(all_rec.index[all_rec])

def fix(cls):
    c = [x for x in cls if x not in ("chimera_interspecies", "chimera_intraspecies")]
    c.append("recurrent_junction")
    return c

R = pd.read_csv(raw_r, sep="\t")
m = R.id.isin(move)
R.loc[m, "cls"] = R.loc[m, "cls"].map(lambda s: ",".join(fix(s.split(",")) if s != "clean" else ["recurrent_junction"]))
R.to_csv(P + ".reads.tsv.gz", sep="\t", index=False, compression={"method": "gzip", "compresslevel": 3})
with gzip.open(raw_a, "rt") as fi, gzip.open(P + ".artefacts.jsonl.gz", "wt", compresslevel=3) as fo:
    for line in fi:
        r = json.loads(line)
        if r["id"] in move:
            r["cls"] = fix(r["cls"])
        fo.write(json.dumps(r, separators=(",", ":")) + "\n")
top = site[site.recurrent].sort_values("n", ascending=False).head(15)
out = {"recurrent_sites": int(site.recurrent.sum()), "reads_moved": int(m.sum()),
       "inter_before": int(pd.read_csv(raw_r, sep="\t", usecols=["cls"]).cls.str.contains("chimera_interspecies", regex=False).sum()),
       "inter_after": int(R.cls.str.contains("chimera_interspecies", regex=False).sum()),
       "intra_after": int(R.cls.str.contains("chimera_intraspecies", regex=False).sum()),
       "top_sites": [{"ctg": c, "pos": int(b * BIN), "reads": int(n)} for c, b, n in top[["ctg", "bin", "n"]].values]}
json.dump(out, open(P + "_recurrent_junctions.json", "w"), indent=1)
print(json.dumps(out, indent=1))
