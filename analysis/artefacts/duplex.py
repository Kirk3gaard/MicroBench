"""Duplex follow-ons: the complement strand of a molecule read right after its template strand in the
same pore. Both reads stay in the simplex data, so the molecule is counted twice.

usage: duplex.py <prefix> out.json
"""
import sys, json
import numpy as np, pandas as pd

P, OUT = sys.argv[1], sys.argv[2]
R = pd.read_csv(P + ".reads.tsv.gz", sep="\t", usecols=["id", "len", "ch", "st", "du", "split", "ctg", "r0", "r1", "strand"])
R = R.dropna(subset=["ctg"]).sort_values(["ch", "st"]).reset_index(drop=True)
nxt = R.groupby("ch").shift(-1)
gap = nxt.st - (R.st + R.du)
ov = (np.minimum(R.r1, nxt.r1) - np.maximum(R.r0, nxt.r0)) / np.minimum(R.r1 - R.r0, nxt.r1 - nxt.r0)
pair = (R.ctg == nxt.ctg) & (R.strand != nxt.strand) & (ov > 0.8) & (gap >= -1) & (gap < 20)
n_mapped_total = len(pd.read_csv(P + ".reads.tsv.gz", sep="\t", usecols=["id"]))
from_split = pair & (R.split == 1) & (nxt.split == 1) & (gap.abs() < 0.05)
out = {"from_one_split_read": int(from_split.sum()), "separate_captures": int((pair & ~from_split).sum()),
       "pairs": int(pair.sum()),
       "reads_in_pairs": int(len(set(R.id[pair]) | set(nxt.id[pair]))),
       "pct_reads_in_pairs": round(100 * len(set(R.id[pair]) | set(nxt.id[pair])) / n_mapped_total, 2),
       "second_reads": int(nxt.id[pair].nunique()),
       "pct_reads_second": round(100 * nxt.id[pair].nunique() / n_mapped_total, 2),
       "bases_in_second_reads_pct": round(100 * nxt.len[pair].sum() / R.len.sum(), 2),
       "gap_s_median": round(float(gap[pair].median()), 2) if pair.any() else None}
json.dump(out, open(OUT, "w")); print(out)
