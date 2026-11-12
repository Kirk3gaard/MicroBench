"""Recompute the headline numbers of the artefact page directly from the per-read table, with
every read counted once, so that totals over overlapping classes are unions, not sums.

usage: count_numbers.py <prefix> out.json [page.json]
If the page JSON is given, its per-class numbers are checked against the recomputed ones.
"""
import sys, json, os
import pandas as pd

P, OUT = sys.argv[1], sys.argv[2]
R = pd.read_csv(P + ".reads.tsv.gz", sep="\t", usecols=["id", "len", "split", "n_int_adapt", "cls"])
N, B = len(R), int(R["len"].sum())
has = lambda c: R.cls.str.contains(c, regex=False)
C = {c: has(c) for c in ["chimera_interspecies", "chimera_intraspecies", "foldback", "internal_adapter", "end_adapter",
                         "low_complexity", "skip_prone", "unaligned_head", "unaligned_tail", "unmapped"]}
gap_ids = set(open(P + ".gap_reads.txt").read().split()) if os.path.exists(P + ".gap_reads.txt") else set()
C["signal_gap"] = R.id.isin(gap_ids)
def stat(m):
    return {"reads": int(m.sum()), "pct_reads": round(100 * float(m.mean()), 3), "pct_bases": round(100 * float(R.loc[m, "len"].sum()) / B, 3)}
chim = C["chimera_interspecies"] | C["chimera_intraspecies"]
chim_fold = chim | C["foldback"]
adapters = C["internal_adapter"] | C["end_adapter"]
core = chim_fold | adapters | C["low_complexity"] | C["skip_prone"]
out = {"n_reads": N, "n_bases": B,
       "classes": {k: stat(v) for k, v in C.items()},
       "chimeric": stat(chim), "chimeric_or_foldback": stat(chim_fold), "adapter_any": stat(adapters),
       "any_core": stat(core), "any_core_or_gap": stat(core | C["signal_gap"]),
       "split_pieces": stat(R.split == 1),
       "chimeric_with_internal_adapter": int((chim & (R.n_int_adapt > 0)).sum()),
       "foldback_with_internal_adapter": int((C["foldback"] & (R.n_int_adapt > 0)).sum()),
       "overlap_note": {"sum_of_class_reads": int(sum(int(C[k].sum()) for k in ["chimera_interspecies", "chimera_intraspecies", "foldback", "internal_adapter", "end_adapter", "low_complexity", "skip_prone"])),
                        "union_reads": int(core.sum())}}
problems = []
if len(sys.argv) > 3:
    page = json.load(open(sys.argv[3]))
    for k, v in out["classes"].items():
        pv = page.get("summary", {}).get(k)
        if pv and (pv["reads"] != v["reads"] or abs(pv["pct_bases"] - v["pct_bases"]) > 0.01):
            problems.append(f"{k}: page {pv} vs recomputed {v}")
    if page.get("stats", {}).get("n_reads") != N: problems.append("n_reads differs")
    if page.get("gaps") and page["gaps"]["reads_hit"] != out["classes"]["signal_gap"]["reads"]:
        problems.append(f"gap reads: page {page['gaps']['reads_hit']} vs {out['classes']['signal_gap']['reads']}")
out["check_problems"] = problems
json.dump(out, open(OUT, "w"), indent=1)
print(json.dumps({k: out[k] for k in ("chimeric_or_foldback", "any_core", "any_core_or_gap", "split_pieces", "overlap_note", "check_problems")}, indent=1))
