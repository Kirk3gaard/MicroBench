"""Assemble docs/data/artefacts_zymohmw.json for the artefact page."""
import sys, json, re, random, collections, statistics, datetime, zlib, subprocess
import pysam

JSONL, BAM, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
SYNC = [json.loads(l) for l in open(sys.argv[4])] if len(sys.argv) > 4 else []
random.seed(7)
R = [json.loads(l) for l in open(JSONL)]
byid = {r["id"]: r for r in R}
seqs = {}
for a in pysam.AlignmentFile(BAM, check_sq=False):
    seqs[a.query_name] = a.query_sequence

LOWC = re.compile(r"(?:GGAAT){6,}|(?:ATTCC){6,}|A{25,}|T{25,}|(?:TA){12,}|(?:AT){12,}|(?:AG){12,}|(?:CT){12,}")
CHIM = {"chimera_interspecies", "chimera_intraspecies", "foldback"}

def lowc(seq):
    return [[m.start(), m.end(), m.group()[:10]] for m in LOWC.finditer(seq)]

for r in R:
    s = seqs[r["id"]]
    r["lowc"] = lowc(s)
    # long low-complexity insertions inside alignments
    r["junk_ins"] = 0
    for sg in r["segs"]:
        for t, q, rp, ln in sg["long"]:
            if t == "ins" and ln >= 50:
                part = s[q:q + ln] if sg["strand"] == "+" else s[q - ln:q]
                if len(zlib.compress(part.encode(), 9)) / len(part) < 0.25:
                    r["junk_ins"] += 1
    r["n_skip"] = sum(1 for sg in r["segs"] for e in sg["long"] if e[0] == "del" and e[3] >= 10)
    if r["lowc"] and r["segs"]:
        r["cls"].append("low_complexity")
    if r["n_skip"] >= 5:
        r["cls"].append("skip_prone")
    r["cls"] = [c for c in r["cls"] if c not in ("lowq_stretch", "split_alignment")] or ["clean"]
    if r["cls"] != ["clean"] and "clean" in r["cls"]:
        r["cls"].remove("clean")

N = len(R); B = sum(r["len"] for r in R)
CATS = ["chimera_interspecies", "chimera_intraspecies", "foldback", "internal_adapter", "end_adapter",
        "low_complexity", "skip_prone", "unaligned_head", "unaligned_tail", "unmapped"]
summary = {}
for c in CATS:
    X = [r for r in R if c in r["cls"]]
    summary[c] = {"reads": len(X), "pct_reads": round(100 * len(X) / N, 2),
                  "pct_bases": round(100 * sum(r["len"] for r in X) / B, 2)}
anyart = [r for r in R if set(r["cls"]) & (CHIM | {"internal_adapter", "end_adapter", "low_complexity", "skip_prone"})]
summary["_any"] = {"reads": len(anyart), "pct_reads": round(100 * len(anyart) / N, 2),
                   "pct_bases": round(100 * sum(r["len"] for r in anyart) / B, 2)}

# would the usual length/quality filters remove these reads? (dorado mean qscore tag as proxy)
filters = {}
for c in CATS[:7] + ["clean"]:
    X = [r for r in R if c in r["cls"]]
    if not X: continue
    f = lambda q, l: round(100 * sum(1 for r in X if r["qs"] >= q and r["len"] >= l) / len(X), 1)
    filters[c] = {"q10_l1000": f(10, 1000), "q15_l1000": f(15, 1000), "q20_l5000": f(20, 5000),
                  "median_qs": round(statistics.median(r["qs"] for r in X), 1)}

# junction statistics for chimeras
gaps, shortest, adapter_at_junction = [], [], 0
for r in R:
    if set(r["cls"]) & {"chimera_interspecies", "chimera_intraspecies"}:
        sg = r["segs"]
        for a, b in zip(sg, sg[1:]):
            gaps.append(b["q0"] - a["q1"])
        shortest.append(min(x["q1"] - x["q0"] for x in sg))
fold_ident = []
for r in R:
    if "foldback" in r["cls"] and len(r["segs"]) == 2:
        a, b = r["segs"]
        ia = a["match"] / (a["match"] + a["mm"] + a["ins"] + a["del"])
        ib = b["match"] / (b["match"] + b["mm"] + b["ins"] + b["del"])
        fold_ident.append([round(100 * ia, 2), round(100 * ib, 2)])

stats = {
    "n_reads": N, "n_bases": B, "n_mapped": sum(1 for r in R if r["segs"]),
    "n_channels": len({r["ch"] for r in R}),
    "dorado_split_reads": sum(1 for r in R if r["pi"]),
    "chimera_gap_hist": collections.Counter(max(-50, min(150, g // 10 * 10)) for g in gaps),
    "chimera_direct_frac": round(sum(1 for g in gaps if -30 <= g <= 15) / len(gaps), 3) if gaps else None,
    "chimera_shortest_seg_median": statistics.median(shortest) if shortest else None,
    "chimera_shortest_lt1k_frac": round(sum(1 for x in shortest if x < 1000) / len(shortest), 3) if shortest else None,
    "chimera_with_adapter": {c: sum(1 for r in R if c in r["cls"] and r["n_int_adapt"]) for c in ("chimera_interspecies", "chimera_intraspecies", "foldback")},
    "foldback_ident": fold_ident,
    "internal_adapter_near_end": sum(1 for r in R for a in r["adapters"] if 100 < a["start"] < r["len"] - 100 and (a["start"] < 600 or a["start"] > r["len"] - 600)),
    "internal_adapter_total": sum(1 for r in R for a in r["adapters"] if 100 < a["start"] < r["len"] - 100),
    "junk_insertions": sum(r["junk_ins"] for r in R),
    "lowc_motifs": collections.Counter(m[2][:5] for r in R for m in r["lowc"]).most_common(8),
}
stats["chimera_gap_hist"] = sorted(stats["chimera_gap_hist"].items())

# error rates through the run
ts = lambda s: datetime.datetime.fromisoformat(s).timestamp()
t0 = min(ts(r["st"]) for r in R)
bins = collections.defaultdict(list)
for r in R:
    if r["aln"] and r["len"] >= 1000:
        bins[int((ts(r["st"]) - t0) // 3600)].append(r)
timeline = []
for b in sorted(bins):
    X = bins[b]
    if len(X) < 50: continue
    timeline.append({"h": b, "n": len(X), "channels": len({r["ch"] for r in X}),
                     "del": round(1000 * statistics.median(r["del"] / r["aln"] for r in X), 3),
                     "ins": round(1000 * statistics.median(r["ins"] / r["aln"] for r in X), 3),
                     "mm": round(1000 * statistics.median(r["mm"] / r["aln"] for r in X), 3),
                     "ident": round(100 * statistics.median(r["ident"] for r in X), 3),
                     "chim": round(100 * sum(1 for r in X if set(r["cls"]) & CHIM) / len(X), 2)})

# raster of long-deletion events (time x channel) for the synchrony figure
raster = []
for r in R:
    st = ts(r["st"]) - t0
    for sg in r["segs"]:
        for t, q, rp, ln in sg["long"]:
            if t == "del" and ln >= 20:
                raster.append([round((st + r["du"] * q / r["len"]) / 3600, 4), r["ch"], ln])

# examples
def pick(cat, k=14, extra=lambda r: True):
    X = [r for r in R if cat in r["cls"] and 1500 <= r["len"] <= 30000 and extra(r)]
    if len(X) < k:
        X = [r for r in R if cat in r["cls"] and extra(r)]
    return random.sample(X, min(k, len(X)))

def thin(qw, maxn=300):
    if len(qw) <= maxn: return qw, 50
    f = -(-len(qw) // maxn)
    return [round(sum(qw[i:i + f]) / len(qw[i:i + f]), 1) for i in range(0, len(qw), f)], 50 * f

def example(r):
    s = seqs[r["id"]]
    qw, w = thin(r["qwin"])
    junctions = []
    for a, b in zip(r["segs"], r["segs"][1:]):
        lo, hi = max(0, a["q1"] - 40), min(len(s), b["q0"] + 40)
        junctions.append({"at": a["q1"], "gap": b["q0"] - a["q1"], "seq": s[lo:hi], "lo": lo,
                          "left_end": a["q1"] - lo, "right_start": b["q0"] - lo})
    ads = [{**a, "ctx": s[max(0, a["start"] - 25):a["end"] + 25], "ctx_lo": max(0, a["start"] - 25)} for a in r["adapters"]]
    segs = []
    for sg in r["segs"]:
        cols = sg["match"] + sg["mm"] + sg["ins"] + sg["del"]
        segs.append({"q0": sg["q0"], "q1": sg["q1"], "sp": sg["sp"], "ctg": sg["ctg"], "strand": sg["strand"],
                     "r0": sg["r0"], "r1": sg["r1"], "mapq": sg["mapq"], "ident": round(100 * sg["match"] / cols, 2),
                     "indels": [[t, q, ln] for t, q, rp, ln in sg["long"] if ln >= 20]})
    lc = [[a, b2, m, s[a:min(b2, a + 60)]] for a, b2, m in r["lowc"]][:12]
    return {"id": r["id"], "len": r["len"], "qs": r["qs"], "ch": r["ch"], "st": r["st"], "du": r["du"],
            "split_from": r["pi"], "cls": r["cls"], "segs": segs, "adapters": ads, "qwin": qw, "qw": w,
            "junctions": junctions, "lowc": lc,
            "head": s[:min(len(s), 120)], "tail": s[-120:]}

examples = {
    "chimera_interspecies": [example(r) for r in pick("chimera_interspecies")],
    "chimera_intraspecies": [example(r) for r in pick("chimera_intraspecies")],
    "foldback": [example(r) for r in pick("foldback", 10, lambda r: r["n_int_adapt"] > 0)] +
                [example(r) for r in pick("foldback", 6, lambda r: r["n_int_adapt"] == 0)],
    "internal_adapter": [example(r) for r in pick("internal_adapter", 14, lambda r: not set(r["cls"]) & CHIM and r["segs"])],
    "end_adapter": [example(r) for r in pick("end_adapter", 10)],
    "low_complexity": [example(r) for r in pick("low_complexity", 16, lambda r: any(b - a >= 100 for a, b, m in r["lowc"]))],
    "skip_prone": [example(r) for r in pick("skip_prone", 12)],
    "clean": [example(r) for r in pick("clean", 6)],
}

out = {"dataset": {"sample": "ZymoHMW", "flowcell": "PAS01578", "run_date": "2023-08-07",
                   "file": "PAS01578.dorado1.0.0.bm5.2.0_sup.sim.bam", "accession": "ERR15459724",
                   "basecaller": "dorado 1.0.0, dna_r10.4.1_e8.2_400bps_sup@v5.2.0 (default adapter trimming)",
                   "reference": "ZymoBIOMICS D6322 reference genomes (Zymo Research)"},
       "summary": summary, "stats": stats, "filters": filters, "timeline": timeline, "sync": SYNC,
       "raster": raster, "examples": examples}
json.dump(out, open(OUT, "w"), separators=(",", ":"))
print(json.dumps(summary, indent=0))
print({k: v for k, v in stats.items() if k not in ("foldback_ident",)})
