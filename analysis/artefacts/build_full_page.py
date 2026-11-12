"""Page data for the whole-run artefact scan: summary stats (fullstats.py) + example reads.

usage: build_full_page.py <prefix> reads.bam stats.json out.json [dataset-json]
"""
import sys, json, gzip, random, re, collections, statistics
import pysam

P, BAM, STATS, OUT = sys.argv[1:5]
random.seed(11)
S = json.load(open(STATS))
CATS = ["chimera_interspecies", "chimera_intraspecies", "foldback", "internal_adapter", "end_adapter", "low_complexity", "skip_prone"]
CHIM = {"chimera_interspecies", "chimera_intraspecies", "foldback"}
K = {"chimera_interspecies": 14, "chimera_intraspecies": 14, "foldback": 16, "internal_adapter": 14, "end_adapter": 10,
     "low_complexity": 16, "skip_prone": 12}

# reservoir samples per category + chimera junction stats over all artefact reads
res = {c: [] for c in CATS}; seen = collections.Counter()
gaps = collections.Counter(); shortest = []; fold = []; with_ad = collections.Counter()
n_int_ad = n_int_ad_near = 0; lowc_motif = collections.Counter()
def ok(c, r):
    if not (1500 <= r["len"] <= 30000): return False
    if c == "internal_adapter": return not (set(r["cls"]) & CHIM) and r["segs"]
    if c == "low_complexity": return any(b - a >= 100 for a, b in r["lowc"])
    return True
POOLIDS = None
if BAM.endswith((".fastq.gz", ".fq.gz")):
    POOLIDS = {r.name for r in pysam.FastxFile(BAM)}
for line in gzip.open(P + ".artefacts.jsonl.gz", "rt"):
    r = json.loads(line)
    cls = set(r["cls"])
    for c in CATS:
        if c in cls and ok(c, r) and (POOLIDS is None or r["id"] in POOLIDS):
            seen[c] += 1
            if len(res[c]) < 400: res[c].append(r)
            else:
                j = random.randrange(seen[c])
                if j < 400: res[c][j] = r
    if cls & {"chimera_interspecies", "chimera_intraspecies"}:
        sg = r["segs"]
        for a, b in zip(sg, sg[1:]):
            gaps[max(-50, min(150, (b["q0"] - a["q1"]) // 10 * 10))] += 1
        shortest.append(min(x["q1"] - x["q0"] for x in sg))
    if "foldback" in cls and len(r["segs"]) == 2 and len(fold) < 200000:
        a, b = r["segs"]
        f = lambda s: s["match"] / max(1, s["match"] + s["mm"] + s["ins"] + s["del"])
        fold.append([round(100 * f(a), 2), round(100 * f(b), 2)])
    ads = [a for a in r["adapters"] if 100 < a["start"] < r["len"] - 100]
    for c in CHIM:
        if c in cls and ads: with_ad[c] += 1
    n_int_ad += len(ads); n_int_ad_near += sum(1 for a in ads if a["start"] < 600 or a["start"] > r["len"] - 600)

# pick examples; prefer adapter-containing foldbacks for half of them
picked = {}
for c in CATS:
    pool = res[c]
    if c == "foldback":
        wa = [r for r in pool if any(100 < a["start"] < r["len"] - 100 for a in r["adapters"])]
        na = [r for r in pool if r not in wa]
        picked[c] = random.sample(wa, min(10, len(wa))) + random.sample(na, min(6, len(na)))
    else:
        picked[c] = random.sample(pool, min(K[c], len(pool)))
want = {r["id"] for v in picked.values() for r in v}

seqs = {}
if BAM.endswith((".fastq.gz", ".fq.gz")):
    for r in pysam.FastxFile(BAM):
        if r.name in want: seqs[r.name] = (r.sequence, [ord(c) - 33 for c in r.quality])
else:
    for a in pysam.AlignmentFile(BAM, check_sq=False, threads=4):
        if a.query_name in want:
            seqs[a.query_name] = (a.query_sequence, list(a.query_qualities))
            if len(seqs) == len(want): break

def thin(q, w=50, maxn=300):
    qw = [sum(q[i:i + w]) / len(q[i:i + w]) for i in range(0, len(q), w)]
    if len(qw) <= maxn: return [round(x, 1) for x in qw], w
    f = -(-len(qw) // maxn)
    return [round(sum(qw[i:i + f]) / len(qw[i:i + f]), 1) for i in range(0, len(qw), f)], w * f

def example(r):
    s, q = seqs[r["id"]]
    qw, w = thin(q)
    junctions = []
    for a, b in zip(r["segs"], r["segs"][1:]):
        lo, hi = max(0, a["q1"] - 40), min(len(s), b["q0"] + 40)
        junctions.append({"at": a["q1"], "gap": b["q0"] - a["q1"], "seq": s[lo:hi], "lo": lo, "left_end": a["q1"] - lo, "right_start": b["q0"] - lo})
    ads = [{**a, "ctx": s[max(0, a["start"] - 25):a["end"] + 25], "ctx_lo": max(0, a["start"] - 25)} for a in r["adapters"]]
    segs = []
    for sg in r["segs"]:
        cols = sg["match"] + sg["mm"] + sg["ins"] + sg["del"]
        segs.append({"q0": sg["q0"], "q1": sg["q1"], "sp": sg["sp"], "ctg": sg["ctg"], "strand": sg["strand"], "r0": sg["r0"], "r1": sg["r1"],
                     "mapq": sg["mapq"], "ident": round(100 * sg["match"] / max(1, cols), 2), "indels": [[t, qq, ln] for t, qq, rp, ln in sg["long"] if ln >= 20]})
    lc = [[a, b, s[a:a + 5], s[a:min(b, a + 60)]] for a, b in r["lowc"]][:12]
    return {"id": r["id"], "len": r["len"], "qs": r["qs"], "ch": r["ch"], "st": r["st"], "du": r["du"], "split_from": r["pi"], "cls": r["cls"],
            "segs": segs, "adapters": ads, "qwin": qw, "qw": w, "junctions": junctions, "lowc": lc}

examples = {c: [example(r) for r in v if r["id"] in seqs] for c, v in picked.items()}
random.shuffle(fold)
stats = {
    "n_reads": S["n_reads"], "n_bases": S["n_bases"], "n_channels": S["n_channels"], "n_mapped": S["n_reads"] - S["summary"]["unmapped"]["reads"],
    "dorado_split_reads": S["dorado_split"],
    "chimera_gap_hist": sorted(gaps.items()),
    "chimera_direct_frac": round(sum(v for k, v in gaps.items() if -30 <= k <= 10) / max(1, sum(gaps.values())), 3),
    "chimera_shortest_seg_median": statistics.median(shortest) if shortest else None,
    "chimera_shortest_lt1k_frac": round(sum(1 for x in shortest if x < 1000) / max(1, len(shortest)), 3),
    "chimera_with_adapter": {c: with_ad[c] for c in CHIM},
    "foldback_ident": fold[:4000], "foldback_n": len(fold),
    "foldback_second_worse": round(sum(1 for a, b in fold if b < a) / max(1, len(fold)), 3),
    "internal_adapter_total": n_int_ad, "internal_adapter_near_end": n_int_ad_near,
}
DATASET = {"sample": "ZymoHMW", "flowcell": "PAS01578", "run_date": "2023-08-07", "scope": "full",
           "file": "PAS01578.dorado1.0.0.bm5.2.0_sup.sim.bam", "accession": "ERR15459724",
           "basecaller": "dorado 1.0.0, dna_r10.4.1_e8.2_400bps_sup@v5.2.0 (default adapter trimming)",
           "reference": "ZymoBIOMICS D6322 reference genomes (Zymo Research)"}
if len(sys.argv) > 5: DATASET.update(json.loads(sys.argv[5]))
out = {"dataset": DATASET,
       "summary": S["summary"], "filters": S["filters"], "stats": stats, "timeline": S["timeline"], "sync": S["sync"],
       "deletions": S["deletions"], "systematic_sites": S["systematic_sites"], "n_systematic_sites": S["n_systematic_sites"],
       "mass_end": S["mass_end_seconds"], "mass_end_deletions": S.get("mass_end_deletions"), "channels": S["channels"],
       "duplicate_groups": S["duplicate_groups"], "species_bases": S["species_bases"], "run_hours": S["run_hours"],
       "examples": examples}
json.dump(out, open(OUT, "w"), separators=(",", ":"))
print({c: len(v) for c, v in examples.items()}, stats["chimera_direct_frac"], stats["foldback_second_worse"])
