"""Classify read-level artefacts in a sampled ONT BAM against reference genomes.

Output: per-read TSV + JSON with segment layouts for the web page.
"""
import sys, json, re, math, collections, statistics
import pysam, mappy, edlib

BAM, REF, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
MAX_READS = int(sys.argv[4]) if len(sys.argv) > 4 else 10**9

# ONT ligation adapter (LA) as listed for SQK-LSK114 by dorado/porechop
ADAPTERS = {
    "LA_top": "CCTGTACTTCGTTCAGTTACGTATTGCT",
    "LA_bottom": "AGCAATACGTAACTGAACGAAGTACAGG",
}

def rc(s):
    return s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))

def find_adapters(seq, max_frac=0.2):
    hits = []
    for name, a in ADAPTERS.items():
        k = int(len(a) * max_frac)
        for strand, q in (("+", a), ("-", rc(a))):
            r = edlib.align(q, seq, mode="HW", task="locations", k=k)
            for (s, e) in r["locations"]:
                hits.append({"name": name, "strand": strand, "start": s, "end": e + 1, "ed": r["editDistance"]})
    # merge overlapping hits, keep best
    hits.sort(key=lambda h: (h["start"], h["ed"]))
    merged = []
    for h in hits:
        if merged and h["start"] < merged[-1]["end"]:
            if h["ed"] < merged[-1]["ed"]: merged[-1] = h
            continue
        merged.append(h)
    return merged

contig_genome = {}
for name, seq, _ in mappy.fastx_read(REF):
    contig_genome[name] = name
aligner = mappy.Aligner(REF, preset="map-ont", best_n=5)

# map contig -> species label from the reference FASTA headers file names
species = {}
for line in open(REF + ".species"):
    c, sp = line.rstrip("\n").split("\t")
    species[c] = sp

def cs_stats(cs):
    """From a long-form/short cs string, count matches, mismatches, ins, del and list long indels."""
    m = mm = ins = dele = 0
    long_ev = []  # (type, query_offset_in_aln, ref_offset, length)
    qo = ro = 0
    for op, val in re.findall(r"(:|\*|\+|\-|~)([0-9a-z]+)", cs):
        if op == ":":
            n = int(val); m += n; qo += n; ro += n
        elif op == "*":
            mm += 1; qo += 1; ro += 1
        elif op == "+":
            n = len(val); ins += n
            if n >= 10: long_ev.append(("ins", qo, ro, n))
            qo += n
        elif op == "-":
            n = len(val); dele += n
            if n >= 10: long_ev.append(("del", qo, ro, n))
            ro += n
    return m, mm, ins, dele, long_ev

ctg_len = {}
for name, seq, _ in mappy.fastx_read(REF):
    ctg_len[name] = len(seq)

def junction_kind(a, b):
    """Classify adjacent query segments a (left) and b (right)."""
    if a["ctg"] != b["ctg"]:
        # different contigs of a draft genome (yeast) or chromosome vs plasmid: not decidable
        return "ambiguous" if a["sp"] == b["sp"] else "chimera"
    qgap = b["q0"] - a["q1"]
    if a["strand"] != b["strand"]:
        ov = min(a["r1"], b["r1"]) - max(a["r0"], b["r0"])
        if ov > 0.3 * min(a["r1"] - a["r0"], b["r1"] - b["r0"]) or abs(a["r1"] - b["r1"]) < 2000 or abs(a["r0"] - b["r0"]) < 2000:
            return "foldback"
        return "chimera"
    # same strand: expected ref gap equals query gap
    clen = ctg_len[a["ctg"]]
    if a["strand"] == "+":
        rgap = b["r0"] - a["r1"]
    else:
        rgap = a["r0"] - b["r1"]
    alt = rgap + clen if rgap < 0 else rgap - clen  # across the origin of a circular contig
    d = min(abs(rgap - qgap), abs(alt - qgap))
    return "collinear" if d < 2000 else "chimera"

def windows_q(quals, w=50):
    return [round(sum(quals[i:i+w]) / len(quals[i:i+w]), 1) for i in range(0, len(quals), w)]

bam = pysam.AlignmentFile(BAM, check_sq=False)
rows, examples = [], []
for n, r in enumerate(bam):
    if n >= MAX_READS: break
    seq = r.query_sequence
    if not seq: continue
    L = len(seq)
    q = r.query_qualities
    tags = dict(r.get_tags())
    hits = []
    for h in aligner.map(seq, cs=True):
        if h.mapq < 1 and not h.is_primary: continue
        m, mm, ins, dele, long_ev = cs_stats(h.cs)
        hits.append({
            "q0": h.q_st, "q1": h.q_en, "ctg": h.ctg, "sp": species.get(h.ctg, h.ctg),
            "r0": h.r_st, "r1": h.r_en, "strand": "+" if h.strand == 1 else "-",
            "mapq": h.mapq, "primary": h.is_primary,
            "match": m, "mm": mm, "ins": ins, "del": dele,
            "long": [(t, (h.q_st + qo) if h.strand == 1 else (h.q_en - qo), h.r_st + ro, ln) for t, qo, ro, ln in long_ev],
        })
    # keep non-redundant segments along the query (greedy by aligned length)
    hits.sort(key=lambda h: -(h["q1"] - h["q0"]))
    segs = []
    for h in hits:
        ov = sum(max(0, min(h["q1"], s["q1"]) - max(h["q0"], s["q0"])) for s in segs)
        if ov < 0.2 * (h["q1"] - h["q0"]) and h["q1"] - h["q0"] >= 300:
            segs.append(h)
    segs.sort(key=lambda h: h["q0"])
    adapters = find_adapters(seq)
    covered = sum(s["q1"] - s["q0"] for s in segs)
    # classification
    cls = []
    if not segs:
        cls.append("unmapped")
    if len(segs) >= 2:
        sps = {s["sp"] for s in segs}
        if len(sps) > 1:
            cls.append("chimera_interspecies")
        else:
            # same species: foldback if same contig, opposite strand, overlapping ref span
            fold = False; intra = False; amb = False
            for a, b in zip(segs, segs[1:]):
                kind = junction_kind(a, b)
                if kind == "foldback": fold = True
                elif kind == "chimera": intra = True
                elif kind == "ambiguous": amb = True
            if fold: cls.append("foldback")
            if intra: cls.append("chimera_intraspecies")
            if amb and not (fold or intra): cls.append("split_alignment")
    internal_ad = [a for a in adapters if 100 < a["start"] < L - 100]
    end_ad = [a for a in adapters if a["start"] <= 100 or a["start"] >= L - 100]
    if internal_ad: cls.append("internal_adapter")
    if end_ad: cls.append("end_adapter")
    if segs:
        head = segs[0]["q0"]; tail = L - segs[-1]["q1"]
        if head >= 200: cls.append("unaligned_head")
        if tail >= 200: cls.append("unaligned_tail")
    longdel = sum(1 for s in segs for e in s["long"] if e[0] == "del" and e[3] >= 20)
    longins = sum(1 for s in segs for e in s["long"] if e[0] == "ins" and e[3] >= 20)
    tot = lambda k: sum(s[k] for s in segs)
    aln_cols = tot("match") + tot("mm") + tot("ins") + tot("del")
    ident = tot("match") / aln_cols if aln_cols else None
    qw = windows_q(q, 50)
    lowq_win = sum(1 for x in qw[1:-1] if x < 8)
    if lowq_win: cls.append("lowq_stretch")
    if not cls: cls.append("clean")
    row = {
        "id": r.query_name, "len": L, "qs": round(tags.get("qs", 0), 2), "ch": tags.get("ch"), "st": tags.get("st"),
        "mx": tags.get("mx"), "du": round(tags.get("du", 0), 2), "pi": tags.get("pi"), "ts": tags.get("ts"),
        "nseg": len(segs), "cov": round(covered / L, 3), "ident": round(ident, 4) if ident else None,
        "mm": tot("mm"), "ins": tot("ins"), "del": tot("del"), "aln": aln_cols,
        "longdel": longdel, "longins": longins, "n_adapt": len(adapters),
        "n_int_adapt": len(internal_ad), "lowq_win": lowq_win, "cls": cls,
    }
    rows.append(row)
    examples.append({**row, "segs": [{k: s[k] for k in ("q0", "q1", "ctg", "sp", "r0", "r1", "strand", "mapq", "match", "mm", "ins", "del", "long")} for s in segs],
                     "adapters": adapters, "qwin": qw,
                     "head": seq[:segs[0]["q0"]][-400:] if segs and segs[0]["q0"] > 0 else "",
                     "tail": seq[segs[-1]["q1"]:][:400] if segs and segs[-1]["q1"] < L else ""})
    if n % 2000 == 0:
        print(n, file=sys.stderr, flush=True)

with open(OUT + ".reads.jsonl", "w") as fh:
    for e in examples:
        fh.write(json.dumps(e) + "\n")
c = collections.Counter(x for r in rows for x in r["cls"])
print(len(rows), "reads"); print(c.most_common())
