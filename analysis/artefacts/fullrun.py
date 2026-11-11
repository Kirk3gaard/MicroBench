"""Whole-run artefact scan: align every read of an unaligned dorado BAM to the mock references and
write compact per-read rows, an indel-event table and artefact details.

usage: fullrun.py reads.bam ref.fa outprefix [threads]
outputs: <out>.reads.tsv.gz  one row per read
         <out>.events.tsv.gz indels >=10 bp inside alignments with estimated wall-clock time
         <out>.artefacts.jsonl.gz  full layout for reads with any artefact (for examples)
"""
import sys, os, re, json, gzip, zlib, datetime, multiprocessing as mp
import pysam, mappy, edlib

BAM, REF, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
THREADS = int(sys.argv[4]) if len(sys.argv) > 4 else max(1, os.cpu_count() - 1)
BATCH = 2000
ADAPTERS = {"LA_top": "CCTGTACTTCGTTCAGTTACGTATTGCT", "LA_bottom": "AGCAATACGTAACTGAACGAAGTACAGG"}
LOWC = re.compile(r"(?:GGAAT){6,}|(?:ATTCC){6,}|A{25,}|T{25,}|(?:TA){12,}|(?:AT){12,}|(?:AG){12,}|(?:CT){12,}|(?:GA){12,}|(?:TC){12,}|C{25,}|G{25,}")
COMP = str.maketrans("ACGTN", "TGCAN")
rc = lambda s: s.translate(COMP)[::-1]

def load_species(ref):
    sp = {}
    for line in open(ref + ".species"):
        c, s = line.rstrip("\n").split("\t"); sp[c] = s
    return sp

A = None; SPECIES = None; CTGLEN = None
def init():
    global A, SPECIES, CTGLEN
    A = mappy.Aligner(REF, preset="map-ont", best_n=5)
    SPECIES = load_species(REF)
    CTGLEN = {n: len(s) for n, s, _ in mappy.fastx_read(REF)}

def find_adapters(seq):
    hits = []
    for name, a in ADAPTERS.items():
        k = int(len(a) * 0.2)
        for strand, q in (("+", a), ("-", rc(a))):
            r = edlib.align(q, seq, mode="HW", task="locations", k=k)
            for s, e in r["locations"]:
                hits.append([name, strand, s, e + 1, r["editDistance"]])
    hits.sort(key=lambda h: (h[2], h[4]))
    out = []
    for h in hits:
        if out and h[2] < out[-1][3]:
            if h[4] < out[-1][4]: out[-1] = h
            continue
        out.append(h)
    return out

CS = re.compile(r"(:|\*|\+|\-|~)([0-9a-z]+)")
def cs_stats(cs):
    m = mm = ins = dele = 0; ev = []; qo = ro = 0
    for op, val in CS.findall(cs):
        if op == ":":
            n = int(val); m += n; qo += n; ro += n
        elif op == "*":
            mm += 1; qo += 1; ro += 1
        elif op == "+":
            n = len(val); ins += n
            if n >= 10: ev.append(("ins", qo, ro, n))
            qo += n
        elif op == "-":
            n = len(val); dele += n
            if n >= 10: ev.append(("del", qo, ro, n))
            ro += n
    return m, mm, ins, dele, ev

def junction_kind(a, b):
    if a["ctg"] != b["ctg"]:
        return "ambiguous" if a["sp"] == b["sp"] else "chimera"
    qgap = b["q0"] - a["q1"]
    if a["strand"] != b["strand"]:
        ov = min(a["r1"], b["r1"]) - max(a["r0"], b["r0"])
        if ov > 0.3 * min(a["r1"] - a["r0"], b["r1"] - b["r0"]) or abs(a["r1"] - b["r1"]) < 2000 or abs(a["r0"] - b["r0"]) < 2000:
            return "foldback"
        return "chimera"
    clen = CTGLEN[a["ctg"]]
    rgap = b["r0"] - a["r1"] if a["strand"] == "+" else a["r0"] - b["r1"]
    alt = rgap + clen if rgap < 0 else rgap - clen
    return "collinear" if min(abs(rgap - qgap), abs(alt - qgap)) < 2000 else "chimera"

def ts_epoch(s):
    return datetime.datetime.fromisoformat(s.replace("Z", "+00:00")).timestamp()

def process(batch):
    rows, events, details = [], [], []
    for name, seq, qual, tags in batch:
        L = len(seq)
        hits = []
        for h in A.map(seq, cs=True):
            if h.mapq < 1 and not h.is_primary: continue
            m, mm, ins, dele, ev = cs_stats(h.cs)
            hits.append({"q0": h.q_st, "q1": h.q_en, "ctg": h.ctg, "sp": SPECIES.get(h.ctg, h.ctg), "r0": h.r_st, "r1": h.r_en,
                         "strand": "+" if h.strand == 1 else "-", "mapq": h.mapq, "match": m, "mm": mm, "ins": ins, "del": dele,
                         "long": [(t, (h.q_st + qo) if h.strand == 1 else (h.q_en - qo), h.r_st + ro, ln) for t, qo, ro, ln in ev]})
        hits.sort(key=lambda h: -(h["q1"] - h["q0"]))
        segs = []
        for h in hits:
            ov = sum(max(0, min(h["q1"], s["q1"]) - max(h["q0"], s["q0"])) for s in segs)
            if ov < 0.2 * (h["q1"] - h["q0"]) and h["q1"] - h["q0"] >= 300:
                segs.append(h)
        segs.sort(key=lambda h: h["q0"])
        ads = find_adapters(seq)
        lowc = [(m_.start(), m_.end()) for m_ in LOWC.finditer(seq)]
        cls = []
        if not segs: cls.append("unmapped")
        if len(segs) >= 2:
            kinds = [junction_kind(a, b) for a, b in zip(segs, segs[1:])]
            if any(a["sp"] != b["sp"] for a, b in zip(segs, segs[1:])): cls.append("chimera_interspecies")
            if "foldback" in kinds: cls.append("foldback")
            if any(k == "chimera" for k, a, b in zip(kinds, segs, segs[1:]) if a["sp"] == b["sp"]): cls.append("chimera_intraspecies")
        int_ad = [a for a in ads if 100 < a[2] < L - 100]
        end_ad = [a for a in ads if a[2] <= 100 or a[2] >= L - 100]
        if int_ad: cls.append("internal_adapter")
        if end_ad: cls.append("end_adapter")
        junk = 0
        for sg in segs:
            for t, q, rp, ln in sg["long"]:
                if t == "ins" and ln >= 50:
                    part = seq[q:q + ln] if sg["strand"] == "+" else seq[q - ln:q]
                    if part and len(zlib.compress(part.encode(), 9)) / len(part) < 0.25: junk += 1
        if lowc and segs: cls.append("low_complexity")
        nskip = sum(1 for sg in segs for e in sg["long"] if e[0] == "del")
        if nskip >= 5: cls.append("skip_prone")
        if segs:
            if segs[0]["q0"] >= 200: cls.append("unaligned_head")
            if L - segs[-1]["q1"] >= 200: cls.append("unaligned_tail")
        tot = lambda k: sum(s[k] for s in segs)
        aln = tot("match") + tot("mm") + tot("ins") + tot("del")
        st = ts_epoch(tags["st"]) if "st" in tags else 0.0
        du = float(tags.get("du", 0)); ns = int(tags.get("ns", 0)); tsamp = int(tags.get("ts", 0))
        sr = ns / du if du else 5000.0
        qs = float(tags.get("qs", 0))
        lowc_bp = sum(b - a for a, b in lowc)
        best = max(segs, key=lambda s: s["q1"] - s["q0"]) if segs else None
        primary_sp = best["sp"] if best else ""
        loc = [best["ctg"], best["r0"], best["r1"], best["strand"], segs[0]["q0"], segs[-1]["q1"]] if best else ["", "", "", "", "", ""]
        rows.append("\t".join(map(str, [name, L, round(qs, 2), tags.get("ch", ""), round(st, 3), round(du, 3), tags.get("mx", ""),
                                         1 if "pi" in tags else 0, len(segs), primary_sp, tot("match"), tot("mm"), tot("ins"), tot("del"), aln,
                                         len(ads), len(int_ad), lowc_bp, junk, nskip, *loc, ",".join(cls) or "clean"])))
        tq = lambda q: st + (tsamp + (ns - tsamp) * q / L) / sr if ns else st + du * q / L
        for sg in segs:
            for t, q, rp, ln in sg["long"]:
                events.append("\t".join(map(str, [name, tags.get("ch", ""), round(tq(q), 3), t, ln, sg["ctg"], rp, q, L])))
        # junctions between adjacent alignments: collinear splits (a stretch of the molecule missing
        # or inserted) and chimera/foldback joins, with their time
        for a_, b_ in zip(segs, segs[1:]):
            kind = junction_kind(a_, b_)
            q = a_["q1"]
            if kind == "collinear":
                qgap = b_["q0"] - a_["q1"]
                rgap = b_["r0"] - a_["r1"] if a_["strand"] == "+" else a_["r0"] - b_["r1"]
                clen = CTGLEN[a_["ctg"]]
                if abs(rgap) > clen / 2: rgap = rgap - clen if rgap > 0 else rgap + clen
                miss = rgap - qgap
                rp = a_["r1"] if a_["strand"] == "+" else a_["r0"]
                events.append("\t".join(map(str, [name, tags.get("ch", ""), round(tq(q), 3), "split", miss, a_["ctg"], rp, q, L])))
            else:
                events.append("\t".join(map(str, [name, tags.get("ch", ""), round(tq(q), 3), "junction_" + kind, b_["q0"] - a_["q1"], a_["ctg"], a_["r1"], q, L])))
        if cls and cls != ["unmapped"]:
            details.append(json.dumps({"id": name, "len": L, "qs": round(qs, 2), "ch": tags.get("ch"), "st": tags.get("st"), "du": du,
                                       "pi": tags.get("pi"), "cls": cls,
                                       "segs": [{k: s[k] for k in ("q0", "q1", "ctg", "sp", "r0", "r1", "strand", "mapq", "match", "mm", "ins", "del", "long")} for s in segs],
                                       "adapters": [{"name": a[0], "strand": a[1], "start": a[2], "end": a[3], "ed": a[4]} for a in ads],
                                       "lowc": lowc[:20]}, separators=(",", ":")))
    return rows, events, details

def batches():
    bam = pysam.AlignmentFile(BAM, check_sq=False, threads=2)
    b = []
    for r in bam:
        seq = r.query_sequence
        if not seq: continue
        tags = {k: v for k, v in r.get_tags() if k in ("st", "du", "ns", "ts", "qs", "ch", "mx", "pi")}
        b.append((r.query_name, seq, None, tags))
        if len(b) >= BATCH:
            yield b; b = []
    if b: yield b

if __name__ == "__main__":
    hdr = "id\tlen\tqs\tch\tst\tdu\tmx\tsplit\tnseg\tsp\tmatch\tmm\tins\tdel\taln\tn_adapt\tn_int_adapt\tlowc_bp\tjunk_ins\tn_del10\tctg\tr0\tr1\tstrand\taln_q0\taln_q1\tcls\n"
    with gzip.open(OUT + ".reads.tsv.gz", "wt", compresslevel=3) as fr, \
         gzip.open(OUT + ".events.tsv.gz", "wt", compresslevel=3) as fe, \
         gzip.open(OUT + ".artefacts.jsonl.gz", "wt", compresslevel=3) as fd, \
         mp.Pool(THREADS, initializer=init) as pool:
        fr.write(hdr)
        fe.write("id\tch\tt\ttype\tlen\tctg\trpos\tqpos\treadlen\n")
        n = 0
        for rows, events, details in pool.imap(process, batches(), chunksize=1):
            fr.write("\n".join(rows) + "\n")
            if events: fe.write("\n".join(events) + "\n")
            if details: fd.write("\n".join(details) + "\n")
            n += len(rows)
            if n % 100000 < BATCH:
                print(datetime.datetime.now().isoformat(timespec="seconds"), n, "reads", flush=True)
    print("done", n, flush=True)
