#!/usr/bin/env python3
"""Infer the ONT library kit from the reads themselves.

dorado trims the sequencing adapter by default, but without --kit-name it leaves barcodes and
kit primers in the reads. For each BAM in ENA_data_links.tsv this samples reads from a few
evenly spaced slices (HTTP range requests) and looks at the first and last 150 bases for:
  native barcode construct   ATTGCTAAGGTTAA + NBxx + CAGCACCT       (SQK-NBD114-24/-96)
  rapid barcode construct    TATTGCT + BCxx + GTTTTCGCATTTATCG...   (SQK-RBK114-24/-96)
  rapid adapter primer       GTTTTCGCATTTATCGTGAAACGCTTTCGCG...     (SQK-RAD114 / ULK114)
  ligation adapter remnant   CCTGTACTTCGTTCAGTTACGTATTGC            (any ligation kit)
Kit sequences are from dorado (dorado/utils/barcode_kits.cpp, dorado/demux/adapter_primer_kits.cpp).

usage: kit_detect.py [--slices 4] [--chunk 2000000] out.json
"""
import sys, os, json, csv, collections, argparse
from concurrent.futures import ThreadPoolExecutor
import edlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bamsample as bs
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "analysis", "artefacts"))

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NB_FRONT, NB_REAR = "ATTGCTAAGGTTAA", "CAGCACCT"
RB_FRONT, RB_REAR = "TATTGCT", "GTTTTCGCATTTATCGTGAAACGCTTTCG"
RAD = "GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA"
LA = "CCTGTACTTCGTTCAGTTACGTATTGC"
COMP = str.maketrans("ACGTN", "TGCAN")
rc = lambda s: s.translate(COMP)[::-1]
BC = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "ont_barcodes.json")))
NB = {k: v for k, v in BC.items() if k.startswith("NB")}
RB = {k: v for k, v in BC.items() if k.startswith("BC") or k.startswith("RBK")}

def best(constructs, seq, k):
    top = (None, 99)
    for name, c in constructs.items():
        r = edlib.align(c, seq, mode="HW", task="distance", k=k)
        if r["editDistance"] != -1 and r["editDistance"] < top[1]:
            top = (name, r["editDistance"])
    return top

NBC = {n: NB_FRONT + s + NB_REAR for n, s in NB.items()}
RBC = {n: RB_FRONT + s + RB_REAR for n, s in RB.items()}

def classify(seq):
    head, tail = seq[:150], rc(seq[-150:])
    res = {}
    for end, s in (("start", head), ("end", tail)):
        nb = best(NBC, s, 9)
        rbk = best(RBC, s, 10) if end == "start" else (None, 99)
        res[end] = {"nb": nb, "rb": rbk,
                    "rad": edlib.align(RAD, s, mode="HW", task="distance", k=10)["editDistance"] != -1,
                    "la": edlib.align(LA, s, mode="HW", task="distance", k=5)["editDistance"] != -1}
    return res

def records_raw(data):
    return list(bs.parse_records(data))

def sample_seqs(url, size, slices, chunk):
    import struct
    seqs = []
    starts = [int(size * (k + 0.5) / slices) for k in range(slices)]
    for st in starts:
        data = bs.inflate_bgzf(bs.get_range(url, st, min(size - 1, st + chunk)))
        # walk records and decode sequence
        for m in bs.UUID.finditer(data):
            s = m.start() - 36
            if s >= 0 and data[s + 12] == 37 and struct.unpack("<i", data[s + 4:s + 8])[0] == -1:
                break
        else:
            continue
        i = s
        while i + 36 <= len(data):
            bsz, ref = struct.unpack("<ii", data[i:i + 8])
            if ref != -1 or bsz <= 0 or i + 4 + bsz > len(data): break
            lrn = data[i + 12]; ncig = struct.unpack("<H", data[i + 16:i + 18])[0]
            lseq = struct.unpack("<i", data[i + 20:i + 24])[0]
            p = i + 36 + lrn + 4 * ncig
            packed = data[p:p + (lseq + 1) // 2]
            seq = "".join("=ACMGRSVTWYHKDBN"[b >> 4] + "=ACMGRSVTWYHKDBN"[b & 15] for b in packed)[:lseq]
            if lseq >= 500: seqs.append(seq)
            i += 4 + bsz
    return seqs

def summarise(seqs):
    n = len(seqs); c = collections.Counter(); bcs = collections.Counter()
    for s in seqs:
        r = classify(s)
        nb_s, nb_e = r["start"]["nb"], r["end"]["nb"]
        if nb_s[1] <= 9: c["nb_start"] += 1; bcs[nb_s[0]] += 1
        if nb_e[1] <= 9: c["nb_end"] += 1; bcs[nb_e[0]] += 1
        if r["start"]["rb"][1] <= 10: c["rb_start"] += 1; bcs[r["start"]["rb"][0]] += 1
        if r["start"]["rad"]: c["rad_start"] += 1
        if r["start"]["la"] or r["end"]["la"]: c["la_end"] += 1
    f = {k: round(v / n, 3) for k, v in c.items()} if n else {}
    top = bcs.most_common(3)
    # a barcode counts if it dominates the barcode hits and was seen in enough reads
    barcode = top[0][0] if top and top[0][1] >= max(5, 0.02 * n) and top[0][1] >= 0.8 * sum(bcs.values()) else None
    nb, rb, rad = f.get("nb_start", 0), f.get("rb_start", 0), f.get("rad_start", 0)
    if max(nb, f.get("nb_end", 0)) > 0.2:
        kit, conf = "Native barcoding (SQK-NBD114-24/96)", "high"
        if barcode and int(barcode[2:]) > 24: kit = "Native barcoding (SQK-NBD114-96)"
        retained = "barcode still on {:.0%} of read starts and {:.0%} of read ends".format(nb, f.get("nb_end", 0))
    elif rb + rad >= 0.08 or (rb + rad >= 0.03 and barcode):
        kit, conf = "Rapid barcoding (SQK-RBK114-24/96)" if barcode else "Rapid (SQK-RAD114/RBK114)", "high" if rb + rad >= 0.08 else "medium"
        retained = "rapid barcode/primer left on {:.0%} of read starts".format(max(rb, rad))
    elif rb + rad >= 0.015:
        kit, conf = "Rapid (SQK-RAD114/RBK114)?", "low"
        retained = "traces of the rapid primer on {:.0%} of read starts".format(max(rb, rad))
    else:
        kit, conf = "Ligation (e.g. SQK-LSK114)", "medium"
        retained = "no barcode or rapid primer found"
    return {"reads": n, "fractions": f, "kit": kit, "confidence": conf, "barcode": barcode, "top_barcodes": top, "evidence": retained}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out"); ap.add_argument("--slices", type=int, default=4); ap.add_argument("--chunk", type=int, default=2_000_000)
    ap.add_argument("--threads", type=int, default=6)
    a = ap.parse_args()
    rows = list(csv.DictReader(open(os.path.join(ROOT, "ENA_data_links.tsv")), delimiter="\t"))
    ena = json.load(open(os.path.join(ROOT, "scripts", ".cache", "ena_runs.json")))
    jobs = {}
    for r in rows:
        for k in ("mods", "sup", "hac", "fast"):
            u = r.get(k, "NA")
            if u.endswith(".bam"):
                acc = u.split("/")[-2]
                jobs.setdefault(r["Sample"], []).append(("https://" + u, int(ena[acc]["submitted_bytes"]), r["dorado version"]))
                break
    def run(sample):
        url, size, dv = jobs[sample][0]
        try:
            res = summarise(sample_seqs(url, size, a.slices, a.chunk))
        except Exception as e:
            res = {"error": str(e)}
        res.update(sample=sample, url=url, dorado=dv)
        print(sample, "|", res.get("kit"), res.get("confidence"), res.get("barcode"), res.get("evidence"), flush=True)
        return sample, res
    with ThreadPoolExecutor(a.threads) as ex:
        out = dict(ex.map(run, list(jobs)))
    json.dump(out, open(a.out, "w"), indent=1)

if __name__ == "__main__":
    main()
