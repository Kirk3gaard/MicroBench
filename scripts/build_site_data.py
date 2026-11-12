#!/usr/bin/env python3
"""Build docs/data/microbench.json for the MicroBench GitHub Pages site.

Sources
  * ENA_data_links.tsv / ena_hic.tsv      – the curated tables in this repo
  * ENA_deep.tsv                           – deeply sequenced samples (one row per flow cell
                                             and basecall; volumes are combined per sample)
  * ENA portal API                         – exact read count, bases and file size per run
  * BAM headers (HTTP range request)       – dorado version, command line, models, GPU,
                                             flow cell, device and run date
  * BAM sampling (HTTP range requests)     – read N50 and mean read Q estimated from
                                             ~30 evenly spaced 1 MB slices of each BAM

Usage
  python3 scripts/build_site_data.py              # uses cache in scripts/.cache
  python3 scripts/build_site_data.py --refresh    # re-query everything

Only the Python standard library is needed.
"""
import argparse, csv, json, os, re, statistics, sys, time, urllib.parse, urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

sys.path.insert(0, os.path.dirname(__file__))
import bamsample as bs  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(os.path.dirname(__file__), ".cache")
OUT = os.path.join(ROOT, "docs", "data", "microbench.json")
TIERS = ["fast", "hac", "hacdup", "sup", "supdup", "mods"]
N_CHUNKS = 30
EXACT_BELOW = 250_000_000  # BAMs smaller than this are read completely -> exact N50


def acc_of(url):
    m = re.search(r"/(ERR\d{6,})/", url)
    return m.group(1) if m else None


def https(u):
    return u if u.startswith("http") else "https://" + u


def ena_runs(accs):
    fields = ("run_accession,study_accession,sample_accession,instrument_model,read_count,"
              "base_count,submitted_bytes,submitted_ftp,submitted_md5,first_public")
    out = {}
    accs = sorted(accs)
    for i in range(0, len(accs), 100):
        q = urllib.parse.urlencode({"result": "read_run", "includeAccessions": ",".join(accs[i:i + 100]),
                                    "fields": fields, "format": "json", "limit": 0})
        with urllib.request.urlopen("https://www.ebi.ac.uk/ena/portal/api/search?" + q, timeout=120) as r:
            for rec in json.load(r):
                out[rec["run_accession"]] = rec
    return out


def parse_header(h):
    info = {}
    for line in h.splitlines():
        f = line.split("\t")
        if f[0] == "@PG" and "PN:dorado" in line:
            for x in f[1:]:
                k, v = x.split(":", 1)
                if k == "VN": info["dorado"] = v
                if k == "CL": info["command"] = v
                if k == "DS": info["hardware"] = v.replace("gpu:", "")
        if f[0] == "@RG" and "basecall_model" not in info:
            for x in f[1:]:
                k, v = x.split(":", 1)
                if k == "PU": info["flowcell_id"] = v
                if k == "PM": info["device"] = v
                if k == "DT": info["run_date"] = v[:10]
                if k == "DS":
                    for kv in v.split(" "):
                        if "=" in kv:
                            a, b = kv.split("=", 1)
                            if a == "basecall_model": info["basecall_model"] = b
                            if a == "modbase_models": info["modbase_models"] = b.split(",")
    return info


def cache_path(url):
    key = re.sub(r"[^A-Za-z0-9_.-]", "_", url.split("/")[-2] + "_" + url.split("/")[-1]) + ".json"
    return os.path.join(CACHE, key)


def sample_bam(url, size, keep=False):
    """Header + sampled stats for one BAM, with on-disk cache. keep=True also stores the sampled
    read lengths and qualities, so several flow cells of one sample can be pooled."""
    path = cache_path(url)
    if os.path.exists(path) and not ARGS.refresh:
        res = json.load(open(path))
        if not keep or "lens" in res: return res
    hdr = parse_header(bs.header(url))
    if size < EXACT_BELOW:
        lens, qs = bs.full_pass(url, size)
        res = {"header": hdr, "sampled_reads": len(lens), "exact": True}
        return finish(res, lens, qs, path)
    lens, qs = [], []
    starts = [int(size * (k + 0.5) / N_CHUNKS) for k in range(N_CHUNKS)]
    for st in starts:
        for attempt in range(3):
            try:
                data = bs.inflate_bgzf(bs.get_range(url, st, min(size - 1, st + 1_000_000)))
                break
            except Exception:
                time.sleep(2 + attempt * 5)
        else:
            continue
        for *_, l, q in bs.parse_records(data):
            if l > 0:
                lens.append(l)
                if q is not None: qs.append(q)
    res = {"header": hdr, "sampled_reads": len(lens), "exact": False}
    if keep: res.update(lens=lens, qs=[round(q, 1) for q in qs])
    return finish(res, lens, qs, path)


def finish(res, lens, qs, path):
    if lens:
        res.update(n50_est=bs.n50(lens), median_len_est=int(statistics.median(lens)))
    if qs:
        res.update(median_q_est=round(statistics.median(qs), 1),
                   q20_frac_est=round(sum(q >= 20 for q in qs) / len(qs), 3))
    json.dump(res, open(path, "w"))
    return res


def clean_cmd(v):
    """Drop local directory paths from a dorado command line, keep the settings."""
    if not v: return v
    v = re.sub(r'"?(/[^\s"]*/)?([^/\s"]+)/?"?(?=\s|$)', lambda m: m.group(2) if m.group(1) else m.group(0), v)
    return v.replace('"', "")


def main():
    os.makedirs(CACHE, exist_ok=True)
    rows = list(csv.DictReader(open(os.path.join(ROOT, "ENA_data_links.tsv")), delimiter="\t"))
    hic = list(csv.reader(open(os.path.join(ROOT, "ena_hic.tsv")), delimiter="\t"))
    hic_head, hic = hic[0], hic[1:]
    deep_path = os.path.join(ROOT, "ENA_deep.tsv")
    deep = list(csv.DictReader(open(deep_path), delimiter="\t")) if os.path.exists(deep_path) else []

    accs = set()
    for r in rows:
        for k in TIERS + ["pod5"]:
            if r.get(k, "NA") not in ("NA", ""): accs.add(acc_of(r[k]))
    for h in hic:
        accs.add(acc_of(h[4])); accs.add(acc_of(h[5]))
    for r in deep:
        for k in TIERS + ["pod5"]:
            if r.get(k, "NA") not in ("NA", ""): accs.add(acc_of(r[k]))
    accs.discard(None)

    cache_ena = os.path.join(CACHE, "ena_runs.json")
    if os.path.exists(cache_ena) and not ARGS.refresh:
        ena = json.load(open(cache_ena))
        missing = accs - set(ena)
        if any("submitted_md5" not in v for v in ena.values()):  # older cache without checksums
            missing = accs
        if missing: ena.update(ena_runs(missing))
    else:
        ena = ena_runs(accs)
    json.dump(ena, open(cache_ena, "w"))

    def md5_of(a, url):
        """MD5 of one submitted file of a run (runs can hold several files)."""
        rec = ena.get(a, {})
        files = (rec.get("submitted_ftp") or "").split(";")
        sums = (rec.get("submitted_md5") or "").split(";")
        name = url.rstrip("/").split("/")[-1]
        for f, m in zip(files, sums):
            if f.split("/")[-1] == name: return m or None
        return None

    def bytes_of(a, url):
        rec = ena.get(a, {})
        files = (rec.get("submitted_ftp") or "").split(";")
        sizes = (rec.get("submitted_bytes") or "").split(";")
        name = url.rstrip("/").split("/")[-1]
        for f, s in zip(files, sizes):
            if f.split("/")[-1] == name and s.isdigit(): return int(s)
        return None

    def num(a, k):
        v = ena.get(a, {}).get(k)
        try: return int(v)
        except (TypeError, ValueError): return None

    # sample all BAMs in parallel
    bams = {}
    for r in rows:
        for k in TIERS:
            u = r.get(k, "NA")
            if u not in ("NA", "") and u.endswith(".bam"):
                a = acc_of(u)
                bams[https(u)] = num(a, "submitted_bytes")
    deep_urls = set()
    for r in deep:
        for k in TIERS:
            u = r.get(k, "NA")
            if u not in ("NA", "") and u.endswith(".bam"):
                bams[https(u)] = bytes_of(acc_of(u), u); deep_urls.add(https(u))
    results = {}
    with ThreadPoolExecutor(max_workers=ARGS.threads) as ex:
        todo = {u: s for u, s in bams.items() if s}
        if ARGS.cached_only:
            todo = {u: s for u, s in todo.items() if os.path.exists(cache_path(u))}
        futs = {ex.submit(sample_bam, u, s, u in deep_urls): u for u, s in todo.items()}
        for i, f in enumerate(as_completed(futs), 1):
            u = futs[f]
            try:
                results[u] = f.result()
            except Exception as e:
                print("WARN", u, e, file=sys.stderr)
            print(f"\r{i}/{len(futs)} BAMs sampled", end="", file=sys.stderr)
    print(file=sys.stderr)

    samples = {}
    for r in rows:
        s = samples.setdefault(r["Sample"], {
            "name": r["Sample"], "type": r["SampleType"], "info": r["SampleInfo"],
            "library_kit": set(), "projects": set(), "datasets": [], "hic": []})
        kit = (r.get("library_kit") or "").strip()
        if kit and kit.lower() not in ("na", "unknown"): s["library_kit"].add(kit)
        ds = {"dorado": r["dorado version"], "model_version": r["basecall model version"],
              "files": [], "pod5": None}
        best = None
        for k in TIERS:
            u = r.get(k, "NA")
            if u in ("NA", ""): continue
            a = acc_of(u); fu = https(u)
            f = {"tier": k, "url": fu, "accession": a, "format": "bam" if u.endswith(".bam") else "tar.gz" if u.endswith(".tar.gz") else "fastq.gz",
                 "reads": num(a, "read_count"), "bases": num(a, "base_count"), "bytes": num(a, "submitted_bytes"),
                 "md5": md5_of(a, fu)}
            if ena.get(a, {}).get("study_accession"): s["projects"].add(ena[a]["study_accession"])
            if fu in results:
                res = results[fu]
                f["basecall_model"] = res["header"].get("basecall_model")
                f["modbase_models"] = res["header"].get("modbase_models")
                f["command"] = clean_cmd(res["header"].get("command"))
                for kk in ("n50_est", "median_q_est", "q20_frac_est", "sampled_reads", "exact"):
                    if kk in res: f[kk] = res[kk]
                # prefer the sup/mods BAM for dataset-level run info
                if best is None or k in ("mods", "sup"):
                    best = res
            ds["files"].append(f)
        if r.get("pod5", "NA") not in ("NA", ""):
            a = acc_of(r["pod5"])
            ds["pod5"] = {"url": https(r["pod5"]), "accession": a, "bytes": num(a, "submitted_bytes"), "md5": md5_of(a, r["pod5"])}
            if ena.get(a, {}).get("study_accession"): s["projects"].add(ena[a]["study_accession"])
        if best:
            h = best["header"]
            ds["run"] = {k: h.get(k) for k in ("flowcell_id", "device", "run_date", "hardware")}
            ds["chemistry"] = None
            bm = h.get("basecall_model") or ""
            m = re.match(r"dna_(r[\d.]+)_(e[\d.]+)_(\d+bps)", bm)
            if m: ds["chemistry"] = f"{m.group(1).upper()} {m.group(2)} {m.group(3)}"
            ds["n50_est"] = best.get("n50_est")
            ds["sampled_reads"] = best.get("sampled_reads")
            ds["exact"] = best.get("exact", False)
        # dataset-level read count: prefer sup simplex, then any
        pick = next((f for t in ("sup", "mods", "hac", "fast") for f in ds["files"] if f["tier"] == t and f["reads"]), None)
        if pick:
            ds["reads"], ds["bases"] = pick["reads"], pick["bases"]
            ds["mean_len"] = round(pick["bases"] / pick["reads"]) if pick["reads"] else None
        s["datasets"].append(ds)

    for h in hic:
        s = samples.get(h[0])
        if not s: continue
        a = acc_of(h[4])
        s["hic"].append({"project": h[1], "kit": h[2], "platform": h[3],
                         "r1": https(h[4]), "r2": https(h[5]), "accession": a,
                         "reads": num(a, "read_count"), "bases": num(a, "base_count"),
                         "r1_md5": md5_of(a, h[4]), "r2_md5": md5_of(acc_of(h[5]), h[5]),
                         "r1_bytes": bytes_of(a, h[4]), "r2_bytes": bytes_of(acc_of(h[5]), h[5])})
        s["projects"].add(h[1])

    # deeply sequenced samples: one dataset per flow cell (and basecall); volumes combined per sample
    for r in deep:
        s = samples.setdefault(r["Sample"], {
            "name": r["Sample"], "type": r["SampleType"], "info": "", "description": r["SampleInfo"],
            "biosample": r.get("BioSample"), "library_kit": set(), "projects": set(), "datasets": [], "hic": []})
        kit = (r.get("library_kit") or "").strip()
        if kit and kit.lower() not in ("na", "unknown"): s["library_kit"].add(kit)
        ds = {"dorado": r["dorado version"], "model_version": r["basecall model version"], "flowcell": r["flowcell"],
              "files": [], "pod5": None}
        for k in TIERS:
            u = r.get(k, "NA")
            if u in ("NA", ""): continue
            a = acc_of(u); fu = https(u)
            f = {"tier": k, "url": fu, "accession": a, "format": "bam" if u.endswith(".bam") else "fastq.gz",
                 "reads": num(a, "read_count"), "bases": num(a, "base_count"), "bytes": bytes_of(a, u), "md5": md5_of(a, fu)}
            if ena.get(a, {}).get("study_accession"): s["projects"].add(ena[a]["study_accession"])
            res = results.get(fu)
            if res:
                f["basecall_model"] = res["header"].get("basecall_model")
                f["modbase_models"] = res["header"].get("modbase_models")
                f["command"] = clean_cmd(res["header"].get("command"))
                for kk in ("n50_est", "median_q_est", "q20_frac_est", "sampled_reads", "exact"):
                    if kk in res: f[kk] = res[kk]
                h = res["header"]
                ds["run"] = {k2: h.get(k2) for k2 in ("flowcell_id", "device", "run_date", "hardware")}
                m = re.match(r"dna_(r[\d.]+)_(e[\d.]+)_(\d+bps)", h.get("basecall_model") or "")
                ds["chemistry"] = f"{m.group(1).upper()} {m.group(2)} {m.group(3)}" if m else None
                ds["n50_est"], ds["sampled_reads"], ds["exact"] = res.get("n50_est"), res.get("sampled_reads"), False
            ds["files"].append(f)
        pick = next((f for f in ds["files"] if f["reads"]), None)
        if pick:
            ds["reads"], ds["bases"] = pick["reads"], pick["bases"]
            ds["mean_len"] = round(pick["bases"] / pick["reads"])
        s["datasets"].append(ds)
    for s in samples.values():
        if s["type"] != "deep": continue
        # each flow cell counts once; where a flow cell was basecalled twice, the newest basecall is used
        ver = lambda d: tuple(int(x) for x in re.findall(r"\d+", d["dorado"]))
        per_fc = {}
        for d in s["datasets"]:
            if d["flowcell"] not in per_fc or ver(d) > ver(per_fc[d["flowcell"]]): per_fc[d["flowcell"]] = d
        use = list(per_fc.values())
        for d in s["datasets"]: d["in_combined"] = d in use
        lens, qs = [], []
        for d in use:
            for f in d["files"]:
                res = results.get(f["url"]) or {}
                lens += res.get("lens", []); qs += res.get("qs", [])
        reads = sum(d.get("reads") or 0 for d in use); bases = sum(d.get("bases") or 0 for d in use)
        s["combined"] = {
            "flowcells": len(use), "basecalls": len(s["datasets"]), "reads": reads, "bases": bases,
            "mean_len": round(bases / reads) if reads else None,
            "bytes": sum(f["bytes"] or 0 for d in use for f in d["files"]),
            "bytes_all": sum(f["bytes"] or 0 for d in s["datasets"] for f in d["files"]),
            "n50_est": bs.n50(lens) if lens else None, "median_q_est": round(statistics.median(qs), 1) if qs else None,
            "sampled_reads": len(lens),
            "dorado": sorted({d["dorado"] for d in use}, key=lambda v: tuple(map(int, v.split(".")))),
            "model_version": sorted({d["model_version"] for d in use})}

    kits_path = os.path.join(CACHE, "kits.json")
    kits = json.load(open(kits_path)) if os.path.exists(kits_path) else {}
    out = []
    for s in samples.values():
        k = kits.get(s["name"])
        if k and "kit" in k:
            s["kit_inferred"] = {x: k.get(x) for x in ("kit", "confidence", "barcode", "evidence", "reads")}
            s["kit_inferred"]["fractions"] = k.get("fractions")
        s["library_kit"] = sorted(s["library_kit"])
        s["projects"] = sorted(s["projects"])
        if s["type"] == "deep": s["datasets"].sort(key=lambda d: (d["flowcell"], d["dorado"]))
        else: s["datasets"].sort(key=lambda d: (d["model_version"], d["dorado"]), reverse=True)
        out.append(s)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump({"generated": time.strftime("%Y-%m-%d"), "main_project": "PRJEB85558",
               "n_chunks": N_CHUNKS, "samples": out}, open(OUT, "w"), indent=1)
    print("wrote", OUT, file=sys.stderr)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--refresh", action="store_true")
    ap.add_argument("--cached-only", action="store_true", help="skip BAMs not sampled yet (quick rebuild)")
    ap.add_argument("--threads", type=int, default=12)
    ARGS = ap.parse_args()
    main()
