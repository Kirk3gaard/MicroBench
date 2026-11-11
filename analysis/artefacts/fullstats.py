"""Summaries of a whole-run artefact scan (output of fullrun.py) for the web page.

usage: fullstats.py <prefix> out.json
"""
import sys, json, gzip, math, collections
import numpy as np, pandas as pd
from scipy.stats import poisson

P, OUT = sys.argv[1], sys.argv[2]
R = pd.read_csv(P + ".reads.tsv.gz", sep="\t", dtype={"ch": "Int64", "mx": "Int64", "ctg": "string", "strand": "string"})
E = pd.read_csv(P + ".events.tsv.gz", sep="\t")
t0 = R.loc[R.st > 0, "st"].min()
R["t"] = R.st - t0
R["end"] = R.t + R.du
E["t"] = E.t - t0
N, B = len(R), int(R["len"].sum())
for c in ["chimera_interspecies", "chimera_intraspecies", "foldback", "internal_adapter", "end_adapter",
          "low_complexity", "skip_prone", "unaligned_head", "unaligned_tail", "unmapped", "clean"]:
    R[c] = R.cls.str.contains(c, regex=False)
CHIM = R.chimera_interspecies | R.chimera_intraspecies | R.foldback
ART = CHIM | R.internal_adapter | R.end_adapter | R.low_complexity | R.skip_prone
out = {"n_reads": N, "n_bases": B, "n_channels": int(R.ch.nunique()), "run_hours": float(R.end.max() / 3600),
       "dorado_split": int(R.split.sum())}

summary = {}
for c in ["chimera_interspecies", "chimera_intraspecies", "foldback", "internal_adapter", "end_adapter",
          "low_complexity", "skip_prone", "unaligned_head", "unaligned_tail", "unmapped"]:
    m = R[c]
    summary[c] = {"reads": int(m.sum()), "pct_reads": round(100 * m.mean(), 3), "pct_bases": round(100 * R.loc[m, "len"].sum() / B, 3)}
summary["_any"] = {"reads": int(ART.sum()), "pct_reads": round(100 * ART.mean(), 3), "pct_bases": round(100 * R.loc[ART, "len"].sum() / B, 3)}
out["summary"] = summary

# filters
filt = {}
for c in ["chimera_interspecies", "chimera_intraspecies", "foldback", "internal_adapter", "end_adapter", "low_complexity", "skip_prone", "clean"]:
    X = R[R[c]]
    if len(X) == 0: continue
    f = lambda q, l: round(100 * float(((X.qs >= q) & (X["len"] >= l)).mean()), 1)
    filt[c] = {"q10_l1000": f(10, 1000), "q15_l1000": f(15, 1000), "q20_l5000": f(20, 5000), "median_qs": round(float(X.qs.median()), 1)}
out["filters"] = filt

# ---------- recurrent (genomic) vs sporadic deletions ----------
# coverage per contig from the longest alignment of each read
cov = {}
HAVE_COORDS = "ctg" in R.columns
if HAVE_COORDS:
    M = R[R.ctg.notna() & (R.ctg != "")]
    for ctg, g in M.groupby("ctg"):
        L = int(g.r1.max()) + 2
        d = np.zeros(L + 1, dtype=np.int32)
        np.add.at(d, g.r0.astype(int).values, 1); np.add.at(d, g.r1.astype(int).values, -1)
        cov[ctg] = np.cumsum(d)
else:
    # mean depth per species (aligned bases / genome size), used for every contig of that species
    import mappy
    ref = sys.argv[3] if len(sys.argv) > 3 else "zymo_hmw.fa"
    sp_of = dict(l.rstrip("\n").split("\t") for l in open(ref + ".species"))
    gsize = collections.Counter()
    clen = {}
    for n, s, _ in mappy.fastx_read(ref):
        clen[n] = len(s); gsize[sp_of[n]] += len(s)
    depth = R.groupby("sp")["aln"].sum()
    for c, L in clen.items():
        cov[c] = np.full(L + 1, int(depth.get(sp_of[c], 0) / gsize[sp_of[c]]), dtype=np.int32)
D = E[E.type == "del"].copy()
D["rbin"] = (D.rpos // 20).astype(int)
site = D.groupby(["ctg", "rbin"]).size().rename("n_site").reset_index()
site["cov"] = [int(cov[c][min(len(cov[c]) - 1, b * 20 + 10)]) if c in cov else 0 for c, b in zip(site.ctg, site.rbin)]
site["frac"] = site.n_site / site["cov"].clip(lower=1)
D = D.merge(site, on=["ctg", "rbin"])
D["kind"] = np.where((D.frac >= 0.3) & (D.n_site >= 3), "genomic",
             np.where((D.frac >= 0.03) & (D.n_site >= 5) & (D["cov"] >= 50), "systematic", "sporadic"))
out["deletions"] = {k: int(v) for k, v in D.kind.value_counts().items()}
# systematic hotspots: sites where 3-30% of reads carry the same >=10 bp deletion
hot = site[(site.frac >= 0.03) & (site.frac < 0.3) & (site["cov"] >= 50) & (site.n_site >= 5)].sort_values("frac", ascending=False)
out["systematic_sites"] = hot.head(40).assign(pos=lambda x: x.rbin * 20)[["ctg", "pos", "n_site", "cov", "frac"]].round(3).to_dict("records")
out["n_systematic_sites"] = int(len(hot))

# ---------- synchrony: events per second vs reads in progress ----------
S = int(math.ceil(R.end.max())) + 1
active = np.zeros(S + 1)
ts = R.t.clip(lower=0).values; te = R.end.clip(lower=0).values
# fractional occupancy per second
a = np.floor(ts).astype(int); b = np.floor(te).astype(int) + 1
np.add.at(active, a, 1); np.add.at(active, b, -1)
active = np.cumsum(active)[:S].astype(float)
sp = D[D.kind == "sporadic"]
res = {}
for minlen in (10, 20, 50):
    ev = sp[sp.len >= minlen]
    cnt = np.bincount(np.clip(ev.t.astype(int), 0, S - 1), minlength=S)
    rate = cnt.sum() / active.sum()
    exp = active * rate
    p = poisson.sf(cnt - 1, np.maximum(exp, 1e-9))
    sig = np.where(p < 0.01 / S)[0]
    chans = [int(ev[(ev.t >= s) & (ev.t < s + 1)].ch.nunique()) for s in sig[:200]]
    # 1-minute series
    m = S // 60
    cm = cnt[:m * 60].reshape(m, 60).sum(1); am = active[:m * 60].reshape(m, 60).sum(1)
    res[str(minlen)] = {"events": int(cnt.sum()), "rate_per_read_hour": round(rate * 3600, 4),
                        "mean_active": round(float(active.mean()), 1),
                        "sig_seconds": int(len(sig)), "sig_examples": [{"t_h": round(s / 3600, 4), "n": int(cnt[s]), "exp": round(float(exp[s]), 3), "channels": c} for s, c in zip(sig[:30], chans[:30])],
                        "max_in_one_second": int(cnt.max()), "dispersion_1min": round(float(np.var(cm[am > 0] / am[am > 0] * am.mean()) / max(1e-9, np.mean(cm[am > 0] / am[am > 0] * am.mean()))), 2)}
out["sync"] = res

# read-end synchrony: seconds where many reads end at once (MinKNOW actions) and deletions near them
ends = np.bincount(np.clip(te.astype(int), 0, S - 1), minlength=S)
exp_end = active * (ends.sum() / active.sum())
pend = poisson.sf(ends - 1, np.maximum(exp_end, 1e-9))
mass = np.where((pend < 0.01 / S) & (ends > 3 * exp_end) & (ends >= 5))[0]
out["mass_end_seconds"] = {"n": int(len(mass)), "examples": [{"t_h": round(s / 3600, 4), "ends": int(ends[s]), "exp": round(float(exp_end[s]), 1)} for s in mass[:40]]}
if len(mass) > 2:
    gaps = np.diff(mass) / 3600
    out["mass_end_seconds"]["median_gap_h"] = round(float(np.median(gaps[gaps > 0.05])), 3) if (gaps > 0.05).any() else None
if len(mass):
    near = np.zeros(S, bool)
    for s in mass:
        near[max(0, s - 5):s + 1] = True
    ev = sp[sp.len >= 10]
    sec = np.clip(ev.t.astype(int), 0, S - 1)
    in_near = near[sec].sum(); act_near = active[near].sum(); act_all = active.sum()
    out["mass_end_deletions"] = {"events_near": int(in_near), "expected_near": round(float(len(ev) * act_near / act_all), 1)}

# ---------- timeline (15 min bins) ----------
R["bin"] = (R.t // 900).astype(int)
L1 = R[(R["len"] >= 1000) & (R.aln > 0)].copy()
L1["delkb"] = 1000 * L1["del"] / L1.aln; L1["inskb"] = 1000 * L1.ins / L1.aln; L1["mmkb"] = 1000 * L1.mm / L1.aln
L1["chim"] = CHIM[L1.index]
tl = L1.groupby("bin").agg(n=("id", "size"), channels=("ch", "nunique"), dels=("delkb", "median"), ins=("inskb", "median"),
                           mm=("mmkb", "median"), chim=("chim", "mean"), qs=("qs", "median")).reset_index()
ev_bin = np.bincount((sp[sp.len >= 20].t // 900).astype(int), minlength=int(R.bin.max()) + 1)
act_bin = np.array([active[i * 900:(i + 1) * 900].sum() for i in range(int(R.bin.max()) + 1)])
tl["longdel_per_read_h"] = [round(3600 * ev_bin[b] / act_bin[b], 3) if act_bin[b] else None for b in tl.bin]
tl["active"] = [round(act_bin[b] / 900, 1) for b in tl.bin]
tl["h"] = tl.bin / 4
out["timeline"] = tl.drop(columns="bin").round(4).to_dict("records")

# ---------- channels ----------
ch = R.groupby("ch").agg(n=("id", "size"), bases=("len", "sum"), chim=("chimera_interspecies", "mean"),
                         lowc=("low_complexity", "mean"), skip=("skip_prone", "mean"), qs=("qs", "median"))
ch = ch[ch.n >= 200]
pc = R.chimera_interspecies.mean()
z = (ch.chim - pc) / np.sqrt(pc * (1 - pc) / ch.n)
pl = R.low_complexity.mean(); zl = (ch.lowc - pl) / np.sqrt(pl * (1 - pl) / ch.n)
out["channels"] = {"n": int(len(ch)), "chim_var_ratio": round(float(z.var()), 2), "lowc_var_ratio": round(float(zl.var()), 2),
                   "top_lowc": ch.assign(z=zl).sort_values("z", ascending=False).head(10).reset_index().round(4).to_dict("records"),
                   "top_chim": ch.assign(z=z).sort_values("z", ascending=False).head(10).reset_index().round(4).to_dict("records"),
                   "hist_lowc": np.histogram(ch.lowc * 100, bins=30)[0].tolist(), "hist_lowc_edges": np.histogram(ch.lowc * 100, bins=30)[1].round(3).tolist()}

# ---------- chimera partners vs abundance ----------
sp_ab = R[R.sp.notna() & (R.sp != "")].groupby("sp")["len"].sum()
out["species_bases"] = {k: int(v) for k, v in sp_ab.items()}

# duplicates: same strand, start and end within 10 bp
if HAVE_COORDS:
    Mm = M[["ctg", "strand", "r0", "r1", "ch", "t"]].copy()
    Mm["k0"] = (Mm.r0 // 10); Mm["k1"] = (Mm.r1 // 10)
    dup = Mm.groupby(["ctg", "strand", "k0", "k1"]).size()
    out["duplicate_groups"] = int((dup >= 2).sum())
else:
    out["duplicate_groups"] = None

json.dump(out, open(OUT, "w"), default=lambda o: o.item() if hasattr(o, "item") else str(o))
print(json.dumps({k: v for k, v in out.items() if k in ("n_reads", "summary", "deletions", "n_systematic_sites", "sync", "mass_end_seconds", "mass_end_deletions", "channels", "duplicate_groups")}, indent=1, default=str)[:6000])
