"""Flow-cell-wide signal gaps: seconds in which many reads, on many channels, lose a stretch of
their molecule at the same moment (seen as a deletion inside an alignment, or as a collinear split
of the alignment).

usage: gaps.py <prefix> out.json [ontrim_report.html]
"""
import sys, json, re, html, datetime, math
import numpy as np, pandas as pd
from scipy.stats import poisson

P, OUT = sys.argv[1], sys.argv[2]
R = pd.read_csv(P + ".reads.tsv.gz", sep="\t", usecols=["id", "ch", "st", "du"])
E = pd.read_csv(P + ".events.tsv.gz", sep="\t")
t0 = R.loc[R.st > 0, "st"].min()
R["t"] = R.st - t0; R["end"] = R.t + R.du
E["t"] = E.t - t0
S = int(math.ceil(R.end.max())) + 2
act = np.zeros(S + 1)
np.add.at(act, R.t.clip(lower=0).astype(int), 1); np.add.at(act, (R.end.clip(lower=0).astype(int) + 1).clip(upper=S), -1)
act = np.cumsum(act)[:S]

G = E[((E.type == "del") & (E.len >= 10)) | ((E.type == "split") & (E.len >= 10))].copy()
G["sec"] = G.t.clip(lower=0).astype(int).clip(upper=S - 1)
# count each read once per second
per = G.drop_duplicates(["id", "sec"])
cnt = np.bincount(per.sec, minlength=S)
rate = cnt.sum() / act.sum()
exp = act * rate
p = poisson.sf(cnt - 1, np.maximum(exp, 1e-9))
sig = np.where((p < 0.01 / S) & (cnt >= 10))[0]
# merge neighbouring seconds (timing of an event within a read is estimated to ~1-2 s)
events = []
for s in sig:
    if events and s - events[-1][-1] <= 3: events[-1].append(s)
    else: events.append([s])
rows = []
for ev in events:
    a, b = ev[0] - 2, ev[-1] + 3
    X = per[(per.sec >= a) & (per.sec < b)]
    act_here = act[ev[0]]
    rows.append({"t_utc": datetime.datetime.fromtimestamp(t0 + ev[0], datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
                 "t_h": round(ev[0] / 3600, 4), "reads": int(X.id.nunique()), "channels": int(X.ch.nunique()),
                 "active_reads": int(act_here), "share_of_active": round(X.id.nunique() / max(1, act_here), 3),
                 "missing_median": float(X.len.median()), "split_share": round(float((X.type == "split").mean()), 3)})
hit_ids = per[per.sec.isin([s for ev in events for s in range(ev[0] - 2, ev[-1] + 3)])].id.unique()
reads_hit = len(hit_ids)
open(P + ".gap_reads.txt", "w").write("\n".join(hit_ids) + ("\n" if len(hit_ids) else ""))
out = {"n_reads": int(len(R)), "run_hours": round(float(R.end.max() / 3600), 2), "mean_active": round(float(act.mean()), 1),
       "background_per_read_hour": round(rate * 3600, 3), "n_events": len(rows), "reads_hit": int(reads_hit),
       "pct_reads_hit": round(100 * reads_hit / len(R), 3), "events": rows,
       "per_minute": [[round(i / 60, 3), int(c), round(float(a), 1)] for i, (c, a) in
                      enumerate(zip(cnt[: S // 60 * 60].reshape(-1, 60).sum(1), act[: S // 60 * 60].reshape(-1, 60).mean(1)))]}

# compare with an ontrim report, if given
if len(sys.argv) > 3:
    t = open(sys.argv[3], encoding="utf-8", errors="replace").read()
    t = re.sub(r"<script.*?</script>|<style.*?</style>", "", t, flags=re.S)
    txt = html.unescape(re.sub(r"<[^>]+>", " ", t))
    m = re.findall(r"(\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ)\s+(\d+)\s+(\d+)\s+(\d+)", txt)
    ont = {}
    for ts, r_, c_, b_ in m:
        ont.setdefault(ts, (int(r_), int(c_), int(b_)))  # first table = whole run
    ot = np.array(sorted(datetime.datetime.strptime(x, "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=datetime.timezone.utc).timestamp() for x in ont))
    mine = np.array([t0 + r["t_h"] * 3600 for r in rows])
    def matched(a, b, tol=5):
        if len(b) == 0: return 0
        idx = np.searchsorted(b, a)
        d = np.minimum(np.abs(a - b[np.clip(idx, 0, len(b) - 1)]), np.abs(a - b[np.clip(idx - 1, 0, len(b) - 1)]))
        return int((d <= tol).sum())
    out["ontrim"] = {"events": len(ot), "matched_by_mine": matched(ot, mine), "mine_matched_by_ontrim": matched(mine, ot)}
json.dump(out, open(OUT, "w"))
print({k: v for k, v in out.items() if k not in ("events", "per_minute")})
for r in rows[:15]: print(r)
