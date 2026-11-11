"""Coincidence test: do indel events in different channels happen at the same wall-clock time?
Compares observed cross-channel coincidences with a null where each read's events keep their
within-read timing but reads are shifted by random offsets (preserving the event-per-read structure)."""
import sys, json, random, bisect, datetime, collections, statistics

reads = [json.loads(l) for l in open(sys.argv[1])]
MINLEN = int(sys.argv[2]) if len(sys.argv) > 2 else 10
WIN = float(sys.argv[3]) if len(sys.argv) > 3 else 2.0
ts = lambda s: datetime.datetime.fromisoformat(s).timestamp()
t0 = min(ts(r["st"]) for r in reads)

# recurrent (real) deletions filter
pos = collections.defaultdict(list)
for r in reads:
    for s in r["segs"]:
        for e in s["long"]:
            if e[0] == "del" and e[3] >= MINLEN: pos[s["ctg"]].append(e[2])
for c in pos: pos[c].sort()
def recurrent(ctg, rp):
    a = pos[ctg]; return bisect.bisect_right(a, rp + 30) - bisect.bisect_left(a, rp - 30) >= 3

per_read = []
for r in reads:
    st = ts(r["st"]) - t0
    ev = []
    for s in r["segs"]:
        for e in s["long"]:
            if e[0] == "del" and e[3] >= MINLEN and not recurrent(s["ctg"], e[2]):
                ev.append(r["du"] * e[1] / r["len"])
    per_read.append((st, r["du"], r["ch"], sorted(ev)))

def coincidences(pr):
    ev = sorted((st + x, ch) for st, du, ch, xs in pr for x in xs)
    times = [t for t, _ in ev]
    n = 0
    for i, (t, ch) in enumerate(ev):
        j = bisect.bisect_left(times, t - WIN)
        k = bisect.bisect_right(times, t + WIN)
        n += sum(1 for _, c in ev[j:k] if c != ch)
    return n // 2, len(ev)

obs, nev = coincidences(per_read)
# concurrency: mean number of sampled reads in progress at a random time
span = max(st + du for st, du, *_ in per_read)
conc = sum(du for _, du, *_ in per_read) / span
random.seed(1)
null = []
for _ in range(50):
    shuffled = [(random.uniform(0, span - du), du, ch, xs) for st, du, ch, xs in per_read]
    null.append(coincidences(shuffled)[0])
print(json.dumps({"min_len": MINLEN, "window_s": WIN, "events": nev, "observed_pairs": obs,
                  "null_mean": statistics.mean(null), "null_max": max(null),
                  "mean_concurrent_reads": round(conc, 2), "run_hours": round(span / 3600, 1)}))
