"""Assemble the artefact page data for one run from the outputs of post.sh, and update the run index.

usage: merge_page.py STAGE_DIR TAG ID LABEL DEVICE out_dir [extra-json]
  TAG   prefix of the files in STAGE_DIR (e.g. hmw -> hmw_page.json, hmw_numbers.json, ...)
  ID    run id used on the page (e.g. zymohmw -> artefacts_zymohmw.json)
extra-json is merged into the dataset block (kit, flowcell_note, barcode_note, ...).
"""
import sys, os, json

S, T, ID, LABEL, DEVICE, OUT = sys.argv[1:7]
extra = json.loads(sys.argv[7]) if len(sys.argv) > 7 else {}
L = lambda name: json.load(open(os.path.join(S, f"{T}_{name}.json")))
d = L("page")
num = L("numbers")
if num.get("check_problems"):
    raise SystemExit(f"number check failed: {num['check_problems']}")
d["numbers"] = num
d["gaps"] = L("gaps"); d["duplex"] = L("duplex"); d["saw"] = L("saw")
sc = L("scan_stats"); sc["junk_time"] = L("junk_time")
se = L("scan_examples"); sc["n_post_scan_reads"] = se["n_pool"]
d["scan"] = sc
d["recurrent_junctions"] = L("recurrent")
# summary = recomputed per-class numbers (each read once), core union as "_any"
for k, v in num["classes"].items():
    d["summary"][k] = v
d["summary"]["_any"] = num["any_core"]
# post-scan examples only where the scan effect exists (junk >2x more common right after a scan)
sn = sc["since"]
mid = [x for x in sn[4:8]]
base = sum(x["lowc"] * x["n"] for x in mid) / max(1, sum(x["n"] for x in mid))
if sn and base and sn[0]["lowc"] > 2 * base and se["examples"]:
    d["examples"]["post_scan_junk"] = se["examples"]
else:
    d["examples"].pop("post_scan_junk", None)
for k in ("barcode_note",):
    if k in extra: d[k] = extra.pop(k)
d["dataset"].update(extra)
d["dataset"]["device"] = DEVICE
os.makedirs(OUT, exist_ok=True)
fn = f"artefacts_{ID}.json"
json.dump(d, open(os.path.join(OUT, fn), "w"), separators=(",", ":"))
ix_path = os.path.join(OUT, "artefacts_index.json")
ix = json.load(open(ix_path)) if os.path.exists(ix_path) else []
ix = [x for x in ix if x["id"] != ID]
ix.append({"id": ID, "label": LABEL, "file": fn, "device": DEVICE, "n_reads": d["stats"]["n_reads"],
           "gaps": d["gaps"]["n_events"], "gaps_pct": d["gaps"]["pct_reads_hit"]})
order = ["zymohmw", "zymogut", "zymooral", "zymofecal"]
ix.sort(key=lambda x: order.index(x["id"]) if x["id"] in order else 99)
json.dump(ix, open(ix_path, "w"), indent=1)
print(fn, "any", num["any_core"]["pct_reads"], "chim/fold", num["chimeric_or_foldback"]["pct_reads"],
      "gaps", d["gaps"]["n_events"], "recurrent moved", d["recurrent_junctions"]["reads_moved"])
