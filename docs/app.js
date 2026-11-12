"use strict";
const TYPES = [
  ["monoculture", "Monocultures"],
  ["mock", "Mock communities"],
  ["metagenome", "Metagenomes"],
];
const TIER_LABEL = { fast: "fast", hac: "hac", hacdup: "hac duplex", sup: "sup", supdup: "sup duplex", mods: "sup + mods" };
const state = { q: "", types: new Set(), model: "", has: "", sort: "type" };
let DATA = null;
const open = new Set();

const $ = (s, el = document) => el.querySelector(s);
const esc = (s) => String(s ?? "").replace(/[&<>"']/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
const slug = (s) => s.replace(/ /g, "-").toLowerCase();
const fmtN = (n) => n == null ? "–" : n >= 1e9 ? (n / 1e9).toFixed(1) + " G" : n >= 1e6 ? (n / 1e6).toFixed(2) + " M" : n >= 1e3 ? (n / 1e3).toFixed(1) + " k" : String(n);
const fmtBases = (n) => n == null ? "–" : n >= 1e9 ? (n / 1e9).toFixed(n >= 1e11 ? 0 : 1) + " Gb" : (n / 1e6).toFixed(0) + " Mb";
const fmtLen = (n) => n == null ? "–" : n >= 1000 ? (n / 1000).toFixed(1) + " kb" : n + " bp";
const fmtBytes = (n) => n == null ? "" : n >= 1e12 ? (n / 1e12).toFixed(1) + " TB" : n >= 1e9 ? (n / 1e9).toFixed(1) + " GB" : (n / 1e6).toFixed(0) + " MB";

function sampleBases(s) { return s.datasets.reduce((m, d) => Math.max(m, d.bases || 0), 0); }
function sampleN50(s) { return s.datasets.reduce((m, d) => Math.max(m, d.n50_est || 0), 0); }
function hasData(s, what) {
  if (what === "hic") return s.hic.length > 0;
  if (what === "pod5") return s.datasets.some((d) => d.pod5);
  if (what === "dup") return s.datasets.some((d) => d.files.some((f) => f.tier.endsWith("dup")));
  return s.datasets.some((d) => d.files.some((f) => f.tier === what));
}
function haystack(s) {
  const parts = [s.name, s.type, s.info, ...s.projects, ...s.library_kit, s.kit_inferred?.kit, s.kit_inferred?.barcode];
  for (const d of s.datasets) {
    parts.push(d.dorado, d.model_version, d.chemistry, d.run?.flowcell_id, d.run?.device, d.run?.hardware);
    for (const f of d.files) parts.push(f.accession, f.basecall_model, f.url.split("/").pop());
    if (d.pod5) parts.push(d.pod5.accession);
  }
  for (const h of s.hic) parts.push(h.accession, h.kit, h.platform, "hi-c hic");
  return parts.filter(Boolean).join(" ").toLowerCase();
}

function filtered() {
  const terms = state.q.toLowerCase().split(/\s+/).filter(Boolean);
  let out = DATA.samples.filter((s) => {
    if (state.types.size && !state.types.has(s.type)) return false;
    if (state.model && !s.datasets.some((d) => d.model_version === state.model)) return false;
    if (state.has && !hasData(s, state.has)) return false;
    if (terms.length) { const h = s._hay; if (!terms.every((t) => h.includes(t))) return false; }
    return true;
  });
  const typeIdx = (t) => TYPES.findIndex((x) => x[0] === t);
  const byName = (a, b) => a.name.localeCompare(b.name);
  if (state.sort === "name") out.sort(byName);
  else if (state.sort === "bases") out.sort((a, b) => sampleBases(b) - sampleBases(a));
  else if (state.sort === "n50") out.sort((a, b) => sampleN50(b) - sampleN50(a));
  else out.sort((a, b) => typeIdx(a.type) - typeIdx(b.type) || byName(a, b));
  return out;
}

function visibleDatasets(s) {
  return state.model ? s.datasets.filter((d) => d.model_version === state.model) : s.datasets;
}

const FMT = { bam: "bam", "fastq.gz": "fq", "tar.gz": "tar.gz" };
const approx = (v, exact) => v == null ? "–" : (exact ? "" : "≈ ") + v;
// several tiers can point at the same archive (e.g. one tar.gz holding sup + mods + pod5)
function mergedFiles(d) {
  const by = new Map();
  for (const f of d.files) {
    if (by.has(f.url)) by.get(f.url).labels.push(TIER_LABEL[f.tier] || f.tier);
    else by.set(f.url, { ...f, labels: [TIER_LABEL[f.tier] || f.tier] });
  }
  const out = [...by.values()];
  for (const f of out) if (f.labels.includes("sup + mods")) f.labels = f.labels.filter((l) => l !== "sup");
  return out;
}
function fileChip(f) {
  const label = (f.labels || [TIER_LABEL[f.tier] || f.tier]).join(" + ");
  const title = `${f.url.split("/").pop()}\n${f.accession}${f.reads ? ` · ${f.reads.toLocaleString()} reads` : ""}`;
  return `<a class="file" href="${esc(f.url)}" title="${esc(title)}" download>${esc(label)} <small>${esc(FMT[f.format] || f.format)} · ${fmtBytes(f.bytes)}</small></a>`;
}

function detailsRow(d, key) {
  const r = d.run || {};
  const cmds = [...new Set(d.files.map((f) => f.command).filter(Boolean))];
  const models = [...new Set(d.files.map((f) => f.basecall_model).filter(Boolean))];
  const mods = [...new Set(d.files.flatMap((f) => f.modbase_models || []))];
  const rows = [
    ["Chemistry", d.chemistry],
    ["Basecall models", models.join("<br>")],
    ["Modification models", mods.join("<br>")],
    ["dorado command", cmds.map((c) => `<code>${esc(c)}</code>`).join("<br>")],
    ["GPU", r.hardware], ["Flow cell", r.flowcell_id], ["Device", r.device], ["Run date", r.run_date],
  ].filter(([, v]) => v);
  const files = d.files.map((f) => `<tr><td>${esc(TIER_LABEL[f.tier] || f.tier)}</td>
      <td><a href="https://www.ebi.ac.uk/ena/browser/view/${esc(f.accession)}">${esc(f.accession)}</a></td>
      <td class="num">${f.reads ? f.reads.toLocaleString() : "–"}</td><td class="num">${fmtBases(f.bases)}</td>
      <td class="num">${approx(f.n50_est ? fmtLen(f.n50_est) : null, f.exact)}</td><td class="num">${approx(f.median_q_est, f.exact)}</td>
      <td class="num">${fmtBytes(f.bytes)}</td></tr>`).join("");
  const pod = d.pod5 ? `<tr><td>raw pod5</td><td><a href="https://www.ebi.ac.uk/ena/browser/view/${esc(d.pod5.accession)}">${esc(d.pod5.accession)}</a></td><td class="num">–</td><td class="num">–</td><td class="num">–</td><td class="num">–</td><td class="num">${fmtBytes(d.pod5.bytes)}</td></tr>` : "";
  return `<tr class="details" data-key="${esc(key)}"><td colspan="8">
    <dl>${rows.map(([k, v]) => `<dt>${k}</dt><dd>${k === "dorado command" || k.endsWith("models") ? v : esc(v)}</dd>`).join("")}</dl>
    <div class="pertier table-wrap"><table><thead><tr><th>File</th><th>Run</th><th class="num">Reads</th><th class="num">Yield</th><th class="num">N50</th><th class="num">Read Q</th><th class="num">Size</th></tr></thead>
    <tbody>${files}${pod}</tbody></table></div></td></tr>`;
}

function card(s) {
  const el = $("#cardTpl").content.firstElementChild.cloneNode(true);
  el.id = slug(s.name);
  $(".name", el).innerHTML = `<a href="#${esc(el.id)}">${esc(s.name)}</a>`;
  const meta = [];
  if (s.info) meta.push(`<a href="${esc(s.info)}" target="_blank" rel="noopener">sample description ↗</a>`);
  for (const p of s.projects) meta.push(`<a href="https://www.ebi.ac.uk/ena/browser/view/${esc(p)}">${esc(p)}</a>`);
  $(".meta", el).innerHTML = meta.join("");
  const b = $(".badge", el); b.textContent = s.type; b.classList.add(s.type);

  const ds = visibleDatasets(s);
  const chem = [...new Set(ds.map((d) => d.chemistry).filter(Boolean))];
  const fc = [...new Set(ds.map((d) => d.run?.flowcell_id).filter(Boolean))];
  const facts = [];
  const ki = s.kit_inferred;
  if (s.library_kit.length) facts.push(`<span class="fact"><b>Kit</b> ${esc(s.library_kit.join(", "))}</span>`);
  if (ki) facts.push(`<span class="fact" title="Inferred from ${ki.reads} sampled reads: ${esc(ki.evidence)} (confidence: ${esc(ki.confidence)})"><b>Kit from reads</b> ${esc(ki.kit)}${ki.barcode ? " · " + esc(ki.barcode) : ""}</span>`);
  else if (!s.library_kit.length) facts.push(`<span class="fact warn">Kit not recorded yet</span>`);
  // leftover barcode/primer warnings are held back with the read-artefact page until verified
  if (chem.length) facts.push(`<span class="fact"><b>Chemistry</b> ${esc(chem.join(", "))}</span>`);
  if (fc.length) facts.push(`<span class="fact"><b>Flow cell</b> ${esc(fc.join(", "))}</span>`);
  if (s.hic.length) facts.push(`<span class="fact"><b>Hi-C</b> ${esc(s.hic[0].kit)}</span>`);
  $(".facts", el).innerHTML = facts.join("");

  const tbody = $("tbody", el);
  tbody.innerHTML = ds.map((d, i) => {
    const key = `${s.name}|${d.dorado}|${d.model_version}`;
    const merged = mergedFiles(d);
    const pod5Inside = d.pod5 && merged.some((f) => f.url === d.pod5.url);
    if (pod5Inside) merged.find((f) => f.url === d.pod5.url).labels.push("pod5");
    const tiers = merged.map(fileChip).join("") + (d.pod5 && !pod5Inside ? `<a class="file raw" href="${esc(d.pod5.url)}" download title="Raw signal (pod5 tar.gz)">pod5 <small>${fmtBytes(d.pod5.bytes)}</small></a>` : "");
    const qf = ["mods", "sup", "hac", "fast"].map((t) => d.files.find((f) => f.tier === t && f.median_q_est != null)).find(Boolean);
    const q = qf?.median_q_est;
    const isOpen = open.has(key);
    return `<tr>
      <td class="bc"><span class="model">v${esc(d.model_version)}</span> <span class="sub">dorado ${esc(d.dorado)}</span></td>
      <td class="num" data-label="Reads">${fmtN(d.reads)}</td><td class="num" data-label="Yield">${fmtBases(d.bases)}</td>
      <td class="num" data-label="Mean len">${d.mean_len ? fmtLen(d.mean_len) : "–"}</td>
      <td class="num" data-label="N50" title="${d.exact ? "exact (whole BAM read)" : `estimated from ${d.sampled_reads?.toLocaleString()} sampled reads`}">${approx(d.n50_est ? fmtLen(d.n50_est) : null, d.exact)}</td>
      <td class="num" data-label="Read Q">${approx(q, d.exact)}</td>
      <td class="fcell"><div class="files">${tiers}</div></td>
      <td class="mcell"><button class="more" data-key="${esc(key)}" aria-expanded="${isOpen}">${isOpen ? "less" : "settings"}</button></td>
    </tr>${isOpen ? detailsRow(d, key) : ""}`;
  }).join("");

  if (s.hic.length) {
    $(".hic", el).innerHTML = `<h4>Hi-C</h4>` + s.hic.map((h) =>
      `<div>${esc(h.kit)} · ${esc(h.platform)} · ${fmtN(h.reads)} read pairs · ${fmtBases(h.bases)} ·
        <a href="https://www.ebi.ac.uk/ena/browser/view/${esc(h.accession)}">${esc(h.accession)}</a> ·
        <a class="file" href="${esc(h.r1)}" download>R1</a> <a class="file" href="${esc(h.r2)}" download>R2</a></div>`).join("");
  }
  return el;
}

function render() {
  const list = $("#list");
  const items = filtered();
  list.innerHTML = "";
  if (!items.length) { list.innerHTML = `<p class="empty">No samples match these filters.</p>`; }
  let lastType = null;
  for (const s of items) {
    if (state.sort === "type" && s.type !== lastType) {
      lastType = s.type;
      const n = items.filter((x) => x.type === s.type).length;
      const h = document.createElement("h2");
      h.className = "group-title";
      h.innerHTML = `${esc(TYPES.find((t) => t[0] === s.type)?.[1] || s.type)} <small>${n}</small>`;
      list.appendChild(h);
    }
    list.appendChild(card(s));
  }
  const nds = items.reduce((a, s) => a + visibleDatasets(s).length, 0);
  $("#count").textContent = `${items.length} of ${DATA.samples.length} samples · ${nds} basecalled datasets`;
  writeHash();
  refreshDl();
}

function renderStats() {
  const s = DATA.samples;
  let bases = 0;
  // every distinct file once (several tiers can point to the same file), Hi-C included
  const uniq = new Map();
  for (const x of s) {
    for (const d of x.datasets) {
      for (const f of d.files) uniq.set(f.url, f.bytes || 0);
      if (d.pod5) uniq.set(d.pod5.url, d.pod5.bytes || 0);
    }
    for (const h of x.hic) { uniq.set(h.r1, h.r1_bytes || 0); uniq.set(h.r2, h.r2_bytes || 0); }
  }
  const files = uniq.size, bytes = [...uniq.values()].reduce((a, b) => a + b, 0);
  const runs = {};
  for (const x of s) for (const d of x.datasets) { const k = d.pod5?.url || x.name; runs[k] = Math.max(runs[k] || 0, d.bases || 0); }
  bases = Object.values(runs).reduce((a, b) => a + b, 0);
  const t = (k) => s.filter((x) => x.type === k).length;
  $("#stats").innerHTML = [
    [s.length, "samples"], [t("monoculture"), "monocultures"], [t("mock"), "mock communities"], [t("metagenome"), "metagenomes"],
    [(bases / 1e12).toFixed(2) + " Tb", "unique sequence"], [files, "downloadable files"], [(bytes / 1e12).toFixed(1) + " TB", "total file size"],
  ].map(([v, l]) => `<div class="stat"><b>${v}</b><span>${l}</span></div>`).join("");
}

function buildControls() {
  $("#typeChips").innerHTML = TYPES.map(([k, l]) =>
    `<button class="chip" data-type="${k}" aria-pressed="${state.types.has(k)}">${l}<small>${DATA.samples.filter((s) => s.type === k).length}</small></button>`).join("");
  const models = [...new Set(DATA.samples.flatMap((s) => s.datasets.map((d) => d.model_version)))].sort().reverse();
  $("#model").insertAdjacentHTML("beforeend", models.map((m) => `<option value="${esc(m)}">v${esc(m)}</option>`).join(""));
  $("#q").value = state.q; $("#model").value = state.model; $("#has").value = state.has; $("#sort").value = state.sort;

  $("#q").addEventListener("input", (e) => { state.q = e.target.value; render(); });
  $("#model").addEventListener("change", (e) => { state.model = e.target.value; render(); });
  $("#has").addEventListener("change", (e) => { state.has = e.target.value; render(); });
  $("#sort").addEventListener("change", (e) => { state.sort = e.target.value; render(); });
  $("#typeChips").addEventListener("click", (e) => {
    const b = e.target.closest(".chip"); if (!b) return;
    const k = b.dataset.type;
    state.types.has(k) ? state.types.delete(k) : state.types.add(k);
    b.setAttribute("aria-pressed", state.types.has(k));
    render();
  });
  $("#list").addEventListener("click", (e) => {
    const b = e.target.closest(".more"); if (!b) return;
    const k = b.dataset.key; open.has(k) ? open.delete(k) : open.add(k);
    const cardEl = b.closest(".card");
    const s = DATA.samples.find((x) => slug(x.name) === cardEl.id);
    cardEl.replaceWith(card(s));
  });
  bindDownloads();
}

// ---------- download script ----------
const dirName = (s) => s.replace(/[^A-Za-z0-9._-]+/g, "_").replace(/^_+|_+$/g, "");
function selectedFiles() {
  const out = new Map();
  const add = (url, dir, md5, bytes) => { if (url && !out.has(url)) out.set(url, { url, dir, md5, bytes, name: url.split("/").pop() }); };
  const fileOk = (f) => !state.has || f.tier === state.has || (state.has === "dup" && f.tier.endsWith("dup"));
  for (const s of filtered()) {
    const dir = dirName(s.name);
    for (const d of visibleDatasets(s)) {
      for (const f of d.files) if (fileOk(f)) add(f.url, dir, f.md5, f.bytes);
      if (d.pod5 && (!state.has || state.has === "pod5")) add(d.pod5.url, dir + "/pod5", d.pod5.md5, d.pod5.bytes);
    }
    if (!state.has || state.has === "hic") for (const h of s.hic) { add(h.r1, dir + "/hic", h.r1_md5, h.r1_bytes); add(h.r2, dir + "/hic", h.r2_md5, h.r2_bytes); }
  }
  return [...out.values()];
}
const fmtSize = (b) => b >= 1e12 ? (b / 1e12).toFixed(2) + " TB" : b >= 1e9 ? (b / 1e9).toFixed(1) + " GB" : (b / 1e6).toFixed(0) + " MB";
const shq = (s) => "'" + String(s).replace(/'/g, "'\\''") + "'";
function filterText() {
  const p = [];
  if (state.q) p.push(`search "${state.q}"`);
  if (state.types.size) p.push("type " + [...state.types].join("+"));
  if (state.model) p.push("model v" + state.model);
  if (state.has) p.push("data " + state.has);
  return p.join(", ") || "no filters (everything)";
}
function buildScript(tool, md5) {
  const files = selectedFiles();
  const total = files.reduce((a, f) => a + (f.bytes || 0), 0);
  const L = [];
  L.push("#!/usr/bin/env bash");
  L.push(`# MicroBench download: ${files.length} files, ${fmtSize(total)}`);
  L.push(`# Selection: ${filterText()}`);
  L.push(`# Page: ${location.href}`);
  L.push(`# Generated ${new Date().toISOString().slice(0, 10)} · data from ENA (https://www.ebi.ac.uk/ena)`);
  L.push("# Usage: bash microbench_download.sh [target folder]   (re-run to resume unfinished downloads)");
  if (total > 1e12) L.push(`# NOTE: this is ${fmtSize(total)}. Check free disk space first.`);
  L.push("set -euo pipefail");
  L.push('OUT="${1:-microbench}"');
  L.push('mkdir -p "$OUT" && cd "$OUT"');
  L.push("");
  if (tool === "aria2c") {
    L.push("# aria2c: 4 files at a time, 8 connections each; resumes with -c");
    L.push("aria2c -c -x 8 -s 8 -j 4 --auto-file-renaming=false -i - <<'EOF'");
    for (const f of files) { L.push(f.url); L.push(`  dir=${f.dir}`); if (md5 && f.md5) L.push(`  checksum=md5=${f.md5}`); }
    L.push("EOF");
  } else {
    L.push("# get <folder> <url> <bytes>: skips files that are already complete, resumes partial ones");
    L.push("get() {");
    L.push('  local f="$1/${2##*/}"; mkdir -p "$1"');
    L.push('  if [ -n "${3:-}" ] && [ -f "$f" ] && [ "$(wc -c < "$f" | tr -d " ")" = "$3" ]; then echo "already complete: $f"; return; fi');
    if (tool === "wget") L.push('  wget -c -q --show-progress -P "$1" "$2"');
    else L.push('  (cd "$1" && curl -fL --retry 5 -C - -O "$2")');
    L.push("}");
    for (const f of files) L.push(`get ${shq(f.dir)} ${shq(f.url)} ${f.bytes || ""}  # ${f.bytes ? fmtSize(f.bytes) : "size unknown"}`);
  }
  if (md5) {
    const withSum = files.filter((f) => f.md5);
    L.push("");
    L.push(`# Check MD5 sums (from ENA) for ${withSum.length} files`);
    L.push("cat > md5sums.txt <<'EOF'");
    for (const f of withSum) L.push(`${f.md5}  ${f.dir}/${f.name}`);
    L.push("EOF");
    L.push("if command -v md5sum >/dev/null; then md5sum -c md5sums.txt");
    L.push('else while read -r sum file; do [ "$(md5 -q "$file")" = "$sum" ] && echo "$file: OK" || echo "$file: FAILED"; done < md5sums.txt; fi');
  }
  return { text: L.join("\n") + "\n", files, total };
}
function dlOptions() {
  return { tool: document.querySelector('input[name="dltool"]:checked').value, md5: $("#dlMd5").checked };
}
function refreshDl() {
  if (!$("#dlmenu").open) return;
  const { tool, md5 } = dlOptions();
  const { text, files, total } = buildScript(tool, md5);
  $("#dlSum").innerHTML = files.length
    ? `<b>${files.length} files · ${fmtSize(total)}</b> ${state.q || state.types.size || state.model || state.has ? `matching the current filters (${esc(filterText())})` : "— everything; set filters above to narrow it down"}. One folder per sample, MD5 sums from ENA.`
    : `<b>No files</b> match the current filters.`;
  $("#dlPreview").textContent = text.split("\n").slice(0, 18).join("\n") + (text.split("\n").length > 18 ? "\n…" : "");
  $("#dlCopy").disabled = $("#dlSave").disabled = !files.length;
}
function toast(msg) {
  const t = Object.assign(document.createElement("div"), { className: "toast", textContent: msg });
  document.body.appendChild(t); setTimeout(() => t.remove(), 2600);
}
async function copyText(text) {
  try { await navigator.clipboard.writeText(text); return true; }
  catch { const ta = Object.assign(document.createElement("textarea"), { value: text }); document.body.appendChild(ta); ta.select();
    const ok = document.execCommand("copy"); ta.remove(); return ok; }
}
function saveFile(text, name, type = "text/plain") {
  const a = Object.assign(document.createElement("a"), { href: URL.createObjectURL(new Blob([text], { type })), download: name });
  document.body.appendChild(a); a.click(); a.remove();
}
function bindDownloads() {
  $("#dlmenu").addEventListener("toggle", refreshDl);
  $("#dlmenu").addEventListener("change", refreshDl);
  document.addEventListener("click", (e) => { if (e.isTrusted && $("#dlmenu").open && !e.target.closest("#dlmenu")) $("#dlmenu").open = false; });
  $("#dlCopy").addEventListener("click", async () => {
    const { tool, md5 } = dlOptions(); const { text, files, total } = buildScript(tool, md5);
    toast((await copyText(text)) ? `Copied ${tool} script for ${files.length} files (${fmtSize(total)})` : "Copy failed: use Save .sh");
  });
  $("#dlSave").addEventListener("click", () => { const { tool, md5 } = dlOptions(); saveFile(buildScript(tool, md5).text, "microbench_download.sh", "text/x-sh"); });
  $("#dl").addEventListener("click", () => saveFile(selectedFiles().map((f) => f.url).join("\n") + "\n", "microbench_urls.txt"));
}

// filters live in the query string so a filtered view can be shared; #anchor jumps to a sample
function readHash() {
  const p = new URLSearchParams(location.search);
  state.q = p.get("q") || ""; state.model = p.get("model") || ""; state.has = p.get("has") || ""; state.sort = p.get("sort") || "type";
  (p.get("type") || "").split(",").filter(Boolean).forEach((t) => state.types.add(t));
}
function writeHash() {
  const p = new URLSearchParams();
  if (state.q) p.set("q", state.q);
  if (state.types.size) p.set("type", [...state.types].join(","));
  if (state.model) p.set("model", state.model);
  if (state.has) p.set("has", state.has);
  if (state.sort !== "type") p.set("sort", state.sort);
  const qs = p.toString();
  history.replaceState(null, "", (qs ? "?" + qs : location.pathname) + location.hash);
}

fetch("data/microbench.json").then((r) => r.json()).then((d) => {
  DATA = d;
  for (const s of DATA.samples) s._hay = haystack(s);
  readHash();
  renderStats();
  buildControls();
  render();
  $("#nchunks").textContent = d.n_chunks;
  $("#generated").textContent = `Data table generated ${d.generated} by scripts/build_site_data.py.`;
  if (location.hash) document.getElementById(decodeURIComponent(location.hash.slice(1)))?.scrollIntoView();
}).catch((e) => { $("#list").innerHTML = `<p class="empty">Could not load data (${esc(e.message)}).</p>`; });
