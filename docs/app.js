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
  const parts = [s.name, s.type, s.info, ...s.projects, ...s.library_kit];
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
  facts.push(s.library_kit.length ? `<span class="fact"><b>Kit</b> ${esc(s.library_kit.join(", "))}</span>`
    : `<span class="fact warn">Kit not recorded yet</span>`);
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
}

function renderStats() {
  const s = DATA.samples;
  let files = 0, bases = 0, bytes = 0;
  const seen = new Set();
  for (const x of s) for (const d of x.datasets) {
    for (const f of d.files) { files++; bytes += f.bytes || 0; }
    if (d.pod5 && !seen.has(d.pod5.url)) { seen.add(d.pod5.url); files++; bytes += d.pod5.bytes || 0; }
    // count raw sequence once per sample + run: take the max sup yield per pod5 archive
  }
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
  $("#dl").addEventListener("click", () => {
    const urls = [];
    for (const s of filtered()) for (const d of visibleDatasets(s)) {
      for (const f of d.files) if (!state.has || f.tier === state.has || (state.has === "dup" && f.tier.endsWith("dup"))) urls.push(f.url);
      if (d.pod5 && (!state.has || state.has === "pod5")) urls.push(d.pod5.url);
    }
    if (!state.has || state.has === "hic") for (const s of filtered()) for (const h of s.hic) urls.push(h.r1, h.r2);
    const blob = new Blob([[...new Set(urls)].join("\n") + "\n"], { type: "text/plain" });
    const a = Object.assign(document.createElement("a"), { href: URL.createObjectURL(blob), download: "microbench_urls.txt" });
    document.body.appendChild(a); a.click(); a.remove();
  });
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
