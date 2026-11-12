"use strict";
const $ = (s, el = document) => el.querySelector(s);
const esc = (s) => String(s ?? "").replace(/[&<>"']/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
const SVGNS = "http://www.w3.org/2000/svg";
const fmtLen = (n) => n >= 1000 ? (n / 1000).toFixed(n >= 10000 ? 0 : 1) + " kb" : n + " bp";
const pct = (x) => x.toFixed(x < 1 ? 2 : 1) + "%";

const SPECIES_COLOR = {
  Bacillus_subtilis: "#4e79a7", Enterococcus_faecalis: "#f28e2b", Escherichia_coli: "#e15759",
  Listeria_monocytogenes: "#76b7b2", Pseudomonas_aeruginosa: "#59a14f", Saccharomyces_cerevisiae: "#c9a227",
  Salmonella_enterica: "#b07aa1", Staphylococcus_aureus: "#ff9da7",
};
const spName = (s) => s.replace(/_/g, " ").replace(/^(\w)\w+ /, "$1. ");

const CATS = [
  { key: "chimera_interspecies", title: "Chimeras between species", short: "Inter-species chimeras",
    what: "One read contains two (or more) molecules from different genomes joined end to end. In a mock community the join is unambiguous: the halves align to different species.",
    detect: "Two or more non-overlapping alignments (≥300 bp each) to different species.",
    matters: "Read-level analyses take a chimera at face value: read-based taxonomy, linking plasmids or mobile elements to their host, structural variant calling and per-read methylation profiles. A chimera between two species looks exactly like evidence of horizontal transfer.",
    tools: "Most junctions have no adapter between the two parts, so adapter-based splitting (dorado, Porechop) has nothing to find. Whole-read length/quality filters keep them because both halves are good sequence." },
  { key: "chimera_intraspecies", title: "Chimeras within a genome", short: "Intra-genome chimeras",
    what: "Two fragments of the same genome joined in an order or orientation that does not exist in the genome.",
    detect: "Adjacent alignments on the same contig that are not collinear (wrong distance or strand), allowing for reads that cross the origin of a circular chromosome.",
    matters: "Fake rearrangements and inversions in structural-variant calls. Harmless for high-coverage assembly, but a risk at low coverage, e.g. for rare organisms in a metagenome.",
    tools: "Same as inter-species chimeras: usually no adapter at the join and high read quality. Hard to recognise without a reference." },
  { key: "foldback", title: "Foldback reads (inverted duplicates)", short: "Foldbacks",
    what: "The read continues from one strand into the reverse complement of the same region, so the read is a near-palindrome. The second half is the complement strand passing through the pore after the template.",
    detect: "Adjacent alignments to the same place in the genome on opposite strands.",
    matters: "Each base is counted twice, which inflates coverage and creates fake inverted repeats. The second half is usually less accurate.",
    tools: "Some foldbacks contain the ligation adapter at the fold (see the chimera section), so adapter splitting can catch those; the rest join without an adapter. In the ZymoHMW run foldbacks became rarer during the run (4.8% of reads ≥1 kb at the start, about 2% later) and more common within each scan cycle, so they depend on the state of the pores, not just the library." },
  { key: "internal_adapter", title: "Adapters inside reads", short: "Internal adapters",
    what: "A copy of the ligation adapter sits inside the read instead of at the very end, typically after a short leading fragment.",
    detect: "Match to the ONT ligation adapter (CCTGTACTTCGTTCAGTTACGTATTGCT, ≤20% edits) more than 100 bp from either end.",
    matters: "Adapter sequence in assemblies and k-mer databases. It usually marks a short junk or chimeric piece in front of the real insert.",
    tools: "dorado trims adapters only at the extreme ends. Porechop/Porechop_ABI search inside reads and split or discard them." },
  { key: "end_adapter", title: "Adapters left at read ends", short: "End adapters",
    what: "Adapter sequence within 100 bp of a read end that dorado's default trimming did not remove.",
    detect: "Adapter match within 100 bp of the read start or end.",
    matters: "Rare after dorado trimming (see the numbers), but a common problem when reads were basecalled with --no-trim or by older software.",
    tools: "End trimming by dorado or Porechop_ABI removes it." },
  { key: "low_complexity", title: "Low-complexity junk", short: "Low-complexity junk",
    what: "Hundreds of bases of simple repeats — (GGAAT)n, poly-A, poly-T, (TA)n — inside or at the end of otherwise good reads. They align to nothing in the mock and appear as long insertions.",
    detect: "Runs of ≥6 GGAAT units, ≥25 A or T, or ≥12 dinucleotide units in a read with at least one alignment.",
    matters: "In the ZymoHMW run, stretches like these cluster right after MinKNOW's periodic pore scan (see the scan section). (GGAAT)n is the motif of human satellite III, so these reads can trigger human-contamination screens. Junk insertions also end up in low-coverage assemblies and inflate error rates. In this run the repeat appears almost only as GGAAT and practically never as its reverse complement ATTCC, which points to the basecaller rather than real DNA.",
    tools: "Mean-quality filters average the junk away inside a long read. Only tools that look at local quality or sequence complexity can remove it." },
  { key: "post_scan_junk", title: "Junk right after a pore scan", short: "Post-scan junk",
    what: "Reads that start within a few minutes after MinKNOW's periodic scan often begin with hundreds of bases of (GGAAT)n or poly-A at very low quality (Q≈2–3), followed by normal sequence.",
    detect: "A low-complexity run of ≥100 bp in the first 5% of a read that started 0–330 s after a scan began.",
    matters: "It happens on the whole flow cell at the same time, every 1.6 hours. The junk sits where adapter trimmers look, so it can push the real adapter or the start of the insert out of their search window.",
    tools: "Local quality trimming at read starts removes it. Mean-quality filters average it away in long reads." },
  { key: "skip_prone", title: "Skip-prone reads", short: "Skip-prone reads",
    what: "Reads with five or more deletions of ≥10 bp: stretches of the molecule that never appear in the basecalls, as if the strand slipped through the pore.",
    detect: "≥5 deletions of ≥10 bp that are not shared by other reads (shared ones are real differences from the reference).",
    matters: "A source of false indels and frameshifts in low-coverage regions.",
    tools: "Their mean Q-score is lower but most still pass a Q10 filter." },
  { key: "signal_gap", title: "Signal gaps", short: "Signal gaps",
    what: "A stretch of the molecule missing from every read in progress at one moment.", detect: "See the signal-gap section.", matters: "", tools: "" },
  { key: "unmapped", title: "Unmapped reads", short: "Unmapped",
    what: "Reads with no alignment of ≥300 bp. Almost all are short.",
    detect: "No alignment to the mock references.",
    matters: "Many reads but very few bases. Mostly harmless.",
    tools: "Any length filter removes them." },
];
const CAT = Object.fromEntries(CATS.map((c) => [c.key, c]));

const CORE = ["chimera_interspecies", "chimera_intraspecies", "foldback", "internal_adapter", "end_adapter", "low_complexity", "skip_prone"];
let D, cur = "chimera_interspecies", idx = 0, unit = "pct_bases";

function svg(tag, attrs = {}, parent) {
  const el = document.createElementNS(SVGNS, tag);
  for (const [k, v] of Object.entries(attrs)) el.setAttribute(k, v);
  if (parent) parent.appendChild(el);
  return el;
}

function headline() {
  const st = D.stats, n = D.numbers;
  const any = n.classes.signal_gap.reads ? n.any_core_or_gap : n.any_core;
  const tiles = [
    [pct(any.pct_reads), "of reads have at least one read artefact", `${pct(any.pct_bases)} of all bases are in these reads`],
    [pct(n.chimeric_or_foldback.pct_reads), "of reads are chimeras or foldbacks", `${pct(n.chimeric_or_foldback.pct_bases)} of all bases`],
    ...(n.classes.signal_gap.reads ? [[pct(n.classes.signal_gap.pct_reads), "of reads lost a stretch of sequence in a signal gap", `${D.gaps.n_events.toLocaleString()} gaps in the run`]] : []),
    [pct(n.split_pieces.pct_reads), "of reads are pieces of longer reads that dorado already split", "dorado's own fix; everything else on this page is left after it"],
    [st.n_reads.toLocaleString(), "reads analysed", "every read of the run"],
  ];
  $("#headline").innerHTML = tiles.map(([v, l, sub]) => `<div class="stat"><b>${v}</b><span>${l}</span><small>${sub}</small></div>`).join("");
  const ds = D.dataset;
  $("#datasetLine").innerHTML = `${esc(ds.sample)} · flow cell ${esc(ds.flowcell)} · ${esc(ds.basecaller)} · all ${st.n_reads.toLocaleString()} reads (${(st.n_bases / 1e9).toFixed(1)} Gb) from ${st.n_channels.toLocaleString()} channels in <a href="https://www.ebi.ac.uk/ena/browser/view/${esc(ds.accession)}">${esc(ds.accession)}</a>.`;
  const bc = D.barcode_note;
  $("#barcodeNote").innerHTML = bc ? `<b>Before anything else:</b> ${esc(bc)}` : "";
  $("#barcodeNote").hidden = !bc;
  const inAny = CATS.filter((c) => CORE.includes(c.key) || (c.key === "signal_gap" && n.classes.signal_gap.reads)).map((c) => c.short.toLowerCase());
  $("#overviewNote").innerHTML = `Each bar is the share of all reads, or of all bases, in reads that have that artefact (the whole read is counted). A read can have several artefacts, so the bars do not add up: ${pct(any.pct_reads)} of reads (${pct(any.pct_bases)} of bases) have at least one of ${inAny.slice(0, -1).join(", ")} or ${inAny.slice(-1)}. Unmapped reads are not counted as an artefact. Click a bar for examples.`;
}

function overview() {
  const host = $("#overviewChart");
  const rows = CATS.filter((c) => D.summary[c.key] && c.key !== "untrimmed_barcode" && !(c.key === "signal_gap" && !D.summary.signal_gap.reads))
    .map((c) => ({ c, v: D.summary[c.key][unit], n: D.summary[c.key].reads }));
  const W = 760, rowH = 26, left = 190, H = rows.length * rowH + 24;
  const max = Math.max(...rows.map((r) => r.v)) * 1.15;
  const s = svg("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": "Artefact prevalence" });
  rows.forEach((r, i) => {
    const g = svg("g", { class: "bar-row", transform: `translate(0 ${i * rowH + 4})` }, s);
    svg("text", { x: left - 8, y: 15, "text-anchor": "end" }, g).textContent = r.c.short;
    const w = Math.max(2, (W - left - 90) * r.v / max);
    svg("rect", { class: "bar", x: left, y: 3, width: w, height: 16, rx: 3, fill: r.c.key === "unmapped" ? "var(--muted)" : "var(--bad)" }, g);
    svg("text", { x: left + w + 6, y: 15 }, g).textContent = `${pct(r.v)}  (${r.n.toLocaleString()} reads)`;
    g.addEventListener("click", () => {
      if (!D.examples[r.c.key]?.length) { (r.c.key === "signal_gap" ? $("#gaps") : $("#gallery")).scrollIntoView({ behavior: "smooth" }); return; }
      select(r.c.key); $("#gallery").scrollIntoView({ behavior: "smooth" });
    });
  });
  host.replaceChildren(s);
}

function tabs() {
  $("#tabs").innerHTML = CATS.filter((c) => D.examples[c.key]?.length).map((c) =>
    `<button class="chip" role="tab" data-k="${c.key}" aria-pressed="${c.key === cur}">${esc(c.short)}${D.summary[c.key] ? `<small>${pct(D.summary[c.key].pct_reads)}</small>` : ""}</button>`).join("");
}

function select(k) { cur = k; idx = 0; tabs(); renderCat(); }

function renderCat() {
  const c = CAT[cur], s = D.summary[cur], f = D.filters?.[cur];
  const head = s ? `<div><span class="big">${pct(s.pct_reads)}</span> of reads · <b>${pct(s.pct_bases)}</b> of bases · ${s.reads.toLocaleString()} reads</div>`
    : cur === "post_scan_junk" ? `<div><span class="big">${D.scan.n_post_scan_reads.toLocaleString()}</span> reads of 1.5–30 kb, across ${D.scan.scan_start_h.length - 1} scans</div>` : "";
  $("#catInfo").innerHTML = `<h3>${esc(c.title)}</h3>
    ${head}
    <dl><dt>What it is</dt><dd>${esc(c.what)}</dd>
    <dt>How it was detected</dt><dd>${esc(c.detect)}</dd>
    <dt>Why it matters</dt><dd>${esc(c.matters)}</dd>
    <dt>Tools</dt><dd>${esc(c.tools)}${f ? ` Here, ${f.q10_l1000}% of these reads pass Q≥10 and ≥1 kb (median read Q ${f.median_qs}).` : ""}</dd></dl>`;
  renderExample();
}

function renderExample() {
  const list = D.examples[cur] || [];
  if (!list.length) { $("#readSvg").textContent = "No examples."; return; }
  idx = (idx + list.length) % list.length;
  const r = list[idx];
  $("#exCount").textContent = `Example ${idx + 1} of ${list.length}`;
  $("#exMeta").textContent = `${r.id.slice(0, 8)}… · ${fmtLen(r.len)} · read Q ${r.qs} · channel ${r.ch}${r.split_from ? " · split by dorado" : ""}`;
  drawRead(r);
  details(r);
}

function drawRead(r) {
  const cw = $("#readSvg").clientWidth || 700;
  const W = Math.round(Math.min(1000, Math.max(480, cw * 1.35))), H = 172, x0 = 10, x1 = W - 10, sx = (q) => x0 + (x1 - x0) * q / r.len;
  const s = svg("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": "Read layout" });
  // axis
  const step = niceStep(r.len / (W < 700 ? 4 : 6));
  for (let p = 0; p <= r.len; p += step) {
    svg("line", { x1: sx(p), x2: sx(p), y1: 14, y2: 20, stroke: "var(--muted)" }, s);
    svg("text", { x: sx(p), y: 11, "text-anchor": p === 0 ? "start" : "middle", fill: "var(--muted)" }, s).textContent = p === 0 ? "0" : fmtLen(p);
  }
  // read backbone
  svg("rect", { x: x0, y: 30, width: x1 - x0, height: 26, rx: 4, fill: "var(--surface-2)", stroke: "var(--border)" }, s);
  // segments
  for (const g of r.segs) {
    const col = SPECIES_COLOR[g.sp] || "#888";
    const a = sx(g.q0), b = sx(g.q1), w = Math.max(2, b - a);
    const rect = svg("rect", { x: a, y: 30, width: w, height: 26, rx: 3, fill: col }, s);
    svg("title", {}, rect).textContent = `${spName(g.sp)} ${g.strand} ${g.ctg}:${g.r0.toLocaleString()}-${g.r1.toLocaleString()} · identity ${g.ident}%`;
    // strand chevrons
    for (let x = a + 12; x < b - 8; x += 26) {
      const d = g.strand === "+" ? `M${x - 4} 37 L${x + 3} 43 L${x - 4} 49` : `M${x + 3} 37 L${x - 4} 43 L${x + 3} 49`;
      svg("path", { d, stroke: "rgba(255,255,255,.75)", "stroke-width": 2, fill: "none" }, s);
    }
    if (w > 150) svg("text", { x: a + 4, y: 74, fill: "var(--text)" }, s).textContent = `${spName(g.sp)} (${g.strand}) ${(g.r0 / 1e6).toFixed(3)} Mb · ${g.ident}%`;
    for (const [t, q, ln] of g.indels) {
      svg("line", { x1: sx(q), x2: sx(q), y1: t === "del" ? 56 : 24, y2: t === "del" ? 62 : 30, stroke: "var(--text)", "stroke-width": 1.5 }, s);
    }
  }
  // low complexity
  for (const [a, b] of r.lowc) {
    const rect = svg("rect", { x: sx(a), y: 80, width: Math.max(2, sx(b) - sx(a)), height: 10, fill: "var(--junk)" }, s);
    svg("title", {}, rect).textContent = `low-complexity ${b - a} bp`;
  }
  // adapters
  for (const ad of r.adapters) {
    const x = sx(ad.start);
    const p = svg("path", { d: `M${x} 78 L${x + 6} 88 L${x - 6} 88 Z`, fill: "var(--adapter)" }, s);
    svg("title", {}, p).textContent = `adapter at ${ad.start} (${ad.ed} edits)`;
  }
  // quality track
  const qy = (q) => 168 - Math.min(q, 45) * 1.4;
  svg("line", { x1: x0, x2: x1, y1: qy(10), y2: qy(10), stroke: "var(--muted)", "stroke-dasharray": "3 3" }, s);
  svg("text", { x: x1, y: qy(10) - 3, "text-anchor": "end", fill: "var(--muted)" }, s).textContent = "Q10";
  svg("text", { x: x0, y: 104, fill: "var(--muted)" }, s).textContent = "read quality along the read";
  const pts = r.qwin.map((q, i) => `${sx(Math.min(r.len, (i + 0.5) * r.qw)).toFixed(1)},${qy(q).toFixed(1)}`).join(" ");
  svg("polyline", { points: pts, fill: "none", stroke: "var(--q)", "stroke-width": 1.3 }, s);
  $("#readSvg").replaceChildren(s);
  const sps = [...new Set(r.segs.map((g) => g.sp))];
  $("#legend").innerHTML = sps.map((sp) => `<span><i style="background:${SPECIES_COLOR[sp] || "#888"}"></i>${esc(spName(sp))}</span>`).join("") +
    `<span><i style="background:var(--adapter)"></i>adapter</span><span><i style="background:var(--junk)"></i>low-complexity</span><span>│ ticks: indels ≥20 bp</span><span>› strand</span>`;
}

function niceStep(x) {
  const p = Math.pow(10, Math.floor(Math.log10(x))); const m = x / p;
  return (m < 1.5 ? 1 : m < 3.5 ? 2 : m < 7.5 ? 5 : 10) * p;
}

function colorSeq(seq, lo, marks) {
  // marks: [{a, b, cls}] in read coordinates
  let out = "";
  for (let i = 0; i < seq.length; i++) {
    const p = lo + i, m = marks.find((x) => p >= x.a && p < x.b);
    out += m ? `<span class="${m.cls}"${m.style ? ` style="${m.style}"` : ""}>${seq[i]}</span>` : seq[i];
  }
  return out.replace(/<\/span><span class="(\w+)">/g, (all, c) => all); // keep simple
}

function details(r) {
  const parts = [];
  parts.push(`<div class="table-wrap"><table class="plain segtable"><thead><tr><th>Read part</th><th>Aligns to</th><th>Strand</th><th class="num">Ref position</th><th class="num">Identity</th></tr></thead><tbody>` +
    r.segs.map((g) => `<tr><td>${g.q0.toLocaleString()}–${g.q1.toLocaleString()}</td><td><i style="display:inline-block;width:9px;height:9px;border-radius:2px;background:${SPECIES_COLOR[g.sp] || "#888"}"></i> ${esc(spName(g.sp))}</td><td>${g.strand}</td><td class="num">${g.r0.toLocaleString()}–${g.r1.toLocaleString()}</td><td class="num">${g.ident}%</td></tr>`).join("") +
    (r.segs.length ? "" : `<tr><td colspan="5" class="muted">no alignment</td></tr>`) + `</tbody></table></div>`);
  for (const j of r.junctions) {
    const ri = r.junctions.indexOf(j);
    const cl = SPECIES_COLOR[r.segs[ri].sp] || "#888", cr = SPECIES_COLOR[r.segs[ri + 1].sp] || "#888";
    const marks = [{ a: j.lo + j.left_end, b: j.lo + j.right_start, cls: "gap" },
      { a: j.lo, b: j.lo + j.left_end, cls: "fl", style: `--c:${cl}` },
      { a: j.lo + j.right_start, b: j.lo + j.seq.length, cls: "fl", style: `--c:${cr}` }];
    for (const ad of r.adapters) if (ad.start < j.lo + j.seq.length && ad.end > j.lo) marks.unshift({ a: ad.start, b: ad.end, cls: "ad" });
    parts.push(`<div class="seqbox"><small>Junction at ${j.at.toLocaleString()} · ${j.gap > 0 ? j.gap + " bp between the aligned parts" : j.gap < 0 ? "aligned parts overlap by " + (-j.gap) + " bp" : "parts joined directly"}</small>\n${colorSeq(j.seq, j.lo, marks)}</div>`);
  }
  const shownAd = r.adapters.filter((ad) => !r.junctions.some((j) => ad.start < j.lo + j.seq.length && ad.end > j.lo)).slice(0, 3);
  for (const ad of shownAd) {
    parts.push(`<div class="seqbox"><small>Adapter at ${ad.start.toLocaleString()} (${ad.ed} edits)</small>\n${colorSeq(ad.ctx, ad.ctx_lo, [{ a: ad.start, b: ad.end, cls: "ad" }])}</div>`);
  }
  for (const [a, b, m, seq] of r.lowc.slice(0, 3)) {
    parts.push(`<div class="seqbox"><small>Low-complexity run at ${a.toLocaleString()}, ${b - a} bp (first 60 shown)</small>\n<span class="lc">${esc(seq)}</span>${b - a > 60 ? "…" : ""}</div>`);
  }
  $("#exDetails").innerHTML = parts.join("");
}

function gapChart() {
  const h = D.stats.chimera_gap_hist; // [[bin, n]]
  const W = 460, H = 200, l = 40, b = 30, t = 10;
  const max = Math.max(...h.map((x) => x[1]));
  const bw = (W - l - 10) / h.length;
  const s = svg("svg", { viewBox: `0 0 ${W} ${H}` });
  h.forEach(([bin, n], i) => {
    const hh = (H - b - t) * n / max;
    const r = svg("rect", { x: l + i * bw + 1, y: H - b - hh, width: bw - 2, height: hh, fill: bin >= -30 && bin <= 10 ? "var(--bad)" : "var(--muted)" }, s);
    svg("title", {}, r).textContent = `${bin} to ${bin + 9} bp: ${n} junctions`;
    if (i % 2 === 0) svg("text", { x: l + i * bw + bw / 2, y: H - b + 14, "text-anchor": "middle", fill: "var(--muted)" }, s).textContent = bin === 150 ? "≥150" : bin === -50 ? "≤-50" : bin;
  });
  svg("text", { x: l + (W - l) / 2, y: H - 2, "text-anchor": "middle", fill: "var(--muted)" }, s).textContent = "bases between the two aligned parts";
  svg("text", { x: 4, y: t + 8, fill: "var(--muted)" }, s).textContent = max.toLocaleString();
  $("#gapChart").replaceChildren(s);
  const st = D.stats, ca = st.chimera_with_adapter;
  const nn = D.numbers;
  $("#gapNote").innerHTML = `${Math.round(100 * st.chimera_direct_frac)}% of chimera junctions join the two molecules directly (less than 20 bases between the aligned parts), leaving no room for an adapter. Only ${nn.chimeric_with_internal_adapter.toLocaleString()} of the ${nn.chimeric.reads.toLocaleString()} chimeric reads (${pct(100 * nn.chimeric_with_internal_adapter / nn.chimeric.reads)}) contain an adapter, against ${pct(100 * nn.foldback_with_internal_adapter / nn.classes.foldback.reads)} of foldbacks. The shorter part is often small (median ${fmtLen(st.chimera_shortest_seg_median)}; ${Math.round(100 * st.chimera_shortest_lt1k_frac)}% under 1 kb).`;
}

function foldChart() {
  const pts = D.stats.foldback_ident;
  const W = 460, H = 220, l = 42, b = 34, lo = 85, hi = 100;
  const sc = (v) => (Math.max(lo, Math.min(hi, v)) - lo) / (hi - lo);
  const s = svg("svg", { viewBox: `0 0 ${W} ${H}` });
  const X = (v) => l + sc(v) * (W - l - 10), Y = (v) => H - b - sc(v) * (H - b - 10);
  svg("line", { x1: X(lo), y1: Y(lo), x2: X(hi), y2: Y(hi), stroke: "var(--muted)", "stroke-dasharray": "4 3" }, s);
  for (const [a, c] of pts) svg("circle", { cx: X(a), cy: Y(c), r: 2.2, fill: "var(--bad)", "fill-opacity": .45 }, s);
  for (const v of [85, 90, 95, 100]) {
    svg("text", { x: X(v), y: H - b + 14, "text-anchor": "middle", fill: "var(--muted)" }, s).textContent = v;
    svg("text", { x: l - 6, y: Y(v) + 4, "text-anchor": "end", fill: "var(--muted)" }, s).textContent = v;
  }
  svg("text", { x: l + (W - l) / 2, y: H - 4, "text-anchor": "middle", fill: "var(--muted)" }, s).textContent = "identity of first part (%)";
  svg("text", { x: 10, y: 12, fill: "var(--muted)" }, s).textContent = "second part (%)";
  $("#foldChart").replaceChildren(s);
  $("#foldNote").textContent = `In ${Math.round(100 * D.stats.foldback_second_worse)}% of ${D.stats.foldback_n.toLocaleString()} two-part foldbacks (${pts.length.toLocaleString()} shown) the second (complement) half is less accurate than the first. Points below the dashed line.`;
}

function filterTable() {
  const rows = CATS.filter((c) => D.filters[c.key]).map((c) => [c.short, D.filters[c.key]]);
  rows.push(["Reads without artefacts (includes many short reads)", D.filters.clean]);
  $("#filterTable").innerHTML = `<thead><tr><th>Read class</th><th class="num">Median read Q</th><th class="num">Q≥10 &amp; ≥1 kb</th><th class="num">Q≥15 &amp; ≥1 kb</th><th class="num">Q≥20 &amp; ≥5 kb</th></tr></thead><tbody>` +
    rows.map(([n, f], i) => `<tr><td>${esc(n)}</td><td class="num">${f.median_qs}</td><td class="num">${f.q10_l1000}%</td><td class="num">${f.q15_l1000}%</td><td class="num">${f.q20_l5000}%</td></tr>`).join("") + `</tbody>`;
}

function scanSection() {
  const sc = D.scan, saw = D.saw;
  const starts = sc.scan_start_h.filter((h, i) => sc.scan_dur_s[i] > 0);
  const W = 900, H = 260, l = 46, b = 34, t = 14, r = 12;
  const maxH = Math.ceil(D.run_hours);
  const [hz0, hz1] = [zoom.h0, zoom.h1 ?? maxH];
  const pts = saw.filter((x) => x.h >= hz0 && x.h <= hz1);
  const X = (h) => l + (W - l - r) * (h - hz0) / (hz1 - hz0);
  const dmin = Math.min(...pts.map((x) => x.dels)) * 0.95, dmax = Math.max(...pts.map((x) => x.dels)) * 1.03;
  const Y = (v) => t + (H - b - t) * 0.55 * (1 - (v - dmin) / (dmax - dmin));
  const lmax = Math.max(...pts.map((x) => x.lowc));
  const Yl = (v) => H - b - (H - b - t) * 0.38 * v / lmax;
  const s = svg("svg", { viewBox: `0 0 ${W} ${H}` });
  for (const h0 of starts) if (h0 >= hz0 && h0 <= hz1) svg("line", { x1: X(h0), x2: X(h0), y1: t, y2: H - b, stroke: "var(--muted)", "stroke-dasharray": "2 3", "stroke-opacity": .7 }, s);
  const bw = Math.max(1, (W - l - r) / (pts.length || 1) - 0.5);
  for (const x of pts) svg("rect", { x: X(x.h), y: Yl(x.lowc), width: bw, height: H - b - Yl(x.lowc), fill: "var(--junk)" }, s);
  svg("polyline", { points: pts.map((x) => `${X(x.h + 1 / 24)},${Y(x.dels)}`).join(" "), fill: "none", stroke: "var(--bad)", "stroke-width": 1.6 }, s);
  const step = (hz1 - hz0) > 12 ? 6 : (hz1 - hz0) > 4 ? 1 : 0.5;
  for (let h = Math.ceil(hz0 / step) * step; h <= hz1; h += step) svg("text", { x: X(h), y: H - b + 14, "text-anchor": "middle", fill: "var(--muted)" }, s).textContent = (+h.toFixed(1)) + " h";
  for (const v of [dmin, (dmin + dmax) / 2, dmax]) svg("text", { x: l - 6, y: Y(v) + 4, "text-anchor": "end", fill: "var(--bad)" }, s).textContent = v.toFixed(2);
  svg("text", { x: l, y: t - 3, fill: "var(--bad)" }, s).textContent = "deletions per kb (median read, 5-min bins)";
  svg("text", { x: l - 6, y: Yl(lmax) + 4, "text-anchor": "end", fill: "var(--junk)" }, s).textContent = (100 * lmax).toFixed(0) + "%";
  svg("text", { x: l + 4, y: Yl(lmax) - 4, fill: "var(--junk)" }, s).textContent = "reads with low-complexity junk";
  $("#sawChart").replaceChildren(s);
  $("#sawLegend").innerHTML = `<span>Dashed lines: the ${starts.length} moments when hundreds of reads end at once (MinKNOW's periodic scan).</span>
    <span>Zoom: <button class="chip" data-z="all" aria-pressed="${zoom.h1 == null}">whole run</button> <button class="chip" data-z="6" aria-pressed="${zoom.h1 === 6}">first 6 h</button> <button class="chip" data-z="c" aria-pressed="${zoom.h1 === 13.5}">one cycle</button></span>`;
  const sn = sc.since, first = sn[0], mid = sn.slice(4, 8), during = sn[sn.length - 1];
  const avg = (k) => mid.reduce((a, x) => a + x[k] * x.n, 0) / mid.reduce((a, x) => a + x.n, 0);
  const gap = (starts[starts.length - 2] - starts[1]) / (starts.length - 3);
  const durs = sc.scan_dur_s.filter((x) => x > 100);
  const med = durs.slice().sort((a, b) => a - b)[Math.floor(durs.length / 2)];
  const effect = first.lowc / avg("lowc");
  if (!(effect > 2) || !isFinite(gap)) {
    $("#scanTitle").textContent = "Pore scans: no effect on the reads in this run";
    $("#scanIntro").innerHTML = `MinKNOW still pauses the flow cell from time to time (${starts.length} moments when hundreds of reads ended together, at irregular intervals), but reads that start right after these pauses look like all other reads: ${pct(100 * first.lowc)} carry low-complexity junk against ${pct(100 * avg("lowc"))} otherwise, and the deletion rate does not change. Compare with the <a href="?ds=zymohmw#scans">ZymoHMW run</a>.`;
    $("#scanNote").innerHTML = "";
  }
  const noEffect = !(effect > 2) || !isFinite(gap);
  if (!noEffect) $("#scanTitle").textContent = `Every ${gap.toFixed(1)} hours the whole flow cell changes`;
  if (!noEffect) $("#scanIntro").innerHTML = `Every <b>${gap.toFixed(2)} hours</b>, hundreds of reads across the flow cell stop within the same few seconds. This is MinKNOW's periodic pore/mux scan, which interrupts sequencing in a few waves over about ${Math.round(med / 60 * 10) / 10} minutes (${med} s). Two things follow, on the whole flow cell at once:
    <b>1. Junk at the start of reads.</b> Of the reads (≥1 kb) that start in the first minute after a scan, ${pct(100 * first.lowc)} carry low-complexity junk, ${(first.lowc / avg("lowc")).toFixed(0)} times as many as reads that start 10–90 minutes after a scan (${pct(100 * avg("lowc"))}), usually as the first few hundred bases. They also carry internal adapters ${(first.ia / avg("ia")).toFixed(1)} times as often.
    <b>2. A sawtooth in the error rate.</b> The deletion rate drops right after each scan and climbs again until the next one: ${first.delkb.toFixed(2)} deletions per kb in reads starting in the first minute vs ${sn[7].delkb.toFixed(2)} an hour later. The foldback rate climbs within each cycle too, from ${pct(100 * sn[1].fold)} to ${pct(100 * sn[7].fold)}.`;
  $("#sinceTable").innerHTML = `<thead><tr><th>Read started</th><th class="num">Reads ≥1 kb</th><th class="num">Low-complexity</th><th class="num">Internal adapter</th><th class="num">Foldback</th><th class="num">Inter-species chimera</th><th class="num">Deletions / kb</th><th class="num">Median length</th></tr></thead><tbody>` +
    sn.map((x, i) => {
      const lab = i === sn.length - 1 ? "during a scan" : x.b.replace(/^\((-?\d+), (\d+)\]$/, (m, a, b2) => `${Math.max(0, +a)}–${b2} min after a scan`);
      return `<tr><td>${lab}</td><td class="num">${x.n.toLocaleString()}</td><td class="num ${x.lowc > 2 * avg("lowc") ? "hi" : ""}">${pct(100 * x.lowc)}</td><td class="num ${x.ia > 1.8 * avg("ia") ? "hi" : ""}">${pct(100 * x.ia)}</td><td class="num">${pct(100 * x.fold)}</td><td class="num">${pct(100 * x.inter)}</td><td class="num">${x.delkb.toFixed(2)}</td><td class="num">${fmtLen(x.medlen)}</td></tr>`;
    }).join("") + `</tbody>`;
  const jt = sc.junk_time;
  if (!noEffect) $("#scanNote").innerHTML = `${pct(100 * jt.frac_in_window)} of all low-complexity stretches ≥100 bp sit in reads passing the pore during or just after a scan, which is only ${pct(100 * jt.window_share)} of run time. These junk stretches sit at the very start of the read (median position ${pct(100 * jt.pos_near)} into the read), whereas junk elsewhere in the run sits towards the end (median ${pct(100 * jt.pos_other)}). Their base qualities are very low (Q≈2–3), so local quality trimming removes them, but mean-quality filters on long reads do not. See the <a href="#gallery" data-goto="post_scan_junk">post-scan examples</a>.`;
}
let zoom = { h0: 0, h1: null };

function timeChart() {
  const T = D.timeline;
  const W = 900, H = 240, l = 46, b = 34, t = 14, r = 120;
  const maxH = Math.ceil(D.run_hours);
  const series = [["dels", "deletions /kb", "var(--bad)", 1], ["mm", "mismatches /kb", "#4e79a7", 1], ["ins", "insertions /kb", "#59a14f", 1]];
  const maxV = Math.max(...T.flatMap((x) => series.map(([k, , , f]) => x[k] * f))) * 1.1;
  const X = (h) => l + (W - l - r) * h / maxH, Y = (v) => H - b - (H - b - t) * v / maxV;
  const s = svg("svg", { viewBox: `0 0 ${W} ${H}` });
  for (let v = 0; v <= maxV; v += 0.5) {
    svg("line", { x1: l, x2: W - r, y1: Y(v), y2: Y(v), stroke: "var(--grid)" }, s);
    svg("text", { x: l - 6, y: Y(v) + 4, "text-anchor": "end", fill: "var(--muted)" }, s).textContent = v;
  }
  for (let h = 0; h <= maxH; h += 6) svg("text", { x: X(h), y: H - b + 14, "text-anchor": "middle", fill: "var(--muted)" }, s).textContent = h + " h";
  for (const [k, label, col, f] of series) {
    svg("polyline", { points: T.map((x) => `${X(x.h)},${Y(x[k] * f)}`).join(" "), fill: "none", stroke: col, "stroke-width": 1.8 }, s);
    const last = T[T.length - 2];
    svg("text", { x: X(last.h) + 6, y: Y(last[k] * f) + 4, fill: col }, s).textContent = label;
  }
  $("#timeChart").replaceChildren(s);
  const m = (a, k) => a.reduce((x, y) => x + y[k], 0) / a.length;
  const first = T.slice(0, 6), lastT = T.slice(-14, -2);
  const up = (k) => Math.round(100 * (m(lastT, k) / m(first, k) - 1));
  const dUp = up("dels"), mUp = up("mm");
  const trend = Math.max(dUp, mUp) >= 10
    ? `error rates drift upwards: deletions by about ${dUp}% and mismatches by about ${mUp}%`
    : `error rates stay roughly flat (deletions ${dUp >= 0 ? "+" : ""}${dUp}%, mismatches ${mUp >= 0 ? "+" : ""}${mUp}%)`;
  $("#timeIntro").innerHTML = `Over the ${Math.round(D.run_hours)}-hour run (15-minute bins, median of reads ≥1 kb), ${trend} between the first 1.5 hours and the last 3. Chimeras and foldbacks together go from ${(100 * m(T.slice(0, 2), "chim")).toFixed(1)}% of reads ≥1 kb at the start to about ${(100 * m(lastT, "chim")).toFixed(1)}% at the end.`;
}

function channelSection() {
  const c = D.channels;
  const h = c.hist_lowc, e = c.hist_lowc_edges;
  const W = 900, H = 200, l = 46, b = 34, t = 10;
  const max = Math.max(...h);
  const bw = (W - l - 10) / h.length;
  const s = svg("svg", { viewBox: `0 0 ${W} ${H}` });
  h.forEach((n, i) => {
    const hh = n ? Math.max(2, (H - b - t) * Math.sqrt(n / max)) : 0;
    const rr = svg("rect", { x: l + i * bw + 1, y: H - b - hh, width: bw - 2, height: hh, fill: e[i] > 3 ? "var(--bad)" : "var(--junk)" }, s);
    svg("title", {}, rr).textContent = `${e[i].toFixed(1)}–${e[i + 1].toFixed(1)}%: ${n} channels`;
    if (i % 3 === 0) svg("text", { x: l + i * bw, y: H - b + 14, fill: "var(--muted)" }, s).textContent = e[i].toFixed(1) + "%";
  });
  svg("text", { x: l + (W - l) / 2, y: H - 4, "text-anchor": "middle", fill: "var(--muted)" }, s).textContent = "share of a channel's reads with low-complexity junk (bar height: √channels)";
  $("#chanChart").replaceChildren(s);
  $("#chanIntro").innerHTML = `Across ${c.n.toLocaleString()} channels with ≥200 reads, the share of reads with low-complexity junk varies ${c.lowc_var_ratio.toFixed(0)} times more between channels than chance would allow${c.chim_var_ratio >= 2 ? `, and the inter-species chimera share ${c.chim_var_ratio.toFixed(0)} times more` : `, while the inter-species chimera share varies no more than chance`}. The worst pores produce junk in up to ${pct(100 * Math.max(...c.top_lowc.map((x) => x.lowc)))} of their reads, against ${pct(D.summary.low_complexity.pct_reads)} overall. Channel number is in every read's <code>ch</code> tag, so these pores can be identified and their reads flagged.`;
  $("#chanTable").innerHTML = `<thead><tr><th>Channel</th><th class="num">Reads</th><th class="num">Low-complexity</th><th class="num">Inter-species chimeras</th><th class="num">Skip-prone</th><th class="num">Median read Q</th></tr></thead><tbody>` +
    c.top_lowc.slice(0, 8).map((x) => `<tr><td>${x.ch}</td><td class="num">${x.n.toLocaleString()}</td><td class="num hi">${pct(100 * x.lowc)}</td><td class="num">${pct(100 * x.chim)}</td><td class="num">${pct(100 * x.skip)}</td><td class="num">${x.qs}</td></tr>`).join("") + `</tbody>`;
}

function toolTable() {
  const Y = '<span class="yes">yes</span>', N = '<span class="no">no</span>', P = '<span class="part">partly</span>';
  const cols = ["End adapters", "Internal adapters", "Chimeras without adapter", "Foldbacks", "Low-complexity junk", "Low-quality reads"];
  const rows = [
    ["dorado (default trimming + read splitting)", [Y, P, N, P, N, N], "Trims adapters/primers at read ends and splits reads when it recognises an adapter inside. In this run it had already split " + pct(100 * D.stats.dorado_split_reads / D.stats.n_reads) + " of reads, and the artefacts above are what remained."],
    ["Porechop / Porechop_ABI", [Y, Y, N, P, N, N], "Searches read ends and middles for adapters (Porechop_ABI also infers unknown ones). Splits or discards reads with middle adapters. Finds nothing at direct joins."],
    ["Filtlong, chopper, NanoFilt, fastplong", [N, N, N, N, P, Y], "Length and mean-quality filters. They score the whole read, so long, high-quality chimeras pass (see the filter table). fastplong also trims adapters."],
    ["Assemblers (Flye, metaMDBG, myloasm…)", [P, P, P, P, P, P], "Consensus over many reads usually outvotes artefacts at high coverage. Rare genomes in metagenomes, plasmids and repeats are the risk zones."],
    ["Reference/overlap-based chimera detection", [N, N, Y, Y, N, N], "Splitting reads where alignments to a reference or to other reads break (e.g. yacrd, or a mock-community reference as here). Needs a reference or enough coverage."],
  ];
  $("#toolTable").innerHTML = `<thead><tr><th>Tool type</th>${cols.map((c) => `<th>${c}</th>`).join("")}<th>Notes</th></tr></thead><tbody>` +
    rows.map(([n, v, note]) => `<tr><td><b>${esc(n)}</b></td>${v.map((x) => `<td>${x}</td>`).join("")}<td>${esc(note)}</td></tr>`).join("") + `</tbody>`;
}

function methods() {
  const ds = D.dataset, st = D.stats;
  const items = [
    `Reads: all ${st.n_reads.toLocaleString()} reads in <code>${esc(ds.file)}</code> (${esc(ds.accession)}, MD5 checked against ENA), ${(st.n_bases / 1e9).toFixed(1)} Gb from ${st.n_channels.toLocaleString()} channels over ${Math.round(D.run_hours)} hours.`,
    `Basecalling: ${esc(ds.basecaller)}. Reads are analysed as published, i.e. after dorado's own adapter trimming and read splitting (${pct(100 * st.dorado_split_reads / st.n_reads)} of reads are products of a split).`,
    `Reference: ${esc(ds.reference)}. Alignment with minimap2 (mappy, map-ont preset). Alignments shorter than 300 bp or overlapping a longer one by >20% are ignored.`,
    `Adapter search: edlib, ≤20% edit distance, both orientations of the SQK-LSK114 ligation adapter.`,
    `Chimeras: adjacent alignments to different species (inter-species), or to the same contig more than 2 kb away from where a continuous molecule would be, allowing for circular chromosomes (intra-genome). Opposite strands at the same place = foldback. Joins between contigs of the draft yeast assembly are not counted.`,
    `Low-complexity junk: runs of ≥6 GGAAT units, ≥25 identical bases, or ≥12 dinucleotide units, in reads that also align to the mock.`,
    `Scan cycle: seconds in which significantly more reads ended than the number of reads in progress predicts (Poisson, Bonferroni-corrected, ≥5 reads and ≥3× expected), grouped into events. Reads are binned by start time relative to the end of the last event.`,
    `Deletion timing: deletions ≥10 bp inside alignments, placed in time from the read start, trimmed samples and the sample position (ns/ts/du tags). Deletions found in ≥30% of reads at the same genome position are real differences from the reference and were excluded.`,
    `Signal gaps: deletions ≥10 bp inside alignments and collinear split alignments, counted once per read per second and compared with the number of reads in progress that second (Poisson, Bonferroni-corrected, ≥10 reads); neighbouring seconds within 3 s are merged into one gap.`,
    `Molecules read twice: consecutive reads in the same channel on opposite strands that overlap by ≥80% on the genome and start within 20 s of each other.`,
    `Caveats: the Zymo reference genomes are not from the same lot as this DNA. Library: ${esc(D.dataset.kit || "not recorded")}. Unmapped reads (${pct(D.summary.unmapped.pct_reads)} of reads, ${pct(D.summary.unmapped.pct_bases)} of bases) are mostly <500 bp. The cause of the scan-cycle effects cannot be read from the BAM; MinKNOW logs or POD5 run information (bias voltage, mux changes) would show it.`,
  ];
  $("#methodsList").innerHTML = items.map((x) => `<li>${x}</li>`).join("");
}

function bind() {
  $("#tabs").addEventListener("click", (e) => { const b = e.target.closest(".chip"); if (b) select(b.dataset.k); });
  $("#prev").addEventListener("click", () => { idx--; renderExample(); });
  $("#next").addEventListener("click", () => { idx++; renderExample(); });
  document.querySelectorAll("[data-unit]").forEach((b) => b.addEventListener("click", () => {
    unit = b.dataset.unit;
    document.querySelectorAll("[data-unit]").forEach((x) => x.setAttribute("aria-pressed", x === b));
    overview();
  }));
  $("#scans").addEventListener("click", (e) => {
    const z = e.target.closest("[data-z]");
    if (z) {
      const s0 = D.scan.scan_start_h;
      zoom = z.dataset.z === "all" ? { h0: 0, h1: null } : z.dataset.z === "6" ? { h0: 0, h1: 6 } : { h0: s0[7] - 0.25, h1: 13.5 };
      if (z.dataset.z === "c") zoom.h1 = 13.5, zoom.h0 = 10.5;
      scanSection();
    }
    const g = e.target.closest("[data-goto]");
    if (g) { e.preventDefault(); select(g.dataset.goto); $("#gallery").scrollIntoView({ behavior: "smooth" }); }
  });
  let rt; window.addEventListener("resize", () => { clearTimeout(rt); rt = setTimeout(() => renderExample(), 200); });
  document.addEventListener("keydown", (e) => {
    if (e.target.closest("input,select,textarea")) return;
    if (e.key === "ArrowRight") { idx++; renderExample(); }
    if (e.key === "ArrowLeft") { idx--; renderExample(); }
  });
}

let INDEX = [], DSID = "zymohmw", bound = false;

function gapSection() {
  const sec = $("#gaps"), g = D.gaps;
  if (!g) { sec.hidden = true; return; }
  sec.hidden = false;
  const pm = g.per_minute.filter((x) => x[2] > 20);
  const W = 900, H = 220, l = 46, b = 34, t = 16, r = 10;
  const maxH = Math.ceil(g.run_hours);
  const X = (h) => l + (W - l - r) * h / maxH;
  const vals = pm.map((x) => 1000 * x[1] / (x[2] * 60));
  const vmax = Math.max(...vals) * 1.05 || 1;
  const Y = (v) => H - b - (H - b - t) * v / vmax;
  const s = svg("svg", { viewBox: `0 0 ${W} ${H}` });
  for (let h = 0; h <= maxH; h += 6) {
    svg("line", { x1: X(h), x2: X(h), y1: t, y2: H - b, stroke: "var(--grid)" }, s);
    svg("text", { x: X(h), y: H - b + 14, "text-anchor": "middle", fill: "var(--muted)" }, s).textContent = h + " h";
  }
  svg("polyline", { points: pm.map((x, i) => `${X(x[0]).toFixed(1)},${Y(vals[i]).toFixed(1)}`).join(" "), fill: "none", stroke: "var(--q)", "stroke-width": 1 }, s);
  for (const e of g.events) {
    const tick = svg("line", { x1: X(e.t_h), x2: X(e.t_h), y1: t - 10, y2: t - 2, stroke: "var(--bad)", "stroke-width": 1.5 }, s);
    svg("title", {}, tick).textContent = `${e.t_utc}: ${e.reads} reads on ${e.channels} channels`;
  }
  svg("text", { x: l, y: H - 2, fill: "var(--muted)" }, s).textContent = "reads losing ≥10 bases per 1,000 reads in progress, per second (1-min means) · red ticks: detected gaps";
  svg("text", { x: l - 6, y: Y(vmax / 1.05) + 4, "text-anchor": "end", fill: "var(--muted)" }, s).textContent = (vmax / 1.05).toFixed(1);
  svg("text", { x: l - 6, y: Y(0) + 4, "text-anchor": "end", fill: "var(--muted)" }, s).textContent = "0";
  $("#gapTimeChart").replaceChildren(s);
  const dev = D.dataset.device ? ` on <b>${esc(D.dataset.device)}</b>` : "";
  $("#gapIntro").innerHTML = g.n_events
    ? `If the flow cell stops delivering signal for a moment, every read being sequenced loses the stretch of its molecule that passed through the pore meanwhile (about 400 bases per second). In the reads this shows up as a deletion, or as a read whose alignment is split in two, in hundreds of channels at the same second. This run${dev} has <b>${g.n_events.toLocaleString()} such gaps</b>, hitting <b>${g.reads_hit.toLocaleString()} reads (${pct(g.pct_reads_hit)})</b>. The missing stretch is too short to break a read and too long to be a normal error, so neither dorado nor the usual filters notice it.`
    : `If the flow cell stops delivering signal for a moment, every read being sequenced loses the stretch of its molecule that passed through the pore meanwhile (about 400 bases per second). Every second of this run${dev} was tested (on average ${Math.round(g.mean_active).toLocaleString()} reads in progress). <b>No signal gaps were found</b>: no second has many reads on many channels losing sequence together.`;
  $("#gapRuns").innerHTML = `<thead><tr><th>Run</th><th>Device</th><th class="num">Reads</th><th class="num">Gaps</th><th class="num">Reads hit</th></tr></thead><tbody>` +
    INDEX.map((x) => `<tr><td>${x.id === DSID ? "<b>" + esc(x.label) + "</b>" : `<a href="?ds=${x.id}#gaps">${esc(x.label)}</a>`}</td><td>${esc(x.device || "")}</td><td class="num">${x.n_reads ? x.n_reads.toLocaleString() : "–"}</td><td class="num ${x.gaps ? "hi" : ""}">${x.gaps ?? "–"}</td><td class="num">${x.gaps_pct != null ? pct(x.gaps_pct) : "–"}</td></tr>`).join("") + `</tbody>`;
  const top = g.events.slice().sort((a, b) => b.reads - a.reads).slice(0, 12);
  $("#gapTable").innerHTML = top.length ? `<thead><tr><th>Time (UTC)</th><th class="num">Run hour</th><th class="num">Reads hit</th><th class="num">Channels</th><th class="num">Share of reads in progress</th><th class="num">Bases missing (median)</th><th class="num">Seen as split alignment</th></tr></thead><tbody>` +
    top.map((e) => `<tr><td>${e.t_utc}</td><td class="num">${e.t_h.toFixed(2)}</td><td class="num">${e.reads}</td><td class="num">${e.channels}</td><td class="num">${pct(100 * e.share_of_active)}</td><td class="num">${e.missing_median}</td><td class="num">${pct(100 * e.split_share)}</td></tr>`).join("") + `</tbody>` : "";
  $("#gapTable").parentElement.previousElementSibling.hidden = !top.length;
  $("#gapNote2").innerHTML = g.ontrim ? `Independent check: ontrim's report for the same run lists ${g.ontrim.events} gaps; ${g.ontrim.matched_by_mine} of them fall within 5 s of a gap found here, and ${g.ontrim.mine_matched_by_ontrim} of the ${g.n_events} gaps found here are in ontrim's list.` : "";
}

function duplexSection() {
  const d = D.duplex, sec = $("#duplex");
  if (!d) { sec.hidden = true; return; }
  sec.hidden = false;
  $("#duplexIntro").innerHTML = `${d.pairs.toLocaleString()} times, the next read in the same pore is the other strand of the same molecule (opposite strand, ≥80% overlap on the genome, starting within 20 s). ${pct(d.pct_reads_in_pairs)} of reads belong to such a pair. The second read of each pair (${pct(d.pct_reads_second)} of all reads, ${pct(d.bases_in_second_reads_pct)} of all bases) is a molecule that is already in the data. ${d.from_one_split_read.toLocaleString()} of the pairs were one read that dorado split in two. ${D.dataset.flowcell_note ? D.dataset.flowcell_note + " " : ""}The simplex BAM keeps both reads, so these molecules count twice in coverage, abundance and allele frequencies. Duplex basecalling, or dropping the second read of each pair, avoids that.`;
}

function picker() {
  $("#dsPicker").innerHTML = INDEX.map((x) => `<button class="chip" data-ds="${x.id}" aria-pressed="${x.id === DSID}">${esc(x.label)}</button>`).join("");
}

function renderAll() {
  headline(); overview(); scanSection(); gapSection(); tabs(); renderCat(); gapChart(); foldChart(); filterTable(); timeChart(); channelSection(); duplexSection(); toolTable(); methods();
  if (!bound) { bind(); bound = true; }
}

function load(id) {
  const meta = INDEX.find((x) => x.id === id) || INDEX[0];
  DSID = meta.id; picker();
  return fetch("data/" + meta.file).then((r) => r.json()).then((d) => {
    D = d; idx = 0; zoom = { h0: 0, h1: null };
    if (!D.examples[cur]?.length) cur = "chimera_interspecies";
    renderAll();
  });
}

fetch("data/artefacts_index.json").then((r) => r.json()).then((ix) => {
  INDEX = ix;
  const p = new URLSearchParams(location.search);
  const hash = location.hash.slice(1);
  if (CAT[hash]) cur = hash;
  $("#dsPicker").addEventListener("click", (e) => {
    const b = e.target.closest("[data-ds]"); if (!b) return;
    history.replaceState(null, "", "?ds=" + b.dataset.ds + location.hash);
    load(b.dataset.ds);
  });
  return load(p.get("ds") || ix[0].id).then(() => { if (location.hash) document.getElementById(location.hash.slice(1))?.scrollIntoView(); });
}).catch((e) => { $("#overview").insertAdjacentHTML("beforeend", `<p class="empty">Could not load data (${esc(e.message)}).</p>`); });
