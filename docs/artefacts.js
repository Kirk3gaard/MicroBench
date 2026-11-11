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
    tools: "In this run about 4 in 10 foldbacks contain the ligation adapter at the fold, so adapter splitting can catch those. The rest join without an adapter." },
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
    matters: "(GGAAT)n is the motif of human satellite III, so these reads can trigger human-contamination screens. Junk insertions also end up in low-coverage assemblies and inflate error rates. In this run the repeat appears almost only as GGAAT and practically never as its reverse complement ATTCC, which points to the basecaller rather than real DNA.",
    tools: "Mean-quality filters average the junk away inside a long read. Only tools that look at local quality or sequence complexity can remove it." },
  { key: "skip_prone", title: "Skip-prone reads", short: "Skip-prone reads",
    what: "Reads with five or more deletions of ≥10 bp: stretches of the molecule that never appear in the basecalls, as if the strand slipped through the pore.",
    detect: "≥5 deletions of ≥10 bp that are not shared by other reads (shared ones are real differences from the reference).",
    matters: "A source of false indels and frameshifts in low-coverage regions.",
    tools: "Their mean Q-score is lower but most still pass a Q10 filter." },
  { key: "unmapped", title: "Unmapped reads", short: "Unmapped",
    what: "Reads with no alignment of ≥300 bp. Almost all are short.",
    detect: "No alignment to the mock references.",
    matters: "Many reads but very few bases. Mostly harmless.",
    tools: "Any length filter removes them." },
];
const CAT = Object.fromEntries(CATS.map((c) => [c.key, c]));

let D, cur = "chimera_interspecies", idx = 0, unit = "pct_bases";

function svg(tag, attrs = {}, parent) {
  const el = document.createElementNS(SVGNS, tag);
  for (const [k, v] of Object.entries(attrs)) el.setAttribute(k, v);
  if (parent) parent.appendChild(el);
  return el;
}

function headline() {
  const s = D.summary, st = D.stats;
  const chimReads = ["chimera_interspecies", "chimera_intraspecies", "foldback"].reduce((a, k) => a + s[k].reads, 0);
  const tiles = [
    [pct(s._any.pct_bases), "of bases are in reads with at least one artefact"],
    [pct(s._any.pct_reads), "of reads"],
    [pct(100 * chimReads / st.n_reads), "of reads are chimeric or foldback"],
    [pct(100 * st.dorado_split_reads / st.n_reads), "of reads were already split by dorado"],
    [st.n_reads.toLocaleString(), "reads analysed"],
  ];
  $("#headline").innerHTML = tiles.map(([v, l]) => `<div class="stat"><b>${v}</b><span>${l}</span></div>`).join("");
  const ds = D.dataset;
  $("#datasetLine").innerHTML = `${esc(ds.sample)} · flow cell ${esc(ds.flowcell)} · ${esc(ds.basecaller)} · ${st.n_reads.toLocaleString()} reads (${(st.n_bases / 1e6).toFixed(0)} Mb) sampled from ${st.n_channels} channels of <a href="https://www.ebi.ac.uk/ena/browser/view/${esc(ds.accession)}">${esc(ds.accession)}</a>.`;
}

function overview() {
  const host = $("#overviewChart");
  const rows = CATS.map((c) => ({ c, v: D.summary[c.key][unit], n: D.summary[c.key].reads }));
  const W = 760, rowH = 26, left = 190, H = rows.length * rowH + 24;
  const max = Math.max(...rows.map((r) => r.v)) * 1.15;
  const s = svg("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": "Artefact prevalence" });
  rows.forEach((r, i) => {
    const g = svg("g", { class: "bar-row", transform: `translate(0 ${i * rowH + 4})` }, s);
    svg("text", { x: left - 8, y: 15, "text-anchor": "end" }, g).textContent = r.c.short;
    const w = Math.max(2, (W - left - 90) * r.v / max);
    svg("rect", { class: "bar", x: left, y: 3, width: w, height: 16, rx: 3, fill: r.c.key === "unmapped" ? "var(--muted)" : "var(--bad)" }, g);
    svg("text", { x: left + w + 6, y: 15 }, g).textContent = `${pct(r.v)}  (${r.n.toLocaleString()} reads)`;
    g.addEventListener("click", () => { select(r.c.key); $("#gallery").scrollIntoView({ behavior: "smooth" }); });
  });
  host.replaceChildren(s);
}

function tabs() {
  $("#tabs").innerHTML = CATS.filter((c) => D.examples[c.key]?.length).map((c) =>
    `<button class="chip" role="tab" data-k="${c.key}" aria-pressed="${c.key === cur}">${esc(c.short)}<small>${pct(D.summary[c.key].pct_reads)}</small></button>`).join("");
}

function select(k) { cur = k; idx = 0; tabs(); renderCat(); }

function renderCat() {
  const c = CAT[cur], s = D.summary[cur], f = D.filters?.[cur];
  $("#catInfo").innerHTML = `<h3>${esc(c.title)}</h3>
    <div><span class="big">${pct(s.pct_reads)}</span> of reads · <b>${pct(s.pct_bases)}</b> of bases · ${s.reads.toLocaleString()} reads</div>
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
  $("#gapNote").innerHTML = `${Math.round(100 * st.chimera_direct_frac)}% of chimera junctions join the two molecules directly (within ±15 bp), with no room for an adapter. Only ${ca.chimera_interspecies + ca.chimera_intraspecies} of ${(D.summary.chimera_interspecies.reads + D.summary.chimera_intraspecies.reads).toLocaleString()} chimeric reads contain an adapter. The shorter part is often small (median ${fmtLen(st.chimera_shortest_seg_median)}; ${Math.round(100 * st.chimera_shortest_lt1k_frac)}% under 1 kb).`;
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
  const worse = pts.filter(([a, c]) => c < a).length;
  $("#foldNote").textContent = `In ${Math.round(100 * worse / pts.length)}% of ${pts.length.toLocaleString()} foldbacks the second (complement) half is less accurate than the first. Points below the dashed line.`;
}

function filterTable() {
  const rows = CATS.filter((c) => D.filters[c.key]).map((c) => [c.short, D.filters[c.key]]);
  rows.push(["Clean reads (reference)", D.filters.clean]);
  $("#filterTable").innerHTML = `<thead><tr><th>Read class</th><th class="num">Median read Q</th><th class="num">Q≥10 &amp; ≥1 kb</th><th class="num">Q≥15 &amp; ≥1 kb</th><th class="num">Q≥20 &amp; ≥5 kb</th></tr></thead><tbody>` +
    rows.map(([n, f], i) => `<tr><td>${esc(n)}</td><td class="num">${f.median_qs}</td><td class="num ${i < 3 && f.q10_l1000 > D.filters.clean.q10_l1000 ? "hi" : ""}">${f.q10_l1000}%</td><td class="num">${f.q15_l1000}%</td><td class="num">${f.q20_l5000}%</td></tr>`).join("") + `</tbody>`;
}

function timeChart() {
  const T = D.timeline;
  const W = 900, H = 240, l = 46, b = 34, t = 12, r = 120;
  const maxH = Math.max(...T.map((x) => x.h)) + 1;
  const series = [["del", "deletions", "var(--bad)"], ["mm", "mismatches", "#4e79a7"], ["ins", "insertions", "#59a14f"]];
  const maxV = Math.max(...T.flatMap((x) => series.map(([k]) => x[k]))) * 1.1;
  const X = (h) => l + (W - l - r) * (h + 0.5) / maxH, Y = (v) => H - b - (H - b - t) * v / maxV;
  const s = svg("svg", { viewBox: `0 0 ${W} ${H}` });
  for (let v = 0; v <= maxV; v += 0.5) {
    svg("line", { x1: l, x2: W - r, y1: Y(v), y2: Y(v), stroke: "var(--grid)" }, s);
    svg("text", { x: l - 6, y: Y(v) + 4, "text-anchor": "end", fill: "var(--muted)" }, s).textContent = v.toFixed(1);
  }
  for (let h = 0; h <= maxH; h += 6) svg("text", { x: X(h) , y: H - b + 14, "text-anchor": "middle", fill: "var(--muted)" }, s).textContent = h + " h";
  svg("text", { x: l, y: t - 2, fill: "var(--muted)" }, s).textContent = "errors per kb (median read)";
  for (const [k, label, col] of series) {
    svg("polyline", { points: T.map((x) => `${X(x.h)},${Y(x[k])}`).join(" "), fill: "none", stroke: col, "stroke-width": 2 }, s);
    const last = T[T.length - 1];
    svg("text", { x: X(last.h) + 6, y: Y(last[k]) + 4, fill: col }, s).textContent = label;
  }
  $("#timeChart").replaceChildren(s);
  const first = T.slice(0, 3), lastT = T.slice(-6, -1);
  const m = (a, k) => a.reduce((x, y) => x + y[k], 0) / a.length;
  const up = (k) => Math.round(100 * (m(lastT, k) / m(first, k) - 1));
  $("#timeIntro").innerHTML = `Median per-read error rates for reads starting in each hour of the ${T.length}-hour run (reads ≥1 kb). They creep up across the flow cell: deletions rise by about ${up("del")}% and mismatches by about ${up("mm")}% between the first hours and the end of the run. The chimera rate fluctuates (${Math.min(...T.map((x) => x.chim)).toFixed(0)}–${Math.max(...T.map((x) => x.chim)).toFixed(0)}% of reads ≥1 kb per hour), with no trend, which fits chimeras forming during library preparation rather than during sequencing.`;
}

function raster() {
  const ev = D.raster;
  const chans = [...new Set(ev.map((e) => e[1]))].sort((a, b) => a - b);
  const W = 900, H = 260, l = 46, b = 30, t = 10;
  const maxH = Math.max(...D.timeline.map((x) => x.h)) + 1;
  const X = (h) => l + (W - l - 10) * h / maxH;
  const cidx = Object.fromEntries(chans.map((c, i) => [c, i]));
  const Y = (c) => t + (H - b - t) * cidx[c] / Math.max(1, chans.length - 1);
  const s = svg("svg", { viewBox: `0 0 ${W} ${H}` });
  for (let h = 0; h <= maxH; h += 6) {
    svg("line", { x1: X(h), x2: X(h), y1: t, y2: H - b, stroke: "var(--grid)" }, s);
    svg("text", { x: X(h), y: H - b + 14, "text-anchor": "middle", fill: "var(--muted)" }, s).textContent = h + " h";
  }
  for (const [h, c, ln] of ev) {
    const p = svg("circle", { cx: X(h), cy: Y(c), r: Math.min(5, 1.5 + Math.log10(ln)), fill: "var(--bad)", "fill-opacity": .55 }, s);
    svg("title", {}, p).textContent = `${ln} bp deletion · channel ${c} · ${h.toFixed(2)} h`;
  }
  svg("text", { x: 12, y: H / 2, fill: "var(--muted)", transform: `rotate(-90 12 ${H / 2})`, "text-anchor": "middle" }, s).textContent = `${chans.length} channels`;
  $("#raster").replaceChildren(s);
  const sy = D.sync;
  const conc = sy[0]?.mean_concurrent_reads;
  $("#syncIntro").innerHTML = `Each dot is a deletion of ≥20 bp (not present in other reads, so not a real difference from the reference), placed at the moment it passed through the pore (read start + position × duration) and on its channel's row. If something hit the whole flow cell at once, such as a voltage step or a pause in data acquisition, you would see vertical stripes. On average ${conc} of the sampled reads were being sequenced at any moment. The table compares the number of deletion pairs in different channels within ±2 s of each other with the same reads placed at random times (50 shuffles).`;
  $("#syncTable").innerHTML = `<thead><tr><th>Deletion size</th><th class="num">Events</th><th class="num">Coincident pairs</th><th class="num">Expected by chance (max)</th><th>Verdict</th></tr></thead><tbody>` +
    sy.map((x) => `<tr><td>≥${x.min_len} bp</td><td class="num">${x.events.toLocaleString()}</td><td class="num">${x.observed_pairs}</td><td class="num">${x.null_mean.toFixed(1)} (${x.null_max})</td><td>${x.observed_pairs > x.null_max ? '<span class="no">synchronised</span>' : '<span class="yes">no excess</span>'}</td></tr>`).join("") + `</tbody>`;
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
    `Reads: complete BAM records from ${st.n_reads.toLocaleString()} reads, taken from evenly spaced 1 MB slices of <code>${esc(ds.file)}</code> (${esc(ds.accession)}). Dorado writes the BAM grouped by channel, so the sample covers ${st.n_channels} channels over the whole run.`,
    `Basecalling: ${esc(ds.basecaller)}. Reads are analysed as published, i.e. after dorado's own adapter trimming and read splitting.`,
    `Reference: ${esc(ds.reference)}. Alignment with minimap2 (mappy, map-ont preset). Alignments shorter than 300 bp or overlapping a longer one by >20% are ignored.`,
    `Adapter search: edlib, ≤20% edit distance, both orientations of the SQK-LSK114 ligation adapter.`,
    `Chimeras: adjacent alignments to different species (inter-species), or to the same contig more than 2 kb away from where a continuous molecule would be, allowing for circular chromosomes (intra-genome). Opposite strands at the same place = foldback. Joins between contigs of the draft yeast assembly are not counted.`,
    `Deletion events: deletions ≥10 bp inside alignments that occur in fewer than 3 reads at the same reference position. Their time is estimated as read start + duration × position / length.`,
    `Caveats: the Zymo reference genomes are not from the same lot as this DNA, so a few real differences remain. Chimera rates depend on library prep (the ligation adapter is present, so this was a ligation library; the exact kit is not recorded yet). The run had about ${(st.n_mapped / st.n_reads * 100).toFixed(0)}% mappable reads; unmapped reads are mostly <500 bp.`,
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
  let rt; window.addEventListener("resize", () => { clearTimeout(rt); rt = setTimeout(() => renderExample(), 200); });
  document.addEventListener("keydown", (e) => {
    if (e.target.closest("input,select,textarea")) return;
    if (e.key === "ArrowRight") { idx++; renderExample(); }
    if (e.key === "ArrowLeft") { idx--; renderExample(); }
  });
}

fetch("data/artefacts_zymohmw.json").then((r) => r.json()).then((d) => {
  D = d;
  const hash = location.hash.slice(1);
  if (CAT[hash]) cur = hash;
  headline(); overview(); tabs(); renderCat(); gapChart(); foldChart(); filterTable(); timeChart(); raster(); toolTable(); methods(); bind();
}).catch((e) => { $("#overview").insertAdjacentHTML("beforeend", `<p class="empty">Could not load data (${esc(e.message)}).</p>`); });
