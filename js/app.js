(async function () {
  const $ = id => document.getElementById(id);
  const { el, clear, scale, extent, axes, fmt, showTip, hideTip, group, esc, COLORS, NIGHT, ORDER, SHORT, LINE_COLORS } = PT;
  let R, G, synthetic = false;

  async function get(url) {
    const r = await fetch(url, { cache: "no-cache" });
    if (!r.ok) throw new Error(url);
    return r.json();
  }
  try {
    [R, G] = await Promise.all([get("data/results.json"), get("data/genes.json")]);
  } catch (e) {
    try {
      [R, G] = await Promise.all([get("data/results_synthetic.json"), get("data/genes_synthetic.json")]);
      synthetic = true;
    } catch (e2) {
      document.querySelectorAll(".plate .fig, .stage").forEach(c => c.innerHTML = '<p class="note">The analysis has not been published yet. The GitHub Actions workflow downloads the GEO data and fills this page.</p>');
      return;
    }
  }
  if (synthetic || R.meta.synthetic) {
    $("banner").textContent = "Synthetic test data. Nothing on this page is a real result.";
    $("banner").classList.add("show");
  }

  const passports = R.passports || [];
  const lines = R.meta.cell_lines;
  const EXT = R.external || {};
  const med = a => { if (!a.length) return null; const b = [...a].sort((x, y) => x - y), m = b.length >> 1; return b.length % 2 ? b[m] : (b[m - 1] + b[m]) / 2; };
  const memOf = tp => R.memory.summary.filter(r => r.timepoint === tp && r.median_memory !== null).map(r => r.median_memory);
  const memLong = memOf("wo_long"), memShort = memOf("wo_short");
  const betas = Object.values(R.betas);
  const hyp = id => R.hypotheses.find(h => h.id === id);
  const halfLives = Object.values(R.memory.decay).map(d => d.median_half_life_days).filter(Number.isFinite);
  const pct = v => v === null || v === undefined || !Number.isFinite(v) ? "n/a" : Math.round(v * 100) + "%";
  const slug = s => String(s).replace(/\s+/g, "-");

  function ledger() {
    const h1 = hyp("H1");
    const sc = EXT.pc9_single_cell;
    const items = [
      [R.meta.n_samples, "public RNA-seq samples, four cell lines"],
      [`${Math.round(Math.min(...betas) * 100)}–${Math.round(Math.max(...betas) * 100)}%`, "of the day-one stress still carried by persisters"],
      [h1 ? `${fmt(h1.naive, 2)} → ${fmt(h1.decomposed, 2)}` : "n/a", "AUROC on an unseen cell line, naive → stress removed", "red"],
      [halfLives.length ? `${Math.round(Math.min(...halfLives))}–${Math.round(Math.max(...halfLives))} d` : "n/a", "half-life of the persister memory"],
      [sc ? sc.cells_after_qc.toLocaleString() : "n/a", "single cells re-analysed from GSE275199"],
      [passports.length, "evidence passports, six dimensions each"],
    ];
    $("ledger").innerHTML = items.map(([v, l, c]) => `<div><b class="${c || ""}">${v}</b><span>${l}</span></div>`).join("");
  }

  const PLATES = [...document.querySelectorAll("section.plate")];
  function spark(id) {
    const W = 140, H = 34;
    if (id === "problem") return `<svg viewBox="0 0 ${W} ${H}">${Array.from({ length: 40 }, (_, i) => `<circle cx="${4 + (i % 20) * 7}" cy="${8 + Math.floor(i / 20) * 14}" r="2.4" fill="${i < 11 ? "#d9707b" : i < 34 ? "#8db4cf" : "#6b655b"}"/>`).join("")}</svg>`;
    if (id === "universe") { const xs = R.samples.map(s => s.pc[0]), ys = R.samples.map(s => s.pc[1]); const [a, b] = extent(xs), [c, d] = extent(ys); return `<svg viewBox="0 0 ${W} ${H}">${R.samples.map(s => `<circle cx="${(s.pc[0] - a) / (b - a) * (W - 6) + 3}" cy="${H - 3 - (s.pc[1] - c) / (d - c) * (H - 6)}" r="1.6" fill="${NIGHT[s.condition]}"/>`).join("")}</svg>`; }
    if (id === "reversibility" || id === "paradox") { const g = G.TACSTD2 || Object.values(G)[0]; const t = (g.t[lines[0]] || [0, 0, 0, 0, 0]).map(v => v ?? 0); const [lo, hi] = extent(t.concat([0])); return `<svg viewBox="0 0 ${W} ${H}"><path d="${t.map((v, i) => `${i ? "L" : "M"}${8 + i * 31},${H - 4 - (v - lo) / (hi - lo) * (H - 8)}`).join("")}" fill="none" stroke="#d9707b" stroke-width="2"/></svg>`; }
    if (id === "arena") { const f = R.benchmark.folds; return `<svg viewBox="0 0 ${W} ${H}">${f.map((d, i) => { const a = d["naive_treated_vs_control|auroc_dtp_vs_acute"] || 0, b = d["persister_specific|auroc_dtp_vs_acute"] || 0; return `<rect x="${6 + i * 34}" y="${H - 30 * a}" width="12" height="${30 * a}" fill="#e0ad55"/><rect x="${19 + i * 34}" y="${H - 30 * b}" width="12" height="${30 * b}" fill="#d9707b"/>`; }).join("")}</svg>`; }
    if (id === "decompose") { const pts = (R.decomposition_scatter[lines[0]] || []).slice(0, 160); return `<svg viewBox="0 0 ${W} ${H}">${pts.map(p => `<circle cx="${W / 2 + (p[0] || 0) * 10}" cy="${H / 2 - (p[1] || 0) * 3}" r="1" fill="${Math.abs(p[2]) > 1.5 ? "#d9707b" : "#6b655b"}"/>`).join("")}</svg>`; }
    const pts = Array.from({ length: 24 }, (_, i) => Math.sin(i / 3 + id.length) * .4 + Math.cos(i / 5 + id.length * 2) * .3);
    return `<svg viewBox="0 0 ${W} ${H}"><path d="${pts.map((v, i) => `${i ? "L" : "M"}${i * (W / 23)},${H / 2 - v * 18}`).join("")}" fill="none" stroke="#6b655b" stroke-width="1.4"/></svg>`;
  }
  function plateIndex() {
    const sc = EXT.pc9_single_cell;
    const live = {
      problem: "2.48 M diagnoses, 1.82 M deaths, 2022",
      drug: "PDB 4ZAU in live 3D",
      design: "4 lines × 5 conditions × 3",
      universe: `PC1 holds ${pct(R.meta.pca_variance[0])} of variance`,
      decompose: `β from ${fmt(Math.min(...betas), 2)} to ${fmt(Math.max(...betas), 2)}`,
      reversibility: memLong.length ? `${pct(med(memLong))} still there after regrowth` : "",
      paradox: "bulk, sorted and single-cell",
      singlecell: sc ? `${sc.cells_after_qc.toLocaleString()} cells and 2 PDXs` : "pending",
      graph: `${R.graph.nodes.length} nodes, ${R.graph.edges.length} edges`,
      passport: `${passports.length} genes`,
      arena: hyp("H1") ? `H1 ${hyp("H1").result}` : "",
      math: "eigenvalues decide relapse",
      methods: `${Object.values(R.gates).filter(g => g === "passed").length} gates passed`,
    };
    $("plateIndex").innerHTML = PLATES.map(p => `<li><a href="#${p.id}"><span class="pn">Plate ${p.dataset.plate}</span><span class="pt">${p.dataset.title}</span><span class="pv">${esc(live[p.id] || "")}</span>${spark(p.id)}</a></li>`).join("");
    $("platenav").innerHTML = PLATES.map(p => `<a href="#${p.id}">${p.dataset.title}</a>`).join("");
  }

  function dotField() {
    const c = $("dotField"), ctx = c.getContext("2d");
    const total = 248, nsclc = Math.round(total * .85), egfr = Math.round(total * .85 * .323), cols = 31;
    let t0 = null;
    function draw(ts) {
      const dpr = window.devicePixelRatio || 1, r = c.getBoundingClientRect();
      c.width = r.width * dpr; c.height = r.height * dpr;
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      const W = r.width, H = r.height, gap = Math.min(W / (cols + 1), (H - 30) / 8.4);
      const p = PT.reduced || t0 === null ? 1 : Math.min(1, (ts - t0) / 2800);
      ctx.clearRect(0, 0, W, H);
      for (let i = 0; i < total; i++) {
        const x = gap + (i % cols) * gap, y = gap * .8 + Math.floor(i / cols) * gap;
        const stage = i < egfr ? 2 : i < nsclc ? 1 : 0;
        const lit = Math.min(1, Math.max(0, p * 3 - stage));
        ctx.beginPath();
        ctx.arc(x, y, gap * .32, 0, Math.PI * 2);
        ctx.fillStyle = stage === 2 && lit > 0 ? `rgba(158,42,58,${.3 + .7 * lit})` : stage >= 1 && p > .34 ? "#3c6e8f" : "#bdb4a5";
        ctx.fill();
      }
      ctx.font = "italic 17px Garamond, serif";
      ctx.fillStyle = "#1c1a16";
      if (p >= 1) ctx.fillText(`about ${Math.round(egfr) * 10},000 EGFR-mutant NSCLC diagnoses a year, a rough estimate`, gap * .5, H - 6);
      if (p < 1) requestAnimationFrame(draw);
    }
    draw(0);
    const io = new IntersectionObserver(es => { if (es[0].isIntersecting) { io.disconnect(); requestAnimationFrame(ts => { t0 = ts; draw(ts); }); } }, { threshold: .3 });
    io.observe(c);
    addEventListener("resize", () => draw(1e12));
    const stats = [
      ["2.48 M", "new lung cancer cases in 2022", "https://gco.iarc.who.int/media/globocan/factsheets/cancers/15-trachea-bronchus-and-lung-fact-sheet.pdf", "GLOBOCAN 2022"],
      ["1.82 M", "lung cancer deaths in 2022, first among cancers", "https://gco.iarc.who.int/media/globocan/factsheets/cancers/15-trachea-bronchus-and-lung-fact-sheet.pdf", "GLOBOCAN 2022"],
      ["32.3%", "of NSCLC carries an EGFR mutation, pooled", "https://doi.org/10.18632/oncotarget.12587", "Zhang et al. 2016"],
      ["18.9 mo", "median progression-free survival on osimertinib", "https://pubmed.ncbi.nlm.nih.gov/29151359/", "FLAURA, NEJM 2018"],
    ];
    $("problemStats").innerHTML = stats.map(([v, l, u, s]) => `<div><b>${v}</b><span>${l}</span><a href="${u}" target="_blank" rel="noopener">${s}</a></div>`).join("");
  }

  function molecule() {
    const host = $("mol3d");
    const go = () => {
      if (!window.$3Dmol) { host.innerHTML = '<p class="note" style="padding:20px">The 3D viewer needs WebGL.</p>'; return; }
      fetch("media/4ZAU.pdb").then(r => r.text()).then(pdb => {
        const v = $3Dmol.createViewer(host, { backgroundColor: "#0d0c0a", antialias: true });
        v.addModel(pdb, "pdb");
        v.setStyle({}, {});
        v.setStyle({ chain: "A", hetflag: false }, { cartoon: { color: "#6b6356", opacity: .95 } });
        v.addStyle({ resn: "YY3", chain: "A" }, { stick: { colorscheme: { prop: "elem", map: { C: "#d9707b", N: "#8db4cf", O: "#efe8da" } }, radius: .24 } });
        v.addStyle({ chain: "A", resi: 797 }, { stick: { color: "#e0ad55", radius: .28 } });
        v.addLabel("Cys797", { position: { chain: "A", resi: 797, atom: "SG" }, backgroundColor: "#13110e", fontColor: "#e0ad55", fontSize: 14, borderThickness: 0, font: "Garamond" });
        v.zoomTo({ resn: "YY3", chain: "A" }); v.zoom(.42); v.render();
        if (!PT.reduced) v.spin("y", .35);
        host.addEventListener("pointerdown", () => v.spin(false), { once: true });
      }).catch(() => host.innerHTML = '<p class="note" style="padding:20px">Structure file missing in this build.</p>');
    };
    const io = new IntersectionObserver(es => { if (es[0].isIntersecting) { io.disconnect(); go(); } }, { rootMargin: "300px" });
    io.observe(host);
  }

  function timeline() {
    const svg = $("timeline");
    const xs = [90, 330, 600, 850, 1100];
    const days = ["day 0", "day 1", "day 21", "day 22", "day 24 to 31"];
    const what = ["untreated", "24 h, 500 nM osimertinib", "21 days on drug", "24 h after removal", "regrowth"];
    el("line", { x1: 40, x2: 1160, y1: 116, y2: 116, stroke: "var(--rule)" }, svg);
    el("rect", { x: xs[1] - 44, y: 108, width: xs[2] - xs[1] + 88, height: 16, fill: COLORS.dtp, opacity: .12 }, svg);
    el("text", { x: (xs[1] + xs[2]) / 2, y: 100, "text-anchor": "middle", class: "lbl" }, svg).textContent = "drug on";
    ORDER.forEach((c, i) => {
      const g = el("g", { class: "tl-step", style: `opacity:${PT.reduced ? 1 : 0};transition:opacity .7s ease ${i * .25}s` }, svg);
      el("circle", { cx: xs[i], cy: 116, r: 38, fill: "var(--paper)", stroke: COLORS[c], "stroke-width": 1.6 }, g);
      const dots = c === "dtp" ? 5 : c === "acute" ? 13 : c === "wo_short" ? 8 : c === "wo_long" ? 18 : 24;
      for (let k = 0; k < dots; k++) { const a = k * 2.399, rr = 6 + Math.sqrt(k) * 5.6; el("circle", { cx: xs[i] + Math.cos(a) * rr, cy: 116 + Math.sin(a) * rr, r: 3.2, fill: COLORS[c], opacity: .85 }, g); }
      el("text", { x: xs[i], y: 182, "text-anchor": "middle", class: "ttl" }, g).textContent = SHORT[c];
      el("text", { x: xs[i], y: 204, "text-anchor": "middle" }, g).textContent = days[i];
      el("text", { x: xs[i], y: 224, "text-anchor": "middle" }, g).textContent = what[i];
    });
    const io = new IntersectionObserver(es => { if (es[0].isIntersecting) { io.disconnect(); svg.querySelectorAll(".tl-step").forEach(g => g.style.opacity = 1); } }, { threshold: .3 });
    io.observe(svg);
  }

  function universe() {
    let mode = "pc", filt = "all";
    const a = document.createElement("span"), b = document.createElement("span");
    $("uniTools").append(a, b);
    group(a, [["pc", "centred per line"], ["pc_raw", "raw"]], mode, m => { mode = m; draw(); });
    group(b, [["all", "all lines"], ...lines.map(l => [l, l])], filt, f => { filt = f; draw(); });
    $("uniLegend").innerHTML = '<span class="legend">' + ORDER.map(c => `<span><i style="background:${COLORS[c]}"></i>${SHORT[c]}</span>`).join("") + "<span>● PC9 ■ H1975 ▲ HCC827 ◆ HCC2935</span></span>";
    const svg = $("uniChart");
    function shape(g, s, x, y, r, fill) {
      if (s.line === "H1975") return el("rect", { x: x - r, y: y - r, width: 2 * r, height: 2 * r, fill }, g);
      if (s.line === "HCC827") return el("path", { d: `M${x},${y - r * 1.2}L${x + r * 1.1},${y + r * .8}L${x - r * 1.1},${y + r * .8}Z`, fill }, g);
      if (s.line === "HCC2935") return el("path", { d: `M${x},${y - r * 1.3}L${x + r * 1.1},${y}L${x},${y + r * 1.3}L${x - r * 1.1},${y}Z`, fill }, g);
      return el("circle", { cx: x, cy: y, r, fill }, g);
    }
    function draw() {
      clear(svg);
      const box = { l: 64, r: 1080, t: 16, b: 470 }, S = R.samples;
      const x = scale(...extent(S.map(s => s[mode][0])), box.l, box.r), y = scale(...extent(S.map(s => s[mode][1])), box.b, box.t);
      const v = mode === "pc" ? R.meta.pca_variance : R.meta.pca_raw_variance;
      axes(svg, x, y, box, { xlab: `PC1, ${Math.round(v[0] * 100)}% of variance`, ylab: `PC2, ${Math.round(v[1] * 100)}%` });
      const g = el("g", {}, svg);
      for (const s of S) {
        const on = filt === "all" || s.line === filt;
        const n = shape(g, s, x(s[mode][0]), y(s[mode][1]), on ? 7.5 : 4, COLORS[s.condition]);
        n.setAttribute("opacity", on ? .92 : .12);
        n.style.cursor = "pointer";
        n.addEventListener("mousemove", e => showTip(e, `<b>${s.gsm}</b><br>${s.line}, ${SHORT[s.condition]}, replicate ${s.rep}<br>${(s.lib / 1e6).toFixed(1)} M reads`));
        n.addEventListener("mouseleave", hideTip);
        n.addEventListener("click", () => window.open(`https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${s.gsm}`, "_blank", "noopener"));
      }
    }
    draw();
  }

  function decomposition() {
    let line = lines[0];
    group($("decTools"), lines.map(l => [l, l]), line, l => { line = l; draw(); });
    const svg = $("decChart");
    function draw() {
      clear(svg);
      const box = { l: 58, r: 620, t: 14, b: 410 };
      const x = scale(-6, 6, box.l, box.r), y = scale(-8, 8, box.b, box.t);
      axes(svg, x, y, box, { xlab: "24-hour change, A (log2)", ylab: "persister change, P (log2)" });
      const g = el("g", {}, svg);
      for (const [a, p, s, sym] of R.decomposition_scatter[line] || []) {
        if (a === null || p === null) continue;
        const strong = Math.abs(s) > 1.5;
        const c = el("circle", { cx: x(Math.max(-6, Math.min(6, a))), cy: y(Math.max(-8, Math.min(8, p))), r: strong ? 2.8 : 1.5, fill: strong ? (s > 0 ? COLORS.dtp : COLORS.control) : "#a79e8f", opacity: strong ? .9 : .45 }, g);
        if (strong) { c.addEventListener("mousemove", e => showTip(e, `<b>${esc(sym)}</b><br>A ${fmt(a)}, P ${fmt(p)}, S ${fmt(s)}`)); c.addEventListener("mouseleave", hideTip); }
      }
      const b = R.betas[line];
      el("line", { x1: x(-6), y1: y(-6 * b), x2: x(6), y2: y(6 * b), stroke: COLORS.acute, "stroke-width": 2.4 }, svg);
      el("line", { x1: x(-6), y1: y(-6), x2: x(6), y2: y(6), stroke: "var(--rule)", "stroke-dasharray": "4 4" }, svg);
      el("text", { x: x(3.6), y: y(3.6 * b) + 22, class: "lbl" }, svg).textContent = `β = ${fmt(b, 2)}`;
      const sens = (R.beta_sensitivity || []).find(r => r.cell_line === line);
      $("betaText").innerHTML = `In <b>${line}</b> a gene's persister change is on average <b>${fmt(b, 2)}</b> times its 24-hour change. That is the stress persisters still carry. Points off the gold line changed for other reasons.${sens ? ` A plain least-squares slope gives ${fmt(sens.beta_ols, 2)}.` : ""}`;
    }
    const top = passports.filter(p => p.dimensions.persistence_enrichment.status !== "not supported")
      .sort((a, b) => Math.abs(b.dimensions.acute_stress_specificity.persister_specific_lfc || 0) - Math.abs(a.dimensions.acute_stress_specificity.persister_specific_lfc || 0)).slice(0, 12);
    $("decTop").innerHTML = `<table class="t"><thead><tr><th>gene</th><th class="num">P</th><th class="num">S</th></tr></thead><tbody>${top.map(p => `<tr><td><a href="#passport" data-gene="${p.symbol}">${p.symbol}</a></td><td class="num">${fmt(p.dimensions.persistence_enrichment.meta_lfc)}</td><td class="num">${fmt(p.dimensions.acute_stress_specificity.persister_specific_lfc)}</td></tr>`).join("")}</tbody></table>`;
    draw();
  }

  const CLASS_TEXT = {
    acute_only: "24-hour stress only", dtp_reversed_by_24h: "fades within 24 h", dtp_persists_24h_reversed_late: "holds a day, fades by regrowth",
    dtp_persists_after_regrowth: "still there after regrowth", dtp_overshoot: "overshoots past baseline", dtp_unresolved: "persister change, unresolved",
    not_changed: "not changed", context_dependent: "differs between lines", single_line: "one line only",
  };

  function reversibility() {
    const names = Object.keys(G);
    $("geneList").innerHTML = names.slice(0, 5000).map(n => `<option value="${esc(n)}">`).join("");
    const svg = $("revChart");
    $("revLegend").innerHTML = '<span class="legend">' + lines.map(l => `<span><i style="background:${LINE_COLORS[l]}"></i>${l}</span>`).join("") + "</span>";
    function draw(name) {
      const g = G[name];
      clear(svg);
      if (!g) { el("text", { x: 320, y: 190, "text-anchor": "middle", class: "lbl" }, svg).textContent = "not in the tested gene set"; return; }
      const box = { l: 58, r: 620, t: 20, b: 320 };
      const vals = Object.values(g.t).flat().filter(v => v !== null);
      const x = scale(0, 4, box.l + 30, box.r - 30), y = scale(...extent(vals.concat([0]), .15), box.b, box.t);
      axes(svg, x, y, box, { noX: true, ylab: "log2 change from control" });
      ORDER.forEach((c, i) => {
        el("rect", { x: x(i) - 46, y: box.t, width: 92, height: box.b - box.t, fill: COLORS[c], opacity: c === "dtp" ? .08 : .03 }, svg);
        el("text", { x: x(i), y: box.b + 20, "text-anchor": "middle", class: "lbl" }, svg).textContent = SHORT[c];
      });
      el("line", { x1: box.l, x2: box.r, y1: y(0), y2: y(0), stroke: "var(--rule)", "stroke-dasharray": "3 4" }, svg);
      for (const l of lines) {
        const t = g.t[l];
        if (!t) continue;
        const p = el("path", { d: t.map((v, i) => `${i ? "L" : "M"}${x(i)},${y(v ?? 0)}`).join(""), fill: "none", stroke: LINE_COLORS[l], "stroke-width": 2.4 }, svg);
        if (!PT.reduced) { const len = p.getTotalLength(); p.style.strokeDasharray = len; p.style.strokeDashoffset = len; p.getBoundingClientRect(); p.style.transition = "stroke-dashoffset 1s ease"; p.style.strokeDashoffset = 0; }
        t.forEach((v, i) => el("circle", { cx: x(i), cy: y(v ?? 0), r: 3.8, fill: LINE_COLORS[l] }, svg));
      }
      el("text", { x: box.l + 10, y: box.t + 18, class: "ttl" }, svg).textContent = name;
      $("geneClass").innerHTML = `<span class="badge">${esc(CLASS_TEXT[g.c] || g.c || "n/a")}</span>`;
    }
    const input = $("geneSearch");
    input.addEventListener("change", () => { const v = input.value.trim(); draw(v in G ? v : v.toUpperCase()); });
    const start = G.TACSTD2 ? "TACSTD2" : names[0];
    input.value = start; draw(start);
    window.PTdrawGene = n => { input.value = n; draw(n); };

    const ms = R.memory.summary, svgM = $("memChart");
    const box = { l: 44, r: 450, t: 10, b: 250 };
    const xs = scale(-1, 2, box.l, box.r), y = scale(0, 1, box.b, box.t);
    axes(svgM, xs, y, box, { xlab: "share of the persister change remaining" });
    const bw = (box.r - box.l) / 30;
    for (const tp of ["wo_short", "wo_long"]) {
      const sum = Array(30).fill(0);
      ms.filter(r => r.timepoint === tp).forEach(r => r.hist.forEach((v, i) => sum[i] += v));
      const m = Math.max(...sum, 1);
      sum.forEach((v, i) => el("rect", { x: box.l + i * bw, y: y(v / m), width: bw - 1, height: box.b - y(v / m), fill: COLORS[tp], opacity: .55 }, svgM));
    }
    el("line", { x1: xs(1), x2: xs(1), y1: box.t, y2: box.b, stroke: COLORS.dtp, "stroke-dasharray": "3 3" }, svgM);
    el("line", { x1: xs(0), x2: xs(0), y1: box.t, y2: box.b, stroke: COLORS.control, "stroke-dasharray": "3 3" }, svgM);
    el("text", { x: xs(1) + 5, y: box.t + 14, class: "lbl" }, svgM).textContent = "persister level";
    el("text", { x: xs(0) + 5, y: box.t + 14, class: "lbl" }, svgM).textContent = "baseline";
    const per = lines.map(l => { const r = ms.find(q => q.cell_line === l && q.timepoint === "wo_long"); const d = R.memory.decay[l]; return `${l} ${pct(r && r.median_memory)}${d && d.median_half_life_days ? ` (half-life ${fmt(d.median_half_life_days, 1)} d)` : ""}`; });
    $("memText").innerHTML = `Median memory 24 h after washout: <b>${pct(med(memShort))}</b>. After regrowth: ${per.join(", ")}. <span class="k" style="color:${COLORS.wo_short}">Violet</span> 24 h, <span class="k" style="color:${COLORS.wo_long}">green</span> regrowth.`;
    const dist = R.reversibility.distribution;
    const keys = Object.keys(dist).filter(k => k !== "not_changed").sort((a, b) => dist[b] - dist[a]);
    const mx = Math.max(1, ...keys.map(k => dist[k]));
    $("classBars").innerHTML = keys.map(k => `<div style="display:grid;grid-template-columns:190px 1fr 56px;gap:8px;align-items:center;font-size:16px;margin:3px 0"><span>${CLASS_TEXT[k] || k}</span><span style="height:7px;background:${k.includes("after_regrowth") ? COLORS.dtp : k.includes("acute") ? COLORS.acute : COLORS.wo_short};width:${Math.max(1, dist[k] / mx * 100)}%"></span><span class="num" style="text-align:right">${dist[k].toLocaleString()}</span></div>`).join("");
  }

  const ACR = { dna: "DNA", e2f: "E2F", g2m: "G2/M", tnfa: "TNFα", nfkb: "NF-κB", myc: "MYC", uv: "UV", mtorc1: "mTORC1", p53: "p53", tgf: "TGF", il6: "IL-6", jak: "JAK", stat3: "STAT3", kras: "KRAS", pi3k: "PI3K", akt: "AKT", mtor: "mTOR", dn: "down" };
  const nice = p => { const t = p.replace("HALLMARK_", "").toLowerCase().split("_").map(w => ACR[w] || w).join(" "); return t.charAt(0).toUpperCase() + t.slice(1); };

  function paradox() {
    const tc = R.programs.time_course || {};
    const progs = Object.keys(tc[lines[0]] || {});
    const svg = $("progChart");
    if (progs.length) {
      const pref = ["HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_DNA_REPAIR", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_FATTY_ACID_METABOLISM"];
      const choices = pref.filter(p => progs.includes(p)).concat(progs.filter(p => !pref.includes(p))).slice(0, 8);
      let cur = choices[0];
      group($("progTools"), choices.map(p => [p, nice(p)]), cur, p => { cur = p; draw(); });
      $("progLegend").innerHTML = '<span class="legend">' + lines.map(l => `<span><i style="background:${LINE_COLORS[l]}"></i>${l}</span>`).join("") + "</span>";
      const draw = () => {
        clear(svg);
        const box = { l: 58, r: 620, t: 20, b: 320 };
        const vals = lines.flatMap(l => ORDER.map(c => tc[l]?.[cur]?.[c])).filter(Number.isFinite);
        const x = scale(0, 4, box.l + 30, box.r - 30), y = scale(...extent(vals.concat([0]), .15), box.b, box.t);
        axes(svg, x, y, box, { noX: true, ylab: "mean gene z-score" });
        ORDER.forEach((c, i) => el("text", { x: x(i), y: box.b + 20, "text-anchor": "middle", class: "lbl" }, svg).textContent = SHORT[c]);
        for (const l of lines) {
          const v = ORDER.map(c => tc[l]?.[cur]?.[c]);
          el("path", { d: v.map((q, i) => `${i ? "L" : "M"}${x(i)},${y(q ?? 0)}`).join(""), fill: "none", stroke: LINE_COLORS[l], "stroke-width": 2.4 }, svg);
          v.forEach((q, i) => el("circle", { cx: x(i), cy: y(q ?? 0), r: 3.8, fill: LINE_COLORS[l] }, svg));
        }
        el("text", { x: box.l + 10, y: box.t + 18, class: "ttl" }, svg).textContent = nice(cur);
      };
      draw();
      const rec = key => { const k = progs.find(q => q.endsWith(key)); if (!k) return null; const f = lines.map(l => { const v = tc[l][k]; return (v.dtp - v.acute) / ((v.control - v.acute) || 1e-9); }).filter(Number.isFinite); return f.length ? f.reduce((a, b) => a + b, 0) / f.length : null; };
      const e2f = rec("E2F_TARGETS"), g2m = rec("G2M_CHECKPOINT");
      const tf = (R.tf_activity || []).find(t => t.source === "E2F1");
      $("progLit").innerHTML = e2f !== null ? `<p class="note"><b>The paradox in numbers.</b> After 24 h of drug, proliferation programmes collapse. In 21-day persisters, still on drug, E2F targets are <b>${pct(e2f)}</b> of the way back to normal and G2/M genes <b>${pct(g2m)}</b>.${tf ? ` Inferred E2F1 activity is higher than stress predicts (score ${fmt(tf.score)}).` : ""}</p>` : "";
    }
    const rows = R.programs.persister_specific || [];
    $("progTable").innerHTML = `<table class="t"><tbody>${rows.slice(0, 9).map(r => `<tr><td>${nice(r.source)}</td><td class="num" style="color:${r.score > 0 ? COLORS.dtp : COLORS.control}">${r.score > 0 ? "▲" : "▼"} ${fmt(Math.abs(r.score), 1)}</td></tr>`).join("")}</tbody></table>`;
    const cc = EXT.cell_cycle_sorted, svg2 = $("ccChart");
    if (!cc) { svg2.outerHTML = '<p class="note">The cell-cycle sorted dataset was not processed in this build.</p>'; return; }
    const box = { l: 60, r: 1080, t: 34, b: 250 };
    const phases = ["G1", "S", "G2"], sigs = ["naive_treated_vs_control", "persister_specific"];
    const x = scale(-.5, 5.5, box.l, box.r), y = scale(...extent(cc.samples.map(s => s.score).concat([0]), .15), box.b, box.t);
    axes(svg2, x, y, box, { noX: true, ylab: "signature score" });
    el("line", { x1: box.l, x2: box.r, y1: y(0), y2: y(0), stroke: "var(--rule)" }, svg2);
    el("line", { x1: x(2.5), x2: x(2.5), y1: box.t, y2: box.b, stroke: "var(--rule)" }, svg2);
    sigs.forEach((sg, si) => phases.forEach((ph, pi) => {
      const i = si * 3 + pi;
      el("text", { x: x(i), y: box.b + 22, "text-anchor": "middle", class: "lbl" }, svg2).textContent = ph;
      cc.samples.filter(s => s.signature === sg && s.phase === ph).forEach((s, k) => {
        const n = el("circle", { cx: x(i) + (s.day === 28 ? 16 : -16) + (k % 2) * 6, cy: y(s.score), r: 6.5, fill: s.day === 28 ? COLORS.dtp : COLORS.control, opacity: .9 }, svg2);
        n.addEventListener("mousemove", e => showTip(e, `${esc(s.sample)}<br>score ${fmt(s.score)}`)); n.addEventListener("mouseleave", hideTip);
      });
    }));
    el("text", { x: x(1), y: box.t - 12, "text-anchor": "middle", class: "ttl" }, svg2).textContent = "naive treated-vs-control signature";
    el("text", { x: x(4), y: box.t - 12, "text-anchor": "middle", class: "ttl" }, svg2).textContent = "persister-specific signature";
    const sh = cc.summary.filter(s => s.signature === "persister_specific");
    $("ccCap").innerHTML = `<b>Fig.</b> PC9 sorted by cell-cycle phase, <span class="k k-blue">day 0</span> against <span class="k k-red">day 28 on osimertinib</span> (GSE335846, May, Martin and Houseley, Babraham Institute, 2026). Signatures frozen from GSE193258. Day-28 shift in the persister-specific score: ${sh.map(s => `${s.phase} ${s.shift > 0 ? "+" : ""}${fmt(s.shift, 2)}`).join(", ")}. If S- and G2-phase persisters shift as much as G1 ones, the cycling persisters carry the programme too.`;
  }

  function singleCell() {
    const sc = EXT.pc9_single_cell, pdx = EXT.pdx_single_cell, canvas = $("scCanvas");
    if (!sc && !pdx) { canvas.outerHTML = `<p class="note">Single-cell validation did not run in this build${EXT.errors ? ": " + esc(Object.values(EXT.errors).join("; ")) : ""}.</p>`; return; }
    let ds = sc ? "pc9" : "pdx", color = "exposure";
    const a = document.createElement("span"), b = document.createElement("span");
    $("scTools").append(a, b);
    const opts = []; if (sc) opts.push(["pc9", "PC9 cells in culture"]); if (pdx) opts.push(["pdx", "PDX tumours"]);
    group(a, opts, ds, v => { ds = v; draw(); });
    group(b, [["exposure", "treatment"], ["cycling", "cycling"], ["score", "persister score"]], color, v => { color = v; draw(); });
    const ctx = canvas.getContext("2d");
    function draw() {
      const dpr = window.devicePixelRatio || 1, r = canvas.getBoundingClientRect();
      canvas.width = r.width * dpr; canvas.height = r.height * dpr; ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      const W = r.width, H = r.height, D = ds === "pc9" ? sc : pdx, pts = D.points;
      ctx.clearRect(0, 0, W, H);
      const xs = extent(pts.map(p => p[0]), .04), ys = extent(pts.map(p => p[1]), .04);
      const X = v => 12 + (v - xs[0]) / (xs[1] - xs[0]) * (W - 24), Y = v => H - 12 - (v - ys[0]) / (ys[1] - ys[0]) * (H - 24);
      const sv = pts.map(p => p[5]).filter(Number.isFinite), lo = Math.min(...sv), hi = Math.max(...sv);
      for (const p of pts) {
        let c;
        if (color === "exposure") c = ds === "pc9" ? (!p[2] ? "#4a443b" : p[3] ? NIGHT.dtp : NIGHT.control) : (p[3] === "residual" ? NIGHT.dtp : NIGHT.control);
        else if (color === "cycling") c = p[4] ? NIGHT.acute : "#4a443b";
        else { const t = (p[5] - lo) / (hi - lo || 1); c = `rgb(${Math.round(55 + 170 * t)},${Math.round(60 + 52 * t)},${Math.round(78 + 20 * t)})`; }
        ctx.fillStyle = c; ctx.globalAlpha = .82;
        ctx.beginPath(); ctx.arc(X(p[0]), Y(p[1]), 2, 0, Math.PI * 2); ctx.fill();
      }
      ctx.globalAlpha = 1;
      const legends = {
        exposure: ds === "pc9" ? `<span class="k k-blue">unexposed PC9</span>, <span class="k k-red">drug-exposed PC9</span> (inferred from MAPK output), grey: the second cell line in the pool` : `<span class="k k-blue">untreated</span> and <span class="k k-red">residual after six weeks</span>, tumour cells from YU-003 and YU-006`,
        cycling: `<span class="k k-gold">cycling</span> cells, with E2F or G2/M programmes high`,
        score: "persister-specific signature score, dark to light",
      };
      $("scCap").innerHTML = `<b>Fig.</b> UMAP of ${pts.length.toLocaleString()} sampled cells. ${legends[color]}.`;
    }
    const notes = [];
    if (sc) {
      const st = (sc.states || []);
      notes.push(`<p class="note"><b>PC9 under osimertinib.</b> ${sc.cells_after_qc.toLocaleString()} cells passed quality control; ${sc.exposed_cells.toLocaleString()} are drug-exposed PC9 cells. Of these, <b>${pct(sc.cycling_fraction_exposed)}</b> are cycling, against ${pct(sc.cycling_fraction_unexposed)} of unexposed PC9.</p>`);
      if (st.length) notes.push(`<p class="note"><b>One state or several?</b> The exposed cells fall into ${st.length} clusters whose cycling share runs from ${pct(Math.min(...st.map(s => s.cycling_fraction)))} to ${pct(Math.max(...st.map(s => s.cycling_fraction)))}. Persister score against cycling score: Spearman ρ ${fmt(sc.persister_vs_cycling_spearman, 2)}. Cluster agreement across resolutions, ARI ${fmt(sc.cluster_stability_ari, 2)}.</p>`);
      notes.push(`<p class="note"><b>Signature transfer.</b> Exposed against unexposed PC9 cells, AUROC: naive ${fmt(sc.auroc_exposed_vs_unexposed.naive_treated_vs_control, 2)}, persister-specific ${fmt(sc.auroc_exposed_vs_unexposed.persister_specific, 2)}. MAPK genes were removed from both signatures first, because they define exposure.</p>`);
      notes.push(`<p class="note" style="font-style:italic">${esc(sc.caveat)}</p>`);
    }
    if (pdx) {
      const au = pdx.auroc_residual_vs_untreated;
      notes.push(`<p class="note"><b>Patient-derived tumours.</b> Residual against untreated tumour cells, AUROC per PDX. Naive: ${Object.entries(au.naive_treated_vs_control || {}).map(([k, v]) => `${k} ${fmt(v, 2)}`).join(", ")}. Persister-specific: ${Object.entries(au.persister_specific || {}).map(([k, v]) => `${k} ${fmt(v, 2)}`).join(", ")}. Cycling share: ${Object.entries(pdx.cycling_fraction).map(([k, v]) => `${k} ${pct(v)}`).join(", ")}.</p>`);
    }
    $("scNotes").innerHTML = notes.join("");
    draw();
    addEventListener("resize", draw);
  }

  function graph() {
    const svg = $("graphSvg");
    const types = Object.keys(R.graph.edge_types), on = new Set(types);
    types.forEach(t => {
      const bt = PT.button(t.replace(/_/g, " "), true, () => { on.has(t) ? on.delete(t) : on.add(t); bt.classList.toggle("on"); render(); });
      bt.title = R.graph.edge_types[t];
      $("graphTools").appendChild(bt);
    });
    const nodes = R.graph.nodes.map(n => ({ ...n }));
    const idx = Object.fromEntries(nodes.map((n, i) => [n.id, i]));
    const edges = R.graph.edges.filter(e => e.source in idx && e.target in idx);
    const W = 1000, H = 640, cx = W / 2, cy = H / 2;
    const order = { contested: 0, dependency: 1, negative: 2, regulator: 3, marker: 4, stress: 5, weak: 6 };
    const genesN = nodes.filter(n => n.kind === "gene").sort((a, b) => (a.direction > b.direction ? 1 : -1) || (order[a.claim] - order[b.claim]));
    const N = genesN.length, ups = genesN.filter(n => n.direction === "up").length;
    const keyFirst = list => { const k = list.filter(n => order[n.claim] < 4), r = list.filter(n => order[n.claim] >= 4), out = []; const step = r.length / (k.length + 1); let j = 0; k.forEach((n, i) => { while (j < Math.round((i + 1) * step)) out.push(r[j++]); out.push(n); }); while (j < r.length) out.push(r[j++]); return out; };
    keyFirst(genesN.slice(0, ups)).concat(keyFirst(genesN.slice(ups))).forEach((n, i) => { n.ang = -Math.PI / 2 + (i + .5) / N * Math.PI * 2; n.x = cx + Math.cos(n.ang) * 215; n.y = cy + Math.sin(n.ang) * 190; });
    const nb = {};
    edges.forEach(e => { (nb[e.source] ??= []).push(e.target); (nb[e.target] ??= []).push(e.source); });
    const outer = nodes.filter(n => n.kind !== "gene" && n.kind !== "state");
    outer.forEach(n => { const g = (nb[n.id] || []).map(id => nodes[idx[id]]).filter(m => m.kind === "gene"); n.ang = g.length ? Math.atan2(g.reduce((a, m) => a + Math.sin(m.ang), 0), g.reduce((a, m) => a + Math.cos(m.ang), 0)) : 0; });
    outer.sort((a, b) => a.ang - b.ang);
    for (let it = 0; it < 200; it++) for (let i = 0; i < outer.length; i++) { const a = outer[i], b = outer[(i + 1) % outer.length]; let d = b.ang - a.ang; if (i === outer.length - 1) d += Math.PI * 2; if (outer.length > 1 && d < .17) { a.ang -= (.17 - d) / 2; b.ang += (.17 - d) / 2; } }
    outer.forEach(n => { n.x = cx + Math.cos(n.ang) * 365; n.y = cy + Math.sin(n.ang) * 285; });
    nodes.filter(n => n.kind === "state").forEach(n => { n.x = cx; n.y = cy; });
    const KC = { state: COLORS.dtp, gene: "#1c1a16", tf: COLORS.acute, program: COLORS.wo_long, source: COLORS.control };
    const ES = { observed_expression: [COLORS.dtp, ""], inferred_regulatory: [COLORS.acute, "5 4"], pathway_membership: [COLORS.wo_long, "2 3"], published_perturbation: [COLORS.control, ""], published_regulatory: [COLORS.control, "8 3"], published_chromatin: [COLORS.wo_short, ""] };
    $("graphLegend").innerHTML = '<span class="legend">' + types.map(t => `<span><svg width="26" height="8"><line x1="0" y1="4" x2="26" y2="4" stroke="${ES[t][0]}" stroke-width="2" stroke-dasharray="${ES[t][1]}"/></svg> ${t.replace(/_/g, " ")}</span>`).join("") + "</span>";
    function render() {
      clear(svg);
      const ge = el("g", {}, svg), gn = el("g", {}, svg);
      for (const e of edges) {
        if (!on.has(e.type)) continue;
        const a = nodes[idx[e.source]], b = nodes[idx[e.target]];
        const ln = el("line", { x1: a.x, y1: a.y, x2: b.x, y2: b.y, stroke: ES[e.type][0], "stroke-dasharray": ES[e.type][1], "stroke-width": e.type === "published_perturbation" ? 2 : 1, opacity: e.type === "observed_expression" ? .28 : .8 }, ge);
        ln.addEventListener("mousemove", ev => showTip(ev, `<b>${e.type.replace(/_/g, " ")}</b><br>${esc(e.method)}${e.effect ? "<br>effect: " + esc(e.effect) : ""}${e.outcome ? "<br>outcome: " + esc(e.outcome) : ""}`));
        ln.addEventListener("mouseleave", hideTip);
      }
      for (const n of nodes) {
        const r = n.kind === "state" ? 18 : n.kind === "gene" ? 4.5 : 6.5;
        const c = el("circle", { cx: n.x, cy: n.y, r, fill: KC[n.kind], stroke: n.claim === "contested" ? COLORS.acute : n.claim === "dependency" ? COLORS.wo_long : "none", "stroke-width": 2.5, tabindex: 0 }, gn);
        if (n.kind === "gene") {
          const deg = n.ang * 180 / Math.PI, flip = Math.cos(n.ang) < 0, tx = n.x + Math.cos(n.ang) * 9, ty = n.y + Math.sin(n.ang) * 9;
          el("text", { x: tx, y: ty, "dominant-baseline": "middle", "text-anchor": flip ? "end" : "start", style: `font-size:11px;fill:${order[n.claim] < 4 ? "#1c1a16" : "#8a8273"};font-weight:${order[n.claim] < 4 ? 700 : 400}`, transform: `rotate(${flip ? deg + 180 : deg} ${tx} ${ty})` }, gn).textContent = n.label;
        } else {
          const left = n.x < cx - 5;
          el("text", { x: n.x + (n.kind === "state" ? 0 : left ? -r - 4 : r + 4), y: n.y + (n.kind === "state" ? 36 : 4), "text-anchor": n.kind === "state" ? "middle" : left ? "end" : "start", class: "lbl" }, gn).textContent = n.label;
        }
        c.style.cursor = "pointer";
        c.addEventListener("mousemove", ev => showTip(ev, `<b>${esc(n.label)}</b><br>${n.kind}${n.doi ? "<br>doi " + esc(n.doi) : ""}`));
        c.addEventListener("mouseleave", hideTip);
        const open = () => { if (n.kind === "gene") { openPassport(n.label); location.hash = "passport"; } else if (n.doi) window.open("https://doi.org/" + n.doi, "_blank", "noopener"); };
        c.addEventListener("click", open);
        c.addEventListener("keydown", ev => { if (ev.key === "Enter") open(); });
      }
    }
    render();
  }

  const DIM_TITLE = { persistence_enrichment: "persister enrichment", acute_stress_specificity: "beyond 24-hour stress", reversibility: "after washout", regulatory_support: "regulators", perturbation_support: "perturbation", reproducibility: "reproducibility" };
  function openPassport(sym) {
    const p = passports.find(q => q.symbol === sym);
    if (!p) return;
    document.querySelectorAll("#plist button").forEach(b => b.classList.toggle("on", b.dataset.sym === sym));
    const d = p.dimensions;
    const dimHtml = Object.entries(d).map(([k, v]) => {
      let body = "";
      if (k === "persistence_enrichment") body = `log2FC ${fmt(v.meta_lfc)}, FDR ${fmt(v.fdr)}, ${v.lines_agreeing} of ${lines.length} lines`;
      if (k === "acute_stress_specificity") body = `specific ${fmt(v.persister_specific_lfc)}, stress ${fmt(v.acute_explained_lfc)}`;
      if (k === "reversibility") body = esc(CLASS_TEXT[v.class] || v.class);
      if (k === "regulatory_support") body = (v.regulators || []).slice(0, 3).map(r => `${r.tf} (${r.weight > 0 ? "activates" : "represses"})`).join(", ") || (v.literature.length ? "published record" : "none found");
      if (k === "perturbation_support") body = v.records.length ? `${v.records.length} published record${v.records.length > 1 ? "s" : ""}` : "none found";
      if (k === "reproducibility") body = `${v.lines} of ${v.of} lines, I² ${fmt(v.I2)}`;
      return `<div class="dim"><h4>${DIM_TITLE[k]}</h4><span class="badge ${slug(v.status)}">${esc(v.status)}</span><div style="margin-top:5px">${body}</div></div>`;
    }).join("");
    const recs = d.perturbation_support.records.concat(d.regulatory_support.literature || []);
    const recHtml = recs.length ? `<p class="smallcaps">Published evidence</p><ul class="records">${recs.map(r => { const s = R.sources[r.source] || {}; return `<li><i>${esc(r.kind.replace(/_/g, " "))}</i>, ${esc(r.model || "")}. Outcome: ${esc(r.outcome || "")}. Effect: <b>${esc(r.effect || "")}</b>. ${esc(r.summary || "")} <a href="https://doi.org/${s.doi}" target="_blank" rel="noopener">${esc(s.journal)} ${s.year}</a></li>`; }).join("")}</ul>` : "";
    const contra = p.contradictions.length ? `<p class="smallcaps">Contradictions and caveats</p><ul class="contra">${p.contradictions.map(c => `<li>${esc(c)}</li>`).join("")}</ul>` : "";
    const perLine = Object.entries(p.per_line).map(([l, v]) => `<tr><td>${l}</td><td class="num">${fmt(v.A)}</td><td class="num">${fmt(v.P)}</td><td class="num">${fmt(v.S)}</td><td>${esc(CLASS_TEXT[v.class] || v.class || "")}</td></tr>`).join("");
    $("pcard").innerHTML = `<div style="display:flex;justify-content:space-between;align-items:baseline;gap:10px;flex-wrap:wrap"><h3>${p.symbol}</h3><button class="ghost" id="toRev">see its trajectory</button></div>
      <div class="verdict"><p class="smallcaps">Strongest claim the evidence allows</p><b>${esc(p.claim.label)}</b><span class="note">${esc(p.claim.meaning)}</span></div>
      <div class="dims">${dimHtml}</div>${contra}${recHtml}
      <p class="smallcaps">Per cell line, log2 fold change</p><table class="t"><thead><tr><th>line</th><th class="num">A, 24 h</th><th class="num">P, persister</th><th class="num">S, specific</th><th>after washout</th></tr></thead><tbody>${perLine}</tbody></table>
      <p class="smallcaps">Evidence-to-claim trace</p><div class="trace"><span>${esc(p.claim.label)}</span><span class="arr">←</span><span>${esc(p.trace.contrasts.join(", "))}</span><span class="arr">←</span><span>${esc(p.trace.test)}</span><span class="arr">←</span><span>${p.trace.samples.length} samples</span><span class="arr">←</span><span><a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${p.trace.dataset}" target="_blank" rel="noopener">${p.trace.dataset}</a></span><span class="arr">←</span><span><a href="https://doi.org/${p.trace.publication.doi}" target="_blank" rel="noopener">${esc(p.trace.publication.journal)} ${p.trace.publication.year}</a></span></div>`;
    $("toRev").addEventListener("click", () => { window.PTdrawGene(p.symbol); location.hash = "reversibility"; });
  }

  function passportPanel() {
    const t = passports.find(q => q.symbol === "TACSTD2");
    if (t) {
      const up = Object.entries(t.per_line).filter(([, v]) => v.P > .5).map(([l]) => l);
      const flat = Object.entries(t.per_line).filter(([, v]) => v.P <= .5).map(([l]) => l);
      const share = Object.values(t.per_line).filter(v => v.P > .5).map(v => v.A / v.P);
      $("trop2Plain").insertAdjacentHTML("beforeend", ` In these 60 samples it rose in persisters of ${up.join(" and ") || "no line"}${flat.length ? `, not in ${flat.join(" or ")}` : ""}${share.length && Math.min(...share) > .7 ? ", and most of that rise was already there after 24 hours" : ""}.`);
    }
    const labels = [...new Set(passports.map(p => p.claim.label))];
    $("claimFilter").innerHTML += labels.map(l => `<option>${esc(l)}</option>`).join("");
    const rank = { contested: 0, dependency: 1, negative: 2, regulator: 3, marker: 4, stress: 5, weak: 6 };
    const sorted = passports.slice().sort((a, b) => rank[a.claim.key] - rank[b.claim.key] || Math.abs(b.dimensions.persistence_enrichment.meta_lfc || 0) - Math.abs(a.dimensions.persistence_enrichment.meta_lfc || 0));
    function list(f) {
      $("plist").innerHTML = sorted.filter(p => !f || p.claim.label === f).map(p => `<button role="option" data-sym="${p.symbol}"><span>${p.symbol} <small style="color:var(--muted)">${p.direction === "up" ? "▲" : "▼"}</small></span><span class="badge ${p.claim.key}">${esc(p.claim.label)}</span></button>`).join("");
      $("plist").querySelectorAll("button").forEach(b => b.addEventListener("click", () => openPassport(b.dataset.sym)));
    }
    $("claimFilter").addEventListener("change", e => list(e.target.value));
    list("");
    $("trapBtn").addEventListener("click", () => openPassport(t ? "TACSTD2" : sorted[0].symbol));
    document.addEventListener("click", e => { const a = e.target.closest("[data-gene]"); if (a) openPassport(a.dataset.gene); });
    openPassport(t ? "TACSTD2" : sorted[0]?.symbol);
  }

  function arena() {
    $("hyps").innerHTML = R.hypotheses.map(h => {
      let detail = "";
      if (h.naive !== undefined && h.naive !== null) detail = `naive ${fmt(h.naive, 3)} against persister-specific ${fmt(h.decomposed, 3)}`;
      if (h.classes) detail = Object.entries(h.classes).filter(([k]) => k !== "not_changed").map(([k, v]) => `${CLASS_TEXT[k] || k} ${v}`).join(", ");
      if (h.value !== undefined && h.value !== null) detail = `${h.value} regulators`;
      if (h.genes) detail = h.genes.join(", ") || "none";
      if (h.ranks) detail = h.ranks.map(r => `${r.gene} ${fmt(r.naive_percentile * 100, 0)}→${fmt(r.persister_specific_percentile * 100, 0)}%`).join(", ");
      return `<div class="hyp"><span class="id">${h.id}</span><span class="txt">${esc(h.text)}<span class="sub">${esc(h.metric || "")}${detail ? ": " + esc(detail) : ""}${h.note ? ". " + esc(h.note) : ""}</span></span><span class="res ${slug(h.result)}">${esc(h.result)}</span></div>`;
    }).join("");
    const svg = $("foldChart"), box = { l: 58, r: 450, t: 14, b: 240 }, folds = R.benchmark.folds;
    const x = scale(-.5, folds.length - .5, box.l, box.r), y = scale(0, 1, box.b, box.t);
    axes(svg, x, y, box, { noX: true, ylab: "AUROC" });
    el("line", { x1: box.l, x2: box.r, y1: y(.5), y2: y(.5), stroke: "var(--rule)", "stroke-dasharray": "3 3" }, svg);
    folds.forEach((f, i) => {
      const a = f["naive_treated_vs_control|auroc_dtp_vs_acute"] ?? 0, b = f["persister_specific|auroc_dtp_vs_acute"] ?? 0;
      el("rect", { x: x(i) - 30, y: y(a), width: 28, height: box.b - y(a), fill: COLORS.acute }, svg);
      el("rect", { x: x(i) + 2, y: y(b), width: 28, height: box.b - y(b), fill: COLORS.dtp }, svg);
      el("text", { x: x(i), y: box.b + 18, "text-anchor": "middle", class: "lbl" }, svg).textContent = f.held_out;
    });
    const rows = R.benchmark.summary;
    $("benchTable").innerHTML = `<p class="smallcaps">Five-condition classifiers on the held-out line</p><table class="t"><thead><tr><th>features</th><th>model</th><th class="num">balanced accuracy</th><th class="num">95% interval</th><th class="num">macro F1</th><th class="num">persister/acute confusion</th></tr></thead><tbody>${rows.map(r => `<tr><td>${r.features.replace(/_/g, " ")}</td><td>${r.model.replace(/_/g, " ")}</td><td class="num">${fmt(r.balanced_accuracy, 3)}</td><td class="num">${fmt(r.ci_low, 2)}–${fmt(r.ci_high, 2)}</td><td class="num">${fmt(r.macro_f1, 3)}</td><td class="num">${fmt(r.dtp_acute_confusion, 3)}</td></tr>`).join("")}</tbody></table><p class="note">Chance is 0.20. Intervals resample the 60 held-out predictions and are optimistic, because only four lines exist.</p>`;
  }

  function mathModel() {
    const P = { on: { r_S: -0.8, r_P: 0, a: 0.05, b: 0 }, off: { r_S: 0.6, r_P: 0.1, a: 0.001, b: 0.3 } };
    const kData = halfLives.length ? Math.LN2 / med(halfLives) : null;
    const sliders = [["on", "r_S", "death rate of sensitive cells on drug", -2, 0, .05], ["on", "a", "switching into the persister state, a", 0, .3, .005], ["off", "b", "switching back after washout, b", 0, 1.5, .01], ["off", "r_S", "regrowth rate after washout", 0, 1.2, .05]];
    $("modelControls").innerHTML = sliders.map(([ph, k, label, lo, hi, st], i) => `<label class="note" style="display:block">${label}: <b id="mv${i}">${P[ph][k]}</b><input type="range" min="${lo}" max="${hi}" step="${st}" value="${P[ph][k]}" id="ms${i}"></label>`).join("") +
      `<p class="note" id="eig"></p>` + (kData ? `<button class="ghost" id="useData">set b from the measured memory, ${fmt(kData, 3)} per day</button><p class="note" style="font-style:italic;margin-top:8px">This treats memory decay as a stand-in for the switch-back rate. It is an assumption, said out loud.</p>` : "");
    sliders.forEach(([ph, k], i) => $("ms" + i).addEventListener("input", e => { P[ph][k] = +e.target.value; $("mv" + i).textContent = e.target.value; draw(); }));
    if (kData) $("useData").addEventListener("click", () => { P.off.b = +kData.toFixed(3); $("ms2").value = P.off.b; $("mv2").textContent = P.off.b; draw(); });
    const lam = p => { const a = p.r_S - p.a, d = p.r_P - p.b, tr = a + d, det = a * d - p.b * p.a; return (tr + Math.sqrt(Math.max(tr * tr - 4 * det, 0))) / 2; };
    const svg = $("modelChart");
    function draw() {
      let s = .999, q = .001;
      const out = [[0, s, q]], dt = .02;
      for (let t = dt; t <= 42 + 1e-9; t += dt) { const p = t <= 21 ? P.on : P.off; const ds = (p.r_S - p.a) * s + p.b * q, dq = p.a * s + (p.r_P - p.b) * q; s += ds * dt; q += dq * dt; if (Math.abs(t * 4 - Math.round(t * 4)) < 1e-6) out.push([t, s, q]); }
      clear(svg);
      const bx = { l: 58, r: 620, t: 14, b: 310 };
      const vals = out.flatMap(o => [o[1], o[2], o[1] + o[2]]).filter(v => v > 0).map(Math.log10);
      const x = scale(0, 42, bx.l, bx.r), y = scale(Math.max(-8, Math.min(...vals)), Math.max(...vals) + .3, bx.b, bx.t);
      axes(svg, x, y, bx, { xlab: "days", ylab: "log10 cell fraction" });
      el("rect", { x: x(0), y: bx.t, width: x(21) - x(0), height: bx.b - bx.t, fill: NIGHT.dtp, opacity: .07 }, svg);
      el("text", { x: x(10.5), y: bx.t + 18, "text-anchor": "middle", class: "lbl" }, svg).textContent = "drug on";
      el("text", { x: x(31.5), y: bx.t + 18, "text-anchor": "middle", class: "lbl" }, svg).textContent = "drug off";
      const line = (j, c, w) => el("path", { d: out.map((o, i) => { const v = j === 3 ? o[1] + o[2] : o[j]; return `${i ? "L" : "M"}${x(o[0])},${y(Math.log10(Math.max(v, 1e-8)))}`; }).join(""), fill: "none", stroke: c, "stroke-width": w }, svg);
      line(1, NIGHT.control, 2.2); line(2, NIGHT.dtp, 2.2); line(3, "#efe8da", 1.4);
      const lo = lam(P.on), lf = lam(P.off), e21 = out.find(o => o[0] >= 21);
      $("eig").innerHTML = `λ<sub>max</sub> on drug <b>${fmt(lo, 3)}</b> per day (${lo < -1e-3 ? "shrinks" : lo > 1e-3 ? "grows" : "holds steady"}); off drug <b>${fmt(lf, 3)}</b>. Persister share at day 21: <b>${fmt(e21[2] / (e21[1] + e21[2]) * 100, 1)}%</b>.`;
    }
    draw();
  }

  function methods() {
    $("gates").innerHTML = Object.entries(R.gates).map(([k, v]) => `<div class="gate"><span>${k.replace(/^(\w)_/, "$1 · ").replace(/_/g, " ")}</span><span class="${v === "passed" ? "passed" : v === "failed" ? "failed" : ""}">${esc(v)}</span></div>`).join("");
    const rb = R.robustness;
    $("robust").innerHTML = `Withdrawal classes over ${rb.threshold_grid.length} threshold settings: <i>${esc(rb.labels.threshold_agreement)}</i>. Signature with one line dropped: <i>${esc(rb.labels.signature_lolo)}</i> (Jaccard ${rb.leave_one_line_out.map(r => r.dropped + " " + fmt(r.jaccard_top, 2)).join(", ")}). Gene lists move; programme-level results hold.`;
    const man = (R.manifest || []).concat(EXT.manifest || []);
    $("manifest").innerHTML = man.length ? `<table class="t"><tbody>${man.map(m => `<tr><td><a href="${esc(m.url)}" target="_blank" rel="noopener">${esc(m.path.split("/").pop())}</a><br><span style="font:12px ui-monospace,monospace;color:var(--muted)">${m.sha256.slice(0, 20)}…</span></td><td class="num">${(m.bytes / 1e6).toFixed(1)} MB</td></tr>`).join("")}</tbody></table>` : "";
    $("runinfo").innerHTML = `Generated ${esc(R.meta.generated_utc)}, Python ${esc(R.meta.python)}, pipeline ${esc(R.meta.version)}. Gene sets: ${esc(R.meta.networks.hallmark)}. Regulons: ${esc(R.meta.networks.collectri)}.`;
    const limits = [
      "No new wet-lab experiment: every functional claim comes from a cited paper.",
      "Bulk RNA cannot separate reversal from outgrowth; there is no lineage tracing here.",
      "Four cell lines give four test folds, so every accuracy has wide uncertainty.",
      "GSE275199 hashtags were not deposited, so drug exposure of single cells is inferred.",
      "Gene-level lists are threshold-sensitive; programme-level results are steadier.",
      "Cell lines and xenografts are not patients. Nothing here is clinical advice.",
    ];
    if (EXT.errors && Object.keys(EXT.errors).length) limits.unshift("This build failed on: " + Object.keys(EXT.errors).join(", ") + ".");
    $("limits").innerHTML = limits.map(l => `<li>${esc(l)}</li>`).join("");
    $("refs").innerHTML = Object.values(R.sources).map(s => `<div>${esc(s.title)}. <i>${esc(s.journal)}</i> ${s.year}. <a href="https://doi.org/${s.doi}" target="_blank" rel="noopener">doi</a></div>`).join("");
    $("trajVideo").src = `media/trajectory${synthetic ? "_synthetic" : ""}.mp4`;
  }

  function chrome() {
    document.querySelectorAll(".plate").forEach(s => s.classList.add("reveal"));
    const io = new IntersectionObserver(es => es.forEach(e => { if (e.isIntersecting) { e.target.classList.add("in"); io.unobserve(e.target); } }), { threshold: .03 });
    document.querySelectorAll(".reveal").forEach(s => io.observe(s));
    const links = [...document.querySelectorAll("#platenav a")], head = $("runhead");
    const so = new IntersectionObserver(es => es.forEach(e => {
      if (!e.isIntersecting) return;
      links.forEach(a => a.classList.toggle("on", a.getAttribute("href") === "#" + e.target.id));
      head.classList.toggle("dark", e.target.classList.contains("night"));
    }), { rootMargin: "-48px 0px -90% 0px" });
    document.querySelectorAll("main > section").forEach(s => so.observe(s));
    const bar = $("progress");
    addEventListener("scroll", () => { const h = document.documentElement; bar.style.width = (h.scrollTop / (h.scrollHeight - h.clientHeight) * 100) + "%"; }, { passive: true });
  }

  const steps = [ledger, plateIndex, () => Hero.init(R), dotField, molecule, timeline, universe, decomposition, reversibility, paradox, singleCell, graph, passportPanel, arena, mathModel, methods, chrome];
  for (const f of steps) { try { f(); } catch (err) { console.error(f.name || "step", err); } }
})();
