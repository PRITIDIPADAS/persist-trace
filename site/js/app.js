(async function () {
  const $ = id => document.getElementById(id);
  const { el, clear, scale, extent, axes, fmt, showTip, hideTip, group, esc, COLORS, ORDER, SHORT, LINE_COLORS } = PT;
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
      document.querySelectorAll("section.mod .card, .stage").forEach(c => c.innerHTML = '<div class="pending">The analysis has not been published yet. The GitHub Actions workflow downloads GSE193258 and fills this page with real results.</div>');
      return;
    }
  }
  if (synthetic || R.meta.synthetic) {
    const b = $("banner");
    b.textContent = "Synthetic test data. Nothing on this page is a real result. The real run replaces it automatically.";
    b.classList.add("show");
  }

  const passports = R.passports || [];
  const byClaim = k => passports.filter(p => p.claim.key === k).length;
  const lines = R.meta.cell_lines;
  const reproducibleDtp = Object.entries(R.reversibility.distribution).filter(([k]) => k.startsWith("dtp_")).reduce((a, [, v]) => a + v, 0);
  const memLong = R.memory.summary.filter(r => r.timepoint === "wo_long" && r.median_memory !== null).map(r => r.median_memory);
  const memShort = R.memory.summary.filter(r => r.timepoint === "wo_short" && r.median_memory !== null).map(r => r.median_memory);
  const med = a => a.length ? [...a].sort((x, y) => x - y)[Math.floor(a.length / 2)] : null;
  const h1 = R.hypotheses.find(h => h.id === "H1");
  const betaMean = Object.values(R.betas).reduce((a, b) => a + b, 0) / lines.length;

  function heroStats() {
    const stats = [
      ["", R.meta.n_samples, "RNA-seq samples"],
      ["", R.meta.n_genes_tested.toLocaleString(), "genes tested in all 4 lines"],
      ["acute", Math.round(betaMean * 100) + "%", "of the acute response still carried by persisters (β)"],
      ["dtp", reproducibleDtp.toLocaleString(), "persister genes reproducible in 3+ lines"],
      ["dtp", med(memShort) !== null ? Math.round(med(memShort) * 100) + "%" : "n/a", "of the persister change left 24 h after washout"],
      ["ok", passports.length, "evidence passports"],
    ];
    $("heroStats").innerHTML = stats.map(([c, v, l]) => `<div class="stat ${c}"><b>${v}</b><span>${l}</span></div>`).join("");
  }

  function tiles() {
    const bestBench = (R.benchmark.summary || [])[0];
    const t = [
      ["01", "question", "The question", "Stress, persistence and recovery, pulled apart.", `${lines.length} lines × 5 conditions`],
      ["02", "universe", "Cell-State Universe", "Every sample on one map.", `PC1 ${Math.round(R.meta.pca_variance[0] * 100)}% of variance`],
      ["03", "decompose", "Decomposition", "P = βA + S, per gene.", `β ${fmt(Math.min(...Object.values(R.betas)), 2)} to ${fmt(Math.max(...Object.values(R.betas)), 2)}`],
      ["04", "reversibility", "Reversibility map", "Memory after the drug is gone.", med(memLong) !== null ? `${Math.round(med(memLong) * 100)}% left after regrowth` : "see panel"],
      ["05", "graph", "Escape graph", "Measured and inferred edges, labelled.", `${R.graph.nodes.length} nodes · ${R.graph.edges.length} edges`],
      ["06", "paradox", "Replication paradox", "Cell cycle and repair through time.", "bulk average, honest limit"],
      ["07", "passport", "Evidence passports", "Marker or mechanism, six answers.", `${byClaim("dependency")} supported · ${byClaim("contested")} contested`],
      ["08", "arena", "Validation arena", "Scored on an unseen cell line.", h1 ? `H1: ${h1.result}` : ""],
      ["09", "math", "Switching math", "Eigenvalues decide relapse.", "interactive"],
      ["10", "methods", "Methods", "Checksums, gates, tests.", Object.values(R.gates).filter(g => g === "passed").length + " gates passed"],
    ];
    $("tiles").innerHTML = t.map(([n, id, h, p, live]) => `<a class="tile" href="#${id}"><div class="n">${n}</div><h3>${h}</h3><p>${p}</p><div class="live">${esc(live)}</div></a>`).join("");
  }

  function video() {
    const v = $("trajVideo");
    const tag = synthetic ? "_synthetic" : "";
    v.src = `media/trajectory${tag}.mp4`;
    v.poster = "";
    v.addEventListener("error", () => v.replaceWith(Object.assign(document.createElement("p"), { className: "sub", textContent: "Video not rendered in this build." })), { once: true });
  }

  function universe() {
    let mode = "pc", filt = "all";
    const tools = $("uniTools");
    const a = document.createElement("span"), b = document.createElement("span");
    tools.append(a, b);
    group(a, [["pc", "Centred per line"], ["pc_raw", "Raw"]], mode, m => { mode = m; draw(); });
    group(b, [["all", "All lines"], ...lines.map(l => [l, l])], filt, f => { filt = f; draw(); });
    $("uniLegend").innerHTML = ORDER.map(c => `<span><i style="background:${COLORS[c]}"></i>${SHORT[c]}</span>`).join("") + '<span>shape: ● PC9 ■ H1975 ▲ HCC827 ◆ HCC2935</span>';
    const svg = $("uniChart");
    function shape(g, s, x, y, r, fill) {
      const line = s.line;
      if (line === "H1975") return el("rect", { x: x - r, y: y - r, width: 2 * r, height: 2 * r, fill }, g);
      if (line === "HCC827") return el("path", { d: `M${x},${y - r * 1.2}L${x + r * 1.1},${y + r * .8}L${x - r * 1.1},${y + r * .8}Z`, fill }, g);
      if (line === "HCC2935") return el("path", { d: `M${x},${y - r * 1.3}L${x + r * 1.1},${y}L${x},${y + r * 1.3}L${x - r * 1.1},${y}Z`, fill }, g);
      return el("circle", { cx: x, cy: y, r, fill }, g);
    }
    function draw() {
      clear(svg);
      const box = { l: 60, r: 880, t: 20, b: 470 };
      const S = R.samples;
      const x = scale(...extent(S.map(s => s[mode][0])), box.l, box.r);
      const y = scale(...extent(S.map(s => s[mode][1])), box.b, box.t);
      const v = mode === "pc" ? R.meta.pca_variance : R.meta.pca_raw_variance;
      axes(svg, x, y, box, { xlab: `PC1 (${Math.round(v[0] * 100)}%)`, ylab: `PC2 (${Math.round(v[1] * 100)}%)` });
      const g = el("g", {}, svg);
      for (const s of S) {
        const on = filt === "all" || s.line === filt;
        const n = shape(g, s, x(s[mode][0]), y(s[mode][1]), on ? 7 : 4, COLORS[s.condition]);
        n.setAttribute("opacity", on ? .92 : .12);
        n.setAttribute("tabindex", on ? 0 : -1);
        n.addEventListener("mousemove", e => showTip(e, `<b>${s.gsm}</b><br>${s.line} · ${SHORT[s.condition]} · rep ${s.rep}<br>${(s.lib / 1e6).toFixed(1)} M reads${s.outlier ? "<br>QC flag: low within-line correlation" : ""}`));
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
      const box = { l: 56, r: 620, t: 16, b: 410 };
      const pts = R.decomposition_scatter[line] || [];
      const x = scale(-6, 6, box.l, box.r), y = scale(-8, 8, box.b, box.t);
      axes(svg, x, y, box, { xlab: "acute log2 fold change (A)", ylab: "persister log2 fold change (P)" });
      const g = el("g", {}, svg);
      for (const [a, p, s, sym] of pts) {
        if (a === null || p === null) continue;
        const strong = Math.abs(s) > 1.5;
        const c = el("circle", { cx: x(Math.max(-6, Math.min(6, a))), cy: y(Math.max(-8, Math.min(8, p))), r: strong ? 2.6 : 1.6, fill: strong ? (s > 0 ? COLORS.dtp : COLORS.control) : "#5b587a", opacity: strong ? .85 : .4 }, g);
        if (strong) c.addEventListener("mousemove", e => showTip(e, `<b>${esc(sym)}</b><br>A ${fmt(a)} · P ${fmt(p)} · S ${fmt(s)}`));
        c.addEventListener("mouseleave", hideTip);
      }
      const b = R.betas[line];
      el("line", { x1: x(-6), y1: y(-6 * b), x2: x(6), y2: y(6 * b), stroke: COLORS.acute, "stroke-width": 2 }, svg);
      el("line", { x1: x(-6), y1: y(-6), x2: x(6), y2: y(6), stroke: "#3a3757", "stroke-dasharray": "4 4" }, svg);
      el("text", { x: x(4.2), y: y(4.2 * b) - 8, class: "lbl" }, svg).textContent = `β = ${fmt(b, 2)}`;
      el("text", { x: x(4.6), y: y(5.4), fill: "#66637f" }, svg).textContent = "P = A";
      const sens = (R.beta_sensitivity || []).find(r => r.cell_line === line);
      $("betaText").innerHTML = `In ${line}, persisters keep about <b>${Math.round(b * 100)}%</b> of the size of the day-one response, on average across genes. Points far from the yellow line (coloured) changed in persisters for reasons the acute response does not explain.` +
        (sens ? ` A plain least-squares slope gives ${fmt(sens.beta_ols, 2)}; the residuals of the two fits correlate at r = ${fmt(sens.residual_correlation, 3)}.` : "");
    }
    const top = passports.filter(p => p.dimensions.persistence_enrichment.status !== "not supported").slice().sort((a, b) => Math.abs(b.dimensions.acute_stress_specificity.persister_specific_lfc || 0) - Math.abs(a.dimensions.acute_stress_specificity.persister_specific_lfc || 0)).slice(0, 16);
    $("decTop").innerHTML = `<table class="t"><thead><tr><th>Gene</th><th>P</th><th>S</th><th>Label</th></tr></thead><tbody>${top.map(p => `<tr><td><a href="#passport" data-gene="${p.symbol}">${p.symbol}</a></td><td class="num">${fmt(p.dimensions.persistence_enrichment.meta_lfc)}</td><td class="num">${fmt(p.dimensions.acute_stress_specificity.persister_specific_lfc)}</td><td><span class="badge ${p.claim.key}">${p.claim.label}</span></td></tr>`).join("")}</tbody></table>`;
    draw();
  }

  const CLASS_TEXT = {
    acute_only: "Acute only", dtp_reversed_by_24h: "Reverses within 24 h", dtp_persists_24h_reversed_late: "Holds 24 h, reverses by regrowth",
    dtp_persists_after_regrowth: "Still present after regrowth", dtp_overshoot: "Overshoots past baseline", dtp_unresolved: "Persister change, unresolved",
    not_changed: "Not changed", context_dependent: "Differs between lines", single_line: "One line only",
  };

  function reversibility() {
    const names = Object.keys(G);
    $("geneList").innerHTML = names.slice(0, 4000).map(n => `<option value="${esc(n)}">`).join("");
    const svg = $("revChart");
    $("revLegend").innerHTML = lines.map(l => `<span><i style="background:${LINE_COLORS[l]}"></i>${l}</span>`).join("");
    function draw(name) {
      const g = G[name];
      clear(svg);
      if (!g) { el("text", { x: 320, y: 190, "text-anchor": "middle" }, svg).textContent = "Gene not in the tested set"; return; }
      const box = { l: 56, r: 620, t: 20, b: 320 };
      const vals = Object.values(g.t).flat().filter(v => v !== null);
      const x = scale(0, 4, box.l + 20, box.r - 20);
      const y = scale(...extent(vals.concat([0]), .15), box.b, box.t);
      axes(svg, x, y, box, { noX: true, ylab: "log2 change vs control" });
      ORDER.forEach((c, i) => {
        el("text", { x: x(i), y: box.b + 18, "text-anchor": "middle", fill: COLORS[c] }, svg).textContent = SHORT[c];
        el("rect", { x: x(i) - 44, y: box.t, width: 88, height: box.b - box.t, fill: COLORS[c], opacity: c === "dtp" ? .06 : .025 }, svg);
      });
      el("line", { x1: box.l, x2: box.r, y1: y(0), y2: y(0), stroke: "#3a3757", "stroke-dasharray": "3 4" }, svg);
      for (const l of lines) {
        const t = g.t[l];
        if (!t) continue;
        const d = t.map((v, i) => `${i ? "L" : "M"}${x(i)},${y(v ?? 0)}`).join("");
        const p = el("path", { d, fill: "none", stroke: LINE_COLORS[l], "stroke-width": 2.2 }, svg);
        if (!PT.reduced) {
          const len = p.getTotalLength();
          p.style.strokeDasharray = len; p.style.strokeDashoffset = len;
          p.getBoundingClientRect();
          p.style.transition = "stroke-dashoffset .9s ease"; p.style.strokeDashoffset = 0;
        }
        t.forEach((v, i) => el("circle", { cx: x(i), cy: y(v ?? 0), r: 3.5, fill: LINE_COLORS[l] }, svg));
      }
      el("text", { x: box.l + 8, y: box.t + 14, class: "lbl" }, svg).textContent = name;
      $("geneClass").innerHTML = `<span class="badge">${esc(CLASS_TEXT[g.c] || g.c || "n/a")}</span> <span class="sub">persister z ${fmt(g.P)} · specific z ${fmt(g.S)}</span>`;
    }
    const input = $("geneSearch");
    input.addEventListener("change", () => draw(input.value.trim().toUpperCase() in G ? input.value.trim().toUpperCase() : input.value.trim()));
    const start = G.TACSTD2 ? "TACSTD2" : (passports[0] && passports[0].symbol) || names[0];
    input.value = start;
    draw(start);
    window.PTdrawGene = n => { input.value = n; draw(n); };

    const ms = R.memory.summary;
    const svgM = $("memChart");
    const box = { l: 40, r: 450, t: 10, b: 260 };
    const xs = scale(-1, 2, box.l, box.r);
    const maxH = Math.max(1, ...ms.map(r => Math.max(...r.hist)));
    const y = scale(0, maxH, box.b, box.t);
    axes(svgM, xs, y, box, { xlab: "fraction of persister change remaining" });
    const bw = (box.r - box.l) / 30;
    for (const tp of ["wo_short", "wo_long"]) {
      const sum = Array(30).fill(0);
      ms.filter(r => r.timepoint === tp).forEach(r => r.hist.forEach((v, i) => sum[i] += v));
      const m = Math.max(...sum, 1);
      sum.forEach((v, i) => el("rect", { x: box.l + i * bw, y: y(v / m * maxH), width: bw - 1, height: box.b - y(v / m * maxH), fill: COLORS[tp], opacity: .5 }, svgM));
    }
    el("line", { x1: xs(1), x2: xs(1), y1: box.t, y2: box.b, stroke: COLORS.dtp, "stroke-dasharray": "3 3" }, svgM);
    el("line", { x1: xs(0), x2: xs(0), y1: box.t, y2: box.b, stroke: COLORS.control, "stroke-dasharray": "3 3" }, svgM);
    el("text", { x: xs(1) + 4, y: box.t + 12, fill: COLORS.dtp }, svgM).textContent = "persister level";
    el("text", { x: xs(0) + 4, y: box.t + 12, fill: COLORS.control }, svgM).textContent = "baseline";
    const dec = R.memory.decay;
    const hl = Object.entries(dec).filter(([, d]) => d.status === "fitted" && d.median_half_life_days).map(([l, d]) => `${l} ${fmt(d.median_half_life_days, 1)} d`);
    $("memText").innerHTML = `Median memory 24 h after washout: <b>${med(memShort) !== null ? Math.round(med(memShort) * 100) + "%" : "n/a"}</b>; after regrowth: <b>${med(memLong) !== null ? Math.round(med(memLong) * 100) + "%" : "n/a"}</b> (medians across lines, genes with a significant persister change). ` +
      (hl.length ? `Median memory half-life: ${hl.join(", ")}.` : "The regrowth day is not in the GEO metadata, so no half-life is fitted; the two washout points are still shown.") +
      ` <span style="color:${COLORS.wo_short}">■</span> 24 h <span style="color:${COLORS.wo_long}">■</span> regrowth`;
    const dist = R.reversibility.distribution;
    const keys = Object.keys(dist).filter(k => k !== "not_changed");
    const mx = Math.max(1, ...keys.map(k => dist[k]));
    $("classBars").innerHTML = keys.sort((a, b) => dist[b] - dist[a]).map(k => `<div style="display:grid;grid-template-columns:200px 1fr 50px;gap:8px;align-items:center;font-size:13px;margin:4px 0"><span>${CLASS_TEXT[k] || k}</span><span style="height:9px;border-radius:4px;background:${k.includes("after_regrowth") ? COLORS.dtp : k.includes("acute") ? COLORS.acute : COLORS.wo_short};width:${dist[k] / mx * 100}%"></span><span class="num">${dist[k]}</span></div>`).join("");
  }

  function graph() {
    const svg = $("graphSvg");
    const types = Object.keys(R.graph.edge_types);
    const on = new Set(types);
    const tools = $("graphTools");
    tools.innerHTML = "";
    types.forEach(t => {
      const b = PT.button(t.replace(/_/g, " "), true, () => { on.has(t) ? on.delete(t) : on.add(t); b.classList.toggle("on"); b.setAttribute("aria-pressed", String(on.has(t))); render(); });
      b.title = R.graph.edge_types[t];
      tools.appendChild(b);
    });
    const nodes = R.graph.nodes.map(n => ({ ...n }));
    const idx = Object.fromEntries(nodes.map((n, i) => [n.id, i]));
    const edges = R.graph.edges.filter(e => e.source in idx && e.target in idx);
    const W = 1000, H = 620, cx = W / 2, cy = H / 2;
    const order = { contested: 0, dependency: 1, negative: 2, regulator: 3, marker: 4, stress: 5, weak: 6 };
    const genesN = nodes.filter(n => n.kind === "gene").sort((a, b) => (a.direction > b.direction ? 1 : -1) || (order[a.claim] - order[b.claim]));
    const N = genesN.length, ups = genesN.filter(n => n.direction === "up").length;
    genesN.forEach((n, i) => {
      n.ang = -Math.PI / 2 + (i + .5) / N * Math.PI * 2;
      n.x = cx + Math.cos(n.ang) * 215; n.y = cy + Math.sin(n.ang) * 185;
    });
    const keyFirst = list => { const k = list.filter(n => order[n.claim] < 4), r = list.filter(n => order[n.claim] >= 4), out = []; const step = r.length / (k.length + 1); let j = 0; k.forEach((n, i) => { while (j < Math.round((i + 1) * step)) out.push(r[j++]); out.push(n); }); while (j < r.length) out.push(r[j++]); return out; };
    const arranged = keyFirst(genesN.slice(0, ups)).concat(keyFirst(genesN.slice(ups)));
    arranged.forEach((n, i) => { n.ang = -Math.PI / 2 + (i + .5) / N * Math.PI * 2; n.x = cx + Math.cos(n.ang) * 215; n.y = cy + Math.sin(n.ang) * 185; });
    const nb = {};
    edges.forEach(e => { (nb[e.source] ??= []).push(e.target); (nb[e.target] ??= []).push(e.source); });
    const outer = nodes.filter(n => n.kind !== "gene" && n.kind !== "state");
    outer.forEach(n => {
      const g = (nb[n.id] || []).map(id => nodes[idx[id]]).filter(m => m.kind === "gene");
      const sx = g.reduce((a, m) => a + Math.cos(m.ang), 0), sy = g.reduce((a, m) => a + Math.sin(m.ang), 0);
      n.ang = g.length ? Math.atan2(sy, sx) : Math.random() * 6.28;
    });
    outer.sort((a, b) => a.ang - b.ang);
    for (let it = 0; it < 200; it++) for (let i = 0; i < outer.length; i++) {
      const a = outer[i], b = outer[(i + 1) % outer.length];
      let d = b.ang - a.ang; if (i === outer.length - 1) d += Math.PI * 2;
      const min = .16;
      if (outer.length > 1 && d < min) { a.ang -= (min - d) / 2; b.ang += (min - d) / 2; }
    }
    outer.forEach(n => { n.x = cx + Math.cos(n.ang) * 360; n.y = cy + Math.sin(n.ang) * 282; });
    nodes.filter(n => n.kind === "state").forEach(n => { n.x = cx; n.y = cy; });
    const KC = { state: COLORS.dtp, gene: "#ecebf5", tf: COLORS.acute, program: COLORS.wo_long, source: COLORS.control };
    const ES = { observed_expression: ["#e879f9", ""], inferred_regulatory: ["#fbbf24", "5 4"], pathway_membership: ["#34d399", "2 3"], published_perturbation: ["#7dd3fc", ""], published_regulatory: ["#7dd3fc", "8 3"], published_chromatin: ["#a78bfa", ""] };
    $("graphLegend").innerHTML = Object.entries(KC).map(([k, c]) => `<span><i style="background:${c}"></i>${k}</span>`).join("") + types.map(t => `<span><svg width="26" height="8"><line x1="0" y1="4" x2="26" y2="4" stroke="${ES[t][0]}" stroke-width="2" stroke-dasharray="${ES[t][1]}"/></svg> ${t.replace(/_/g, " ")}</span>`).join("");
    function render() {
      clear(svg);
      const ge = el("g", {}, svg), gn = el("g", {}, svg);
      for (const e of edges) {
        if (!on.has(e.type)) continue;
        const a = nodes[idx[e.source]], b = nodes[idx[e.target]];
        const ln = el("line", { x1: a.x, y1: a.y, x2: b.x, y2: b.y, stroke: ES[e.type][0], "stroke-dasharray": ES[e.type][1], "stroke-width": e.type === "published_perturbation" ? 2 : 1.1, opacity: e.type === "observed_expression" ? .35 : .8 }, ge);
        ln.addEventListener("mousemove", ev => showTip(ev, `<b>${e.type.replace(/_/g, " ")}</b><br>${esc(e.method)}${e.effect ? "<br>effect: " + esc(e.effect) : ""}${e.outcome ? "<br>outcome: " + esc(e.outcome) : ""}${e.fdr !== undefined ? "<br>FDR " + fmt(e.fdr) : ""}`));
        ln.addEventListener("mouseleave", hideTip);
      }
      for (const n of nodes) {
        const r = n.kind === "state" ? 16 : n.kind === "gene" ? 5.5 : 7;
        const c = el("circle", { cx: n.x, cy: n.y, r, fill: KC[n.kind], stroke: n.claim === "contested" ? COLORS.acute : n.claim === "dependency" ? "#34d399" : "none", "stroke-width": 2.5, tabindex: 0 }, gn);
        if (n.kind === "gene") {
          const deg = n.ang * 180 / Math.PI, flip = Math.cos(n.ang) < 0;
          const tx = n.x + Math.cos(n.ang) * 9, ty = n.y + Math.sin(n.ang) * 9;
          el("text", { x: tx, y: ty, "font-size": 9.5, "dominant-baseline": "middle", "text-anchor": flip ? "end" : "start", fill: order[n.claim] < 4 ? "#ecebf5" : "#8d8aa6", "font-weight": order[n.claim] < 4 ? 700 : 400, transform: `rotate(${flip ? deg + 180 : deg} ${tx} ${ty})` }, gn).textContent = n.label;
        } else {
          const left = n.x < cx - 5;
          el("text", { x: n.x + (left ? -r - 4 : r + 4), y: n.y + 4, fill: "#cfcde0", "text-anchor": n.kind === "state" ? "middle" : left ? "end" : "start", dy: n.kind === "state" ? 30 : 0 }, gn).textContent = n.label;
        }
        c.addEventListener("mousemove", ev => showTip(ev, `<b>${esc(n.label)}</b><br>${n.kind}${n.claim ? "<br>" + esc((passports.find(p => p.symbol === n.label) || {}).claim?.label || "") : ""}${n.doi ? "<br>doi " + esc(n.doi) : ""}`));
        c.addEventListener("mouseleave", hideTip);
        const open = () => { if (n.kind === "gene") { openPassport(n.label); location.hash = "passport"; } else if (n.doi) window.open("https://doi.org/" + n.doi, "_blank", "noopener"); };
        c.addEventListener("click", open);
        c.addEventListener("keydown", ev => { if (ev.key === "Enter") open(); });
      }
    }
    render();
  }

  function paradox() {
    const tc = R.programs.time_course || {};
    const progs = Object.keys(tc[lines[0]] || {});
    const svg = $("progChart");
    if (!progs.length) { svg.outerHTML = '<div class="pending">Gene-set download failed in this build, so programme scores are not available.</div>'; return; }
    const pref = ["HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_DNA_REPAIR", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"];
    const choices = pref.filter(p => progs.includes(p)).concat(progs.filter(p => !pref.includes(p))).slice(0, 9);
    let cur = choices[0];
    const nice = p => p.replace("HALLMARK_", "").replace(/_/g, " ").toLowerCase().replace(/^\w/, c => c.toUpperCase());
    group($("progTools"), choices.map(p => [p, nice(p)]), cur, p => { cur = p; draw(); });
    $("progLegend").innerHTML = lines.map(l => `<span><i style="background:${LINE_COLORS[l]}"></i>${l}</span>`).join("");
    function draw() {
      clear(svg);
      const box = { l: 56, r: 620, t: 20, b: 320 };
      const vals = lines.flatMap(l => ORDER.map(c => tc[l]?.[cur]?.[c])).filter(Number.isFinite);
      const x = scale(0, 4, box.l + 20, box.r - 20), y = scale(...extent(vals.concat([0]), .15), box.b, box.t);
      axes(svg, x, y, box, { noX: true, ylab: "mean gene z-score" });
      ORDER.forEach((c, i) => el("text", { x: x(i), y: box.b + 18, "text-anchor": "middle", fill: COLORS[c] }, svg).textContent = SHORT[c]);
      for (const l of lines) {
        const v = ORDER.map(c => tc[l]?.[cur]?.[c]);
        el("path", { d: v.map((q, i) => `${i ? "L" : "M"}${x(i)},${y(q ?? 0)}`).join(""), fill: "none", stroke: LINE_COLORS[l], "stroke-width": 2.2 }, svg);
        v.forEach((q, i) => el("circle", { cx: x(i), cy: y(q ?? 0), r: 3.5, fill: LINE_COLORS[l] }, svg));
      }
      el("text", { x: box.l + 8, y: box.t + 14, class: "lbl" }, svg).textContent = nice(cur);
      const lit = (R.programs.literature || []).filter(p => p.program === cur);
      $("progLit").innerHTML = lit.length ? lit.map(p => `<p><b>Published context.</b> ${esc(p.summary)} <a href="https://doi.org/${R.sources[p.source].doi}" target="_blank" rel="noopener">${esc(R.sources[p.source].journal)} ${R.sources[p.source].year}</a></p>`).join("") : "";
    }
    const rows = R.programs.persister_specific || [];
    $("progTable").innerHTML = `<table class="t"><thead><tr><th>Programme</th><th>score</th><th>FDR</th></tr></thead><tbody>${rows.slice(0, 14).map(r => `<tr><td>${nice(r.source)}</td><td class="num">${fmt(r.score)}</td><td class="num">${fmt(r.fdr)}</td></tr>`).join("")}</tbody></table><p class="sub">Score: t-statistic of a linear model of persister-specific z on programme membership. Positive means the programme rises in persisters beyond acute stress.</p>`;
    draw();
  }

  const DIM_TITLE = { persistence_enrichment: "Persister enrichment", acute_stress_specificity: "Beyond acute stress", reversibility: "Reversibility", regulatory_support: "Regulatory support", perturbation_support: "Perturbation support", reproducibility: "Reproducibility" };
  let current = null;

  function openPassport(sym) {
    const p = passports.find(q => q.symbol === sym);
    if (!p) return;
    current = sym;
    document.querySelectorAll("#plist button").forEach(b => b.classList.toggle("on", b.dataset.sym === sym));
    const d = p.dimensions;
    const slug = s => String(s).replace(/\s+/g, "-");
    const dimHtml = Object.entries(d).map(([k, v]) => {
      let body = "";
      if (k === "persistence_enrichment") body = `meta log2FC ${fmt(v.meta_lfc)}, FDR ${fmt(v.fdr)}, ${v.lines_agreeing} of ${lines.length} lines`;
      if (k === "acute_stress_specificity") body = `specific part ${fmt(v.persister_specific_lfc)} vs acute part ${fmt(v.acute_explained_lfc)}`;
      if (k === "reversibility") body = esc(CLASS_TEXT[v.class] || v.class);
      if (k === "regulatory_support") body = (v.regulators || []).slice(0, 3).map(r => `${r.tf} (${r.weight > 0 ? "activates" : "represses"}, FDR ${fmt(r.tf_fdr)})`).join(", ") || (v.literature.length ? "published record" : "none found");
      if (k === "perturbation_support") body = v.records.length ? `${v.records.length} published record${v.records.length > 1 ? "s" : ""}` : "no published perturbation found";
      if (k === "reproducibility") body = `${v.lines} of ${v.of} lines, I² ${fmt(v.I2)}`;
      return `<div class="dim"><h4>${DIM_TITLE[k]}</h4><span class="badge ${slug(v.status)}">${esc(v.status)}</span><div style="margin-top:6px">${body}</div></div>`;
    }).join("");
    const recs = d.perturbation_support.records.concat(d.regulatory_support.literature || []);
    const recHtml = recs.length ? `<h3 style="margin-top:14px">Published evidence</h3><ul class="records">${recs.map(r => { const s = R.sources[r.source] || {}; return `<li><b>${esc(r.kind.replace(/_/g, " "))}</b> · ${esc(r.model || "")} · outcome: ${esc(r.outcome || "")} · effect: <b>${esc(r.effect || "")}</b><br>${esc(r.summary || "")} <a href="https://doi.org/${s.doi}" target="_blank" rel="noopener">${esc(s.journal)} ${s.year}</a></li>`; }).join("")}</ul>` : "";
    const contra = p.contradictions.length ? `<h3>Contradictions and caveats</h3><ul class="contra">${p.contradictions.map(c => `<li>${esc(c)}</li>`).join("")}</ul>` : "";
    const perLine = Object.entries(p.per_line).map(([l, v]) => `<tr><td>${l}</td><td class="num">${fmt(v.A)}</td><td class="num">${fmt(v.P)}</td><td class="num">${fmt(v.S)}</td><td>${esc(CLASS_TEXT[v.class] || v.class || "")}</td></tr>`).join("");
    $("pcard").innerHTML = `<div style="display:flex;justify-content:space-between;align-items:baseline;gap:10px;flex-wrap:wrap"><h3 style="font-size:24px;margin:0">${p.symbol}</h3><button class="btn small" id="toRev">See its trajectory</button></div>
      <div class="verdict"><span class="sub">Strongest claim the evidence allows</span><br><b>${esc(p.claim.label)}</b><div class="sub">${esc(p.claim.meaning)}</div></div>
      <div class="dims">${dimHtml}</div>${contra}${recHtml}
      <h3 style="margin-top:14px">Per cell line (log2FC)</h3><table class="t"><thead><tr><th>Line</th><th>A acute</th><th>P persister</th><th>S specific</th><th>Withdrawal</th></tr></thead><tbody>${perLine}</tbody></table>
      <h3 style="margin-top:14px">Evidence-to-claim trace</h3><div class="trace"><span>${esc(p.claim.label)}</span><span class="arr">←</span><span>${esc(p.trace.contrasts.join(" · "))}</span><span class="arr">←</span><span>${esc(p.trace.test)}</span><span class="arr">←</span><span>${p.trace.samples.length} samples</span><span class="arr">←</span><span><a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${p.trace.dataset}" target="_blank" rel="noopener">${p.trace.dataset}</a></span><span class="arr">←</span><span><a href="https://doi.org/${p.trace.publication.doi}" target="_blank" rel="noopener">${esc(p.trace.publication.journal)} ${p.trace.publication.year}</a></span></div>`;
    $("toRev").addEventListener("click", () => { window.PTdrawGene(p.symbol); location.hash = "reversibility"; });
  }

  function passportPanel() {
    const labels = [...new Set(passports.map(p => p.claim.label))];
    $("claimFilter").innerHTML += labels.map(l => `<option>${esc(l)}</option>`).join("");
    const rank = { contested: 0, dependency: 1, negative: 2, regulator: 3, marker: 4, stress: 5, weak: 6 };
    const sorted = passports.slice().sort((a, b) => rank[a.claim.key] - rank[b.claim.key] || Math.abs(b.dimensions.persistence_enrichment.meta_lfc || 0) - Math.abs(a.dimensions.persistence_enrichment.meta_lfc || 0));
    function list(filter) {
      $("plist").innerHTML = sorted.filter(p => !filter || p.claim.label === filter).map(p => `<button role="option" data-sym="${p.symbol}"><span>${p.symbol} <small style="color:var(--dim)">${p.direction === "up" ? "▲" : "▼"}</small></span><span class="badge ${p.claim.key}">${esc(p.claim.label)}</span></button>`).join("");
      $("plist").querySelectorAll("button").forEach(b => b.addEventListener("click", () => openPassport(b.dataset.sym)));
    }
    $("claimFilter").addEventListener("change", e => list(e.target.value));
    list("");
    $("trapBtn").addEventListener("click", () => openPassport(passports.some(p => p.symbol === "TACSTD2") ? "TACSTD2" : sorted[0].symbol));
    document.addEventListener("click", e => { const a = e.target.closest("[data-gene]"); if (a) openPassport(a.dataset.gene); });
    openPassport(passports.some(p => p.symbol === "TACSTD2") ? "TACSTD2" : sorted[0]?.symbol);
  }

  function arena() {
    const slug = s => String(s).replace(/\s+/g, "-");
    $("hyps").innerHTML = R.hypotheses.map(h => {
      let detail = "";
      if (h.naive !== undefined) detail = `naive ${fmt(h.naive, 3)} vs persister-specific ${fmt(h.decomposed, 3)}`;
      if (h.classes) detail = Object.entries(h.classes).filter(([k]) => k !== "not_changed").map(([k, v]) => `${CLASS_TEXT[k] || k}: ${v}`).join(" · ");
      if (h.value !== undefined && h.value !== null) detail = `${h.value} regulators`;
      if (h.genes) detail = h.genes.join(", ") || "none";
      if (h.ranks) detail = h.ranks.map(r => `${r.gene} ${fmt(r.naive_percentile * 100, 0)}→${fmt(r.persister_specific_percentile * 100, 0)}%`).join(" · ");
      return `<div class="hyp"><span class="id">${h.id}</span><span>${esc(h.text)}<br><span class="sub">${esc(h.metric || "")}${detail ? ": " + esc(detail) : ""}${h.note ? " " + esc(h.note) : ""}</span></span><span class="res ${slug(h.result)}">${esc(h.result)}</span></div>`;
    }).join("");
    const svg = $("foldChart");
    const box = { l: 46, r: 450, t: 14, b: 240 };
    const folds = R.benchmark.folds;
    const x = scale(-.5, folds.length - .5, box.l, box.r), y = scale(0, 1, box.b, box.t);
    axes(svg, x, y, box, { noX: true, ylab: "AUROC" });
    el("line", { x1: box.l, x2: box.r, y1: y(.5), y2: y(.5), stroke: "#3a3757", "stroke-dasharray": "3 3" }, svg);
    folds.forEach((f, i) => {
      const a = f["naive_treated_vs_control|auroc_dtp_vs_acute"], b = f["persister_specific|auroc_dtp_vs_acute"];
      el("rect", { x: x(i) - 30, y: y(a ?? 0), width: 28, height: box.b - y(a ?? 0), fill: COLORS.acute }, svg);
      el("rect", { x: x(i) + 2, y: y(b ?? 0), width: 28, height: box.b - y(b ?? 0), fill: COLORS.dtp }, svg);
      el("text", { x: x(i), y: box.b + 16, "text-anchor": "middle" }, svg).textContent = "held out " + f.held_out;
    });
    const rows = R.benchmark.summary;
    $("benchTable").innerHTML = `<table class="t"><thead><tr><th>Features</th><th>Model</th><th>Balanced acc.</th><th>95% CI</th><th>Macro F1</th><th>Persister/acute confusion</th></tr></thead><tbody>${rows.map(r => `<tr><td>${r.features.replace(/_/g, " ")}</td><td>${r.model.replace(/_/g, " ")}</td><td class="num">${fmt(r.balanced_accuracy, 3)}</td><td class="num">${fmt(r.ci_low, 2)} to ${fmt(r.ci_high, 2)}</td><td class="num">${fmt(r.macro_f1, 3)}</td><td class="num">${fmt(r.dtp_acute_confusion, 3)}</td></tr>`).join("")}</tbody></table><p class="sub">Chance is 0.20 for five balanced classes. Intervals come from resampling the 60 held-out predictions and are optimistic, because only four lines exist.</p>`;
  }

  function mathModel() {
    const P = { on: { r_S: -0.8, r_P: 0, a: 0.05, b: 0 }, off: { r_S: 0.6, r_P: 0.1, a: 0.001, b: 0.3 } };
    const dec = Object.values(R.memory.decay).filter(d => d.median_half_life_days);
    const kData = dec.length ? Math.LN2 / med(dec.map(d => d.median_half_life_days)) : null;
    const sliders = [
      ["on", "r_S", "Death rate of sensitive cells under drug", -2, 0, .05],
      ["on", "a", "Switch into persister state under drug (a)", 0, .3, .005],
      ["off", "b", "Switch back after washout (b)", 0, 1.5, .01],
      ["off", "r_S", "Regrowth rate after washout", 0, 1.2, .05],
    ];
    const box = $("modelControls");
    box.innerHTML = sliders.map(([ph, k, label, lo, hi, st], i) => `<label style="display:block;margin-bottom:10px;font-size:13px">${label}: <b id="mv${i}">${P[ph][k]}</b><input type="range" min="${lo}" max="${hi}" step="${st}" value="${P[ph][k]}" id="ms${i}"></label>`).join("") +
      `<div id="eig" class="sub"></div>` + (kData ? `<button class="btn small" id="useData" style="margin-top:8px">Set b from the measured memory decay (${fmt(kData, 2)} per day)</button><p class="sub">This treats transcriptional memory decay as a stand-in for the switch-back rate. That is an assumption, stated here so nobody mistakes it for a measurement.</p>` : "");
    sliders.forEach(([ph, k], i) => $("ms" + i).addEventListener("input", e => { P[ph][k] = +e.target.value; $("mv" + i).textContent = e.target.value; draw(); }));
    if (kData) $("useData").addEventListener("click", () => { P.off.b = +kData.toFixed(3); $("ms2").value = P.off.b; $("mv2").textContent = P.off.b; draw(); });
    const lam = p => { const a = p.r_S - p.a, d = p.r_P - p.b, tr = a + d, det = a * d - p.b * p.a; return (tr + Math.sqrt(Math.max(tr * tr - 4 * det, 0))) / 2; };
    const svg = $("modelChart");
    function draw() {
      let s = .999, q = .001;
      const out = [[0, s, q]];
      const dt = .02;
      for (let t = dt; t <= 42 + 1e-9; t += dt) {
        const p = t <= 21 ? P.on : P.off;
        const ds = (p.r_S - p.a) * s + p.b * q, dq = p.a * s + (p.r_P - p.b) * q;
        s += ds * dt; q += dq * dt;
        if (Math.abs(t * 4 - Math.round(t * 4)) < 1e-6) out.push([t, s, q]);
      }
      clear(svg);
      const bx = { l: 56, r: 620, t: 14, b: 310 };
      const vals = out.flatMap(o => [o[1], o[2], o[1] + o[2]]).filter(v => v > 0).map(Math.log10);
      const x = scale(0, 42, bx.l, bx.r), y = scale(Math.max(-8, Math.min(...vals)), Math.max(...vals) + .3, bx.b, bx.t);
      axes(svg, x, y, bx, { xlab: "days", ylab: "log10 cell fraction" });
      el("rect", { x: x(0), y: bx.t, width: x(21) - x(0), height: bx.b - bx.t, fill: COLORS.dtp, opacity: .05 }, svg);
      el("text", { x: x(10.5), y: bx.t + 14, "text-anchor": "middle", fill: COLORS.dtp }, svg).textContent = "drug on";
      el("text", { x: x(31.5), y: bx.t + 14, "text-anchor": "middle", fill: COLORS.wo_long }, svg).textContent = "drug off";
      const line = (j, c, w) => el("path", { d: out.map((o, i) => { const v = j === 3 ? o[1] + o[2] : o[j]; return `${i ? "L" : "M"}${x(o[0])},${y(Math.log10(Math.max(v, 1e-8)))}`; }).join(""), fill: "none", stroke: c, "stroke-width": w }, svg);
      line(1, COLORS.control, 2); line(2, COLORS.dtp, 2); line(3, "#ecebf5", 1.4);
      const lo = lam(P.on), lf = lam(P.off);
      const end21 = out.find(o => o[0] >= 21);
      $("eig").innerHTML = `λ<sub>max</sub> under drug: <b>${fmt(lo, 3)}</b> per day (${lo < -1e-3 ? "population shrinks" : lo > 1e-3 ? "population grows" : "population holds steady"}). After washout: <b>${fmt(lf, 3)}</b> per day. Persister share at day 21: <b>${fmt(end21[2] / (end21[1] + end21[2]) * 100, 1)}%</b>.`;
    }
    draw();
  }

  function methods() {
    $("gates").innerHTML = Object.entries(R.gates).map(([k, v]) => `<div class="gate"><span>${k.replace(/_/g, " ").replace(/^(\w) /, "$1 · ")}</span><span class="${v === "passed" ? "passed" : v === "failed" ? "failed" : ""}">${esc(v)}</span></div>`).join("");
    const rb = R.robustness;
    $("robust").innerHTML = `Withdrawal classes across ${rb.threshold_grid.length} threshold settings: <b>${esc(rb.labels.threshold_agreement)}</b> (lowest agreement ${fmt(Math.min(...rb.threshold_grid.map(r => r.agreement_with_reference)), 2)}). Signature when one line is dropped: <b>${esc(rb.labels.signature_lolo)}</b> (Jaccard ${rb.leave_one_line_out.map(r => r.dropped + " " + fmt(r.jaccard_top, 2)).join(", ")}).`;
    $("manifest").innerHTML = R.manifest.length ? `<table class="t"><thead><tr><th>File</th><th>SHA-256</th><th>Size</th></tr></thead><tbody>${R.manifest.map(m => `<tr><td><a href="${esc(m.url)}" target="_blank" rel="noopener">${esc(m.path.split("/").pop())}</a></td><td style="font-family:var(--mono);font-size:11px">${m.sha256.slice(0, 16)}…</td><td class="num">${(m.bytes / 1e6).toFixed(1)} MB</td></tr>`).join("")}</tbody></table>` : "Synthetic run: no downloaded files.";
    $("runinfo").innerHTML = `Generated ${esc(R.meta.generated_utc)} · Python ${esc(R.meta.python)} · pipeline ${esc(R.meta.version)}<br>Hallmark sets: ${esc(R.meta.networks.hallmark)}<br>Regulons: ${esc(R.meta.networks.collectri)}`;
    $("refs").innerHTML = Object.values(R.sources).map(s => `<div style="margin-bottom:6px">${esc(s.title)}. <i>${esc(s.journal)}</i> ${s.year}. <a href="https://doi.org/${s.doi}" target="_blank" rel="noopener">doi</a></div>`).join("");
  }

  function reveal() {
    const io = new IntersectionObserver(es => es.forEach(e => { if (e.isIntersecting) { e.target.classList.add("in"); io.unobserve(e.target); } }), { threshold: .08 });
    document.querySelectorAll(".reveal").forEach(s => io.observe(s));
    const links = [...document.querySelectorAll("nav.sections a")];
    const so = new IntersectionObserver(es => es.forEach(e => { if (e.isIntersecting) links.forEach(a => a.classList.toggle("on", a.getAttribute("href") === "#" + e.target.id)); }), { rootMargin: "-45% 0px -50% 0px" });
    document.querySelectorAll("section.mod").forEach(s => so.observe(s));
  }

  const steps = [heroStats, tiles, () => Hero.init(R), video, universe, decomposition, reversibility, graph, paradox, passportPanel, arena, mathModel, methods, reveal];
  for (const f of steps) {
    try { f(); } catch (err) { console.error(f.name, err); }
  }
})();
