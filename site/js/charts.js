const PT = (() => {
  const NS = "http://www.w3.org/2000/svg";
  const COLORS = { control: "#7dd3fc", acute: "#fbbf24", dtp: "#e879f9", wo_short: "#a78bfa", wo_long: "#34d399" };
  const ORDER = ["control", "acute", "dtp", "wo_short", "wo_long"];
  const SHORT = { control: "Control", acute: "Acute 24 h", dtp: "Persister 21 d", wo_short: "Washout 24 h", wo_long: "Regrowth" };
  const LINE_COLORS = { PC9: "#7dd3fc", H1975: "#fbbf24", HCC827: "#e879f9", HCC2935: "#34d399" };
  const reduced = window.matchMedia("(prefers-reduced-motion: reduce)").matches;

  function el(tag, attrs = {}, parent) {
    const n = document.createElementNS(NS, tag);
    for (const [k, v] of Object.entries(attrs)) n.setAttribute(k, v);
    if (parent) parent.appendChild(n);
    return n;
  }

  function clear(svg) { while (svg.firstChild) svg.removeChild(svg.firstChild); }

  function scale(d0, d1, r0, r1) {
    const f = v => r0 + (v - d0) / ((d1 - d0) || 1) * (r1 - r0);
    f.domain = [d0, d1];
    return f;
  }

  function extent(vals, pad = 0.06) {
    const v = vals.filter(Number.isFinite);
    let lo = Math.min(...v), hi = Math.max(...v);
    if (!v.length) { lo = 0; hi = 1; }
    const p = (hi - lo || 1) * pad;
    return [lo - p, hi + p];
  }

  function ticks(lo, hi, n = 5) {
    const step = Math.pow(10, Math.floor(Math.log10((hi - lo) / n || 1)));
    const err = (hi - lo) / n / step;
    const m = err >= 7.5 ? 10 : err >= 3.5 ? 5 : err >= 1.5 ? 2 : 1;
    const s = m * step, out = [];
    for (let t = Math.ceil(lo / s) * s; t <= hi + 1e-9; t += s) out.push(+t.toFixed(10));
    return out;
  }

  function axes(svg, x, y, box, opts = {}) {
    const g = el("g", { class: "axis" }, svg);
    el("line", { x1: box.l, x2: box.r, y1: box.b, y2: box.b }, g);
    el("line", { x1: box.l, x2: box.l, y1: box.t, y2: box.b }, g);
    if (!opts.noX) for (const t of ticks(...x.domain)) {
      el("line", { x1: x(t), x2: x(t), y1: box.b, y2: box.b + 4 }, g);
      el("text", { x: x(t), y: box.b + 16, "text-anchor": "middle" }, g).textContent = fmt(t);
    }
    for (const t of ticks(...y.domain)) {
      el("line", { x1: box.l - 4, x2: box.l, y1: y(t), y2: y(t) }, g);
      el("text", { x: box.l - 7, y: y(t) + 4, "text-anchor": "end" }, g).textContent = fmt(t);
    }
    if (opts.xlab) el("text", { x: (box.l + box.r) / 2, y: box.b + 34, "text-anchor": "middle", class: "lbl" }, g).textContent = opts.xlab;
    if (opts.ylab) el("text", { x: -(box.t + box.b) / 2, y: box.l - 40, transform: "rotate(-90)", "text-anchor": "middle", class: "lbl" }, g).textContent = opts.ylab;
    return g;
  }

  function fmt(v, d) {
    if (v === null || v === undefined || !Number.isFinite(+v)) return "n/a";
    v = +v;
    if (d !== undefined) return v.toFixed(d);
    const a = Math.abs(v);
    if (Number.isInteger(v)) return String(v);
    if (a !== 0 && (a < 0.001 || a >= 1e5)) return v.toExponential(1);
    return a >= 100 ? v.toFixed(0) : a >= 10 ? v.toFixed(1) : +v.toFixed(2) + "";
  }

  const tip = () => document.getElementById("tip");
  function showTip(evt, html) {
    const t = tip();
    t.innerHTML = html;
    t.classList.add("show");
    const x = Math.min(evt.clientX + 14, window.innerWidth - 300);
    t.style.left = x + "px";
    t.style.top = (evt.clientY + 14) + "px";
  }
  function hideTip() { tip().classList.remove("show"); }

  function button(label, on, onClick) {
    const b = document.createElement("button");
    b.className = "btn small" + (on ? " on" : "");
    b.type = "button";
    b.textContent = label;
    b.setAttribute("aria-pressed", on ? "true" : "false");
    b.addEventListener("click", onClick);
    return b;
  }

  function group(container, labels, active, onPick) {
    container.innerHTML = "";
    const btns = labels.map(([key, text]) => {
      const b = button(text, key === active, () => {
        btns.forEach(x => { x.classList.remove("on"); x.setAttribute("aria-pressed", "false"); });
        b.classList.add("on");
        b.setAttribute("aria-pressed", "true");
        onPick(key);
      });
      container.appendChild(b);
      return b;
    });
    return btns;
  }

  function esc(s) { return String(s ?? "").replace(/[&<>"]/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c])); }

  return { el, clear, scale, extent, ticks, axes, fmt, showTip, hideTip, button, group, esc, COLORS, ORDER, SHORT, LINE_COLORS, reduced };
})();
