const Hero = (() => {
  let ctx, canvas, W, H, dpr, pts = [], cents = {}, box, playing = !PT.reduced, t0 = 0, last = 0, tElapsed = 0;
  const STEP = 2600, HOLD = 1600;
  const order = PT.ORDER;

  function init(res) {
    canvas = document.getElementById("heroCanvas");
    ctx = canvas.getContext("2d");
    resize();
    window.addEventListener("resize", resize);
    const xs = res.samples.map(s => s.pc[0]), ys = res.samples.map(s => s.pc[1]);
    box = { x: PT.extent(xs, .12), y: PT.extent(ys, .12) };
    const byKey = {};
    res.samples.forEach((s, i) => {
      const seed = Math.sin(i * 12.9898) * 43758.5453;
      pts.push({ s, jx: (seed - Math.floor(seed) - .5), jy: (Math.cos(i * 7.31) * .5), born: order.indexOf(s.condition) });
      (byKey[s.line] ??= {})[s.condition] ??= [];
      byKey[s.line][s.condition].push(s.pc);
    });
    for (const [l, d] of Object.entries(byKey)) {
      cents[l] = {};
      for (const [c, arr] of Object.entries(d)) cents[l][c] = [mean(arr.map(a => a[0])), mean(arr.map(a => a[1]))];
    }
    const btn = document.getElementById("heroPlay");
    const sync = () => { btn.textContent = playing ? "pause" : "play"; btn.setAttribute("aria-pressed", String(playing)); };
    btn.addEventListener("click", () => { playing = !playing; sync(); if (playing) requestAnimationFrame(frame); });
    sync();
    document.getElementById("heroCap").textContent = `${res.meta.dataset} · ${res.samples.length} real RNA-seq samples · each dot is one experiment`;
    if (!playing) tElapsed = (order.length - 1) * STEP;
    requestAnimationFrame(frame);
  }

  function mean(a) { return a.reduce((x, y) => x + y, 0) / a.length; }

  function resize() {
    dpr = window.devicePixelRatio || 1;
    const r = canvas.getBoundingClientRect();
    W = r.width; H = r.height;
    canvas.width = W * dpr; canvas.height = H * dpr;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    if (!playing) draw();
  }

  function X(v) { return 30 + (v - box.x[0]) / (box.x[1] - box.x[0]) * (W - 130); }
  function Y(v) { return H - 46 - (v - box.y[0]) / (box.y[1] - box.y[0]) * (H - 100); }
  const ease = u => u * u * (3 - 2 * u);

  function frame(ts) {
    if (!last) last = ts;
    if (playing) tElapsed += Math.min(ts - last, 50);
    last = ts;
    draw();
    if (playing) requestAnimationFrame(frame);
  }

  function draw() {
    if (!ctx || !box) return;
    const total = (order.length - 1) * STEP + HOLD;
    const t = tElapsed % total;
    const pos = Math.min(t / STEP, order.length - 1);
    const k = Math.min(Math.floor(pos), order.length - 2);
    const u = ease(Math.min(pos - k, 1));
    const current = order[pos >= order.length - 1 ? order.length - 1 : (u > .5 ? k + 1 : k)];
    ctx.clearRect(0, 0, W, H);
    for (const p of pts) {
      const vis = p.born <= k ? 1 : p.born === k + 1 ? u : 0;
      if (vis <= 0) continue;
      const pulse = p.s.condition === current ? 1 : .55;
      ctx.globalAlpha = vis * pulse;
      ctx.fillStyle = PT.NIGHT[p.s.condition];
      ctx.beginPath();
      ctx.arc(X(p.s.pc[0]), Y(p.s.pc[1]), p.s.condition === current ? 5.2 : 4, 0, Math.PI * 2);
      ctx.fill();
    }
    ctx.globalAlpha = 1;
    const heads = [];
    Object.entries(cents).forEach(([l, c]) => {
      const path = order.slice(0, k + 1).map(o => c[o]).filter(Boolean);
      const a = c[order[k]], b = c[order[k + 1]];
      if (a && b) path.push([a[0] + u * (b[0] - a[0]), a[1] + u * (b[1] - a[1])]);
      ctx.strokeStyle = "rgba(239,232,218,.5)";
      ctx.lineWidth = 1.4;
      ctx.beginPath();
      path.forEach((q, i) => i ? ctx.lineTo(X(q[0]), Y(q[1])) : ctx.moveTo(X(q[0]), Y(q[1])));
      ctx.stroke();
      const h = path[path.length - 1];
      ctx.fillStyle = PT.NIGHT[current];
      ctx.strokeStyle = "#efe8da";
      ctx.lineWidth = 1.5;
      ctx.beginPath();
      ctx.arc(X(h[0]), Y(h[1]), 8, 0, Math.PI * 2);
      ctx.fill();
      ctx.stroke();
      heads.push({ l, x: X(h[0]), y: Y(h[1]), ly: Y(h[1]) });
    });
    heads.sort((p, q) => p.ly - q.ly);
    for (let i = 1; i < heads.length; i++) if (heads[i].ly - heads[i - 1].ly < 19) heads[i].ly = heads[i - 1].ly + 19;
    ctx.font = "italic 16px Garamond, serif";
    for (const hd of heads) {
      const lx = hd.x + 26;
      ctx.strokeStyle = "rgba(239,232,218,.35)";
      ctx.lineWidth = 1;
      ctx.beginPath(); ctx.moveTo(hd.x + 9, hd.y); ctx.lineTo(lx - 3, hd.ly - 4); ctx.stroke();
      ctx.fillStyle = "#efe8da";
      ctx.fillText(hd.l, lx, hd.ly);
    }
    const ph = document.getElementById("heroPhase");
    ph.textContent = PT.SHORT[current];
    ph.style.color = PT.NIGHT[current];
  }

  return { init };
})();
