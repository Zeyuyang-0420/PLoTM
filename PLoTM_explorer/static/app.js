"use strict";
let DATA = null, DID = null, COND = "rest", SORT = {k: "meta_FDR", asc: true}, TRAIT = "";
let NET = {group: null, focus: "", annot: false, sig: false};
let REG = {shared: false, prog: "", search: ""}, REGSORT = {k: "n_clean_traits", asc: false};
let netDrag = null;
const $ = s => document.querySelector(s);
const $$ = s => Array.from(document.querySelectorAll(s));

// ---------- formatting ----------
const isNum = x => typeof x === "number" && !isNaN(x);
const fP = x => !isNum(x) ? "—" : (x < 0.001 ? x.toExponential(1) : x.toFixed(3));
const fB = x => !isNum(x) ? "—" : (x >= 0 ? "+" : "") + x.toFixed(3);
const fI2 = x => !isNum(x) ? "—" : Math.round(x) + "%";
function el(tag, attrs, kids){
  const e = document.createElement(tag);
  if (attrs) for (const k in attrs){
    if (k === "class") e.className = attrs[k];
    else if (k === "html") e.innerHTML = attrs[k];
    else if (k.startsWith("on")) e.addEventListener(k.slice(2), attrs[k]);
    else if (attrs[k] != null) e.setAttribute(k, attrs[k]);
  }
  (kids || []).forEach(c => e.appendChild(typeof c === "string" ? document.createTextNode(c) : c));
  return e;
}
const dz = () => DATA.diseases[DID];
const progs = () => dz().programs[COND];

// ---------- init ----------
async function init(){
  const r = await fetch("/api/data");
  if (!r.ok){ $("#subtitle").textContent = "run build_demo_data.py first"; return; }
  DATA = await r.json();
  const ids = Object.keys(DATA.diseases);
  DID = ids[0];
  const sel = $("#diseaseSel");
  ids.forEach(id => sel.appendChild(el("option", {value: id}, [DATA.diseases[id].name])));
  sel.addEventListener("change", () => { DID = sel.value; onDisease(); });
  $$("#condSeg button").forEach(b => b.addEventListener("click", () => {
    COND = b.dataset.cond; $$("#condSeg button").forEach(x => x.classList.toggle("on", x === b));
    renderAll();
    netPost(); // drive the network iframe's condition from the main explorer
  }));
  // the network iframe is controlled entirely from here: push condition + theme when it loads,
  // and reply whenever it announces it is ready
  const nf = $("#netframe");
  if (nf) nf.addEventListener("load", netPost);
  window.addEventListener("message", e => { if (e.data && e.data.type === "proglof-net-ready") netPost(); });
  $$(".tabs button").forEach(b => b.addEventListener("click", () => {
    $$(".tabs button").forEach(x => x.classList.toggle("on", x === b));
    $$(".tabpane").forEach(p => p.classList.toggle("on", p.id === "tab-" + b.dataset.tab));
    if (b.dataset.tab === "network") renderNetwork();   // needs visible container to size the SVG
  }));
  $("#fAnnot").addEventListener("change", renderTable);
  $("#fSig").addEventListener("change", renderTable);
  $("#traitSel").addEventListener("change", () => { TRAIT = $("#traitSel").value; renderTable(); });
  $$("#progTable thead th").forEach(th => { if (th.dataset.k) th.addEventListener("click", () => {
    if (SORT.k === th.dataset.k) SORT.asc = !SORT.asc; else { SORT.k = th.dataset.k; SORT.asc = true; }
    renderTable();
  });});
  // network controls (legacy vanilla network; the tab now hosts the React app via an iframe, so these
  // elements are absent — guard so init still runs)
  if ($("#netFocus")) {
    $("#netFocus").addEventListener("change", () => { NET.focus = $("#netFocus").value.trim(); renderNetwork(); });
    $("#netFocus").addEventListener("input", () => { if ($("#netFocus").value === "") { NET.focus = ""; renderNetwork(); } });
    $("#netClear").addEventListener("click", () => { NET.focus = ""; $("#netFocus").value = ""; renderNetwork(); });
    $("#netAnnot").addEventListener("change", () => { NET.annot = $("#netAnnot").checked; renderNetwork(); });
    $("#netSig").addEventListener("change", () => { NET.sig = $("#netSig").checked; renderNetwork(); });
  }
  // regulator controls
  $("#regShared").addEventListener("change", () => { REG.shared = $("#regShared").checked; renderRegulators(); });
  $("#regProg").addEventListener("change", () => { REG.prog = $("#regProg").value; renderRegulators(); });
  $("#regSearch").addEventListener("input", () => { REG.search = $("#regSearch").value.trim().toUpperCase(); renderRegulators(); });
  $$("#regTable thead th").forEach(th => { if (th.dataset.k) th.addEventListener("click", () => {
    if (REGSORT.k === th.dataset.k) REGSORT.asc = !REGSORT.asc; else { REGSORT.k = th.dataset.k; REGSORT.asc = true; }
    renderRegulators();
  });});
  window.addEventListener("resize", () => { if ($("#tab-network").classList.contains("on")) renderNetwork(); });
  initNetInteract();
  initTheme();
  $("#chatForm").addEventListener("submit", onSend);
  $$(".examples li").forEach(li => li.addEventListener("click", () => {
    $("#chatInput").value = li.textContent; $("#chatInput").focus();
  }));
  $("#drawer").addEventListener("click", e => { if (e.target.classList.contains("close")) closeDrawer(); });
  $("#lightbox").addEventListener("click", () => $("#lightbox").classList.add("hidden"));
  onDisease();
}

function onDisease(){
  const d = dz();
  $("#subtitle").textContent = `${d.name} · ${d.cell_type}`;
  const ts = $("#traitSel");
  ts.innerHTML = "<option value=''>— cross-source meta —</option>";
  d.traits.forEach(t => ts.appendChild(el("option", {value: t.short}, [`${t.short} (${t.source})`])));
  // network group buttons (legacy vanilla network; absent when the iframe hosts the React app)
  const groups = Object.keys((d.networks && d.networks.rest) || {});
  NET.group = groups[0] || null;
  const gseg = $("#netGroup");
  if (gseg) {
    gseg.innerHTML = "";
    groups.forEach((g, i) => gseg.appendChild(el("button", {class: i === 0 ? "on" : "",
      onclick: () => { NET.group = g; NET.focus = ""; $("#netFocus").value = "";
        $$("#netGroup button").forEach(b => b.classList.toggle("on", b.textContent === netLabel(g))); renderNetwork(); }},
      [netLabel(g)])));
  }
  // regulator program filter
  const rp = $("#regProg"); rp.innerHTML = "<option value=''>— any —</option>";
  const progsWithReg = Array.from(new Set((d.regulators || []).flatMap(r => r.targets.map(t => t.program))))
    .sort((a, b) => (+a.slice(1)) - (+b.slice(1)));
  progsWithReg.forEach(p => rp.appendChild(el("option", {value: p}, [p])));
  renderAll();
}
const netLabel = g => (dz().networks.rest[g] || {}).label || g;
function renderAll(){ renderTable(); renderTraits(); renderRegulators(); renderNotes();
  if ($("#tab-network").classList.contains("on")) renderNetwork(); }

// ---------- theme ----------
function initTheme(){
  const saved = localStorage.getItem("proglof-theme") || "dark";
  applyTheme(saved);
  const btn = $("#themeBtn");
  if (btn) btn.addEventListener("click", () => {
    applyTheme(document.documentElement.dataset.theme === "light" ? "dark" : "light");
    if ($("#tab-network").classList.contains("on")) renderNetwork();
  });
}
function applyTheme(t){
  document.documentElement.dataset.theme = t;
  localStorage.setItem("proglof-theme", t);
  const btn = $("#themeBtn"); if (btn) btn.textContent = t === "light" ? "🌙" : "☀️";
  netPost(); // keep the network iframe's theme in sync
}
// push the explorer's condition + theme into the network iframe (single source of truth)
function netPost(){
  const f = $("#netframe");
  if (f && f.contentWindow) f.contentWindow.postMessage(
    { type: "proglof-net", condition: COND, theme: document.documentElement.dataset.theme }, "*");
}

// ---------- program table ----------
function rowsForTable(){
  let rows = progs().slice();
  if ($("#fAnnot").checked) rows = rows.filter(p => p.annotation);
  if ($("#fSig").checked) rows = rows.filter(p => p.meta && isNum(p.meta.FDR) && p.meta.FDR < 0.05);
  const val = p => {
    const m = p.meta || {};
    switch (SORT.k){
      case "num": return p.num;
      case "annotation": return p.annotation ? 0 : 1;   // annotated first
      case "meta_pooled": return isNum(m.pooled) ? -Math.abs(m.pooled) : 1e9;
      case "meta_FDR": return isNum(m.FDR) ? m.FDR : 1e9;
      case "meta_I2": return isNum(m.I2) ? m.I2 : 1e9;
      case "best_program_P": return isNum(traitP(p)) ? traitP(p) : 1e9;
      case "xsrc": return m.cross_source_concordant ? 0 : 1;
      default: return p.num;
    }
  };
  rows.sort((a, b) => { const x = val(a), y = val(b); return SORT.asc ? (x > y ? 1 : x < y ? -1 : a.num - b.num) : (x < y ? 1 : x > y ? -1 : a.num - b.num); });
  return rows;
}
const traitP = p => TRAIT ? (p.traits[TRAIT] ? p.traits[TRAIT].prog_P : null) : p.best_prog_P;

function renderTable(){
  const tb = $("#progTable tbody"); tb.innerHTML = "";
  const rows = rowsForTable();
  const lastHead = $$("#progTable thead th").find(th => th.dataset.k === "best_program_P");
  lastHead.textContent = TRAIT ? `${TRAIT} prog P` : "best prog P";
  $$("#progTable thead th").forEach(th => {
    th.classList.toggle("sorted", th.dataset.k === SORT.k);
    th.classList.toggle("asc", th.dataset.k === SORT.k && SORT.asc);
  });
  rows.forEach(p => {
    const m = p.meta || {};
    const sigFDR = isNum(m.FDR) && m.FDR < 0.05;
    const tr = el("tr", {"data-prog": p.program, onclick: () => openDrawer(p.program)}, [
      el("td", {class: "num dim"}, [String(p.num)]),
      el("td", {}, [el("span", {class: "progname"}, [p.program + " "]),
                    el("span", {class: "anno" + (p.annotation ? "" : " none")}, [p.annotation || "unlabeled"])]),
      el("td", {class: "num"}, [fB(m.pooled)]),
      el("td", {class: "num" + (sigFDR ? " sig" : "")}, [fP(m.FDR)]),
      el("td", {class: "num"}, [fI2(m.I2)]),
      el("td", {}, [m.cross_source_concordant ? el("span", {class: "pill ok"}, ["both"]) : el("span", {class: "dim"}, ["—"])]),
      el("td", {class: "num"}, [fP(traitP(p))]),
      el("td", {class: "genes"}, [(p.top_genes || []).slice(0, 6).join(", ")]),
    ]);
    tb.appendChild(tr);
  });
  $("#progCount").textContent = `${rows.length} programs · ${COND}` + (TRAIT ? ` · trait ${TRAIT}` : "");
}

// ---------- trait table ----------
function renderTraits(){
  const tb = $("#traitTable tbody"); tb.innerHTML = "";
  dz().traits.forEach(t => {
    ["rest", "stim48hr"].forEach((c, i) => {
      const pn = (t.panels || {})[c] || {};
      const flagCls = t.flag === "ok" ? "ok" : (t.flag.includes("SIGN") || t.flag.includes("PRIOR") ? "bad" : "warn");
      const permSig = isNum(pn.fig5a_perm_P) && pn.fig5a_perm_P < 0.05;
      const tr = el("tr", {}, [
        el("td", {}, i === 0 ? [el("b", {}, [t.short]), el("div", {class: "dim", html: t.desc})] : [el("span", {class: "dim"}, ["↳"])]),
        el("td", {}, i === 0 ? [t.source] : [""]),
        el("td", {}, i === 0 ? [el("span", {class: "pill " + flagCls}, [t.flag])] : [""]),
        el("td", {class: "dim"}, [c]),
        el("td", {}, [pn.fig4c_prog_hit && pn.fig4c_prog_hit !== "-" ? el("a", {onclick: () => { COND = c; syncCond(); openDrawer(pn.fig4c_prog_hit); }}, [pn.fig4c_prog_hit]) : el("span", {class: "dim"}, ["—"])]),
        el("td", {class: "num" + (permSig ? " sig" : "")}, [fP(pn.fig5a_perm_P)]),
        el("td", {class: "genes dim"}, [selectedProgs(pn.fig5a_selected)]),
      ]);
      tb.appendChild(tr);
    });
  });
}
function selectedProgs(s){
  if (!s) return "—";
  const parts = s.split(" ").filter(x => /^P\d/.test(x));
  return parts.join("  ");
}
function syncCond(){ $$("#condSeg button").forEach(x => x.classList.toggle("on", x.dataset.cond === COND)); renderTable(); }

// ---------- network module (interactive combined Fig 5a) ----------
const NS = "http://www.w3.org/2000/svg";
const svgEl = (t, a) => { const e = document.createElementNS(NS, t); for (const k in a) if (a[k] != null) e.setAttribute(k, a[k]); return e; };
const CONS = ["#3d4a5c", "#4f6d9a", "#5a86c9", "#6ea8fe", "#8fc0ff"];  // conservation ramp 1..5+
const consColor = n => CONS[Math.min((n || 1) - 1, CONS.length - 1)];

function curNet(){ const n = dz().networks; return n && NET.group ? n[COND][NET.group] : null; }
function canonOf(f){ const net = curNet(); if (!net || !f) return null; f = f.toUpperCase();
  const m = net.nodes.find(n => n.id.toUpperCase() === f); return m ? m.id : null; }

function filteredGraph(){
  const net = curNet(); if (!net) return {nodes: [], edges: []};
  let nodes = net.nodes.map(n => Object.assign({}, n));
  const byId = {}; nodes.forEach(n => byId[n.id] = n);
  let edges = net.edges.map(e => Object.assign({}, e));  // gene→program (regulates/member) + program→trait (selects)
  // program drop from the filter toggles
  const drop = new Set();
  nodes.forEach(n => { if (n.type === "program") {
    if (NET.annot && !n.annotated) drop.add(n.id);
    if (NET.sig && !n.sig) drop.add(n.id);
  }});
  edges = edges.filter(e => !drop.has(e.source) && !drop.has(e.target));
  nodes = nodes.filter(n => !drop.has(n.id));
  if (NET.sig) edges = edges.filter(e => e.cls === "selects"
      ? true : ((byId[e.source] && byId[e.source].n_traits >= 2) || (byId[e.target] && byId[e.target].sig)));
  // focus: keep the path through the focused node (gene→program→trait)
  const focusId = canonOf(NET.focus);
  if (focusId) {
    const F = byId[focusId];
    const nb = (id, t) => { const o = []; edges.forEach(e => {
      if (e.source === id && (!t || (byId[e.target] || {}).type === t)) o.push(e.target);
      if (e.target === id && (!t || (byId[e.source] || {}).type === t)) o.push(e.source); }); return o; };
    const keep = new Set([focusId]);
    if (F && F.type === "gene") nb(focusId, "program").forEach(p => { keep.add(p); nb(p, "trait").forEach(t => keep.add(t)); });
    else if (F && F.type === "trait") nb(focusId, "program").forEach(p => { keep.add(p); nb(p, "gene").forEach(g => keep.add(g)); });
    else nb(focusId).forEach(x => keep.add(x));
    edges = edges.filter(e => keep.has(e.source) && keep.has(e.target));
    nodes = nodes.filter(n => keep.has(n.id));
  }
  // drop orphans (but keep the focus node itself)
  const deg = {}; edges.forEach(e => { deg[e.source] = (deg[e.source] || 0) + 1; deg[e.target] = (deg[e.target] || 0) + 1; });
  nodes = nodes.filter(n => deg[n.id] || n.id === focusId);
  const ids = new Set(nodes.map(n => n.id));
  edges = edges.filter(e => ids.has(e.source) && ids.has(e.target));
  return {nodes, edges};
}

function layout(nodes, edges, W, H){
  const pad = 26, colX = {gene: W * 0.16, program: W * 0.52, trait: W * 0.86};
  const col = t => colX[t] != null ? colX[t] : W * 0.5;
  const groups = {gene: [], program: [], trait: []};
  nodes.forEach(n => (groups[n.type] || (groups[n.type] = [])).push(n));
  Object.values(groups).forEach(arr => {
    arr.sort((a, b) => (b.n_traits || 1) - (a.n_traits || 1) || String(a.id).localeCompare(b.id));
    const n = arr.length;
    arr.forEach((nd, i) => { nd.x = col(nd.type); nd.y = pad + (n <= 1 ? (H - 2 * pad) / 2 : (H - 2 * pad) * i / (n - 1)); nd.vy = 0; });
  });
  const pos = {}; nodes.forEach(n => pos[n.id] = n);
  for (let it = 0; it < 240; it++) {
    edges.forEach(e => { const a = pos[e.source], b = pos[e.target]; if (!a || !b) return;
      const m = (a.y + b.y) / 2; a.vy += (m - a.y) * 0.02; b.vy += (m - b.y) * 0.02; });
    Object.values(groups).forEach(arr => { for (let i = 0; i < arr.length; i++) for (let j = i + 1; j < arr.length; j++) {
      const a = arr[i], b = arr[j], d = a.y - b.y, min = 18; if (Math.abs(d) < min) { const f = (min - Math.abs(d)) * 0.5 * Math.sign(d || 1); a.vy += f; b.vy -= f; } } });
    nodes.forEach(n => { n.y += n.vy; n.vy *= 0.55; n.y = Math.max(pad, Math.min(H - pad, n.y)); });
  }
}

function renderNetwork(){
  const net = curNet(); const svg = $("#netsvg");
  if (!svg) return; // network tab now hosts the React app (iframe); legacy renderer disabled
  $$("#netGroup button").forEach(b => b.classList.toggle("on", b.textContent === netLabel(NET.group)));
  buildLegend();
  while (svg.firstChild) svg.removeChild(svg.firstChild);
  if (!net) { $("#netCount").textContent = "no network"; return; }
  const focusList = $("#netFocusList"); focusList.innerHTML = "";
  net.nodes.forEach(n => focusList.appendChild(el("option", {value: n.id})));
  const wrap = $("#netwrap"), W = wrap.clientWidth || 700, H = wrap.clientHeight || 460;
  svg.setAttribute("viewBox", `0 0 ${W} ${H}`);
  const {nodes, edges} = filteredGraph();
  $("#netCount").textContent = `${nodes.filter(n => n.type === "gene").length} genes · ${nodes.filter(n => n.type === "program").length} programs · ${nodes.filter(n => n.type === "trait").length} traits · ${edges.length} edges` + (NET.focus ? ` · focus ${NET.focus}` : "");
  layout(nodes, edges, W, H);
  const pos = {}; nodes.forEach(n => pos[n.id] = n);
  ZOOM = {k: 1, tx: 0, ty: 0};
  const vp = svgEl("g", {id: "netvp"}); svg.appendChild(vp);
  const gEdges = svgEl("g"), gNodes = svgEl("g"); vp.appendChild(gEdges); vp.appendChild(gNodes);
  edges.forEach(e => { const a = pos[e.source], b = pos[e.target]; if (!a || !b) return;
    const w = e.cls === "regulates" ? Math.max(0.6, Math.min(4, Math.abs(e.beta || 1) * 0.5)) : (e.cls === "selects" ? 0.8 : 0.7);
    const ln = svgEl("line", {class: "edge " + e.cls, x1: a.x, y1: a.y, x2: b.x, y2: b.y,
      "stroke-width": w, "stroke-opacity": e.cls === "regulates" ? 0.6 : 0.35});
    ln.dataset.s = e.source; ln.dataset.t = e.target; gEdges.appendChild(ln); });
  nodes.forEach(n => {
    const g = svgEl("g", {class: "node" + (n.type === "gene" && n.role !== "regulator" ? " member" : "")});
    g.dataset.id = n.id;
    if (n.type === "program") {
      const r = 7 + (n.n_traits || 1) * 1.6;
      g.appendChild(svgEl("circle", {class: "pn" + (n.sig ? " sig" : ""), cx: n.x, cy: n.y, r, fill: consColor(n.n_traits)}));
      const lab = svgEl("text", {x: n.x + r + 3, y: n.y + 3.5}); lab.textContent = n.id + (n.annotation ? " " + n.annotation : "");
      g.appendChild(lab);
    } else if (n.type === "trait") {
      const w = 10 + n.id.length * 6.3, h = 15;
      g.appendChild(svgEl("rect", {class: "tn", x: n.x, y: n.y - h / 2, width: w, height: h, rx: 7}));
      const lab = svgEl("text", {class: "tlabel", x: n.x + w / 2, y: n.y + 3.5, "text-anchor": "middle"}); lab.textContent = n.id;
      g.appendChild(lab);
    } else {  // gene — always labelled; members faded via the group's .member class
      const r = 3.2 + (n.n_traits || 1) * 1.05, reg = n.role === "regulator";
      g.appendChild(svgEl(reg ? "circle" : "rect", reg
        ? {class: "gn " + (n.sign > 0 ? "pos" : "neg"), cx: n.x, cy: n.y, r, fill: consColor(n.n_traits)}
        : {class: "gn mem", x: n.x - r, y: n.y - r, width: 2 * r, height: 2 * r}));
      const lab = svgEl("text", {x: n.x - r - 3, y: n.y + 3.5, "text-anchor": "end"}); lab.textContent = n.id;
      g.appendChild(lab);
    }
    g._oy = n.y;
    g.addEventListener("mouseenter", ev => netTip(ev, n, edges));
    g.addEventListener("mouseleave", () => { $("#nettip").classList.add("hidden"); undim(); });
    g.addEventListener("mousedown", e => { netDrag = {g, n, moved: false}; e.preventDefault(); e.stopPropagation(); });
    g.addEventListener("click", () => { if (netDrag && netDrag.moved) return;
      if (n.type === "program") openDrawer(n.id);
      else { $("#netFocus").value = n.id; NET.focus = n.id; renderNetwork(); } });
    gNodes.appendChild(g);
  });
  applyVP();
  vp.style.opacity = 0; requestAnimationFrame(() => { vp.style.opacity = 1; });
}
function undim(){ $$("#netsvg .dim").forEach(e => e.classList.remove("dim")); }
function netTip(ev, n, edges){
  const tip = $("#nettip"); let h = "";
  if (n.type === "program") h = `<b>${n.id}</b> ${n.annotation || "(unlabeled)"}<br><span class="mut">selected in ${n.n_traits} phenotype(s)${n.sig ? " · meta-significant" : ""} · role ${n.role} · ${(n.traits || []).join(", ")}</span>`;
  else if (n.type === "trait") h = `<b>${n.id}</b><br><span class="mut">phenotype — the programs/genes linked to it were selected by it</span>`;
  else h = `<b>${n.id}</b> ${n.role} gene<br><span class="mut">predicted γ-sign ${n.sign > 0 ? "+" : n.sign < 0 ? "−" : "?"} · in ${n.n_traits} phenotype(s)</span>`;
  tip.innerHTML = h;
  const wrap = $("#netwrap").getBoundingClientRect();
  tip.style.left = Math.min(wrap.width - 240, ev.clientX - wrap.left + 12) + "px";
  tip.style.top = (ev.clientY - wrap.top + 12) + "px"; tip.classList.remove("hidden");
  const nb = new Set([n.id]); edges.forEach(e => { if (e.source === n.id) nb.add(e.target); if (e.target === n.id) nb.add(e.source); });
  $$("#netsvg .node").forEach(g => g.classList.toggle("dim", !nb.has(g.dataset.id)));
  $$("#netsvg .edge").forEach(l => l.classList.toggle("dim", !(l.dataset.s === n.id || l.dataset.t === n.id)));
}
let ZOOM = {k: 1, tx: 0, ty: 0}, netPan = null;
function applyVP(){ const vp = document.getElementById("netvp"); if (vp) vp.setAttribute("transform", `translate(${ZOOM.tx},${ZOOM.ty}) scale(${ZOOM.k})`); }
function svgPt(e){ const svg = $("#netsvg"), r = svg.getBoundingClientRect(), vb = svg.viewBox.baseVal;
  return {x: (e.clientX - r.left) * (vb.width / r.width), y: (e.clientY - r.top) * (vb.height / r.height)}; }
function initNetInteract(){
  const svg = $("#netsvg");
  if (!svg) return; // legacy vanilla network disabled (React app runs in the iframe)
  svg.addEventListener("wheel", e => { e.preventDefault(); const p = svgPt(e);
    const f = e.deltaY < 0 ? 1.12 : 1 / 1.12, nk = Math.max(0.3, Math.min(6, ZOOM.k * f)), r = nk / ZOOM.k;
    ZOOM.tx = p.x - (p.x - ZOOM.tx) * r; ZOOM.ty = p.y - (p.y - ZOOM.ty) * r; ZOOM.k = nk; applyVP(); }, {passive: false});
  svg.addEventListener("mousedown", e => { if (netDrag) return; const p = svgPt(e); netPan = {x: p.x, y: p.y, tx: ZOOM.tx, ty: ZOOM.ty}; svg.style.cursor = "grabbing"; });
  window.addEventListener("mousemove", e => {
    if (netDrag) { const p = svgPt(e), y = (p.y - ZOOM.ty) / ZOOM.k, n = netDrag.n; netDrag.moved = true;
      netDrag.g.setAttribute("transform", `translate(0,${y - netDrag.g._oy})`);
      const vp = document.getElementById("netvp"); if (vp) vp.querySelectorAll("line").forEach(l => {
        if (l.dataset.s === n.id) l.setAttribute("y1", y); if (l.dataset.t === n.id) l.setAttribute("y2", y); }); return; }
    if (netPan) { const p = svgPt(e); ZOOM.tx = netPan.tx + (p.x - netPan.x); ZOOM.ty = netPan.ty + (p.y - netPan.y); applyVP(); }
  });
  window.addEventListener("mouseup", () => { netPan = null; const s = $("#netsvg"); if (s) s.style.cursor = "grab"; setTimeout(() => { netDrag = null; }, 0); });
}
function buildLegend(){
  $("#netLegend").innerHTML =
    '<span><i style="background:#6ea8fe"></i>program (size = #phenotypes)</span>' +
    '<span><i style="background:#5a86c9;border:1.6px solid var(--accent2);border-radius:50%"></i>regulator +γ</span>' +
    '<span><i style="background:#5a86c9;border:1.6px solid var(--bad);border-radius:50%"></i>−γ</span>' +
    '<span><i style="background:var(--memfill);border-radius:2px"></i>member gene (faded)</span>' +
    '<span><i style="background:var(--tnode);border-radius:6px"></i>phenotype</span>' +
    '<span><span class="line"></span>regulates (|β|)</span>' +
    '<span><span class="line dash"></span>member</span>' +
    '<span class="mut">scroll = zoom · drag bg = pan · drag node = move</span>';
}

// ---------- regulator module ----------
function renderRegulators(){
  const tb = $("#regTable tbody"); tb.innerHTML = "";
  let regs = (dz().regulators || []).slice();
  if (REG.shared) regs = regs.filter(r => r.shared);
  if (REG.prog) regs = regs.filter(r => r.targets.some(t => t.program === REG.prog));
  if (REG.search) regs = regs.filter(r => r.gene.toUpperCase().includes(REG.search));
  const dir = REGSORT.asc ? 1 : -1;
  regs.sort((a, b) => { const x = a[REGSORT.k], y = b[REGSORT.k];
    return (x > y ? 1 : x < y ? -1 : (a.gene > b.gene ? 1 : -1)) * dir; });
  $$("#regTable thead th").forEach(th => th.classList.toggle("sorted", th.dataset.k === REGSORT.k));
  regs.forEach(r => {
    const tr = el("tr", {onclick: () => openReg(r.gene)}, [
      el("td", {}, [el("b", {}, [r.gene]), r.shared ? el("span", {class: "pill ok", html: "&nbsp;shared"}) : ""]),
      el("td", {}, [r.sign > 0 ? el("span", {class: "sig"}, ["+"]) : el("span", {class: "dim"}, ["−"])]),
      el("td", {class: "num"}, [String(r.n_programs)]),
      el("td", {class: "num"}, [`${r.n_clean_traits}/${r.n_traits}`]),
      el("td", {class: "dim"}, [r.conditions.join(", ")]),
      el("td", {class: "genes"}, [r.targets.slice(0, 5).map(t =>
        `${t.program}(${t.dir > 0 ? "+" : "−"}β${Math.abs(t.beta || 0).toFixed(1)})`).join("  ")]),
    ]);
    tb.appendChild(tr);
  });
  $("#regCount").textContent = `${regs.length} regulator genes` + (REG.shared ? " · shared" : "") + (REG.prog ? ` · →${REG.prog}` : "");
}
function openReg(gene){
  const r = (dz().regulators || []).find(x => x.gene === gene); if (!r) return;
  const box = $("#drawer .drawer-inner"); box.innerHTML = "";
  box.appendChild(el("span", {class: "close"}, ["×"]));
  box.appendChild(el("h2", {}, [`${r.gene}  `, el("span", {class: r.sign > 0 ? "sig" : "dim"}, [`γ-sign ${r.sign > 0 ? "+" : "−"}`])]));
  box.appendChild(el("p", {class: "sub"}, [`regulator gene · ${r.n_programs} programs · ${r.n_clean_traits}/${r.n_traits} phenotypes (clean/all) · ${r.conditions.join(", ")}`]));
  if (r.shared) box.appendChild(el("p", {}, [el("span", {class: "pill ok"}, ["shared regulator — recurs across QC-clean phenotypes"])]));
  box.appendChild(el("h3", {}, ["Regulates"]));
  const tbl = el("table"); tbl.appendChild(el("thead", {html:
    "<tr><th>program</th><th class='num'>β</th><th>dir</th><th>trait</th><th>cond</th></tr>"}));
  const tb = el("tbody");
  r.targets.forEach(t => tb.appendChild(el("tr", {}, [
    el("td", {}, [el("a", {onclick: () => openDrawer(t.program)}, [t.program]), " ", el("span", {class: "dim"}, [t.annotation || ""])]),
    el("td", {class: "num"}, [(t.beta || 0).toFixed(2)]),
    el("td", {class: t.dir > 0 ? "sig" : "dim"}, [t.dir > 0 ? "↑" : "↓"]),
    el("td", {}, [t.trait, t.flag !== "ok" ? el("span", {class: "pill warn", html: "&nbsp;" + t.flag.split(";")[0]}) : ""]),
    el("td", {class: "dim"}, [t.condition]),
  ])));
  tbl.appendChild(tb); box.appendChild(tbl);
  box.appendChild(el("p", {class: "regfocus"}, [el("a", {onclick: () => {
    $$(".tabs button").forEach(x => x.classList.toggle("on", x.dataset.tab === "network"));
    $$(".tabpane").forEach(p => p.classList.toggle("on", p.id === "tab-network"));
    $("#netFocus").value = r.gene; NET.focus = r.gene; renderNetwork(); closeDrawer();
  }}, ["→ focus this gene in the network"])]));
  $("#drawer").classList.remove("hidden");
}

function lightbox(src, cap){
  const lb = $("#lightbox"); lb.querySelector("img").src = src; lb.querySelector("figcaption").textContent = cap || "";
  lb.classList.remove("hidden");
}

// ---------- conclusions ----------
function renderNotes(){
  const txt = (dz().narrative && dz().narrative.conclusions) || "no conclusions available";
  $("#notes").innerHTML = mdBlock(txt);
}

// ---------- program drawer ----------
function openDrawer(program){
  const p = progs().find(x => x.program === program);
  if (!p){ return; }
  $$("#progTable tbody tr").forEach(tr => tr.classList.toggle("hl", tr.dataset.prog === program));
  const tr = $(`#progTable tbody tr[data-prog='${program}']`);
  if (tr) tr.scrollIntoView({block: "nearest"});
  const m = p.meta || {};
  const box = $("#drawer .drawer-inner"); box.innerHTML = "";
  box.appendChild(el("span", {class: "close"}, ["×"]));
  box.appendChild(el("h2", {}, [`${p.program} — ${p.annotation || "unlabeled"}`]));
  box.appendChild(el("p", {class: "sub"}, [`${dz().name} · ${COND} · cNMF program`]));
  box.appendChild(el("p", {class: "dim", html: "<b>Top-loaded genes:</b> " + (p.top_genes || []).join(", ")}));

  box.appendChild(el("h3", {}, ["Cross-source meta-analysis"]));
  const kv = el("div", {class: "kv"});
  const add = (k, v) => { kv.appendChild(el("b", {}, [k])); kv.appendChild(el("span", {}, [v])); };
  add("pooled β", `${fB(m.pooled)} ± ${isNum(m.pooled_se) ? m.pooled_se.toFixed(3) : "—"}`);
  add("meta P / FDR", `${fP(m.P)} / ${fP(m.FDR)}${isNum(m.FDR) && m.FDR < 0.05 ? "  ✓" : ""}`);
  add("I² (heterogeneity)", fI2(m.I2));
  add("phenotypes (k)", `${m.k || "—"}  (pos ${m.n_pos ?? "—"} / neg ${m.n_neg ?? "—"})`);
  add("LOO sign-stable", m.loo_sign_stable ? "yes" : "no");
  add("cross-source concordant", m.cross_source_concordant ? "yes (Backman & GeneBass agree)" : "no");
  add("Backman / GeneBass β", `${fB(m.backman_pooled)} / ${fB(m.genebass_pooled)}`);
  box.appendChild(kv);

  box.appendChild(el("h3", {}, ["Per-phenotype burden"]));
  const tbl = el("table"); tbl.appendChild(el("thead", {html:
    "<tr><th>trait</th><th>flag</th><th class='num'>prog P</th><th class='num'>reg P</th><th class='num'>reg β</th></tr>"}));
  const tb = el("tbody");
  Object.keys(p.traits).forEach(short => {
    const t = p.traits[short];
    const fl = t.flag === "ok" ? "ok" : (t.flag.includes("SIGN") || t.flag.includes("PRIOR") ? "bad" : "warn");
    tb.appendChild(el("tr", {}, [
      el("td", {}, [short]),
      el("td", {}, [el("span", {class: "pill " + fl}, [t.flag])]),
      el("td", {class: "num" + (isNum(t.prog_P) && t.prog_P < 0.05 ? " sig" : "")}, [fP(t.prog_P)]),
      el("td", {class: "num" + (isNum(t.reg_P) && t.reg_P < 0.05 ? " sig" : "")}, [fP(t.reg_P)]),
      el("td", {class: "num"}, [fB(t.reg_beta)]),
    ]));
  });
  tbl.appendChild(tb); box.appendChild(tbl);

  const forest = dz().figures.find(f => f.kind === "condition_comparison" && f.file.includes("3_forest"));
  if (forest) box.appendChild(el("p", {}, [el("a", {onclick: () => lightbox("/static/" + forest.file, forest.caption)}, ["→ open meta-analysis forest plot"])]));
  $("#drawer").classList.remove("hidden");
}
function closeDrawer(){ $("#drawer").classList.add("hidden"); $$("#progTable tbody tr").forEach(tr => tr.classList.remove("hl")); }

// ---------- chat ----------
const convo = [];
async function onSend(e){
  e.preventDefault();
  const q = $("#chatInput").value.trim(); if (!q) return;
  $("#chatInput").value = ""; $("#sendBtn").disabled = true;
  addMsg("user", q); convo.push({role: "user", content: q});
  const bot = addMsg("bot", ""); bot.innerHTML = "<span class='spinner'></span> thinking…";
  try{
    const r = await fetch("/api/chat", {method: "POST", headers: {"Content-Type": "application/json"},
      body: JSON.stringify({messages: convo, disease: DID})});
    const j = await r.json();
    if (j.error){ bot.innerHTML = "<span class='dim'>⚠ " + esc(j.error) + "</span>"; }
    else { renderAnswer(bot, j); convo.push({role: "assistant", content: j.answer || ""}); }
  } catch (err){ bot.innerHTML = "<span class='dim'>⚠ " + esc(String(err)) + "</span>"; }
  $("#sendBtn").disabled = false; $("#messages").scrollTop = 1e9;
}
function addMsg(role, text){
  const m = el("div", {class: "msg " + role}); m.innerHTML = role === "user" ? esc(text) : text;
  $("#messages").appendChild(m); $("#messages").scrollTop = 1e9; return m;
}
function renderAnswer(node, j){
  node.innerHTML = mdBlock(j.answer || "(no answer)");
  if (j.trace && j.trace.length){
    const names = j.trace.map(t => t.tool);
    const det = el("details", {class: "trace"});
    det.appendChild(el("summary", {}, [`grounded via ${j.trace.length} tool call(s)`]));
    j.trace.forEach(t => det.appendChild(el("div", {html: `<code>${esc(t.tool)}</code> ${esc(JSON.stringify(t.input || {}))}`})));
    node.appendChild(det);
  }
  if (j.figures && j.figures.length){
    const row = el("div", {class: "figrow"});
    j.figures.slice(0, 6).forEach(f => row.appendChild(
      el("img", {src: "/static/" + f.file, title: f.caption, onclick: () => lightbox("/static/" + f.file, f.caption)})));
    node.appendChild(row);
  }
  wireProgLinks(node);
}
function wireProgLinks(node){
  node.querySelectorAll("a[data-prog]").forEach(a =>
    a.addEventListener("click", () => openDrawer(a.dataset.prog)));
}

// ---------- tiny markdown ----------
const esc = s => s.replace(/[&<>]/g, c => ({"&": "&amp;", "<": "&lt;", ">": "&gt;"}[c]));
function inline(s){
  s = esc(s)
    .replace(/\*\*(.+?)\*\*/g, "<b>$1</b>")
    .replace(/`([^`]+)`/g, "<code>$1</code>")
    .replace(/\bhttps?:\/\/[^\s)]+/g, u => `<a href="${u}" target="_blank" rel="noopener">${u}</a>`)
    .replace(/\b(P\d{1,2})\b/g, '<a data-prog="$1">$1</a>');
  return s;
}
function mdBlock(md){
  const lines = md.split("\n"); let html = "", inList = false, inTable = false;
  const closeList = () => { if (inList){ html += "</ul>"; inList = false; } };
  const closeTable = () => { if (inTable){ html += "</table>"; inTable = false; } };
  for (let raw of lines){
    const line = raw.replace(/\s+$/, "");
    if (/^#{1,6}\s/.test(line)){ closeList(); closeTable(); const lvl = line.match(/^#+/)[0].length; html += `<h${Math.min(lvl+1,4)}>${inline(line.replace(/^#+\s/, ""))}</h${Math.min(lvl+1,4)}>`; }
    else if (/^\s*[-*]\s+/.test(line)){ closeTable(); if (!inList){ html += "<ul>"; inList = true; } html += `<li>${inline(line.replace(/^\s*[-*]\s+/, ""))}</li>`; }
    else if (/^\|.*\|/.test(line)){ closeList(); if (/^\|[\s:|-]+\|?$/.test(line)) continue; if (!inTable){ html += "<table>"; inTable = true; } const cells = line.split("|").slice(1, -1); html += "<tr>" + cells.map(c => `<td>${inline(c.trim())}</td>`).join("") + "</tr>"; }
    else if (line.trim() === ""){ closeList(); closeTable(); }
    else { closeList(); closeTable(); html += `<p>${inline(line)}</p>`; }
  }
  closeList(); closeTable(); return html;
}

document.addEventListener("DOMContentLoaded", init);
