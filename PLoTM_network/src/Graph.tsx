import {
  forwardRef,
  useCallback,
  useEffect,
  useImperativeHandle,
  useMemo,
  useRef,
  useState,
} from "react";
import type React from "react";
import type { Simulation, ForceLink } from "d3-force";
import { forceX, forceY, forceManyBody, forceCollide } from "d3-force";
import type { NetworkGroup, SimNode, SimEdge, Spacing } from "./types";
import { buildGraph, buildSimulation, seedPositions, seedRandom, layerX } from "./sim";
import {
  nodeRadius,
  collidePad,
  edgeColor,
  edgeWidth,
  edgePath,
  nid,
  UP,
  DOWN,
  MEM,
  GFILL,
  gLabel,
} from "./viz";

export interface GraphHandle {
  regenerate: () => void;
  reset: () => void;
  start: (animate: boolean) => void;
  fit: () => void;
  focus: (id: string) => boolean;
  exportSVG: () => void;
  exportPNG: () => void;
  ids: () => string[];
}

interface Props {
  group: NetworkGroup;
  groupKey: string; // condition|group — rebuild when this changes
  spacing: Spacing;
  theme: "dark" | "light";
  reduced: boolean;
  started: boolean; // false until the user presses Start — sim stays parked, graph hidden
  animate: boolean; // whether Start plays the physics reveal
  onSelect?: (n: SimNode | null) => void;
}

interface Palette {
  ink: string;
  tnode: string;
  memFill: string;
  memStroke: string;
  ring: string;
  pin: string;
  focus: string;
  pbox: string; // program box fill
  tfill: string; // trait node fill
}
const palettes: Record<string, Palette> = {
  dark: { ink: "#e6edf3", tnode: "#33415a", memFill: "#243046", memStroke: "#39475c", ring: "#0b0e12", pin: "#e9a23b", focus: "#8fc0ff", pbox: "#28313e", tfill: "#5a2e34" },
  light: { ink: "#1a2230", tnode: "#d7e2f1", memFill: "#e4eaf3", memStroke: "#c4cfdf", ring: "#f0f4fa", pin: "#b9760f", focus: "#2f6fe0", pbox: "#eef1f6", tfill: "#ffe3e3" },
};

// wrap a program annotation into up to 2 short lines for the box label
function wrapLabel(s: string, width = 20): string[] {
  if (!s) return [];
  const words = s.split(/\s+/);
  const lines: string[] = [];
  let cur = "";
  for (const w of words) {
    if ((cur + " " + w).trim().length > width && cur) { lines.push(cur); cur = w; }
    else cur = (cur + " " + w).trim();
    if (lines.length === 2) break;
  }
  if (cur && lines.length < 2) lines.push(cur);
  if (lines.length === 2 && cur !== lines[1]) lines[1] = lines[1].slice(0, width - 1) + "…";
  return lines.slice(0, 2);
}

export const Graph = forwardRef<GraphHandle, Props>(function Graph(
  { group, groupKey, spacing, theme, reduced, started, animate, onSelect },
  ref,
) {
  const wrapRef = useRef<HTMLDivElement>(null);
  const svgRef = useRef<SVGSVGElement>(null);
  const simRef = useRef<Simulation<SimNode, SimEdge> | null>(null);
  const nodesRef = useRef<SimNode[]>([]);
  const edgesRef = useRef<SimEdge[]>([]);
  const builtKey = useRef<string>("");
  const rafRef = useRef(0);
  const drawRaf = useRef(0);
  const dragRef = useRef<{ n: SimNode; moved: boolean; lx: number; ly: number; vx: number; vy: number } | null>(null);
  const panRef = useRef<{ x: number; y: number; tx: number; ty: number } | null>(null);

  const [dims, setDims] = useState({ w: 0, h: 0 });
  const [, setVersion] = useState(0);
  const [builtN, setBuiltN] = useState(0); // bumped after each (re)build so adjacency recomputes on populated refs
  const [tf, setTf] = useState({ k: 1, x: 0, y: 0 });
  const [hover, setHover] = useState<string | null>(null);
  const [hoverEdge, setHoverEdge] = useState<SimEdge | null>(null);
  const [edgeSel, setEdgeSel] = useState<SimEdge | null>(null);
  const [pinned, setPinned] = useState<Set<string>>(() => new Set<string>()); // persistent highlight set
  const [draw, setDraw] = useState<number | null>(null); // connection-animation progress 0..1 (null = off)
  const pal = palettes[theme];
  // live mirrors so the (dependency-light) build effect can read current values without re-running
  const startedRef = useRef(started); startedRef.current = started;
  const animateRef = useRef(animate); animateRef.current = animate;
  const startRef = useRef<(a: boolean) => void>(() => {});
  const clearSel = useCallback(() => { setEdgeSel(null); setHoverEdge(null); }, []);
  const unpinAll = useCallback(() => {
    nodesRef.current.forEach((n) => { if (n.pinned) { n.pinned = false; n.fx = null; n.fy = null; } });
    setPinned(new Set());
  }, []);

  const bump = useCallback(() => setVersion((v) => (v + 1) % 1e9), []);
  const scheduleTick = useCallback(() => {
    if (rafRef.current) return;
    rafRef.current = requestAnimationFrame(() => {
      rafRef.current = 0;
      bump();
    });
  }, [bump]);

  const settle = useCallback(
    (alpha = 0.6) => {
      const sim = simRef.current;
      if (!sim) return;
      if (reduced) {
        sim.stop();
        for (let i = 0; i < 220; i++) sim.tick();
        bump();
      } else {
        sim.alpha(alpha).restart();
      }
    },
    [reduced, bump],
  );

  // ResizeObserver → responsive canvas
  useEffect(() => {
    const wrap = wrapRef.current;
    if (!wrap) return;
    const ro = new ResizeObserver(() => {
      const r = wrap.getBoundingClientRect();
      setDims({ w: Math.max(320, r.width), h: Math.max(360, r.height) });
    });
    ro.observe(wrap);
    return () => ro.disconnect();
  }, []);

  const applyForces = useCallback(() => {
    const sim = simRef.current;
    if (!sim) return;
    const { w, h } = dims;
    sim.force("x", forceX<SimNode>((d) => layerX(d.type, w, spacing.spread)).strength(0.92));
    sim.force("y", forceY<SimNode>(h / 2).strength(0.045));
    sim.force("charge", forceManyBody<SimNode>().strength((d) => spacing.charge * (1 + (d.degree || 0) * 0.35)));
    sim.force("collide", forceCollide<SimNode>((d) => nodeRadius(d) + collidePad(d) + spacing.collide).iterations(2));
    (sim.force("link") as ForceLink<SimNode, SimEdge> | undefined)?.distance(spacing.linkDist);
  }, [dims, spacing]);

  // world-space bounding box of the current node positions (finite nodes only)
  const bounds = useCallback(() => {
    const ns = nodesRef.current;
    let minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
    ns.forEach((n) => {
      if (!Number.isFinite(n.x) || !Number.isFinite(n.y)) return;
      const r = nodeRadius(n) + 20;
      minX = Math.min(minX, n.x! - r); maxX = Math.max(maxX, n.x! + r);
      minY = Math.min(minY, n.y! - r); maxY = Math.max(maxY, n.y! + r);
    });
    if (!Number.isFinite(minX) || !Number.isFinite(maxX)) return null;
    return { minX, minY, maxX, maxY };
  }, []);

  // keep the graph from being panned/zoomed entirely off-screen (a lost panel reads as a "crash")
  const clampTf = useCallback((t: { k: number; x: number; y: number }) => {
    if (!Number.isFinite(t.x) || !Number.isFinite(t.y) || !Number.isFinite(t.k)) return t;
    const b = bounds();
    if (!b || !dims.w) return t;
    const m = 90; // keep at least this many px of the graph on screen
    let { x, y } = t;
    const loX = m - t.k * b.maxX, hiX = dims.w - m - t.k * b.minX;
    if (loX <= hiX) x = Math.max(loX, Math.min(hiX, x));
    const loY = m - t.k * b.maxY, hiY = dims.h - m - t.k * b.minY;
    if (loY <= hiY) y = Math.max(loY, Math.min(hiY, y));
    return { k: t.k, x, y };
  }, [bounds, dims]);

  const fit = useCallback(() => {
    const b = bounds();
    if (!b || !dims.w) return;
    const gw = b.maxX - b.minX + 80, gh = b.maxY - b.minY + 80;
    if (gw <= 0 || gh <= 0) return;
    const k = Math.min(1.2, Math.min(dims.w / gw, dims.h / gh));
    if (!Number.isFinite(k) || k <= 0) return;
    const cx = (b.minX + b.maxX) / 2, cy = (b.minY + b.maxY) / 2;
    setTf({ k, x: dims.w / 2 - k * cx, y: dims.h / 2 - k * cy });
  }, [bounds, dims]);

  // build / rebuild
  useEffect(() => {
    if (!dims.w) return;
    if (builtKey.current !== groupKey) {
      simRef.current?.stop();
      const g = buildGraph(group, dims.w, dims.h);
      nodesRef.current = g.nodes;
      edgesRef.current = g.edges;
      const sim = buildSimulation(g.nodes, g.edges, dims.w, dims.h, spacing, scheduleTick);
      simRef.current = sim;
      builtKey.current = groupKey;
      setHover(null);
      setEdgeSel(null);
      setHoverEdge(null);
      setPinned(new Set());
      setBuiltN((v) => v + 1); // refs are now populated → force adjacency recompute this render cycle
      setTf({ k: 1, x: 0, y: 0 });
      if (startedRef.current) {
        // already past the intro gate (e.g. condition switch) → play the reveal for the new layout
        sim.stop();
        startRef.current(animateRef.current);
      } else {
        sim.stop(); // park the simulation until the Start button is pressed
      }
    } else if (startedRef.current) {
      applyForces();
      settle(0.5);
    }
  }, [groupKey, dims, spacing, group, scheduleTick, applyForces, settle]);

  useEffect(() => () => { simRef.current?.stop(); cancelAnimationFrame(drawRaf.current); }, []);

  // The "cool" reveal: scatter every node to a random spot, then let the physics simulation pull
  // them back into the tri-partite layout while the edge groups fade in in order (gene→program,
  // then program→trait). With animation off (or reduced-motion) it snaps straight to the layout.
  const start = useCallback(
    (anim: boolean) => {
      const sim = simRef.current;
      if (!sim) return;
      cancelAnimationFrame(drawRaf.current);
      applyForces();
      if (!anim || reduced) {
        seedPositions(nodesRef.current, dims.w, dims.h, () => true);
        sim.stop();
        for (let i = 0; i < 320; i++) sim.tick();
        setDraw(null);
        bump();
        setTimeout(fit, 0);
        return;
      }
      seedRandom(nodesRef.current, dims.w, dims.h);
      setDraw(0);
      const t0 = performance.now();
      const D = 1700;
      const step = () => {
        const p = Math.min(1, (performance.now() - t0) / D);
        setDraw(p >= 1 ? null : p); // null when done → edges at full opacity
        if (p < 1) drawRaf.current = requestAnimationFrame(step);
        else { drawRaf.current = 0; setTimeout(fit, 250); }
      };
      drawRaf.current = requestAnimationFrame(step);
      sim.alpha(1).restart(); // physics animates the nodes home from the random scatter
    },
    [applyForces, reduced, dims, fit, bump],
  );
  startRef.current = start;

  // ---- coordinate helpers ----
  const toLocal = (clientX: number, clientY: number) => {
    const r = svgRef.current!.getBoundingClientRect();
    return { x: (clientX - r.left - tf.x) / tf.k, y: (clientY - r.top - tf.y) / tf.k };
  };

  // ---- node drag ----
  const onNodePointerDown = (e: React.PointerEvent, n: SimNode) => {
    e.stopPropagation();
    (e.target as Element).setPointerCapture?.(e.pointerId);
    dragRef.current = { n, moved: false, lx: n.x ?? 0, ly: n.y ?? 0, vx: 0, vy: 0 };
    n.fx = n.x; n.fy = n.y;
    // reheat the simulation so the whole graph reacts elastically while you drag
    if (!reduced) simRef.current?.alphaTarget(0.3).restart();
  };
  const onNodePointerMove = (e: React.PointerEvent, n: SimNode) => {
    const drag = dragRef.current;
    if (!drag || drag.n !== n) return;
    const p = toLocal(e.clientX, e.clientY);
    // track the pointer velocity (smoothed) so the release can throw the node with momentum
    drag.vx = 0.6 * drag.vx + 0.4 * (p.x - drag.lx);
    drag.vy = 0.6 * drag.vy + 0.4 * (p.y - drag.ly);
    drag.lx = p.x; drag.ly = p.y;
    n.fx = p.x; n.fy = p.y; drag.moved = true;
    if (reduced) { simRef.current?.tick(); bump(); }
  };
  const onNodePointerUp = (e: React.PointerEvent, n: SimNode) => {
    const drag = dragRef.current;
    dragRef.current = null;
    (e.target as Element).releasePointerCapture?.(e.pointerId);
    if (!drag) return;
    if (!drag.moved) { togglePin(n); return; }
    if (n.pinned) return; // pinned nodes stay where you drop them
    // release: unpin, launch with the drag momentum, and let forceX spring it back with overshoot
    n.fx = null; n.fy = null;
    n.vx = Math.max(-40, Math.min(40, drag.vx * 2.2));
    n.vy = Math.max(-40, Math.min(40, drag.vy * 2.2));
    const sim = simRef.current;
    if (sim) {
      sim.alphaTarget(0);
      if (reduced) { for (let i = 0; i < 220; i++) sim.tick(); bump(); }
      else sim.alpha(0.9).restart(); // high alpha → the spring-back is clearly animated
    }
  };

  const togglePin = (n: SimNode) => {
    clearSel();
    n.pinned = !n.pinned;
    if (n.pinned) { n.fx = n.x; n.fy = n.y; } else { n.fx = null; n.fy = null; settle(0.4); }
    setPinned((prev) => { const s = new Set(prev); if (n.pinned) s.add(n.id); else s.delete(n.id); return s; });
    onSelect?.(n.pinned ? n : null);
    bump();
  };

  // ---- keyboard ----
  const onNodeKey = (e: React.KeyboardEvent, n: SimNode) => {
    const step = e.shiftKey ? 24 : 8;
    if (e.key === "Enter" || e.key === " ") { e.preventDefault(); togglePin(n); return; }
    const d: Record<string, [number, number]> = { ArrowLeft: [-step, 0], ArrowRight: [step, 0], ArrowUp: [0, -step], ArrowDown: [0, step] };
    if (d[e.key]) {
      e.preventDefault();
      n.fx = (n.x || 0) + d[e.key][0]; n.fy = (n.y || 0) + d[e.key][1]; n.pinned = true;
      if (reduced) { simRef.current?.tick(); } else simRef.current?.alpha(0.2).restart();
      bump();
    }
  };

  // ---- pan / zoom ----
  const onBgPointerDown = (e: React.PointerEvent) => {
    if (dragRef.current) return;
    clearSel(); // pressing empty space deselects any edge/path focus
    onSelect?.(null);
    panRef.current = { x: e.clientX, y: e.clientY, tx: tf.x, ty: tf.y };
  };
  const onBgPointerMove = (e: React.PointerEvent) => {
    const p = panRef.current;
    if (!p) return; // capture up front: the pan may end before the state updater runs
    const nx = p.tx + (e.clientX - p.x), ny = p.ty + (e.clientY - p.y);
    setTf((t) => clampTf({ ...t, x: nx, y: ny }));
  };
  const endPan = () => (panRef.current = null);
  // keep the latest clamp available to the (empty-deps) native wheel listener
  const clampRef = useRef(clampTf);
  clampRef.current = clampTf;
  // native, non-passive wheel listener (React's onWheel is passive → preventDefault would no-op)
  useEffect(() => {
    const svg = svgRef.current;
    if (!svg) return;
    const onWheel = (e: WheelEvent) => {
      e.preventDefault();
      const r = svg.getBoundingClientRect();
      const mx = e.clientX - r.left, my = e.clientY - r.top;
      const f = e.deltaY < 0 ? 1.12 : 1 / 1.12;
      setTf((t) => {
        const k = Math.max(0.25, Math.min(6, t.k * f));
        const s = k / t.k;
        return clampRef.current({ k, x: mx - (mx - t.x) * s, y: my - (my - t.y) * s });
      });
    };
    svg.addEventListener("wheel", onWheel, { passive: false });
    return () => svg.removeEventListener("wheel", onWheel);
  }, []);

  useEffect(() => {
    const onKey = (e: KeyboardEvent) => { if (e.key === "Escape") { clearSel(); unpinAll(); onSelect?.(null); } };
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [clearSel, unpinAll, onSelect]);

  // ---- adjacency + hover highlight (complete gene→program→trait paths) ----
  const adj = useMemo(() => {
    const m = new Map<string, { e: SimEdge; other: string; otype: string }[]>();
    edgesRef.current.forEach((e) => {
      const s = nid(e.source), t = nid(e.target);
      const sn = nodesRef.current.find((n) => n.id === s);
      const tn = nodesRef.current.find((n) => n.id === t);
      if (!sn || !tn) return;
      (m.get(s) ?? m.set(s, []).get(s)!).push({ e, other: t, otype: tn.type });
      (m.get(t) ?? m.set(t, []).get(t)!).push({ e, other: s, otype: sn.type });
    });
    return m;
    // builtN bumps once the node/edge refs are populated, so this recomputes on real data (not empty refs)
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [groupKey, builtN]);

  // Highlight priority: an edge selection locks a hide-focus; otherwise the persistent union of the
  // PINNED nodes' paths (stays after the pointer leaves) plus a transient hover overlay on top.
  const hi = useMemo(() => {
    const byId = (id: string) => nodesRef.current.find((n) => n.id === id);
    const nb = (id: string, otype?: string) => (adj.get(id) || []).filter((a) => !otype || a.otype === otype);

    // a node → the complete gene→program→trait paths through it
    const pathOf = (id: string): { nodes: Set<string>; edges: Set<SimEdge> } => {
      const nodes = new Set<string>([id]);
      const edges = new Set<SimEdge>();
      const self = byId(id);
      if (self?.type === "program") nb(id).forEach((a) => { nodes.add(a.other); edges.add(a.e); });
      else if (self?.type === "gene")
        nb(id, "program").forEach((a) => { nodes.add(a.other); edges.add(a.e); nb(a.other, "trait").forEach((b) => { nodes.add(b.other); edges.add(b.e); }); });
      else
        nb(id, "program").forEach((a) => { nodes.add(a.other); edges.add(a.e); nb(a.other, "gene").forEach((b) => { nodes.add(b.other); edges.add(b.e); }); });
      return { nodes, edges };
    };
    // an edge → its two endpoints and everything directly connected to either of them
    const edgeSetOf = (e: SimEdge) => {
      const nodes = new Set<string>();
      const edges = new Set<SimEdge>([e]);
      [nid(e.source), nid(e.target)].forEach((id) => {
        nodes.add(id);
        (adj.get(id) || []).forEach((a) => { nodes.add(a.other); edges.add(a.e); });
      });
      return { nodes, edges };
    };

    if (edgeSel) return { ...edgeSetOf(edgeSel), hide: true as const };

    const N = new Set<string>();
    const E = new Set<SimEdge>();
    let any = false;
    pinned.forEach((id) => { if (byId(id)) { any = true; const s = pathOf(id); s.nodes.forEach((x) => N.add(x)); s.edges.forEach((x) => E.add(x)); } });
    if (hoverEdge) { any = true; const s = edgeSetOf(hoverEdge); s.nodes.forEach((x) => N.add(x)); s.edges.forEach((x) => E.add(x)); }
    else if (hover) { any = true; const s = pathOf(hover); s.nodes.forEach((x) => N.add(x)); s.edges.forEach((x) => E.add(x)); }
    return any ? { nodes: N, edges: E, hide: false as const } : null;
  }, [hover, hoverEdge, edgeSel, pinned, adj]);

  // ---- imperative API ----
  useImperativeHandle(ref, () => ({
    regenerate() { seedPositions(nodesRef.current, dims.w, dims.h, () => true); settle(0.8); setTimeout(fit, reduced ? 0 : 700); },
    reset() { unpinAll(); seedPositions(nodesRef.current, dims.w, dims.h, () => true); setTf({ k: 1, x: 0, y: 0 }); settle(0.9); setHover(null); clearSel(); onSelect?.(null); setTimeout(fit, reduced ? 0 : 700); },
    start,
    fit,
    focus(id: string) {
      const n = nodesRef.current.find((x) => x.id.toUpperCase() === id.toUpperCase());
      if (!n) return false;
      setHover(n.id); onSelect?.(n);
      const k = Math.max(tf.k, 1.4);
      setTf({ k, x: dims.w / 2 - k * (n.x || 0), y: dims.h / 2 - k * (n.y || 0) });
      return true;
    },
    exportSVG() { downloadSVG(svgRef.current!, `proglof-network-${groupKey}.svg`); },
    exportPNG() { downloadPNG(svgRef.current!, dims, `proglof-network-${groupKey}.png`); },
    ids() { return nodesRef.current.map((n) => n.id); },
  }));

  // ---- render ----
  const nodes = nodesRef.current;
  const edges = edgesRef.current;
  const edgeOpacity = (e: SimEdge) => {
    const base = e.cls === "member" ? 0.4 : 0.7;
    if (draw != null) {
      // intro reveal: gene→program edges connect first, program→trait after
      const ph = e.cls === "selects"
        ? Math.max(0, Math.min(1, (draw - 0.4) / 0.55))
        : Math.max(0, Math.min(1, draw / 0.6));
      return base * ph;
    }
    if (!hi) return base;
    if (hi.edges.has(e)) return e.cls === "member" ? 0.55 : 0.9;
    return hi.hide ? 0 : 0.06; // click-select hides others; hover just dims
  };
  const nodeOpacity = (n: SimNode, faded: boolean) => {
    const base = faded ? 0.5 : 1;
    if (draw != null) return base * Math.min(1, draw * 4); // pop in quickly, then fly home
    if (!hi) return base;
    if (hi.nodes.has(n.id)) return base;
    return hi.hide ? 0.05 : 0.12;
  };

  return (
    <div className="graphwrap" ref={wrapRef}>
      <svg
        ref={svgRef}
        width={dims.w}
        height={dims.h}
        role="application"
        aria-label={`Gene to program to trait network, ${group.label}. ${nodes.length} nodes.`}
        onPointerDown={onBgPointerDown}
        onPointerMove={onBgPointerMove}
        onPointerUp={endPan}
        onPointerLeave={endPan}
        style={{ cursor: panRef.current ? "grabbing" : "grab" }}
      >
        <g transform={`translate(${tf.x},${tf.y}) scale(${tf.k})`}>
          <g className="edges">
            {edges.map((e, i) => {
              const s = e.source as SimNode, t = e.target as SimNode;
              if (s.x == null || t.x == null) return null;
              const d = edgePath(s.x, s.y!, t.x, t.y!, e.poff || 0);
              const op = edgeOpacity(e);
              return (
                <g key={i} style={{ cursor: "pointer" }}
                  onMouseEnter={() => setHoverEdge(e)}
                  onMouseLeave={() => setHoverEdge((h) => (h === e ? null : h))}
                  onPointerDown={(ev) => ev.stopPropagation()}
                  onClick={(ev) => { ev.stopPropagation(); onSelect?.(null); setEdgeSel(e); }}
                >
                  {/* fat transparent hit area so thin edges are easy to hover/click */}
                  <path d={d} fill="none" stroke="transparent" strokeWidth={Math.max(10, edgeWidth(e) * 3)} pointerEvents={op === 0 ? "none" : "stroke"} />
                  <path
                    d={d}
                    fill="none"
                    stroke={edgeColor(e)}
                    strokeWidth={edgeWidth(e) + (edgeSel === e ? 1.4 : 0)}
                    strokeOpacity={op}
                    strokeDasharray={e.cls === "member" ? "2 3" : undefined}
                    strokeLinecap="round"
                    pointerEvents="none"
                  />
                  {e.cls === "selects" && e.label && op > 0.25 && draw == null &&
                    (() => {
                      const mx = (s.x! + t.x!) / 2, my = (s.y! + t.y!) / 2;
                      return (
                        <text x={mx} y={my} textAnchor="middle" fontSize={8.5} fill={pal.ink}
                          stroke={pal.ring} strokeWidth={2.5} paintOrder="stroke" pointerEvents="none">
                          {e.label.split("\n").map((ln, li) => (
                            <tspan key={li} x={mx} dy={li === 0 ? -3 : 10}>{ln}</tspan>
                          ))}
                        </text>
                      );
                    })()}
                </g>
              );
            })}
          </g>
          <g className="nodes">
            {nodes.map((n) => {
              if (n.x == null) return null;
              const r = nodeRadius(n);
              const faded = n.type === "gene" && n.role !== "regulator";
              const op = nodeOpacity(n, faded);
              const aria =
                n.type === "program"
                  ? `Regulator program ${n.id}${n.annotation ? ", " + n.annotation : ""}${n.wP != null ? `, weight ${n.wP.toFixed(3)}` : ""}`
                  : n.type === "trait"
                  ? `Trait ${n.id}`
                  : `regulator gene ${n.id}, own gamma sign ${n.sign! > 0 ? "positive" : n.sign! < 0 ? "negative" : "small"}`;
              return (
                <g
                  key={n.id}
                  transform={`translate(${n.x},${n.y})`}
                  opacity={op}
                  tabIndex={0}
                  role="button"
                  aria-label={aria}
                  style={{ cursor: "pointer", outline: "none" }}
                  onPointerDown={(e) => onNodePointerDown(e, n)}
                  onPointerMove={(e) => onNodePointerMove(e, n)}
                  onPointerUp={(e) => onNodePointerUp(e, n)}
                  onKeyDown={(e) => onNodeKey(e, n)}
                  onMouseEnter={() => setHover(n.id)}
                  onMouseLeave={() => setHover((h) => (h === n.id ? null : h))}
                  onFocus={() => setHover(n.id)}
                  onBlur={() => setHover((h) => (h === n.id ? null : h))}
                >
                  {n.type === "program" ? (() => {
                    const lines = [n.id, ...wrapLabel(n.annotation || "")];
                    const bw = Math.max(64, ...lines.map((l) => l.length * 5.7)) + 12;
                    const bh = lines.length * 11 + 8;
                    return (
                      <>
                        {n.pinned && <rect x={-bw / 2 - 3} y={-bh / 2 - 3} width={bw + 6} height={bh + 6} rx={7} fill="none" stroke={pal.pin} strokeWidth={1.5} strokeDasharray="2 2" />}
                        <rect x={-bw / 2} y={-bh / 2} width={bw} height={bh} rx={6} fill={pal.pbox}
                          stroke={n.sig ? UP : pal.memStroke} strokeWidth={n.sig ? 2 : 1} />
                        <text textAnchor="middle" fill={pal.ink} fontSize={9.5} y={-bh / 2 + 12}>
                          {lines.map((ln, li) => (
                            <tspan key={li} x={0} dy={li === 0 ? 0 : 11} fontWeight={li === 0 ? 700 : 400}>{ln}</tspan>
                          ))}
                        </text>
                      </>
                    );
                  })() : n.type === "trait" ? (() => {
                    const bw = 16 + n.id.length * 7;
                    return (
                      <>
                        {n.pinned && <rect x={-3} y={-11} width={bw + 6} height={22} rx={9} fill="none" stroke={pal.pin} strokeWidth={1.5} strokeDasharray="2 2" />}
                        <rect x={0} y={-9} width={bw} height={18} rx={9} fill={pal.tfill} stroke={pal.memStroke} strokeWidth={1} />
                        <text x={bw / 2} y={4} textAnchor="middle" fill={pal.ink} fontSize={11} fontWeight={700}>{n.id}</text>
                      </>
                    );
                  })() : (
                    <>
                      {n.pinned && <circle r={r + 4} fill="none" stroke={pal.pin} strokeWidth={1.5} strokeDasharray="2 2" />}
                      <circle r={r} fill={GFILL[n.gclass || "small"]} stroke="#5b6472" strokeWidth={1} />
                      <text x={-r - 3} y={3} textAnchor="end" fill={gLabel(n.sign)} fontSize={10} fontWeight={600}>{n.id}</text>
                    </>
                  )}
                </g>
              );
            })}
          </g>
        </g>
      </svg>
      {edgeSel &&
        (() => {
          const s = edgeSummary(edgeSel);
          return (
            <aside className="info edge" style={{ borderLeft: `3px solid ${s.accent}` }}>
              <button className="x" aria-label="close" onClick={clearSel}>×</button>
              <b>{s.title}</b>
              {s.sub && <span className="muted"> · {s.sub}</span>}
              <div className="rows">
                {s.rows.map(([k, v]) => (
                  <div key={k}><span className="muted">{k}</span><span>{v}</span></div>
                ))}
              </div>
              <div className="muted small">showing only this edge's connections · Esc / click background to clear</div>
            </aside>
          );
        })()}
    </div>
  );
});

const fmt = (x?: number | null) => (x == null ? "—" : Number(x).toPrecision(3));

function edgeSummary(e: SimEdge): { title: string; sub: string; accent: string; rows: [string, string][] } {
  const a = e.source as SimNode, b = e.target as SimNode;
  const pick = (t: string) => (a.type === t ? a : b);
  if (e.cls === "selects") {
    const prog = pick("program"), trait = a.type === "trait" ? a : b;
    const inc = (e.dir || 0) > 0;
    return {
      title: `${prog.id} → ${trait.id}`, sub: prog.annotation || "", accent: inc ? UP : DOWN,
      rows: [
        ["relation", "program → trait (regulator model)"],
        ["effect on trait γ", inc ? "increases (+)" : "decreases (−)"],
        ["w_P (joint weight)", fmt(e.wP)],
        ["P", fmt(e.wP_P)],
      ],
    };
  }
  if (e.cls === "regulates") {
    const gene = pick("gene"), prog = a.type === "program" ? a : b;
    const up = (e.beta || 0) > 0;
    return {
      title: `${gene.id} → ${prog.id}`, sub: prog.annotation || "", accent: up ? UP : DOWN,
      rows: [
        ["relation", "gene regulates program"],
        ["β_x(P) (Stage-B)", fmt(e.beta)],
        ["direction", up ? "up-regulates (↑)" : "down-regulates (↓)"],
      ],
    };
  }
  const gene = pick("gene"), prog = a.type === "program" ? a : b;
  return {
    title: `${gene.id} ∈ ${prog.id}`, sub: prog.annotation || "", accent: MEM,
    rows: [["relation", "gene is a member of program"]],
  };
}

// ---- export helpers ----
function serialize(svg: SVGSVGElement): string {
  const clone = svg.cloneNode(true) as SVGSVGElement;
  clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
  return '<?xml version="1.0" encoding="UTF-8"?>\n' + new XMLSerializer().serializeToString(clone);
}
function triggerDownload(href: string, name: string) {
  const a = document.createElement("a");
  a.href = href; a.download = name; a.click();
}
function downloadSVG(svg: SVGSVGElement, name: string) {
  const blob = new Blob([serialize(svg)], { type: "image/svg+xml" });
  const url = URL.createObjectURL(blob);
  triggerDownload(url, name);
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}
function downloadPNG(svg: SVGSVGElement, dims: { w: number; h: number }, name: string) {
  const scale = 2;
  const img = new Image();
  const svgUrl = "data:image/svg+xml;charset=utf-8," + encodeURIComponent(serialize(svg));
  img.onload = () => {
    const c = document.createElement("canvas");
    c.width = dims.w * scale; c.height = dims.h * scale;
    const ctx = c.getContext("2d")!;
    ctx.scale(scale, scale);
    ctx.drawImage(img, 0, 0);
    c.toBlob((b) => { if (b) { const u = URL.createObjectURL(b); triggerDownload(u, name); setTimeout(() => URL.revokeObjectURL(u), 1000); } }, "image/png");
  };
  img.src = svgUrl;
}
