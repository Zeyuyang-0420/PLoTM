import type React from "react";
import { useEffect, useMemo, useRef, useState } from "react";
import networks from "./networks.json";
import type { NetworksFile, NetworkGroup, Spacing, SimNode } from "./types";
import { Graph, type GraphHandle } from "./Graph";
import { Legend } from "./Legend";
import { ErrorBoundary } from "./ErrorBoundary";

const DATA = networks as unknown as NetworksFile;
const DEFAULT_SPACING: Spacing = { spread: 1.15, charge: -200, linkDist: 95, collide: 8 };

// every trait that appears anywhere in the data, in stable order
const ALL_TRAITS: string[] = (() => {
  const s: string[] = [];
  Object.values(DATA.networks).forEach((byG) =>
    Object.values(byG).forEach((gr) =>
      gr.nodes.forEach((n) => { if (n.type === "trait" && !s.includes(n.id)) s.push(n.id); }),
    ),
  );
  return s;
})();
const traitShort = (id: string) => id.replace(/^GB_/, "");
const DEFAULT_TRAITS = (() => {
  const c = ALL_TRAITS.filter((t) => /custom/i.test(t));
  return c.length ? c : ALL_TRAITS.slice(0, 1);
})();

// prune the combined graph to the subgraph reachable from the selected trait(s):
// selected traits → the programs they select → those programs' regulator genes.
function filterGroup(g: NetworkGroup, sel: Set<string>): NetworkGroup {
  if (!g.nodes.some((n) => n.type === "trait" && !sel.has(n.id))) return g; // nothing hidden
  const type = new Map(g.nodes.map((n) => [n.id, n.type]));
  const keepProg = new Set<string>();
  g.edges.forEach((e) => {
    if (e.cls !== "selects") return;
    const tr = type.get(e.source) === "trait" ? e.source : e.target;
    const pr = tr === e.source ? e.target : e.source;
    if (sel.has(tr)) keepProg.add(pr);
  });
  const keepGene = new Set<string>();
  g.edges.forEach((e) => {
    if (e.cls === "selects") return;
    const pr = type.get(e.source) === "program" ? e.source : type.get(e.target) === "program" ? e.target : null;
    if (!pr || !keepProg.has(pr)) return;
    keepGene.add(pr === e.source ? e.target : e.source);
  });
  const keep = new Set<string>([...sel, ...keepProg, ...keepGene]);
  return {
    ...g,
    nodes: g.nodes.filter((n) => keep.has(n.id)),
    edges: g.edges.filter((e) => keep.has(e.source) && keep.has(e.target)),
  };
}

export default function App() {
  const graph = useRef<GraphHandle>(null);
  // condition + theme are driven by the parent explorer (single source of truth) via postMessage;
  // they fall back to sensible defaults when the app is opened standalone.
  const [condition, setCondition] = useState(DATA.conditions[0]);
  const [theme, setTheme] = useState<"dark" | "light">(
    () => (localStorage.getItem("proglof-net-theme") as "dark" | "light") || "dark",
  );
  const group = DATA.groups[0].id;
  const spacing: Spacing = DEFAULT_SPACING;
  const [traitSel, setTraitSel] = useState<Set<string>>(() => new Set(DEFAULT_TRAITS));
  const [query, setQuery] = useState("");
  const [sel, setSel] = useState<SimNode | null>(null);
  const [animateStart, setAnimateStart] = useState(true);
  const [started, setStarted] = useState(false);

  // listen for condition/theme pushed from the host explorer; announce readiness on mount
  useEffect(() => {
    const onMsg = (e: MessageEvent) => {
      const d = e.data;
      if (!d || d.type !== "proglof-net") return;
      if (d.condition && DATA.conditions.includes(d.condition)) setCondition(d.condition);
      if (d.theme === "dark" || d.theme === "light") setTheme(d.theme);
    };
    window.addEventListener("message", onMsg);
    if (window.parent && window.parent !== window) window.parent.postMessage({ type: "proglof-net-ready" }, "*");
    return () => window.removeEventListener("message", onMsg);
  }, []);

  const begin = () => {
    graph.current?.start(animateStart); // scatter + restart physics while still hidden…
    setStarted(true); // …then reveal, so the first painted frame is already the assembling graph
  };

  const reduced = useMemo(
    () => window.matchMedia?.("(prefers-reduced-motion: reduce)").matches ?? false,
    [],
  );
  const gRaw = DATA.networks[condition][group];
  const g = useMemo(() => filterGroup(gRaw, traitSel), [gRaw, traitSel]);
  const groupKey = `${condition}|${group}|${[...traitSel].sort().join(",")}`;
  const ids = useMemo(() => g.nodes.map((n) => n.id), [g]);
  const toggleTrait = (t: string) =>
    setTraitSel((prev) => {
      const s = new Set(prev);
      if (s.has(t)) { if (s.size > 1) s.delete(t); } else s.add(t);
      return s;
    });

  document.documentElement.dataset.theme = theme;
  const doSearch = (e?: React.FormEvent) => {
    e?.preventDefault();
    if (query.trim()) graph.current?.focus(query.trim());
  };

  return (
    <div className="app">
      <header className="bar">
        <div className="brand">
          <span className="dot" />
          <div>
            <h1>PLoTM network</h1>
            <p>gene → program → trait · {DATA.name}</p>
          </div>
        </div>

        {ALL_TRAITS.length > 1 && (
          <div className="seg" role="group" aria-label="traits" title="toggle traits (multi-select)">
            {ALL_TRAITS.map((t) => (
              <button key={t} className={traitSel.has(t) ? "on" : ""} aria-pressed={traitSel.has(t)} onClick={() => toggleTrait(t)}>
                {traitShort(t)}
              </button>
            ))}
          </div>
        )}

        <form className="search" onSubmit={doSearch}>
          <input
            list="node-ids"
            placeholder="find gene / program…"
            value={query}
            onChange={(e) => setQuery(e.target.value)}
            aria-label="search node"
          />
          <datalist id="node-ids">{ids.map((i) => <option key={i} value={i} />)}</datalist>
          <button type="submit">go</button>
        </form>

        <div className="actions">
          <button onClick={begin} className="go" title="assemble the network from a random scatter">▶ {started ? "replay" : "start"}</button>
          <label className="chk" title="play the physics reveal on start">
            <input type="checkbox" checked={animateStart} onChange={(e) => setAnimateStart(e.target.checked)} /> animate
          </label>
          <button onClick={() => graph.current?.regenerate()} title="redistribute unpinned nodes">↻ regenerate</button>
          <button onClick={() => graph.current?.fit()} title="fit to view">⤢ fit</button>
          <button onClick={() => graph.current?.reset()} title="unpin all + reset view">reset</button>
        </div>
      </header>

      <main>
        <ErrorBoundary onReset={() => graph.current?.reset()}>
          <Graph ref={graph} group={g} groupKey={groupKey} spacing={spacing} theme={theme} reduced={reduced} started={started} animate={animateStart} onSelect={setSel} />
        </ErrorBoundary>
        {!started && (
          <div className="startgate">
            <button className="startbtn" onClick={begin} autoFocus>▶ Start</button>
            <label className="chk">
              <input type="checkbox" checked={animateStart} onChange={(e) => setAnimateStart(e.target.checked)} /> animate the physics assembly
            </label>
            <p className="ghint">assemble the gene → program → trait network</p>
          </div>
        )}
        {started && <Legend />}
        {sel && (
          <aside className="info">
            <button className="x" aria-label="close" onClick={() => setSel(null)}>×</button>
            <b>{sel.id}</b> <span className="muted">{sel.type === "gene" ? `${sel.role} gene` : sel.type}</span>
            {sel.annotation && <div>{sel.annotation}</div>}
            {sel.type !== "trait" && <div className="muted">in {sel.n_traits} phenotype(s){sel.type === "gene" ? ` · predicted γ ${sel.sign! > 0 ? "+" : "−"}` : sel.sig ? " · meta-significant" : ""}</div>}
            {sel.pinned && <div className="muted">📌 pinned — click again to release</div>}
          </aside>
        )}
      </main>
    </div>
  );
}
