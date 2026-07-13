# PLoTM network — interactive gene → program → trait graph

A physics-based, three-layer network view of the PLoTM results, built with **React + TypeScript +
D3-force + SVG**. It replaces the static combined Fig-5a maps with a spacious, draggable, zoomable graph:

```
regulator / member genes   →   cNMF programs   →   LoF phenotypes (traits)
        (left)                    (centre)                (right)
```

Data comes from the explorer's `demo_data.json` (exported to `src/networks.json`), so this app is fully
self-contained — no server needed.

## Features

- **Three preserved layers** (genes L / programs C / traits R) via per-layer `forceX`, with generous
  horizontal (`spread`) and vertical spacing.
- **Physics**: `d3-force` with link springs, **hub repulsion** (charge scales with node degree),
  **label-aware collision** (`forceCollide` padded by node role so labels don't overlap), and a
  responsive canvas (`ResizeObserver`).
- **Drag** any node; connected Bézier edges update every frame. On release an unpinned node **springs
  back with a slight overshoot** (its drop velocity + the layer spring damp out) instead of snapping.
- **Edges** are smooth cubic Béziers, parallel edges are bowed apart, and width scales with effect
  strength (|β| for *regulates*, −log10 P for *selects*).
- **Direction-coded trait→program edges**: **orange = increases** disease risk, **teal = decreases**
  (sign of the program-burden effect for that phenotype).
- **Hover** any node to highlight the complete gene→program→trait path(s) through it; everything else dims.
- **Zoom / pan** (wheel + drag background), **pin / unpin** (click a node), **search**, **reset**,
  **regenerate layout** (redistributes unpinned nodes, keeps layers + pins), **spacing controls**, and
  **SVG / PNG export**.
- **Stabilises and stops** (alpha decays to rest) — no permanent jitter.
- **Accessibility**: honours `prefers-reduced-motion` (runs the layout synchronously, no animation),
  every node is keyboard-focusable (`Tab`), arrow keys nudge, `Enter`/`Space` pins, and nodes carry ARIA
  labels. Day / night theme toggle.

## Run

```bash
cd PLoTM_network
npm install
npm run export     # (re)generate src/networks.json from ../PLoTM_explorer/data/demo_data.json
npm run dev        # http://127.0.0.1:5178   (hot reload)
# or a production build:
npm run build && npm run preview
```

`npm run build` runs `tsc -b` (full type-check) then Vite. Output in `dist/` is static and works behind
any server or from `file://` (the Vite `base` is `./`).

## Layout

```
scripts/export_network_data.py   demo_data.json → src/networks.json (bundled at build)
src/
  types.ts        graph + simulation types
  viz.ts          layers, colours (incl. increase/decrease), radii, Bézier path, edge width
  sim.ts          d3-force graph build + simulation (layers, hub repulsion, collision)
  Graph.tsx       SVG render, zoom/pan, drag + spring-back, hover path highlight, keyboard, export
  App.tsx         toolbar (condition/group/search/reset/regenerate/spacing/export/theme) + state
  Legend.tsx      colour/interaction key
  styles.css      theming (day/night) + layout
```

## Regenerating the data

The graph reflects whatever is in `../PLoTM_explorer/data/demo_data.json`. After re-running the
explorer's `build_demo_data.py`, run `npm run export` here to refresh `src/networks.json`, then rebuild.
