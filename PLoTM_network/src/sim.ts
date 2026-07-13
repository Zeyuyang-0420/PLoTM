import {
  forceSimulation,
  forceLink,
  forceManyBody,
  forceX,
  forceY,
  forceCollide,
  type Simulation,
} from "d3-force";
import type { SimNode, SimEdge, Spacing, NetworkGroup } from "./types";
import { LAYER_FRAC, nodeRadius, collidePad, nid } from "./viz";

export interface BuiltGraph {
  nodes: SimNode[];
  edges: SimEdge[];
}

/** Clone raw graph data into fresh mutable sim objects and seed positions per layer. */
export function buildGraph(g: NetworkGroup, width: number, height: number): BuiltGraph {
  const nodes: SimNode[] = g.nodes.map((n) => ({ ...n }));
  const byId = new Map(nodes.map((n) => [n.id, n]));
  // parallel-edge offsets per unordered endpoint pair
  const pairIdx = new Map<string, number>();
  const edges: SimEdge[] = g.edges
    .filter((e) => byId.has(e.source) && byId.has(e.target))
    .map((e) => {
      const key = [e.source, e.target].sort().join("|");
      const k = pairIdx.get(key) ?? 0;
      pairIdx.set(key, k + 1);
      return { ...e, poff: k };
    });
  seedPositions(nodes, width, height, () => true);
  return { nodes, edges };
}

/** Place nodes into their layer column with a spread of y; `reseed(n)` decides which get new y. */
export function seedPositions(
  nodes: SimNode[],
  width: number,
  height: number,
  reseed: (n: SimNode) => boolean,
) {
  const cols: Record<string, SimNode[]> = { gene: [], program: [], trait: [] };
  nodes.forEach((n) => cols[n.type]?.push(n));
  for (const type of Object.keys(cols)) {
    const arr = cols[type];
    arr.sort((a, b) => (b.n_traits || 1) - (a.n_traits || 1) || a.id.localeCompare(b.id));
    arr.forEach((n, i) => {
      if (!reseed(n) || n.pinned) return;
      n.x = LAYER_FRAC[n.type] * width + (Math.random() - 0.5) * 10;
      const t = arr.length <= 1 ? 0.5 : i / (arr.length - 1);
      n.y = 40 + t * (height - 80) + (Math.random() - 0.5) * 12;
      n.vx = 0;
      n.vy = 0;
    });
  }
}

/** Scatter every (unpinned) node to a random position with a little kick — the "explode"
 *  starting state the physics simulation then pulls back into the tri-partite layout. */
export function seedRandom(nodes: SimNode[], width: number, height: number) {
  nodes.forEach((n) => {
    if (n.pinned) return;
    n.x = width / 2 + (Math.random() - 0.5) * width * 0.92;
    n.y = height / 2 + (Math.random() - 0.5) * height * 0.92;
    n.fx = null;
    n.fy = null;
    n.vx = (Math.random() - 0.5) * 3;
    n.vy = (Math.random() - 0.5) * 3;
  });
}

export function layerX(type: string, width: number, spread: number): number {
  const f = LAYER_FRAC[type as keyof typeof LAYER_FRAC] ?? 0.5;
  return (0.5 + (f - 0.5) * spread) * width;
}

export function buildSimulation(
  nodes: SimNode[],
  edges: SimEdge[],
  width: number,
  height: number,
  sp: Spacing,
  onTick: () => void,
): Simulation<SimNode, SimEdge> {
  const deg: Record<string, number> = {};
  edges.forEach((e) => {
    deg[nid(e.source)] = (deg[nid(e.source)] || 0) + 1;
    deg[nid(e.target)] = (deg[nid(e.target)] || 0) + 1;
  });
  nodes.forEach((n) => (n.degree = deg[n.id] || 0));

  return forceSimulation<SimNode>(nodes)
    .force(
      "link",
      forceLink<SimNode, SimEdge>(edges)
        .id((d) => d.id)
        .distance(sp.linkDist)
        .strength(0.22),
    )
    // hub repulsion: nodes with more connections push harder
    .force(
      "charge",
      forceManyBody<SimNode>().strength((d) => sp.charge * (1 + (d.degree || 0) * 0.35)),
    )
    .force("x", forceX<SimNode>((d) => layerX(d.type, width, sp.spread)).strength(0.92))
    .force("y", forceY<SimNode>(height / 2).strength(0.045))
    .force(
      "collide",
      forceCollide<SimNode>((d) => nodeRadius(d) + collidePad(d) + sp.collide).iterations(2),
    )
    .alphaDecay(0.03)
    .velocityDecay(0.42)
    .on("tick", onTick);
}
