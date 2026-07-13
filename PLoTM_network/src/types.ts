import type { SimulationNodeDatum, SimulationLinkDatum } from "d3-force";

export type NodeType = "gene" | "program" | "trait";
export type EdgeClass = "regulates" | "member" | "selects";

export interface RawNode {
  id: string;
  type: NodeType;
  role?: string; // gene: "regulator" | "member"
  sign?: number; // gene: sign of its own LoF γ (+1 / −1) — sets the label colour
  gamma?: number | null; // gene: its LoF γ for this trait
  gclass?: "pos" | "neg" | "small"; // gene node fill (γ>0.03 / γ<−0.03 / |γ| small)
  n_traits?: number;
  annotation?: string; // program label
  annotated?: boolean;
  sig?: boolean; // program: w_P significant (P<0.05)
  wP?: number | null; // program: joint-regression weight to the trait
  wP_P?: number | null;
  traits?: string[];
}

export interface RawEdge {
  source: string;
  target: string;
  cls: EdgeClass;
  beta?: number | null; // regulates: Stage-B β (gene→program)
  dir?: number; // sign of the effect: regulates=sign(β), selects=sign(w_P)
  wP?: number | null; // selects: program→trait weight
  wP_P?: number | null;
  label?: string; // selects: "w=…\nP=…"
  effect?: number | null;
  weight?: number;
  n_traits?: number;
  traits?: string[];
}

export interface NetworkGroup {
  label: string;
  nodes: RawNode[];
  edges: RawEdge[];
}

export interface NetworksFile {
  disease: string;
  name: string;
  cell_type: string;
  conditions: string[];
  groups: { id: string; label: string }[];
  networks: Record<string, Record<string, NetworkGroup>>; // [condition][group]
}

export interface SimNode extends RawNode, SimulationNodeDatum {
  x?: number;
  y?: number;
  vx?: number;
  vy?: number;
  fx?: number | null;
  fy?: number | null;
  pinned?: boolean;
  degree?: number;
}

export interface SimEdge extends Omit<RawEdge, "source" | "target">, SimulationLinkDatum<SimNode> {
  source: string | SimNode;
  target: string | SimNode;
  /** parallel-edge offset index, assigned per unordered endpoint pair */
  poff?: number;
}

export interface Spacing {
  spread: number; // horizontal layer separation multiplier
  charge: number; // base repulsion (negative)
  linkDist: number; // link target length
  collide: number; // extra collision padding
}
