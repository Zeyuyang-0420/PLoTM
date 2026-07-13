import type { SimNode, SimEdge, NodeType } from "./types";

// layer positions as fractions of width (before spread); genes left, programs centre, traits right
export const LAYER_FRAC: Record<NodeType, number> = { gene: 0.13, program: 0.5, trait: 0.88 };

// conservation colour ramp (kept for datasets that colour by #phenotypes)
export const CONS = ["#7f93b0", "#5f86c9", "#4f79cf", "#3f6fda", "#2f66e6"];
export const consColor = (n = 1) => CONS[Math.min(Math.max(n, 1) - 1, CONS.length - 1)];

// "Regulator model" palette (matches fig5a_regulators_RAprog.R):
//   gene→program  : sign(beta)  up-regulates (red) / down-regulates (blue)
//   program→trait : sign(w_P)   increases (red) / decreases (blue) trait gamma
export const UP = "#c0392b"; // positive: up-regulates / increases trait gamma
export const DOWN = "#2c6fbb"; // negative: down-regulates / decreases trait gamma
export const MEM = "#7c8798"; // member (dashed; other datasets)
// regulator-gene node fill by its own LoF gamma; label colour by sign(gamma)
export const GFILL: Record<string, string> = { pos: "#e8746a", neg: "#6a9fe8", small: "#c7ccd4" };
export const gLabel = (sign?: number) => ((sign || 0) > 0 ? "#c0392b" : (sign || 0) < 0 ? "#2c6fbb" : "#8a94a3");

export function nodeRadius(n: SimNode): number {
  if (n.type === "program") return 10; // rendered as a box; used for collision + edge anchoring
  if (n.type === "trait") return 9;
  return 4;
}

// label-aware collision padding: labelled nodes need more vertical room so text doesn't collide
export function collidePad(n: SimNode): number {
  if (n.type === "program") return 20; // two-line boxes
  if (n.type === "trait") return 10;
  return 7;
}

export function edgeColor(e: SimEdge): string {
  if (e.cls === "member") return MEM;
  return (e.dir || 0) > 0 ? UP : DOWN; // regulates (sign beta) & selects (sign w_P) share the scheme
}

export function edgeWidth(e: SimEdge): number {
  if (e.cls === "regulates") return clamp(0.4 + Math.abs(e.beta || 1) * 0.4, 0.7, 4.5);
  if (e.cls === "selects") return clamp(0.9 + Math.abs(e.wP || 0) * 45, 1, 5); // ~ |w_P|
  return 0.8;
}

export const clamp = (v: number, lo: number, hi: number) => Math.max(lo, Math.min(hi, v));

/** smooth horizontal cubic Bézier between two points, with a parallel-offset bow */
export function edgePath(ax: number, ay: number, bx: number, by: number, off = 0): string {
  const mx = (ax + bx) / 2;
  const bow = off * 16;
  return `M${ax.toFixed(1)},${ay.toFixed(1)} C${mx.toFixed(1)},${(ay + bow).toFixed(1)} ${mx.toFixed(1)},${(by + bow).toFixed(1)} ${bx.toFixed(1)},${by.toFixed(1)}`;
}

export const nid = (x: string | SimNode): string => (typeof x === "string" ? x : x.id);
