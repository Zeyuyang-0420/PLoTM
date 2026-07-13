import { UP, DOWN, GFILL } from "./viz";

export function Legend() {
  return (
    <div className="legend" aria-label="legend">
      <div className="row"><b>nodes</b>
        <span><i className="c box" />regulator program (w/P to trait)</span>
        <span><i className="c" style={{ background: GFILL.pos }} />regulator γ &gt; 0.03</span>
        <span><i className="c" style={{ background: GFILL.neg }} />γ &lt; −0.03</span>
        <span><i className="c" style={{ background: GFILL.small }} />|γ| small</span>
        <span><i className="c pill" />trait</span>
        <span style={{ color: "var(--muted)" }}>gene label colour = sign(own γ)</span>
      </div>
      <div className="row"><b>edges</b>
        <span><i className="l" style={{ background: UP }} />gene→program: up-regulates (β&gt;0)</span>
        <span><i className="l" style={{ background: DOWN }} />down-regulates (β&lt;0)</span>
        <span><i className="l" style={{ background: UP }} />program→trait: increases γ (w&gt;0)</span>
        <span><i className="l" style={{ background: DOWN }} />decreases γ (w&lt;0)</span>
      </div>
      <div className="hint">width: gene→program ~ |β| · program→trait ~ |w_P| (from γ ~ Σ w_P·β + λ·shet)</div>
    </div>
  );
}
