# Pipeline phases, outputs, and verify gates

Driven by `proglof-mcp` `run_phase(config, <phase>)` → `drivers/run_disease.sh`. Order matters; each
phase depends on the previous. Heavy phases (`cnmf`, `stageB`, `burden`, `figures`, `posterior`) route to
the SSH cluster if `PROGLOF_REMOTE` is set. All are resumable.

| phase | produces | verify (via `read_results` / artifacts) before continuing |
|---|---|---|
| `download` | perturb-seq h5ad + LoF burden summary stats (both sources) | files present; burden tables non-empty |
| `posterior` | GeneBayes γ per trait + QC flags | **`qc_report`**: how many traits flagged? scales sane (GeneBass ~0.005, Backman ~0.3)? |
| `input` | cNMF input h5ad + covariate metadata per condition | one CT per condition; raw counts preserved |
| `cnmf` | 60 programs + usages per condition | consensus files exist; k-selection/consensus plots stable |
| `stageB` | 60 `K60_program*_perturb_effects.txt` per condition | 60 files/condition; self-KDs sane |
| `thresholds` | per-trait q99 |γ| cutoff + flags | `read_results thresholds`: q99 per trait; flags marked |
| `burden` | program & regulator burden per (trait, condition) → Fig 4C inputs | both enrichment tables finite |
| `figures` | Fig 4C + Fig 5a (obs + perm) per (trait, condition) | Fig 4C Bonferroni hits; Fig 5a **permutation** P |
| `annotate` | program labels (GO/Hallmark/markers) per condition | labels attach to the top programs |
| `analysis` | tidy table, cross-source meta, condition interaction, networks, permutation | `read_results meta`: I²/LOO stability; genes-vs-programs overlap permutation |

## The key numbers to read at the end

- `read_results(kind="meta")` — pooled β, P, **FDR**, **I²**, `loo_sign_stable`, `cross_source_concordant`
  per program×condition. Robust = FDR<0.05, I²≈0, LOO-stable, both sources agree.
- `read_results(kind="tidy")` — the full 60×traits×conditions grid (program- and regulator-burden).
- `read_results(kind="summary")` — per (trait, condition): Fig 4C hit + Fig 5a **permutation P** + QC flag.
- `qc_report` — which traits to exclude and why.

## Parameters (config `params`)

`K: 60`, `lof_threshold: q99`, `stageB_min_profiles: 4`, `fig5: {PN:5, RN:3, TOP:200, NPERM:5000,
NPERM_OBS:10000}`. Change K only if the substrate has few pseudobulk samples (SLE used K=30 for 261
donors — size K to the sample count, and matching is by projection, not by number).
