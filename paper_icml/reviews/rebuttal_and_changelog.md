# Reviewer Response & Revision Changelog

This document records every critique raised by the three independent submodel reviewers (`reviewer_1_pro.md`, `reviewer_2_flash.md`, `reviewer_3_flash_lite.md`), our upgraded statistical re-analyses, manuscript revisions in `paper_icml/main.tex`, and codebase fixes/ablation tooling prepared in the repository.

---

## 1. Mixed-Effects Statistical Model Upgrades (Addressing Reviewer 1 W2 & Reviewer 3 Major Issue 1)

| Critique | Action Taken | Verified Result |
| :--- | :--- | :--- |
| **Omitted `(1\|System:Method)` and `(1\|Dataset:Method)` random interactions (pseudo-replication)** | Updated `R/08_mixed_effects_eval.R` to include `(1\|System:Method)` in synthetic LMMs (including `Method * SystemGroup` and GLMM) and `(1\|Dataset:Method)` in Experiment 2 LMMs. Re-ran `scripts/run_mixed_models.R` on all 109,192 + 2,600 rows. | Zero singular fits across all 10 LMMs (`Singular = No` in `table_mixed_effects_icc.tex`). On Synth F1 (`df=132.0`, SE=`0.067`), `BGNLM_Full_Tuned` is statistically indistinguishable from `BGLM_Oracle_Basis` ($p=0.865, p_{\text{adj}}=1.000$) while beating non-oracle baselines. On Exp. 2 (`df=48.0`), `Dataset:Fold_UID` variance recovers cleanly. |
| **Claimed Dunnett vs. printed unadjusted `lmerTest` $p$-values** | Merged Satterthwaite Dunnett-adjusted $p$-values from `emmeans(..., adjust="dunnett")` directly into `fixed_effects` tables (`p_{\text{adj}}`) alongside unadjusted $p$-values; Holm correction added for regime interactions and GLMM. | Updated `table_mixed_effects_synth.tex`, `table_mixed_effects_real.tex`, `table_mixed_effects_regime.tex`, `table_mixed_effects_glmm.tex`, and all text citations in `main.tex`. |
| **Complete ICC decomposition & omitted rows** | Removed `suppressWarnings` hiding singularity, added `Singular` column, and exported all random-effect groups cleanly in `table_mixed_effects_icc.tex` and Table 5 (`tab:icc`). | `System:Method` accounts for 25.1% (F1) and 39.6% (FDR) of variance; `Dataset:Method` accounts for 33.5% (Exp. 2 F1). |
| **Bounded $[0,1]$ outcomes & omitted binomial GLMM** | Exported binomial GLMM (`glmmTMB` logit link on $I_{\text{exact}} \in \{0,1\}$ with `(1\|System) + (1\|System:Method)`) to `table_mixed_effects_glmm.tex` and cited it in Section 4.2 / Appendix B. Explicit disclosure of bounded outcome caveats added. | Confirms statistically significant log-odds superiority of `BGNLM_Full_Tuned` over `SINDy_STLSQ_Poly`, `BGNLM_IWSM_Untuned`, `SINDy_RatPoly2`, `SR_Pool_*`, and `Bag_STLSQ_*`. |
| **Appendix table horizontal & vertical overflow** | Wrapped all exported LaTeX tables in `\resizebox{\textwidth}{!}{...}` with `\scriptsize`, `arraystretch=0.86`, `$t$` header, and stripped `SystemGroupGroup_` prefixes in `R/09_plotting_tables.R`. | Eliminated the 89pt overfull hbox in `table_mixed_effects_regime.tex` and the 13pt float overflow in `table_mixed_effects_synth.tex`. |

---

## 2. Baseline & Dataset Provenance Transparency (Addressing Reviewer 1 W1 & Reviewer 2 W1–W4)

| Critique | Action Taken | Where Addressed |
| :--- | :--- | :--- |
| **Experiment 2 simulated ODEs presented as empirical real-world measurements** | Reframed Experiment 2 throughout the manuscript and figure titles as a **Quasi-Empirical Literature ODE Benchmark (Noisy State Observations + Numerical Differentiation & Temporal CV)**, distinguishing it from Experiment 1's clean-state additive derivative-noise design and noting whole-series pre-split differentiation leakage. Built `prepare_authentic_lynx_hare()` option in `scripts/download_real_data.R` + targeted server ablation. | `paper_icml/main.tex` Abstract, Sec. 1, Sec. 3.1, Sec. 5, `scripts/download_real_data.R`, `scripts/run_targeted_ablation.R` |
| **`SR_PySR` & `SR_Operon` executed Python static-pool OMP fallback (`PYSR_USE_JULIA=0`, `PYOPERON_NATIVE=0`)** | Relabeled across all publication tables, figures, and prose to `SR_Pool_Pareto` and `SR_Pool_MDL` (documenting raw checkpoint IDs and the OMP+mutation heuristic); fixed native `PySR`/`PyOperon` OLS coefficient extraction in `python/run_sr_pysindy.py` + `R/04b_sr_baselines.R` for the targeted server ablation. | `R/08_mixed_effects_eval.R`, `scripts/generate_paper_figures.R`, `paper_icml/main.tex`, `python/run_sr_pysindy.py`, `R/04b_sr_baselines.R` |
| **`UQ_SINDy` unstandardized threshold bug & `SINDy_PI_Rational` scope** | Relabeled benchmark runs transparently to `Bag_STLSQ_Oracle` / `Bag_STLSQ_Poly` and `SINDy_RatPoly2`, explicitly explaining the adaptive threshold scale issue and explicit `Poly(2)` + 6 rational template design. Fixed `R/04_sindy_baselines.R` line 213 (`lam_base <- max(0.01, fit_base$lambda)`) for future/ablation runs. | `paper_icml/main.tex` Sec. 3.1, `R/04_sindy_baselines.R` |
| **Experiment 1 noise model clarity** | Explicitly documented that Experiment 1 applies additive Gaussian noise strictly to finite-difference derivatives with exact state observations (`noise_regime = "derivative_snr"`), and added a `state_noise` + native `PySR`/`PyOperon` + authentic Lynx--Hare server ablation script (`scripts/run_targeted_ablation.R`). | `paper_icml/main.tex` Sec. 2.1, Sec. 4.1, Sec. 6.2, `scripts/run_targeted_ablation.R` |

---

## 3. Theoretical, Algorithmic & Diagnostic Corrections (Addressing Reviewer 1 W3–W5 & Reviewer 2 W5)

1. **Jeffreys-BIC & Complexity Prior Precision (`Sec. 2.2`)**: Corrected marginal likelihood statement to unit-information Laplace/BIC approximation (Kass & Wasserman, 1995), included explicit $(p_\gamma + 1)/2 \log N$ intercept count and $\text{RSS}/N$ scaling, and noted combined $(c(g_k)+0.5)\log N$ penalty.
2. **Expressivity Trade-off, Post-hoc Pruning & Hybrid Intervals (`Sec. 2.3`, `Sec. 3.2`)**: Explicitly noted unit-scale restriction under $p_{\text{proj}}=0, D \le 2$; clarified that $L_2$-energy pruning ($\tau=0.025$ vs. $0.02$ in static-pool SR) operates post-MPM while `PIP_Brier` evaluates unpruned PIPs (and explained denominator dilution for `BGLM_Overcomplete`); replaced inaccurate "Student-$t$" interval wording with the exact normal ($1.96$) + unpruned BMA half-width quadrature formula around the post-pruning OLS fit.
3. **Hyperparameter & Calibration Audit (`Sec. 2.3`, `Sec. 2.4`, `Tab. 1`)**: Corrected GMJMCMC generations $P=15 \to 12$ (`config/best_bgnlm_params.yaml`); corrected calibration slice systems to **Lorenz-63, Lotka--Volterra, and Rössler--Lorenz Hybrid** at $\sigma \in \{0.01, 0.25\}$ ranked by `0.5*F1 + 0.2*Exact + 0.2*R2_ood - 0.1*Brier`; added `Exact` column to Table 1; added **Held-Out 9-System LMM Validation** (Sec. 4.2: `F1=0.808` vs. Oracle `0.832`, $p=0.763, p_{\text{adj}}=0.999$).
4. **Contradicted Explanations Fixed (`Sec. 4.3`, `Sec. 4.4`, `Sec. 5.2`)**:
   - **$N=2500$ F1 drop ($0.883 \to 0.760$)**: Re-diagnosed accurately from data as a **$3.4\times$ FDR surge ($0.063 \to 0.214$)** caused by $O(N)$ likelihood sharpening trapping fixed-budget chains in multi-term Taylor/rational surrogate modes.
   - **`sys01` ($\text{F1} \le 0.330$)**: Corrected upper bound (`0.267 -> 0.330` for `BGNLM_Full_Tuned`) and explained `[0.05T, 0.99T]` holdout removing the initial decay transient.
   - **`Yeast Glycolysis` ($R^2_{\text{traj}} \approx 0$)**: Corrected from "stiff phase-drift" to damped-spiral convergence to equilibrium before the 70% cutoff (`R2_in=0.947`, held-out amplitude $\approx \sigma$).
   - **Added `Bag_STLSQ_Oracle` row to Table 3**, disclosed $[0,1]$ $R^2$ clipping, intercept exclusion, and conditional-on-OK vs. zero-imputed `BGLM_Overcomplete` averages.
   - **Table 6 (`tab:systems_spec`)**: Fixed Thomas indexing to $\sin(x_{(i \bmod 3)+1})$, added Goodwin $p=2$ and $K_m+|x_3|$ notes, and Lorenz-96 cyclic boundary indexing.
   - **Fixed `R/05_symbolic_canonicalize.R`**: Evaluated parenthesized/unary-minus constant exponents so `x1^(-1)` and `x1^(0.5)` no longer emit malformed `^(x1)` tokens.
