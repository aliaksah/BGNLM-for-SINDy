# Bayesian Generalized Nonlinear Models Offer Basis-Free SINDy With Model Uncertainty (ICML Benchmark Suite)

This repository contains the modular, checkpointed, seed-deterministic experimental pipeline expanding *"Bayesian Generalized Nonlinear Models Offer Basis Free SINDy With Model Uncertainty"* (`arXiv:2507.06776`) and the `FBMS` R package (`arXiv:2509.00753`) for ICML.

## Core Contributions & Experimental Architecture

1. **Expanded 12-System Synthetic Benchmark Suite (`R/01_ode_simulators.R`, `config/systems_config.yaml`)**:
   - **Group A (Polynomial ODEs)**: Linear-Quadratic Coupled 3D, Lorenz-63, Rössler Attractor, Lotka–Volterra Predator-Prey, Van der Pol Oscillator, Duffing Oscillator, Brusselator System.
   - **Group B (Non-Polynomial / Transcendental / Rational ODEs)**: Rössler–Lorenz Hybrid 3D, Damped Nonlinear Pendulum, Thomas' Cyclically Symmetric Attractor, Michaelis–Menten / Goodwin Enzymatic Oscillator.
   - **Group C (Higher-Dimensional Coupled ODEs)**: 5D Coupled Lorenz-96 System.
   - Evaluated across sample sizes $N \in \{100, 250, 500, 1000, 2500\}$, standardized derivative noise levels $\eta \in \{0.001, 0.01, 0.05, 0.10, 0.25, 0.50, 1.0\}$, state noise with Savitzky–Golay / TVRegDiff differentiation, and $R = 20$ deterministic replications per condition.

2. **Real-World Empirical Datasets Suite (`R/10_real_data_experiments.R`, `scripts/download_real_data.R`)**:
   - **Hudson Bay Lynx–Hare Pelts** (Ecology),
   - **Gause–Veilleux Paramecium–Didinium Laboratory Microcosm** (Experimental Protozoan Ecology),
   - **Schmidt & Lipson Physical Nonlinear Pendulum** (Lab Mechanics),
   - **Cascaded Water Tanks Benchmark** (Nonlinear System Identification),
   - **Yeast Glycolytic NADH Oscillations** (Biochemistry).

3. **13 Competing Method Configurations (`R/03_fbms_wrappers.R`, `R/04_sindy_baselines.R`, `R/04b_sr_baselines.R`)**:
   - **Proposed Basis-Free BGNLM & Ablations (`FBMS`)**: `BGNLM_Full_Tuned` (`gmjmcmc.parallel`), `BGNLM_IWSM_Untuned`, `BGNLM_Hybrid_Poly`.
   - **Predefined-Basis BGLM (`FBMS`)**: `BGLM_Oracle_Basis` (`mjmcmc.parallel`), `BGLM_Poly_Basis`, `BGLM_Overcomplete` ($q > 60$ candidate columns).
   - **Classical, Ensemble/Bayesian, & Implicit SINDy**: `SINDy_STLSQ_Oracle`, `SINDy_STLSQ_Poly`, `UQ_SINDy_Oracle`, `UQ_SINDy_Poly`, `SINDy_PI_Rational`.
   - **Symbolic Regression**: `SR_PySR` (`PySRRegressor` / Pareto knee selection), `SR_Operon` (`PyOperon` / Levenberg–Marquardt + MDL selection).

4. **Symbolic Canonicalization, Selection Consistency, & Mixed-Effects Evaluation (`R/05_symbolic_canonicalize.R`, `R/06_metrics_consistency.R`, `R/08_mixed_effects_eval.R`)**:
   - Algebraic AST canonicalization for exact structural comparison ($\text{TPR}$, $\text{FDR}$, $\text{F1}$, $I_{\text{exact}}$, $\text{SHD}$, Operation Count `OC`, Jaccard & Nogueira stability $\hat{\Phi}$, PIP Brier calibration, $95\%$ interval coverage & Winkler score, $R^2_{\text{in}}$, $R^2_{\text{ood}}$, and forward ODE trajectory rollout $R^2_{\text{traj}}, T_{\text{valid}}$).
   - Linear and Generalized Linear Mixed-Effects Models (`lme4`, `lmerTest`, `glmmTMB`, `emmeans`) with `BGNLM_Full_Tuned` as the reference baseline level.

---

## Quick Start & Execution Commands

```bash
# 1. Install R packages & Python virtual environment (.venv)
make install

# 2. Run unit test suite (ODE verification, symbolic canonicalizer, method wrappers)
make test

# 3. Run systematic BGNLM hyperparameter & prior tuning
make tune

# 4. Run local end-to-end smoke test (2 synthetic systems + 2 real datasets x 13 methods)
make smoke

# 5. Generate full ICML manifests & launch on multi-core remote Linux server
Rscript scripts/generate_manifest.R
./scripts/run_remote_server.sh both 16 4

# 6. Or submit SLURM job arrays on HPC cluster
sbatch slurm/submit_tuning_array.sh
sbatch slurm/submit_benchmark_array.sh
sbatch slurm/submit_real_data_array.sh
```
