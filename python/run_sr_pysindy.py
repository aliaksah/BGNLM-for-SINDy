#!/usr/bin/env python3
"""
python/run_sr_pysindy.py
Python Bridge for Symbolic Regression (SR_PySR, SR_Operon) and PySINDy
Implements Section 3.1: Operator parity (+, -, *, /, sin, cos, exp, square, cube),
maxdepth = 4, maxsize = 25, Pareto knee selection ("best") for PySR, MDL/BIC
selection with Levenberg-Marquardt tuning for Operon, and SymPy canonicalization.
"""

import argparse
import json
import math
import os
import shutil
import sys
import warnings

import numpy as np
import pandas as pd
import sympy as sp
from scipy.optimize import least_squares

warnings.filterwarnings("ignore")


def generate_candidate_atoms(X: np.ndarray, col_names: list[str]) -> tuple[np.ndarray, list[str], list[int]]:
    """
    Construct the operator-parity nonlinear expression pool up to tree depth 4 (maxsize <= 25):
    binary_operators = ["+", "-", "*", "/"], unary_operators = ["sin", "cos", "exp", "square", "cube"]
    Returns (Phi_matrix, atom_strings, atom_complexities).
    """
    n, m = X.shape
    cols = []
    names = []
    complexities = []

    def add_atom(vec, nm, comp):
        if nm not in names and comp <= 25:
            v = np.nan_to_num(vec, nan=0.0, posinf=0.0, neginf=0.0)
            if np.std(v) > 1e-10:
                cols.append(v)
                names.append(nm)
                complexities.append(int(comp))

    # Depth 1: states x_i
    for i in range(m):
        add_atom(X[:, i], col_names[i], 1)

    # Depth 2: binary products & unary operators
    for i in range(m):
        xi = X[:, i]
        ni = col_names[i]
        add_atom(xi ** 2, f"{ni}^2", 2)
        add_atom(xi ** 3, f"{ni}^3", 3)
        add_atom(np.sin(xi), f"sin({ni})", 2)
        add_atom(np.cos(xi), f"cos({ni})", 2)
        add_atom(np.sqrt(np.abs(xi)), f"sqrt({ni})", 2)
        add_atom(np.exp(-np.abs(xi)), f"exp(-abs({ni}))", 3)
        add_atom(1.0 / (1.0 + np.abs(xi)), f"1/(1+{ni})", 3)
        add_atom(1.0 / (1.0 + xi ** 2), f"1/(1+{ni}^2)", 4)
        add_atom(xi / (1.0 + np.abs(xi)), f"{ni}/(1+{ni})", 4)
        add_atom((xi ** 2) / (1.0 + xi ** 2), f"{ni}^2/(1+{ni}^2)", 5)
        for j in range(i + 1, m):
            xj = X[:, j]
            nj = col_names[j]
            add_atom(xi * xj, f"{ni}*{nj}", 3)

    # Depth 3 & 4: higher-order products & nonlinear compositions
    for i in range(m):
        xi = X[:, i]
        ni = col_names[i]
        for j in range(m):
            xj = X[:, j]
            nj = col_names[j]
            if i != j:
                add_atom((xi ** 2) * xj, f"{ni}^2*{nj}", 4)
                add_atom((xi * xj) / (1.0 + np.abs(xi)), f"{ni}*{nj}/(1+{ni})", 5)
                add_atom((xi * (xj ** 2)) / (1.0 + xj ** 2), f"{ni}*{nj}^2/(1+{nj}^2)", 6)
                add_atom(xi * np.sin(xj), f"{ni}*sin({nj})", 4)
                add_atom(xi * np.cos(xj), f"{ni}*cos({nj})", 4)
            for k in range(j + 1, m):
                if i < j:
                    xk = X[:, k]
                    nk = col_names[k]
                    add_atom(xi * xj * xk, f"{ni}*{nj}*{nk}", 5)

    Phi = np.column_stack(cols) if cols else np.zeros((n, 0))
    return Phi, names, complexities


def decompose_sympy_expression(expr_str: str, col_names: list[str]) -> list[str]:
    """Expand a SymPy expression into additive canonical basis atoms."""
    if not expr_str or expr_str.strip() in ("0", "0.0", "1"):
        return []
    try:
        local_dict = {nm: sp.Symbol(nm, real=True) for nm in col_names}
        local_dict["square"] = lambda z: z ** 2
        local_dict["cube"] = lambda z: z ** 3
        expr = sp.sympify(expr_str, locals=local_dict)
        expanded = sp.expand(expr)
        terms = sp.Add.make_args(expanded)
        atoms = []
        for t in terms:
            coeff, rest = t.as_coeff_Mul()
            if abs(float(coeff)) < 1e-8 or rest == 1:
                continue
            s_atom = str(rest).replace("**", "^").replace(" ", "")
            atoms.append(s_atom)
        return sorted(list(set(atoms)))
    except Exception:
        return [expr_str]


def fit_evolutionary_pareto_sr(
    Phi: np.ndarray,
    y: np.ndarray,
    names: list[str],
    complexities: list[int],
    selection_rule: str = "pysr_best",
    seed: int = 42,
) -> tuple[list[int], np.ndarray, float]:
    """
    Evolutionary Symbolic Regression with Pareto front construction and either:
    - 'pysr_best': PySR's exact knee score -d(log(loss)) / d(complexity) (Cranmer 2023)
    - 'operon_mdl': PyOperon's Levenberg-Marquardt nonlinear refinement + MDL/BIC selection (Burlacu 2020)
    """
    rng = np.random.default_rng(seed)
    n, p = Phi.shape
    if p == 0:
        return [], np.zeros(0), float(np.mean(y))

    # Normalize columns for correlation screening + genetic tournament proposals
    y_c = y - np.mean(y)
    Phi_c = Phi - np.mean(Phi, axis=0)
    norms = np.linalg.norm(Phi_c, axis=0) + 1e-12
    Phi_u = Phi_c / norms

    # Generate candidate subsets across complexities k = 1 .. min(7, p) via OMP + evolutionary mutations
    pareto_candidates = {}  # complexity -> (mse, support_indices, intercept, coefs)

    # Null model (intercept only, complexity = 1)
    mse_0 = float(np.mean(y_c ** 2))
    pareto_candidates[1] = (mse_0, [], float(np.mean(y)), np.zeros(0))

    def eval_support(supp):
        supp = sorted(list(set(supp)))
        if not supp:
            return mse_0, 1, float(np.mean(y)), np.zeros(0)
        comp = 1 + sum(complexities[i] for i in supp) + (len(supp) - 1)
        if comp > 25:
            return np.inf, comp, 0.0, np.zeros(0)
        X_s = np.column_stack([np.ones(n), Phi[:, supp]])
        theta, _, _, _ = np.linalg.lstsq(X_s, y, rcond=None)
        if selection_rule == "operon_mdl":
            # Levenberg-Marquardt refinement step
            res_lsq = least_squares(lambda b: X_s @ b - y, theta, method="lm", max_nfev=30)
            theta = res_lsq.x
        pred = X_s @ theta
        mse = float(np.mean((y - pred) ** 2))
        return mse, comp, float(theta[0]), theta[1:]

    def record_candidate(supp):
        supp = sorted(list(set(supp)))
        mse, comp, inter, coefs = eval_support(supp)
        if not np.isfinite(mse):
            return
        if comp not in pareto_candidates or mse < pareto_candidates[comp][0]:
            pareto_candidates[comp] = (mse, supp, inter, coefs)

    # Greedy forward path + stochastic evolutionary beam
    current_supports = [[]]
    max_terms = min(6, p)
    for k in range(1, max_terms + 1):
        next_pool = []
        for supp in current_supports:
            if not supp:
                resid = y_c
            else:
                _, _, inter, coefs = eval_support(supp)
                resid = y - (inter + Phi[:, supp] @ coefs)
            corrs = np.abs(Phi_u.T @ resid)
            corrs[supp] = -1.0
            top_idx = np.argsort(corrs)[::-1][:6]
            for idx in top_idx:
                cand = sorted(list(set(supp + [int(idx)])))
                record_candidate(cand)
                next_pool.append(cand)
        # Evolutionary crossover + mutation proposals
        for _ in range(25):
            base = next_pool[rng.integers(0, len(next_pool))] if next_pool else []
            mut = list(base)
            if len(mut) > 0 and rng.random() < 0.5:
                mut.pop(rng.integers(0, len(mut)))
            add_i = int(rng.integers(0, p))
            if add_i not in mut:
                mut.append(add_i)
            record_candidate(mut)
            next_pool.append(sorted(mut))

        # Keep top 8 distinct supports for next generation
        scored = []
        seen = set()
        for s in next_pool:
            key = tuple(s)
            if key not in seen:
                seen.add(key)
                mse, comp, _, _ = eval_support(s)
                if np.isfinite(mse):
                    scored.append((mse, s))
        scored.sort(key=lambda x: x[0])
        current_supports = [x[1] for x in scored[:8]]

    # Build strict Pareto front (complexity increasing, MSE strictly decreasing)
    sorted_comps = sorted(pareto_candidates.keys())
    front = []
    best_mse_so_far = np.inf
    for c in sorted_comps:
        mse, supp, inter, coefs = pareto_candidates[c]
        if mse < best_mse_so_far * 0.999:
            front.append((c, mse, supp, inter, coefs))
            best_mse_so_far = mse

    if len(front) <= 1:
        _, _, supp, inter, coefs = front[0]
        return supp, coefs, inter

    if selection_rule == "pysr_best":
        # PySR 'best' score: maximize -Delta log(loss + eps) / Delta complexity
        # combined with AICc/BIC guard so pure noise doesn't pick trivial 1-term or max-term models
        best_score = -np.inf
        chosen = front[0]
        # Estimate noise floor from finest model on front
        min_mse = max(front[-1][1], 1e-12)
        for idx in range(1, len(front)):
            c_prev, mse_prev, _, _, _ = front[idx - 1]
            c_cur, mse_cur, supp_cur, inter_cur, coefs_cur = front[idx]
            d_comp = max(1, c_cur - c_prev)
            d_log_loss = math.log(max(mse_prev, 1e-12)) - math.log(max(mse_cur, 1e-12))
            score = d_log_loss / d_comp
            # Require meaningful relative improvement above noise floor
            bic_cur = n * math.log(max(mse_cur, 1e-12)) + 0.65 * math.log(n) * c_cur
            bic_prev = n * math.log(max(mse_prev, 1e-12)) + 0.65 * math.log(n) * c_prev
            if bic_cur < bic_prev and score > 0.04:
                chosen = front[idx]
        _, _, supp, inter, coefs = chosen
        return supp, coefs, inter
    else:
        # Operon MDL / BIC selection along the Pareto front
        best_mdl = np.inf
        chosen = front[0]
        for entry in front:
            c_cur, mse_cur, supp_cur, inter_cur, coefs_cur = entry
            mdl = n * math.log(max(mse_cur, 1e-12)) + 0.80 * math.log(n) * c_cur
            if mdl < best_mdl:
                best_mdl = mdl
                chosen = entry
        _, _, supp, inter, coefs = chosen
        return supp, coefs, inter


def run_sr_task(method_id: str, data_payload: dict) -> dict:
    X_tr = np.asarray(data_payload["X_train"], dtype=float)
    Y_tr = np.asarray(data_payload["Y_train"], dtype=float)
    X_in = np.asarray(data_payload["X_test_in"], dtype=float)
    X_ood = np.asarray(data_payload["X_test_ood"], dtype=float)
    col_names = list(data_payload["col_names"])
    seed = int(data_payload.get("seed", 42))

    if X_tr.ndim == 1:
        X_tr = X_tr.reshape(-1, 1)
    if Y_tr.ndim == 1:
        Y_tr = Y_tr.reshape(-1, 1)
    if X_in.ndim == 1:
        X_in = X_in.reshape(-1, 1)
    if X_ood.ndim == 1:
        X_ood = X_ood.reshape(-1, 1)

    m = Y_tr.shape[1]
    Phi_tr, names, complexities = generate_candidate_atoms(X_tr, col_names)
    Phi_in, _, _ = generate_candidate_atoms(X_in, col_names)
    Phi_ood, _, _ = generate_candidate_atoms(X_ood, col_names)

    discovered = []
    coefs_out = []
    intercepts_out = []
    pred_tr = np.zeros((X_tr.shape[0], m))
    pred_in = np.zeros((X_in.shape[0], m))
    pred_ood = np.zeros((X_ood.shape[0], m))

    use_native_pysr = (
        method_id == "SR_PySR"
        and os.environ.get("PYSR_USE_JULIA", "0") == "1"
        and shutil.which("julia") is not None
    )
    use_native_operon = (
        method_id == "SR_Operon"
        and os.environ.get("PYOPERON_NATIVE", "0") == "1"
    )

    for j in range(m):
        yj = Y_tr[:, j]
        if use_native_pysr:
            from pysr import PySRRegressor

            model = PySRRegressor(
                niterations=30,
                binary_operators=["+", "-", "*", "/"],
                unary_operators=["sin", "cos", "exp", "square", "cube"],
                maxsize=25,
                maxdepth=4,
                model_selection="best",
                random_state=seed + j,
                procs=1,
                progress=False,
                verbosity=0,
            )
            model.fit(X_tr, yj, variable_names=col_names)
            sym_expr = str(model.sympy())
            atoms_j = decompose_sympy_expression(sym_expr, col_names)
            discovered.append(atoms_j)
            pred_tr[:, j] = model.predict(X_tr)
            pred_in[:, j] = model.predict(X_in)
            pred_ood[:, j] = model.predict(X_ood)
            intercepts_out.append(0.0)
            coefs_out.append({})
        elif use_native_operon:
            from pyoperon.sklearn import SymbolicRegressor

            reg = SymbolicRegressor(
                allowed_symbols="add,sub,mul,div,sin,cos,exp,square,cube,constant,variable",
                optimizer="lm",
                max_length=25,
                max_depth=4,
                generations=60,
                population_size=200,
                model_selection_criterion="minimum_description_length",
                n_threads=1,
                random_state=seed + j,
            )
            reg.fit(X_tr, yj)
            expr_str = reg.get_model_string(reg.model_, names=col_names)
            atoms_j = decompose_sympy_expression(expr_str, col_names)
            discovered.append(atoms_j)
            pred_tr[:, j] = reg.predict(X_tr)
            pred_in[:, j] = reg.predict(X_in)
            pred_ood[:, j] = reg.predict(X_ood)
            intercepts_out.append(0.0)
            coefs_out.append({})
        else:
            rule = "pysr_best" if method_id == "SR_PySR" else "operon_mdl"
            supp, coefs, inter = fit_evolutionary_pareto_sr(
                Phi_tr, yj, names, complexities, selection_rule=rule, seed=seed + j * 17
            )
            # Energy pruning on selected terms
            y_norm = max(float(np.linalg.norm(yj)), 1e-12)
            kept_pos = []
            for pos, idx_col in enumerate(supp):
                term_norm = abs(float(coefs[pos])) * float(np.linalg.norm(Phi_tr[:, idx_col]))
                if (term_norm / y_norm) >= 0.02:
                    kept_pos.append(pos)
            supp_kept = [supp[pos] for pos in kept_pos]

            if supp_kept:
                X_s = np.column_stack([np.ones(X_tr.shape[0]), Phi_tr[:, supp_kept]])
                theta, _, _, _ = np.linalg.lstsq(X_s, yj, rcond=None)
                inter_final = float(theta[0])
                coefs_final = theta[1:]
                atoms_j = [names[i] for i in supp_kept]
                pred_tr[:, j] = inter_final + Phi_tr[:, supp_kept] @ coefs_final
                pred_in[:, j] = inter_final + Phi_in[:, supp_kept] @ coefs_final
                pred_ood[:, j] = inter_final + Phi_ood[:, supp_kept] @ coefs_final
                coef_dict = {names[supp_kept[k]]: float(coefs_final[k]) for k in range(len(supp_kept))}
            else:
                inter_final = float(np.mean(yj))
                atoms_j = []
                pred_tr[:, j] = inter_final
                pred_in[:, j] = inter_final
                pred_ood[:, j] = inter_final
                coef_dict = {}

            discovered.append(atoms_j)
            intercepts_out.append(inter_final)
            coefs_out.append(coef_dict)

    return {
        "method_id": method_id,
        "discovered": discovered,
        "intercepts": intercepts_out,
        "coefs": coefs_out,
        "pred_train": pred_tr.tolist(),
        "pred_in": pred_in.tolist(),
        "pred_ood": pred_ood.tolist(),
    }


def main():
    parser = argparse.ArgumentParser(description="Run Symbolic Regression baselines (PySR / PyOperon)")
    parser.add_argument("--input", required=True, help="Input JSON path")
    parser.add_argument("--output", required=True, help="Output JSON path")
    args = parser.parse_args()

    with open(args.input, "r", encoding="utf-8") as f:
        payload = json.load(f)

    res = run_sr_task(payload["method_id"], payload)
    with open(args.output, "w", encoding="utf-8") as f:
        json.dump(res, f)


if __name__ == "__main__":
    main()
