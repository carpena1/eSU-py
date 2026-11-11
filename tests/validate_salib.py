#!/usr/bin/env python3
"""
validate_salib.py — ballpark cross-check of the eSU screening against SALib's
Method of Morris (optimal trajectories) on an analytic function with KNOWN
importance (the Sobol g-function).

g(x) = prod_i (|4 x_i - 2| + a_i) / (1 + a_i),  x_i in [0,1]
Small a_i  => important factor; large a_i => unimportant. We check that eSU and
SALib rank/important-classify the factors the same way, and report the run counts.

This is a development sanity check. You do not need it to use esu.py.
Requires the dev deps:  pip install -r requirements-dev.txt

Usage:  python tests/validate_salib.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import esu

A = np.array([0.0, 0.5, 3.0, 9.0, 99.0, 99.0, 99.0, 99.0])     # g-function coeffs
D = len(A)
TRUE_IMPORTANT = set(np.where(A <= 9.0)[0].tolist())           # P1..P4 matter


def gfun(X):
    return np.prod((np.abs(4 * X - 2) + A) / (1 + A), axis=1)


def esu_screen(p=8, r=16, seed=0):
    rng = np.random.default_rng(seed)
    unit = esu.esu_sample(D, r, p, rng)
    Y = gfun(unit)                                             # bounds are [0,1]
    EE = esu.elementary_effects(unit, Y, p)
    mu_star, _, _ = esu.ee_measures(EE)
    return mu_star, unit.shape[0]


def salib_screen(p=8, r=16, seed=0):
    from SALib.sample.morris import sample
    from SALib.analyze.morris import analyze
    prob = {"num_vars": D, "names": [f"P{i+1}" for i in range(D)],
            "bounds": [[0.0, 1.0]] * D}
    X = sample(prob, N=r, num_levels=p, optimal_trajectories=None, seed=seed)
    Y = gfun(X)
    res = analyze(prob, X, Y, num_levels=p, print_to_console=False, seed=seed)
    return np.array(res["mu_star"]), X.shape[0]


def top_set(mu_star, k):
    return set(np.argsort(mu_star)[::-1][:k].tolist())


def main():
    p, r = 8, 16
    esu_mus, esu_runs = esu_screen(p, r)
    print(f"g-function: {D} factors, truly important = "
          f"{sorted('P%d' % (i+1) for i in TRUE_IMPORTANT)}\n")
    print(f"eSU   (p={p}, r={r}): {esu_runs} runs")
    print("  mu* =", np.array2string(esu_mus, precision=3))
    k = len(TRUE_IMPORTANT)
    esu_top = top_set(esu_mus, k)
    print(f"  top-{k} by mu* = {sorted('P%d' % (i+1) for i in esu_top)} "
          f"-> {'MATCH' if esu_top == TRUE_IMPORTANT else 'DIFFER'} vs truth")

    try:
        sal_mus, sal_runs = salib_screen(p, r)
    except ImportError:
        print("\nSALib not installed (pip install -r requirements-dev.txt) — skipping benchmark.")
        return
    sal_top = top_set(sal_mus, k)
    print(f"\nSALib (Morris, num_levels={p}, N={r}): {sal_runs} runs")
    print("  mu* =", np.array2string(sal_mus, precision=3))
    print(f"  top-{k} by mu* = {sorted('P%d' % (i+1) for i in sal_top)} "
          f"-> {'MATCH' if sal_top == TRUE_IMPORTANT else 'DIFFER'} vs truth")

    overlap = len(esu_top & sal_top)
    print(f"\neSU vs SALib important-set overlap: {overlap}/{k}  "
          f"(eSU runs {esu_runs} vs SALib runs {sal_runs})")
    print("ballpark agreement OK" if esu_top == TRUE_IMPORTANT and overlap >= k - 1
          else "review: rankings diverge")


if __name__ == "__main__":
    main()
