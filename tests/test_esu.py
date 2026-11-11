#!/usr/bin/env python3
"""
test_esu.py — structural tests for the eSU sampler in esu.py.

eSU is randomized, so we can't reproduce MATLAB's exact samples; instead we check
the *properties* the method guarantees: uniform level coverage, the grid step
Delta, exactly one factor moved per trajectory step, and that the EE estimator
recovers a known linear model's effects (mu* = |slope|, sigma = 0).

Usage:  python tests/test_esu.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import esu


def test_pgen_uniformity():
    rng = np.random.default_rng(1)
    for p, r in [(4, 20), (8, 16), (12, 12), (6, 12)]:
        col = esu.esu_pgen(r, p, rng)
        levels, counts = np.unique(np.round(col, 12), return_counts=True)
        assert len(levels) == p, f"p={p}: got {len(levels)} distinct levels"
        assert np.all(counts == r // p), f"p={p} r={r}: level counts {counts} != {r//p}"
        assert np.allclose(levels, esu.esu_levels(p)), "levels off the grid"


def test_sample_shape_and_bounds():
    rng = np.random.default_rng(2)
    d, r, p = 16, 16, 8
    s = esu.esu_sample(d, r, p, rng)
    assert s.shape == ((d + 1) * r, d), s.shape
    assert s.min() >= -1e-12 and s.max() <= 1 + 1e-12, (s.min(), s.max())


def test_one_factor_per_step():
    rng = np.random.default_rng(3)
    d, r, p = 10, 12, 4
    delta = esu.esu_delta(p)
    s = esu.esu_sample(d, r, p, rng)
    for k in range(r):
        traj = s[k * (d + 1):(k + 1) * (d + 1)]
        moved = []
        for i in range(d):
            diff = traj[i + 1] - traj[i]
            nz = np.where(np.abs(diff) > 1e-12)[0]
            assert len(nz) == 1, f"step {i}: {len(nz)} factors moved"
            assert abs(abs(diff[nz[0]]) - delta) < 1e-12, "step != Delta"
            moved.append(nz[0])
        assert sorted(moved) == list(range(d)), "each factor must move exactly once"


def test_map_uniform_endpoints():
    lo = np.array([0.0, -5.0, 10.0])
    hi = np.array([1.0, 5.0, 20.0])
    assert np.allclose(esu.map_uniform(np.zeros(3), lo, hi), lo)
    assert np.allclose(esu.map_uniform(np.ones(3), lo, hi), hi)


def test_linear_model_effects():
    """Y = sum(a_i * real_i): unit-space EE_i = a_i*(hi-lo), constant -> sigma=0."""
    rng = np.random.default_rng(4)
    d, r, p = 6, 16, 8
    a = np.array([3.0, 0.0, -2.0, 0.5, 0.0, 1.0])
    lo = np.zeros(d); hi = np.array([2.0, 1.0, 4.0, 1.0, 1.0, 3.0])
    unit = esu.esu_sample(d, r, p, rng)
    real = esu.map_uniform(unit, lo, hi)
    Y = real @ a
    EE = esu.elementary_effects(unit, Y, p)
    mu_star, mu, sigma = esu.ee_measures(EE)
    expect = np.abs(a) * (hi - lo)
    assert np.allclose(mu_star, expect, atol=1e-9), (mu_star, expect)
    assert np.allclose(sigma, 0.0, atol=1e-9), sigma
    assert np.allclose(np.abs(mu), expect, atol=1e-9), "linear -> perfectly monotonic"


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    fail = 0
    for t in tests:
        try:
            t()
            print(f"PASS  {t.__name__}")
        except AssertionError as e:
            fail += 1
            print(f"FAIL  {t.__name__}: {e}")
    print(f"\n{len(tests)-fail}/{len(tests)} tests passed")
    sys.exit(1 if fail else 0)


if __name__ == "__main__":
    main()
