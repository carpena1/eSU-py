#!/usr/bin/env python3
"""
validate_esu.py — check the Python EE measures (esu.py) against the group's
MATLAB example runs.

For each MATLAB example we read its factor sample (*_FacSample.sam/.xlsx), one
model-output column (*_Mod_Out*), and *_FacSamChar.txt (method, p, r, d, EEMul),
recompute mu*/mu/sigma with esu.elementary_effects + esu.ee_measures on MATLAB's
OWN samples+outputs (deterministic — independent of the randomized sampler), and
compare to MATLAB's *_EE_Measures.txt. We also render the Python ee_plot for an
eyeball comparison against *_Fig1.jpg.

Goal: consistent measures + identical factor ranking, not bit-perfect.

Usage:
    python tests/validate_esu.py            # all known examples it can read
    python tests/validate_esu.py Example3   # one example

Set ESU_MATLAB_DIR to the folder of the MATLAB package
EE_Sensitivity_Measures_Package_2022Feb18_ForSharing_R1 (it holds the example
files). Get the package from https://abe.ufl.edu/carpena/software/morris-su-sampling/
"""
import os
import sys
import csv
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import esu

MAT = Path(os.environ.get("ESU_MATLAB_DIR", ""))

# example -> candidate model-output files (first that exists & matches wins)
EXAMPLES = {
    "Example1": ["Example1_Mod_Out.xlsx"],                       # OT, p6 r10 d6
    "Example2": ["Example2_Mod_Output.txt", "Example2_Mod_Out.xlsx"],
    "Example3": ["Example3_Mod_Out.txt", "Example3_Mod_Out.out"],
    "SpCase":   ["SpCase_Mod_Out.xlsx"],
    # Example4 is RadialeSU (different EE estimator) -> out of scope here
}


def read_char(path):
    """method, p, r, d, EEMul(d x r) from *_FacSamChar.txt."""
    toks = path.read_text().split()
    method = toks[0]
    ovr, p, r, d = (int(toks[i]) for i in range(1, 5))
    mul = np.array([float(t) for t in "".join(toks[5:]).split(",") if t != ""])
    eemul = mul.reshape(d, r) if mul.size == d * r else np.ones((d, r))
    return method, p, r, d, eemul


def read_table(path):
    """(names, values[n x c]) from a .sam/.txt/.out/.csv (header row) or .xlsx."""
    if path.suffix == ".xlsx":
        import openpyxl
        wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
        ws = wb.active
        rows = [[c for c in row] for row in ws.iter_rows(values_only=True)]
        rows = [r for r in rows if any(c is not None for c in r)]
        names = [str(c) for c in rows[0]]
        vals = np.array([[float(c) for c in r] for r in rows[1:]], float)
        return names, vals
    text = [ln for ln in path.read_text().splitlines() if ln.strip()]
    delim = "," if "," in text[0] else None
    names = text[0].split(delim) if delim else text[0].split()
    data = [(ln.split(delim) if delim else ln.split()) for ln in text[1:]]
    vals = np.array([[float(x) for x in row] for row in data], float)
    return [n.strip() for n in names], vals


def read_measures(path):
    """{output_name: {label: (mu_star, mu, sigma)}} from a MATLAB *_EE_Measures.txt.
    The file has one section per model output: an output-name line, a 'facname...'
    header, the per-factor rows, then a dashed separator."""
    sections, cur = {}, None
    for ln in path.read_text().splitlines():
        if not ln.strip() or ln.startswith("-"):
            continue
        t = ln.split()
        if t[0] == "facname":
            continue
        if len(t) >= 4:                       # data row: label mu* mu sigma ...
            try:
                cur[t[0]] = (float(t[1]), float(t[2]), float(t[3]))
                continue
            except (ValueError, TypeError):
                pass
        cur = sections.setdefault(t[0], {})   # otherwise it's an output-name line
    return sections


def _rel(a, b):
    """max relative diff of a vs reference b (abs ref as denom, guarded)."""
    return float(np.max(np.abs(a - b) / np.where(np.abs(b) > 0, np.abs(b), 1.0)))


def check(example):
    method, p, r, d, eemul = read_char(MAT / f"{example}_FacSamChar.txt")
    sam_path = next((MAT / f"{example}_FacSample{e}" for e in (".sam", ".xlsx")
                     if (MAT / f"{example}_FacSample{e}").exists()), None)
    fac_names, sample = read_table(sam_path)
    sections = read_measures(MAT / f"{example}_EE_Measures.txt")

    # Some examples ship several alternative output files (different model runs);
    # only one reproduces the committed EE_Measures.txt. Try each, keep the best.
    def eval_out(outp):
        out_names, outv = read_table(outp)
        res, total = [], 0.0
        for ci, oname in enumerate(out_names):
            if oname not in sections:
                continue
            ref = sections[oname]
            ref_mus = np.array([ref[n][0] for n in fac_names])
            ref_sig = np.array([ref[n][2] for n in fac_names])
            EE = esu.elementary_effects(sample, outv[:, ci], p) * eemul
            mu_star, mu, sigma = esu.ee_measures(EE)
            rel_mus, rel_sig = _rel(mu_star, ref_mus), _rel(sigma, ref_sig)
            rank_ok = np.allclose(np.sort(mu_star)[::-1], np.sort(ref_mus)[::-1], rtol=2e-3, atol=1e-9)
            ok = rel_mus < 2e-3 and rel_sig < 2e-3 and rank_ok
            total += rel_mus + rel_sig
            res.append((oname, rel_mus, rel_sig, rank_ok, ok, mu_star, mu, sigma, EE))
        return total, res

    best = None
    for o in EXAMPLES[example]:
        if not (MAT / o).exists():
            continue
        try:
            total, res = eval_out(MAT / o)
        except ImportError:
            print(f"{example}: openpyxl needed for {o} — SKIP")
            return None
        if res and (best is None or total < best[0]):
            best = (total, o, res)
    if best is None:
        print(f"{example}: no usable output file — SKIP")
        return None

    _, oname_file, res = best
    print(f"{example} [{method} p={p} r={r} d={d}] sample={sam_path.name} out={oname_file}")
    all_ok = True
    for k, (oname, rel_mus, rel_sig, rank_ok, ok, mu_star, mu, sigma, EE) in enumerate(res):
        all_ok &= ok
        print(f"  {oname}: rel mu*={rel_mus:.1e} sigma={rel_sig:.1e} "
              f"rank={'ok' if rank_ok else 'DIFFER'} -> {'PASS' if ok else 'FAIL'}")
        if k == 0:
            png = Path(__file__).resolve().parent / f"{example}_esu_py.png"
            esu.ee_plot(mu_star, mu, sigma, EE, fac_names, f"{example} {oname} (py)",
                        png, n_label=max(1, int(0.4 * d)), rng=np.random.default_rng(0))
    return all_ok


def main():
    if not os.environ.get("ESU_MATLAB_DIR") or not MAT.is_dir():
        sys.exit("Set ESU_MATLAB_DIR to the MATLAB EE_Sensitivity_Measures package folder.")
    which = sys.argv[1:] or list(EXAMPLES)
    results = {e: check(e) for e in which}
    done = {e: v for e, v in results.items() if v is not None}
    print(f"\n{sum(done.values())}/{len(done)} examples PASS")
    sys.exit(0 if all(done.values()) else 1)


if __name__ == "__main__":
    main()
