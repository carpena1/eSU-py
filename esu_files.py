#!/usr/bin/env python3
"""
esu_files.py — eSU Elementary-Effects analysis through CSV files, for any model.

Use it when your model is NOT written in Python (Fortran, R, Excel, a GUI ...).

Step 1  make the sample:
        python esu_files.py sample factors.csv --p 8 --r 16 --seed 1
        -> esu_sample.csv   (real units, one row per model run)
        -> esu_design.json  (p, r, factor names and ranges; used by step 3)
Step 2  run your model once for each row of esu_sample.csv, in the same order,
        and save esu_outputs.csv (header row, one column per model output).
Step 3  compute the sensitivity measures and the plots:
        python esu_files.py analyze
        -> EE_Measures.txt, EE_raw_<output>.csv, EE_<output>.png

factors.csv has a header row and one row per factor:  name,lo,hi
Needs esu.py in the same folder, in the parent folder, or on PYTHONPATH.
See docs/eSU_Python_Guide.md for the full description.
"""
import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

try:
    import esu
except ImportError:                                   # esu.py one folder up
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    import esu


def read_factors(path):
    with open(path, newline="") as fh:
        rows = [r for r in csv.DictReader(fh) if r.get("name", "").strip()]
    names = [r["name"].strip() for r in rows]
    lo = [float(r["lo"]) for r in rows]
    hi = [float(r["hi"]) for r in rows]
    for n, a, b in zip(names, lo, hi):
        if not b > a:
            sys.exit(f"factor {n}: hi ({b}) must be larger than lo ({a})")
    return names, lo, hi


def read_table(path):
    """(header, values[n x c]) from a CSV file with a header row.
    Empty cells and text such as 'nan' or 'NA' become NaN (a failed run)."""
    with open(path, newline="") as fh:
        rows = [r for r in csv.reader(fh) if r and any(c.strip() for c in r)]
    header = [h.strip() for h in rows[0]]

    def num(c):
        try:
            return float(c)
        except ValueError:
            return np.nan
    return header, np.array([[num(c) for c in r] for r in rows[1:]], float)


def cmd_sample(a):
    names, lo, hi = read_factors(a.factors)
    d = len(names)
    rng = np.random.default_rng(a.seed)
    unit = esu.esu_sample(d, a.r, a.p, rng)
    real = esu.map_uniform(unit, lo, hi)
    out = Path(a.out)
    np.savetxt(out, real, delimiter=",", header=",".join(names), comments="",
               fmt="%.10g")
    design = dict(p=a.p, r=a.r, d=d, seed=a.seed, names=names, lo=lo, hi=hi,
                  sample=out.name)
    Path(a.design).write_text(json.dumps(design, indent=2))
    print(f"{a.r} trajectories x {d + 1} points = {real.shape[0]} model runs")
    print(f"wrote {out} and {a.design}")


def cmd_analyze(a):
    design = json.loads(Path(a.design).read_text())
    p, names = design["p"], design["names"]
    s_names, sample = read_table(a.sample or design["sample"])
    o_names, outs = read_table(a.outputs)
    if s_names != names:
        sys.exit("sample columns do not match the factor names in the design file")
    if outs.shape[0] != sample.shape[0]:
        sys.exit(f"{a.outputs} has {outs.shape[0]} rows but the sample has "
                 f"{sample.shape[0]}: you need one output row per sample row")
    rng = np.random.default_rng(0)                    # fixed bootstrap -> same plot
    report = []
    for k, oname in enumerate(o_names):
        Y = outs[:, k]
        nbad = int(np.sum(~np.isfinite(Y)))
        EE = esu.elementary_effects(sample, Y, p)
        mu_star, mu, sigma = esu.ee_measures(EE)
        report += [oname, esu.format_measures(names, mu_star, mu, sigma), "-" * 49]
        with open(f"EE_raw_{oname}.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["factor"] + [f"traj{i + 1}" for i in range(EE.shape[1])])
            for j, n in enumerate(names):
                w.writerow([n] + [f"{v:.6g}" for v in EE[j]])
        esu.ee_plot(mu_star, mu, sigma, EE, names, oname, f"EE_{oname}.png", rng=rng)
        print(f"\nOutput '{oname}'" + (f"  ({nbad} failed runs ignored)" if nbad else ""))
        print(f"  {'factor':<12}{'mu*':>12}{'mu':>12}{'sigma':>12}")
        for j in np.argsort(-np.nan_to_num(mu_star, nan=-1)):
            print(f"  {names[j]:<12}{mu_star[j]:>12.4g}{mu[j]:>12.4g}{sigma[j]:>12.4g}")
    Path(a.report).write_text("\n".join(report) + "\n")
    print(f"\nwrote {a.report}, EE_raw_<output>.csv and EE_<output>.png")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("sample", help="make the eSU sample")
    s.add_argument("factors", help="CSV with columns name,lo,hi")
    s.add_argument("--p", type=int, default=8, help="grid levels, even (default 8)")
    s.add_argument("--r", type=int, default=16,
                   help="trajectories, even and a multiple of p (default 16)")
    s.add_argument("--seed", type=int, default=None, help="random seed (repeatable sample)")
    s.add_argument("--out", default="esu_sample.csv")
    s.add_argument("--design", default="esu_design.json")
    s.set_defaults(func=cmd_sample)
    z = sub.add_parser("analyze", help="compute EE measures and plots")
    z.add_argument("--design", default="esu_design.json")
    z.add_argument("--sample", default=None, help="default: the file named in the design")
    z.add_argument("--outputs", default="esu_outputs.csv")
    z.add_argument("--report", default="EE_Measures.txt")
    z.set_defaults(func=cmd_analyze)
    a = ap.parse_args()
    if a.cmd == "sample" and (a.p % 2 or a.r % a.p):
        sys.exit(f"--p must be even and --r a multiple of p (got p={a.p}, r={a.r})")
    a.func(a)


if __name__ == "__main__":
    main()
