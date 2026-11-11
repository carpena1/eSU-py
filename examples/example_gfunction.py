#!/usr/bin/env python3
"""
example_gfunction.py — complete eSU Elementary-Effects analysis of a Python model.

Model: the Sobol g-function, a standard test with a KNOWN answer,
    Y = prod_i (|4 x_i - 2| + a_i) / (1 + a_i),   x_i in [0, 1]
A small a_i makes factor i important; a large a_i makes it unimportant.
Here a = [0, 1, 4.5, 9, 99, 99]: x1 > x2 > x3 > x4 >> x5, x6.

Run:  python example_gfunction.py
Out:  gfunction_EE_Measures.txt, gfunction_EE.png
"""
import sys
from pathlib import Path

import numpy as np

try:
    import esu
except ImportError:                                   # esu.py one folder up
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    import esu

# 1. Factors and their ranges (uniform distributions)
names = ["x1", "x2", "x3", "x4", "x5", "x6"]
lo = [0.0] * 6
hi = [1.0] * 6
d = len(names)
a = np.array([0.0, 1.0, 4.5, 9.0, 99.0, 99.0])


# 2. The model: one output value for one row of factor values
def model(x):
    return np.prod((np.abs(4 * x - 2) + a) / (1 + a))


# 3. Sampling design
p, r = 8, 32                                          # levels, trajectories
rng = np.random.default_rng(12345)                    # seed -> repeatable sample
unit = esu.esu_sample(d, r, p, rng)                   # (r*(d+1)) x d, in [0, 1]
real = esu.map_uniform(unit, lo, hi)                  # real factor values
print(f"{real.shape[0]} model runs")

# 4. Run the model for every row
Y = np.array([model(row) for row in real])

# 5. Elementary effects and the sensitivity measures
EE = esu.elementary_effects(unit, Y, p)               # d x r
mu_star, mu, sigma = esu.ee_measures(EE)
table = esu.format_measures(names, mu_star, mu, sigma)
print(table)
Path("gfunction_EE_Measures.txt").write_text("Y\n" + table + "\n")

# 6. Screening plot
esu.ee_plot(mu_star, mu, sigma, EE, names, "g-function", "gfunction_EE.png",
            rng=np.random.default_rng(0))
print("wrote gfunction_EE_Measures.txt and gfunction_EE.png")
