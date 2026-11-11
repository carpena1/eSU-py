#!/usr/bin/env python3
"""
demo_external_model.py — a stand-in for YOUR model in Example B (esu_files.py).

It reads esu_sample.csv (made by `esu_files.py sample factors.csv`), runs a toy
runoff model once per row, and writes esu_outputs.csv with two outputs
(runoff and peak), one row per sample row, in the same order.

Replace this script with your own model runs. Only the file formats matter.

Run (from the examples/ folder):
    python ../esu_files.py sample factors.csv --p 8 --r 16 --seed 1
    python demo_external_model.py
    python ../esu_files.py analyze
"""
import csv
import math

with open("esu_sample.csv", newline="") as fh:
    rows = list(csv.DictReader(fh))

with open("esu_outputs.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["runoff", "peak"])
    for r in rows:
        Ks, CN, slope = float(r["Ks"]), float(r["CN"]), float(r["slope"])
        rain = 80.0                                   # mm, one storm
        S = 25400.0 / CN - 254.0                      # SCS-CN retention, mm
        q = (rain - 0.2 * S) ** 2 / (rain + 0.8 * S) if rain > 0.2 * S else 0.0
        q *= math.exp(-0.01 * Ks)                     # toy infiltration loss
        peak = 10.0 * q * math.sqrt(slope)            # toy peak flow
        w.writerow([f"{q:.6g}", f"{peak:.6g}"])       # porosity, dummy: no effect

print(f"wrote esu_outputs.csv ({len(rows)} runs)")
