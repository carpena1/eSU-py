#!/usr/bin/env python3
"""
esu.py — Enhanced Sampling for Uniformity (eSU) Morris screening.

A self-contained Python port of the eSU sampler and the Elementary-Effects (EE)
sensitivity measures/plots from the group's MATLAB package (Chitale et al. 2017;
Khare et al. 2019, J. Hydrol. Eng. 24(1)). Uniform factor ranges only.

Two stages:
  1. esu_sample(d, r, p) -> r Morris trajectories in the unit cube [0,1]^d, with
     each factor's starting level spread uniformly over the p-level grid (the
     "uniformity" property). map_uniform() then maps them to real [lo, hi].
  2. elementary_effects(sample, Y, p) + ee_measures() -> per-factor mu*, mu, sigma
     from one scalar model output Y per sample row. ee_plot() draws the mu*-sigma
     and mu-sigma screening figure (Khare 2019 layout).

The sampler is randomized; this reproduces the MATLAB *properties* (uniform level
coverage, the grid step Delta, the EE formula), not byte-identical samples. Pass a
numpy Generator for reproducibility.

Source: the Morris SU (Sampling Uniformity) MATLAB packages (EE Sampler/Mapper and
        EE Sensitivity Measures, 2022), https://abe.ufl.edu/carpena/software/morris-su-sampling/
User guide: docs/eSU_Python_Guide.md (or .pdf).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# --------------------------------------------------------------------------- #
# eSU sampling (unit hypercube)
# --------------------------------------------------------------------------- #
def esu_delta(p):
    """Morris grid step in unit space (= p/(2(p-1)))."""
    return (p / 2.0) / (p - 1.0)


def esu_levels(p):
    """The p grid levels {0, 1/(p-1), ..., 1}."""
    return np.arange(p) / (p - 1.0)


def esu_pgen(r, p, rng):
    """Starting levels for ONE factor: each of the p grid levels appears exactly
    r/p times, in random order (the eSU uniformity property; generalizes the
    hard-coded p4GEN..p16GEN tables)."""
    if r % p != 0:
        raise ValueError(f"r ({r}) must be a multiple of p ({p})")
    col = np.repeat(esu_levels(p), r // p)
    rng.shuffle(col)
    return col


def esu_base_points(d, r, p, rng):
    """r x d matrix of trajectory start points; one independent esu_pgen column
    per factor (port of UNIQ_TRAJ_GEN_CHECK)."""
    return np.column_stack([esu_pgen(r, p, rng) for _ in range(d)])


def esu_perturb_order(d, rng):
    """Random order in which the d factors are perturbed along a trajectory
    (port of UNI_DEL)."""
    order = np.arange(d)
    rng.shuffle(order)
    return order


def esu_traj(base_row, order, p):
    """One (d+1) x d trajectory from a start point and a perturbation order
    (port of TRAJ_GEN): perturb each factor once, in `order`, by +Delta if its
    current value < 0.5 else -Delta (keeps points on the [0,1] grid)."""
    d = base_row.size
    delta = esu_delta(p)
    traj = np.empty((d + 1, d))
    traj[0] = base_row
    for i, j in enumerate(order):
        traj[i + 1] = traj[i]
        traj[i + 1, j] += delta if traj[i, j] < 0.5 else -delta
    return traj


def esu_sample(d, r, p, rng=None):
    """eSU sample: r trajectories stacked into a ((d+1)*r) x d unit-cube matrix.
    r must be even and a multiple of p."""
    rng = np.random.default_rng() if rng is None else rng
    if r % 2 != 0:
        raise ValueError(f"r ({r}) must be even")
    if r % p != 0:
        raise ValueError(f"r ({r}) must be a multiple of p ({p})")
    base = esu_base_points(d, r, p, rng)
    blocks = [esu_traj(base[k], esu_perturb_order(d, rng), p) for k in range(r)]
    return np.vstack(blocks)


def map_uniform(unit, lo, hi):
    """Map unit-cube samples to real ranges, per factor (port of Uniform.m)."""
    lo = np.asarray(lo, float)
    hi = np.asarray(hi, float)
    return lo + (hi - lo) * unit


# --------------------------------------------------------------------------- #
# Elementary-effects measures
# --------------------------------------------------------------------------- #
def elementary_effects(sample, Y, p):
    """EE matrix (d x r) from a sample (unit OR real — only each step's sign is
    used) and one scalar output Y per sample row (port of Morris_Measures.m).
    Non-finite outputs propagate to NaN EEs, excluded by the nan-aware measures.
    """
    sample = np.asarray(sample, float)
    Y = np.asarray(Y, float)
    n, d = sample.shape
    r = n // (d + 1)
    delta = esu_delta(p)
    EE = np.full((d, r), np.nan)
    for i in range(r):
        s = sample[i * (d + 1):(i + 1) * (d + 1)]
        y = Y[i * (d + 1):(i + 1) * (d + 1)]
        for ii in range(d):
            diff = s[ii + 1] - s[ii]
            j = int(np.argmax(np.abs(diff)))          # the one factor that moved
            EE[j, i] = np.sign(diff[j]) * (y[ii + 1] - y[ii]) / delta
    return EE


def ee_measures(EE):
    """Per-factor (mu_star, mu, sigma): mean|EE|, mean EE, std EE (sample std,
    ddof=1, matching MATLAB std; nan-aware)."""
    mu_star = np.nanmean(np.abs(EE), axis=1)
    mu = np.nanmean(EE, axis=1)
    sigma = np.nanstd(EE, axis=1, ddof=1)
    return mu_star, mu, sigma


def monotonic_mask(mu_star, mu):
    """True where mu* == |mu| exactly, i.e. all EEs share one sign (perfectly
    monotonic / additive factor) — same test as the MATLAB 'Monotnic' flag."""
    return mu_star == np.abs(mu)


def bootstrap_ci(x, stat, n=100, alpha=0.05, rng=None):
    """Percentile bootstrap CI for `stat` over finite values of x. Approximates
    MATLAB bootci's default (BCa) for the plot error bars; the point measures
    themselves are exact."""
    rng = np.random.default_rng() if rng is None else rng
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return (np.nan, np.nan)
    boots = np.array([stat(rng.choice(x, x.size, replace=True)) for _ in range(n)])
    return tuple(np.percentile(boots, [100 * alpha / 2, 100 * (1 - alpha / 2)]))


def format_measures(labels, mu_star, mu, sigma):
    """Text table with the MATLAB EE_Measures.txt columns. The columns line up
    for any label length."""
    labels = [str(s) for s in labels]
    w = max([len("facname")] + [len(s) for s in labels])
    mono = monotonic_mask(mu_star, mu)
    out = [f"{'facname':<{w}}  {'mu*':>11}  {'mu':>11}  {'sigma':>11}  Monotonic  Interactions"]
    for i, lab in enumerate(labels):
        inter = "Y" if mu_star[i] <= sigma[i] else "N"
        out.append(f"{lab:<{w}}  {mu_star[i]:11.4e}  {mu[i]:11.4e}  {sigma[i]:11.4e}  "
                   f"{'Y' if mono[i] else 'N':^9}  {inter:^12}".rstrip())
    return "\n".join(out)


# --------------------------------------------------------------------------- #
# Screening plot (Khare 2019 mu*-sigma / mu-sigma)
# --------------------------------------------------------------------------- #
def ee_plot(mu_star, mu, sigma, EE, labels, title, out_png, n_label=None,
            label_mask=None, rng=None):
    """Two-panel EE screening figure (port of EE_plots.m, R1):
      left  : mu* vs sigma, with the mu*=sigma (1:1) threshold
      right : mu  vs sigma, with the mu = +/-2 sigma/sqrt(r) (2*SEM) lines
    Monotonic factors (mu*=|mu|) are filled circles, others asterisks; golden
    bars show bootstrap 95% CIs on both axes. Saves `out_png`.

    Which factors are annotated: if `label_mask` (a boolean array over factors)
    is given, only those are labelled (e.g. the important factors kept for
    calibration) — the rest are drawn as unlabelled symbols; otherwise the
    top-`n_label` factors by mu* are labelled."""
    rng = np.random.default_rng() if rng is None else rng
    mu_star = np.asarray(mu_star, float)
    mu = np.asarray(mu, float)
    sigma = np.asarray(sigma, float)
    labels = list(labels)
    d = len(labels)
    r = EE.shape[1]
    mono = monotonic_mask(mu_star, mu)

    order = np.argsort(mu_star)[::-1]            # most important first
    if label_mask is not None:
        label_mask = np.asarray(label_mask, bool)
        lbl_idx = [j for j in order if label_mask[j]]
    else:
        lbl_idx = list(order[:(d if n_label is None else n_label)])

    # bootstrap CIs (asymmetric, as half-widths) per factor
    musL = np.full(d, np.nan); musU = np.full(d, np.nan)
    muL = np.full(d, np.nan); muU = np.full(d, np.nan)
    sigL = np.full(d, np.nan); sigU = np.full(d, np.nan)
    for j in range(d):
        lo, hi = bootstrap_ci(np.abs(EE[j]), np.mean, rng=rng)
        musL[j], musU[j] = mu_star[j] - lo, hi - mu_star[j]
        lo, hi = bootstrap_ci(EE[j], np.mean, rng=rng)
        muL[j], muU[j] = mu[j] - lo, hi - mu[j]
        lo, hi = bootstrap_ci(EE[j], lambda v: np.std(v, ddof=1), rng=rng)
        sigL[j], sigU[j] = sigma[j] - lo, hi - sigma[j]

    gold = "#EDB120"
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 6.2))

    # ---- left: mu* - sigma ----
    dmax = np.nanmax([np.nanmax(mu_star + np.nan_to_num(musU)),
                      np.nanmax(sigma + np.nan_to_num(sigU))]) * 1.3
    dmax = dmax if np.isfinite(dmax) and dmax > 0 else 1.0
    axL.errorbar(mu_star, sigma, xerr=[np.nan_to_num(musL), np.nan_to_num(musU)],
                 yerr=[np.nan_to_num(sigL), np.nan_to_num(sigU)],
                 fmt="none", ecolor=gold, lw=1, zorder=1)
    axL.scatter(mu_star[~mono], sigma[~mono], marker="*", c="b", s=55,
                label="non-monotonic", zorder=2)
    axL.scatter(mu_star[mono], sigma[mono], marker="o", c="r", s=30,
                label="monotonic", zorder=2)
    g = np.linspace(0, dmax / 1.2, 50)
    axL.plot(g, g, "--k", lw=1, label=r"threshold $\mu^*=\sigma$")
    for j in lbl_idx:
        axL.annotate(labels[j], (mu_star[j], sigma[j]), fontsize=7,
                     ha="left", va="bottom")
    axL.set_xlim(0, dmax); axL.set_ylim(0, dmax)
    axL.set_xlabel(r"$\mu^*$", fontsize=13); axL.set_ylabel(r"$\sigma$", fontsize=13)
    axL.set_title(f"{title}  —  $\\mu^*$ vs $\\sigma$")
    axL.grid(True, color="#eee"); axL.legend(loc="upper right", fontsize=8)

    # ---- right: mu - sigma ----
    dr = np.nanmax([np.nanmax(np.abs(mu) + np.nan_to_num(muU)),
                    np.nanmax(sigma + np.nan_to_num(sigU))]) * 1.3
    dr = dr if np.isfinite(dr) and dr > 0 else 1.0
    axR.errorbar(mu, sigma, xerr=[np.nan_to_num(muL), np.nan_to_num(muU)],
                 yerr=[np.nan_to_num(sigL), np.nan_to_num(sigU)],
                 fmt="none", ecolor=gold, lw=1, zorder=1)
    axR.scatter(mu[~mono], sigma[~mono], marker="*", c="b", s=55,
                label="non-monotonic", zorder=2)
    axR.scatter(mu[mono], sigma[mono], marker="o", c="r", s=30,
                label="monotonic", zorder=2)
    sem = 2.0 / np.sqrt(r)
    f = np.linspace(0, dr, 50)
    axR.plot(sem * f, f, "--k", lw=1, label=r"$\mu=\pm2\,\sigma/\sqrt{r}$")
    axR.plot(-sem * f, f, "--k", lw=1)
    for j in lbl_idx:
        axR.annotate(labels[j], (mu[j], sigma[j]), fontsize=7, ha="left", va="bottom")
    axR.set_xlim(-dr, dr); axR.set_ylim(0, dr)
    axR.set_xlabel(r"$\mu$", fontsize=13); axR.set_ylabel(r"$\sigma$", fontsize=13)
    axR.set_title(f"{title}  —  $\\mu$ vs $\\sigma$")
    axR.grid(True, color="#eee"); axR.legend(loc="upper right", fontsize=8)

    fig.tight_layout()
    fig.savefig(out_png, dpi=130)
    plt.close(fig)
    return out_png
