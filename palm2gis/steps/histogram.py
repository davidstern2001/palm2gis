"""Value distribution of the voxel variables over ALL written records.

An ArcGIS voxel layer computes its statistics from the data it has loaded
(the current time step), and the full-run min/max are single extreme
cells - a stretch over either hides most of the field. This module bins
every value of every record, reads robust percentiles off the result and
writes, next to the voxel file:

  <stem>_ranges.csv     per variable: min, p1, p2, p5, p50, p95, p98, p99,
                        max - p2 .. p98 is the suggested display Min/Max
  <stem>_histogram.csv  the binned counts (variable, bin_lo, bin_hi, count)
  <stem>_histogram.png  one histogram per variable with p2 / p98 marked
                        (only if matplotlib is installed)

The percentiles come from a fine histogram (FINE_BINS bins between the
full-run min and max), so each is exact to (max - min) / FINE_BINS.
"""

import csv

import numpy as np

from ..log import debug, info

FINE_BINS = 4000         # percentile resolution: (max - min) / FINE_BINS
PLOT_BINS = 100          # bins in the CSV and the plot
PERCENTILES = (1, 2, 5, 50, 95, 98, 99)
DISPLAY = (2, 98)        # suggested display range

# Reference chart palette (light surface): one series hue, text inks.
_SURFACE = "#fcfcfb"
_BAR = "#2a78d6"
_INK = "#0b0b0b"
_INK_2 = "#52514e"
_GRID = "#e4e3df"


# ------------------------------
# 1. BINNING OVER ALL RECORDS
# ------------------------------
def distribution(ncvar, n_rec, lo, hi, fill):
    """Fine histogram of one (time, ...) netCDF variable over all records.

    Read back one record at a time, so memory stays at one 3D field.
    Returns (edges, counts)."""
    if not hi > lo:                       # constant field
        hi = lo + max(abs(lo) * 1e-6, 1e-6)
    edges = np.linspace(lo, hi, FINE_BINS + 1)
    counts = np.zeros(FINE_BINS, dtype=np.int64)
    for q in range(n_rec):
        a = np.ma.filled(ncvar[q], np.nan).astype("f8", copy=False)
        a = a[np.isfinite(a) & (a != fill)]
        if a.size:
            counts += np.histogram(a, bins=edges)[0]
    return edges, counts


def percentiles(edges, counts, ps=PERCENTILES):
    """Percentiles from a histogram, interpolated linearly inside a bin."""
    total = counts.sum()
    if total == 0:
        return {p: float("nan") for p in ps}
    cdf = np.concatenate([[0.0], np.cumsum(counts) / total])
    return {p: float(np.interp(p / 100.0, cdf, edges)) for p in ps}


def coarse(edges, counts, n=PLOT_BINS):
    """Merge the fine bins into n bins for the CSV and the plot."""
    k = len(counts) // n
    return edges[::k][:n + 1], counts[:k * n].reshape(n, k).sum(axis=1)


# ------------------------------
# 2. STATISTICS FOR EVERY VARIABLE
# ------------------------------
def voxel_statistics(fout, names, n_rec, ranges, fill):
    """{var: dict(units, min, max, p.., edges, counts)} for the written
    variables that hold data; also sets `palm2gis_display_range` = [p2,
    p98] on each variable."""
    stats = {}
    for var in names:
        if var not in ranges:
            continue
        lo, hi = ranges[var]
        edges, counts = distribution(fout.variables[var], n_rec, lo, hi,
                                     fill)
        pct = percentiles(edges, counts)
        units = getattr(fout.variables[var], "units", "")
        fout.variables[var].palm2gis_display_range = np.array(
            [pct[DISPLAY[0]], pct[DISPLAY[1]]], dtype="f4")
        fout.variables[var].palm2gis_display_range_comment = (
            "p{} .. p{} of all values over all {} records: a Min/Max for "
            "the voxel layer stretch that is not set by a few extreme "
            "cells".format(DISPLAY[0], DISPLAY[1], n_rec))
        stats[var] = dict(units=units, min=lo, max=hi, pct=pct,
                          edges=edges, counts=counts)
        info("[voxel] {}: suggested Min/Max {:.4g} .. {:.4g} {}", var,
             pct[DISPLAY[0]], pct[DISPLAY[1]], units)
        debug("[voxel] {}: Min/Max = p{} .. p{}; full range {:.4g} .. {:.4g}",
              var, DISPLAY[0], DISPLAY[1], lo, hi)
    return stats


# ------------------------------
# 3. SIDECAR FILES
# ------------------------------
def write_sidecars(stats, voxel_path, n_rec):
    """<stem>_ranges.csv, <stem>_histogram.csv and, with matplotlib,
    <stem>_histogram.png next to the voxel file. Returns the paths."""
    if not stats:
        return []
    stem = voxel_path.with_suffix("")
    out = []

    p_ranges = stem.with_name(stem.name + "_ranges.csv")
    with open(p_ranges, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["variable", "units", "min"]
                   + ["p{}".format(p) for p in PERCENTILES]
                   + ["max", "suggested_min", "suggested_max"])
        for var, s in stats.items():
            w.writerow([var, s["units"], _g(s["min"])]
                       + [_g(s["pct"][p]) for p in PERCENTILES]
                       + [_g(s["max"]), _g(s["pct"][DISPLAY[0]]),
                          _g(s["pct"][DISPLAY[1]])])
    out.append(p_ranges)

    p_hist = stem.with_name(stem.name + "_histogram.csv")
    with open(p_hist, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["variable", "units", "bin_lo", "bin_hi", "count"])
        for var, s in stats.items():
            e, c = coarse(s["edges"], s["counts"])
            for a, b, n in zip(e[:-1], e[1:], c):
                w.writerow([var, s["units"], _g(a), _g(b), int(n)])
    out.append(p_hist)

    png = _plot(stats, stem.with_name(stem.name + "_histogram.png"), n_rec)
    if png:
        out.append(png)
    return out


def _g(x):
    return "{:.6g}".format(x)


def _plot(stats, path, n_rec):
    """Small multiples: one histogram per variable, each on its own value
    axis, the suggested p2 / p98 Min/Max marked and labelled."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        debug("[voxel] matplotlib not installed - no histogram plot, the "
              "CSV files hold the same numbers.")
        return None

    n = len(stats)
    ncol = min(3, n)
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.0 * nrow),
                             squeeze=False, facecolor=_SURFACE)
    for ax in axes.ravel()[n:]:
        ax.set_visible(False)

    for ax, (var, s) in zip(axes.ravel(), stats.items()):
        e, c = coarse(s["edges"], s["counts"])
        ax.set_facecolor(_SURFACE)
        width = (e[1] - e[0]) * 0.9           # surface gap between bars
        ax.bar(e[:-1], c, width=width, align="edge", color=_BAR,
               linewidth=0)
        lo, hi = s["pct"][DISPLAY[0]], s["pct"][DISPLAY[1]]
        top = c.max() if c.size else 1
        # Enough decimals to tell p2 from p98 (theta spans 299.9 .. 300.3).
        span = abs(hi - lo) or abs(hi) or 1.0
        nd = max(0, 2 - int(np.floor(np.log10(span))))
        for x, lab, ha in ((lo, "p{} {:.{}f}".format(DISPLAY[0], lo, nd),
                            "right"),
                           (hi, "p{} {:.{}f}".format(DISPLAY[1], hi, nd),
                            "left")):
            ax.axvline(x, color=_INK_2, linewidth=1.2, linestyle=(0, (4, 3)))
            ax.text(x, top * 1.02, " " + lab + " ", color=_INK_2,
                    fontsize=8, ha=ha, va="bottom")
        ax.set_ylim(0, top * 1.18)
        ax.set_title("{} [{}]".format(var, s["units"]) if s["units"]
                     else var, color=_INK, fontsize=10, loc="left")
        ax.tick_params(colors=_INK_2, labelsize=8, length=0)
        ax.grid(axis="y", color=_GRID, linewidth=0.8)
        ax.set_axisbelow(True)
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)
        ax.spines["bottom"].set_color(_INK_2)
        ax.yaxis.set_major_formatter(matplotlib.ticker.EngFormatter())
        ax.set_ylabel("cells", color=_INK_2, fontsize=8)

    fig.suptitle("Voxel value distributions over all {} records - dashed: "
                 "suggested Min/Max (p{}..p{})".format(n_rec, *DISPLAY),
                 color=_INK, fontsize=10, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(path, dpi=130, facecolor=_SURFACE)
    plt.close(fig)
    return path
