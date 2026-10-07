#!/usr/bin/env python3
# Copyright (C) 2025 Jakob Linder
# SPDX-License-Identifier: GPL-2.0-only
"""Plot the MINT grid files (xg*.top, topdrawer format) as multi-page PDFs.

One PDF page is produced per dimension. Two styles are available:

  default    Cumulative distribution of the integrand, the old grid, the new
             grid, a rough reconstruction of the integrand ("~ PDF") and,
             optionally, the points where the folding blocks start.
  --classic  Exactly what the topdrawer file describes (the plot topdrawer
             itself would draw).

Note: MINT integrates |f|, so the reconstructed function is always the
absolute value of the integrand.
"""
import argparse
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

NINTERVALS = 50  # Number of grid intervals per dimension, hardcoded in MINT.


@dataclass
class Grid:
    title: str = "Grid Plot"
    x: list = field(default_factory=list)
    y: list = field(default_factory=list)
    joins: list = field(default_factory=list)  # [((x0, y0), (x1, y1)), ...]
    xlim: tuple = None
    ylim: tuple = None
    scatter: bool = False


def _two_floats(line):
    """Return (a, b) if the line consists of exactly two numbers, else None."""
    parts = line.split()
    if len(parts) != 2:
        return None
    try:
        return float(parts[0]), float(parts[1])
    except ValueError:
        return None


def parse_topdrawer_grid(filename):
    """Parse a MINT topdrawer grid file into a list of Grid objects (one per dimension)."""
    grids = []
    current = Grid()
    limits = None
    in_join = False
    pending = []  # first point of a join segment

    def flush():
        nonlocal current
        if current.x and current.y:
            current.xlim, current.ylim = limits if limits else (None, None)
            grids.append(current)
        current = Grid()

    for raw in Path(filename).read_text().splitlines():
        line = raw.strip()
        if line.startswith("set limits"):
            tok = line.split()  # set limits x xmin xmax y ymin ymax
            try:
                limits = ((float(tok[3]), float(tok[4])), (float(tok[6]), float(tok[7])))
            except (IndexError, ValueError):
                pass
        elif line.startswith("title top"):
            flush()
            in_join = False
            if "dim=" in line:
                current.title = f"Dimension {line.split('dim=')[1].split()[0].strip(chr(34))}"
            elif line.count('"') >= 2:
                current.title = line.split('"')[1]
        elif line.startswith("plot symbol"):
            current.scatter = True
        elif line.startswith("join"):
            in_join = True
            pending = []
        elif line.startswith("newplot"):
            flush()
            in_join = False
        else:
            point = _two_floats(line)
            if point is None:
                continue
            if in_join:
                pending.append(point)
                if len(pending) == 2:
                    current.joins.append(tuple(pending))
                    pending = []
            else:
                current.x.append(point[0])
                current.y.append(point[1])
    flush()
    return grids


def folding_per_dimension(folding, ndim):
    """Translate the command line folding list into one factor per dimension.

    The first entry belongs to the LAST dimension, the second to the second to
    last, etc.; missing entries (the first dimensions) are 1. This mirrors
    POWHEG, where only the last three dimensions can be folded: `-f 2 5 10` is
    ifold = [1, 10, 5, 2] for 4 dimensions.
    """
    folding = list(folding or [])
    if len(folding) > ndim:
        raise ValueError(f"{len(folding)} folding factors given, but the grid has only {ndim} dimensions.")
    if any(f < 1 or NINTERVALS % f for f in folding):
        raise ValueError(f"Folding factors must be divisors of {NINTERVALS}.")
    return (folding + [1] * ndim)[:ndim][::-1]


def _plot_classic(ax, grid):
    if grid.scatter:
        ax.scatter(grid.x, grid.y, color="blue")
    else:
        ax.plot(grid.x, grid.y, color="blue")
    for (x0, y0), (x1, y1) in grid.joins:
        ax.plot([x0, x1], [y0, y1], color="black", alpha=0.45, linestyle="-", linewidth=0.8)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(True)


def _plot_useful(ax, grid, fold):
    if grid.scatter:
        ax.scatter(grid.x, grid.y, color="blue")
        ax.grid(True)
        return

    old_points, cdf = grid.x, grid.y
    ax.plot(old_points, cdf, color="blue", label="CDF")

    # The first NINTERVALS join segments redraw the old grid; the rest give the new grid.
    new_points = [seg[0][0] for seg in grid.joins[NINTERVALS:]]

    for i, xval in enumerate(old_points):
        ax.axvline(xval, color="black", linewidth=0.8, alpha=0.45, label="Old grid" if i == 0 else None)

    # New grid in red. With folding, points where a new folding block starts are
    # highlighted: e.g. fold=5 marks the 10th, 20th, ... grid point.
    fold_label_used = False
    for i, xval in enumerate(new_points):
        if 1 < fold < NINTERVALS and i % (NINTERVALS // fold) == 0 and i > 0:
            ax.axvline(xval, color="#EE00FF", linewidth=2.0,
                       label=None if fold_label_used else f"Fold block start (x{fold})")
            fold_label_used = True
        else:
            ax.axvline(xval, color="red", linewidth=0.8, alpha=0.6, label="New grid" if i == 0 else None)

    # Reconstruct the (absolute) integrand as the derivative of the CDF, normalised to 1.
    slope = [(cdf[i + 1] - cdf[i]) / (old_points[i + 1] - old_points[i]) for i in range(len(cdf) - 1)]
    peak = max(slope)
    midpoints = [(old_points[i] + old_points[i + 1]) / 2 for i in range(len(old_points) - 1)]
    ax.plot(midpoints, [s / peak for s in slope], color="green", label="~ PDF")
    ax.set_xlabel("Grid points $x_i$")


def plot_grids_to_pdf(topdrawer_file, pdf_file, classic=False, folding=None):
    grids = parse_topdrawer_grid(topdrawer_file)
    if not grids:
        raise ValueError(f"No grid found in {topdrawer_file}.")
    folds = folding_per_dimension(folding, len(grids))
    with PdfPages(pdf_file) as pdf:
        for grid, fold in zip(grids, folds):
            fig, ax = plt.subplots()
            if classic:
                _plot_classic(ax, grid)
            else:
                _plot_useful(ax, grid, fold)
            ax.set_title(grid.title)
            if grid.xlim and grid.ylim:
                ax.set_xlim(*grid.xlim)
                ax.set_ylim(*grid.ylim)
            if not classic:
                ax.legend()
            pdf.savefig(fig)
            plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Plot MINT topdrawer grid files (xg*.top) to PDF, one page per dimension.")
    parser.add_argument("input", nargs="+", type=Path,
                        help="topdrawer grid file(s), e.g. xg*.top")
    parser.add_argument("--classic", action="store_true",
                        help="plot exactly what the topdrawer file describes")
    parser.add_argument("-f", "--folding", nargs="+", type=int, default=None, metavar="N",
                        help=("highlight where the folding blocks start. Give one divisor of 50 per "
                              "dimension, in REVERSE order: the first number is the folding factor of "
                              "the last dimension, the second of the second to last, etc.; dimensions "
                              "not listed are unfolded. Example: '-f 2 5 10' for a 4D grid means "
                              "ifold = (1, 10, 5, 2). Ignored with --classic."))
    parser.add_argument("-o", "--outdir", type=Path, default=None,
                        help="output directory (default: next to each input file)")
    args = parser.parse_args()

    if args.outdir:
        args.outdir.mkdir(parents=True, exist_ok=True)
    for infile in args.input:
        outfile = (args.outdir or infile.parent) / infile.with_suffix(".pdf").name
        try:
            plot_grids_to_pdf(infile, outfile, classic=args.classic, folding=args.folding)
        except (OSError, ValueError) as err:
            parser.error(f"{infile}: {err}")
        print(f"Wrote {outfile}")


if __name__ == "__main__":
    main()
