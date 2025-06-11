import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pathlib import Path
import re
import argparse

def parse_topdrawer_grid(filename):
    grids = []
    dim_titles = []
    gridlines = []
    limits = []
    scatter_flags = []
    with open(filename) as f:
        lines = f.readlines()
    x, y = [], []
    join_lines = []
    in_grid = False
    join_points = []
    xlim, ylim = None, None
    scatter = False
    for line in lines:
        if line.strip().startswith('set limits'):
            m = re.search(
                r'set limits x (?P<xmin>[\d\.Ee\+\-]+) (?P<xmax>[\d\.Ee\+\-]+) y (?P<ymin>[\d\.Ee\+\-]+) (?P<ymax>[\d\.Ee\+\-]+)',
                line
            )
            if m:
                xlim = (float(m.group('xmin')), float(m.group('xmax')))
                ylim = (float(m.group('ymin')), float(m.group('ymax')))
        elif line.strip().startswith('title top'):
            if x and y:
                grids.append((x, y))
                gridlines.append(join_lines)
                limits.append((xlim, ylim))
                scatter_flags.append(scatter)
                x, y = [], []
                join_lines = []
                scatter = False
            match = re.search(r'dim=\s*(?P<dim>\d+)', line)
            if match:
                dim_titles.append(f"Dimension {match.group('dim')}")
            else:
                match2 = re.match(r'.*\"(?P<title>.*)\"', line)
                dim_titles.append(f"{match2.group('title')}" if match2 else "Grid Plot")
        elif line.strip().startswith('plot symbol'):
            scatter = True
        elif line.strip().startswith('join'):
            in_grid = True
            join_points = []
            continue
        elif in_grid and re.match(r'^\s*[\d\.\-E\+]+', line):
            parts = line.split()
            if len(parts) == 2:
                try:
                    join_points.append((float(parts[0]), float(parts[1])))
                except ValueError:
                    continue
            if len(join_points) == 2:
                join_lines.append(tuple(join_points))
                join_points = []
        elif re.match(r'^\s*[\d\.\-E\+]+', line):
            parts = line.split()
            if len(parts) == 2:
                try:
                    x.append(float(parts[0]))
                    y.append(float(parts[1]))
                except ValueError:
                    continue
        elif line.strip().startswith('newplot'):
            if x and y:
                grids.append((x, y))
                gridlines.append(join_lines)
                limits.append((xlim, ylim))
                scatter_flags.append(scatter)
                x, y = [], []
                join_lines = []
                scatter = False
            in_grid = False
    if x and y:
        grids.append((x, y))
        gridlines.append(join_lines)
        limits.append((xlim, ylim))
        scatter_flags.append(scatter)
    return grids, dim_titles, gridlines, limits, scatter_flags

def plot_grids_to_pdf(topdrawer_file, pdf_file):
    nintervals = 50  # Default number of intervals, hardcoded in MINT.
    grids, titles, gridlines, limits, scatter_flags = parse_topdrawer_grid(topdrawer_file)
    with PdfPages(pdf_file) as pdf:
        for (x, y), title, join_lines, lim, scatter in zip(grids, titles, gridlines, limits, scatter_flags):
            plt.figure()
            # Plot as scatter if plot symbol is present, else as line
            if scatter:
                plt.scatter(x, y, color='blue')
                plt.grid(True)
            else:
                old_gridpoints = x
                xacc = y
                # Plot cumulative distribution of the function:
                plt.plot(old_gridpoints, xacc, color='blue')

                # Throw away the first nintervals points, since they are not directly part of the grid.
                join_lines = join_lines[nintervals:]
                new_gridpoints = [join_lines[i][0][0] for i in range(len(join_lines))]
                # # Plot each grid line segment in black, alpha=0.5, no marker
                # for seg in join_lines:
                #     (x0, y0), (x1, y1) = seg
                #     plt.plot([x0, x1], [y0, y1], color='black', alpha=0.45, linestyle='-', linewidth=0.8)
                # plt.grid(True)
                # Plot old_gridpoints as black vertical lines
                for xval in old_gridpoints:
                    plt.axvline(x=xval, color='black', linestyle='-', linewidth=0.8, alpha=0.45)
                for xval in new_gridpoints:
                    plt.axvline(x=xval, color='red', linestyle='-', linewidth=0.8, alpha=0.6)
                # Try to restore the original function:
                dxaccdx = []
                for i in range(len(xacc)-1):
                    dxaccdx.append( (xacc[i+1] - xacc[i])/ (old_gridpoints[i+1] - old_gridpoints[i]) )

                maxdxaccdx = max(dxaccdx)
                dxaccdx = [dxaccdx[i]/ maxdxaccdx for i in range(len(dxaccdx))]

                old_midths = [(old_gridpoints[i] + old_gridpoints[i+1]) / 2 for i in range(len(old_gridpoints)-1)]
                plt.plot(old_midths, dxaccdx, color='green')
            plt.xlabel('x')
            plt.ylabel('y')
            plt.title(title)
            # Set axis limits if available
            if lim and lim[0] and lim[1]:
                plt.xlim(*lim[0])
                plt.ylim(*lim[1])
            pdf.savefig()
            plt.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot MINT topdrawer grid file to PDF (one page per dimension).")
    parser.add_argument(
        "input",
        nargs="?",
        default="pwg-xg2-xgrid-btl-0050.top",
        type=Path,
        help="Input topdrawer file for grid, e.g.: pwg-xg2-xgrid-btl-0001.top"
    )
    parser.add_argument(
        "output",
        nargs="?",
        default=None,
        type=Path,
        help="Output PDF file"
    )
    args = parser.parse_args()

    if not args.output:
        args.output = args.input.with_suffix('.pdf')

    plot_grids_to_pdf(args.input, args.output)
