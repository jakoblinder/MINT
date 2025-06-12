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

def plot_grids_to_pdf(topdrawer_file, pdf_file, folding=None):
    nintervals = 50  # Default number of intervals, hardcoded in MINT.
    grids, titles, gridlines, limits, scatter_flags = parse_topdrawer_grid(topdrawer_file)
    with PdfPages(pdf_file) as pdf:
        # Adjust the list specifying the folding for each dimension:
        # If is only has a length of k but there are n dimensions,
        # it is assumed that the last n - k entries are not folded.
        if folding:
            fold_tmp = folding.copy()
            folding = [int(fold_tmp[i]) if i < len(fold_tmp) else 1 for i in range(len(grids))]
            assert all([nintervals%f == 0 for f in folding]), "The folding factors must be a divisor of the number of intervals."
        for idim, ((x, y), title, join_lines, lim, scatter) in enumerate(zip(grids, titles, gridlines, limits, scatter_flags)):
            plt.figure()
            # Plot as scatter if plot symbol is present, else as line
            if scatter:
                plt.scatter(x, y, color='blue')
                plt.grid(True)
            else:
                old_gridpoints = x
                xacc = y
                # Plot cumulative distribution of the function:
                plt.plot(old_gridpoints, xacc, color='blue', label='CDF')

                # Throw away the first nintervals points, since they are not directly part of the grid.
                join_lines = join_lines[nintervals:]
                # Get the new grid points from the join_lines array.
                new_gridpoints = [join_lines[i][0][0] for i in range(len(join_lines))]

                # Plot old_gridpoints as black vertical lines
                for i, xval in enumerate(old_gridpoints):
                    if i == 0:
                        plt.axvline(x=xval, color='black', linestyle='-', linewidth=0.8, alpha=0.45, label='Old grid')
                    else:
                        plt.axvline(x=xval, color='black', linestyle='-', linewidth=0.8, alpha=0.45)

                # Plot new_gridpoints as red vertical lines
                # If folding is specified, mark the grid points being on the edge where the new fold starts in orange.
                # E.g.: folding=[1, 2, 5, 10] -> Mark the 50/2=25th, 50/5=10th, and 50/10=5th grid points
                #       of the 2nd, 3rd and 4th dimension in purple.
                for i, xval in enumerate(new_gridpoints):
                    if i == 0:
                        plt.axvline(x=xval, color='red', linestyle='-', linewidth=0.8, alpha=0.6, label='New grid')
                    else:
                        if folding[idim] > 1 and folding[idim] < nintervals and i % (nintervals / folding[idim]) == 0:
                            plt.axvline(x=xval, color="#EE00FFFF", linestyle='-', linewidth=2.0, alpha=1.0)
                        else:
                            plt.axvline(x=xval, color='red', linestyle='-', linewidth=0.8, alpha=0.6)

                # Try to restore the original function:
                dxaccdx = []
                for i in range(len(xacc)-1):
                    dxaccdx.append( (xacc[i+1] - xacc[i])/ (old_gridpoints[i+1] - old_gridpoints[i]) )

                # Normalise dxaccdx to the maximum value, so that it fits in the plot.
                maxdxaccdx = max(dxaccdx)
                dxaccdx = [dxaccdx[i]/ maxdxaccdx for i in range(len(dxaccdx))]

                old_midths = [(old_gridpoints[i] + old_gridpoints[i+1]) / 2 for i in range(len(old_gridpoints)-1)]
                plt.plot(old_midths, dxaccdx, color='green', label='~ PDF')
            plt.xlabel('Grid points $x_i$')
            # plt.ylabel('y')
            plt.title(title)
            # Set axis limits if available
            if lim and lim[0] and lim[1]:
                plt.xlim(*lim[0])
                plt.ylim(*lim[1])
            plt.legend()
            pdf.savefig()
            plt.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot MINT topdrawer grid file to PDF (one page per dimension).")

    parser.add_argument(
        "input",
        nargs="+",
        type=Path,
        help="Input topdrawer file(s) for grid, supports glob patterns (e.g. pwg-xg2-xgrid-btl-*.top)"
    )
    parser.add_argument(
        "-f", "--folding",
        nargs="*",
        type=float,
        default=None,
        help="Optional list of numbers (one per dimension or less)."
    )

    args = parser.parse_args()

    for infile in args.input:
        output = infile.with_suffix('.pdf')
        plot_grids_to_pdf(infile, output, args.folding)
