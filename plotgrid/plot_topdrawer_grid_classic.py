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
    grids, titles, gridlines, limits, scatter_flags = parse_topdrawer_grid(topdrawer_file)
    with PdfPages(pdf_file) as pdf:
        for (x, y), title, join_lines, lim, scatter in zip(grids, titles, gridlines, limits, scatter_flags):
            plt.figure()
            # Plot as scatter if plot symbol is present, else as line
            if scatter:
                plt.scatter(x, y, color='blue')
            else:
                plt.plot(x, y, color='blue')
            # Plot each grid line segment in black, alpha=0.5, no marker
            for seg in join_lines:
                (x0, y0), (x1, y1) = seg
                plt.plot([x0, x1], [y0, y1], color='black', alpha=0.45, linestyle='-', linewidth=0.8)
            plt.xlabel('x')
            plt.ylabel('y')
            plt.title(title)
            plt.grid(True)
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
        nargs="+",
        type=Path,
        help="Input topdrawer file(s) for grid, supports glob patterns (e.g. pwg-xg2-xgrid-btl-*.top)"
    )

    args = parser.parse_args()

    for infile in args.input:
        output = infile.with_suffix('.pdf')
        plot_grids_to_pdf(infile, output)
