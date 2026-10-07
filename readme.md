# MINT Integrator example

A small, self-contained example of how to use the Fortran implementation of the MINT (Monte Carlo INTegrator) algorithm
for adaptive multidimensional integration and unweighted event generation, together with a script to visualise the
adapted integration grids.

MINT is described in [arXiv:0709.2085](https://arxiv.org/abs/0709.2085) and is a variant of the Vegas algorithm
([Lepage 1978](https://doi.org/10.1016/0021-9991(78)90004-9), [manual](https://inspirehep.net/literature/153221)).
For the hit-and-miss event generation see [arXiv:0709.2092](https://arxiv.org/abs/0709.2092); for a gentler
introduction to it see [hep-ph/0606275](https://arxiv.org/abs/hep-ph/0606275).
Notes of a journal club talk on the topic are in `docs/JC_Notes.pdf`.

## Structure

| Path | Content |
| --- | --- |
| `src/integrate_functions.f` | Main program: grid setup, integration, 2D integration, event generation (select with the flags at the top). |
| `src/functions.f` | The integrand and the variable mappings; edit the preset block in `func_wrap` to change the integrand. |
| `src/ndimmax.inc` | Maximum number of dimensions, shared by `functions.f` and `mint-integrator.f`. |
| `third_party/powheg/` | MINT integrator, random number generator and CERNLIB routines (not written by me, see below). |
| `plotgrid/` | Script to plot the grid files `xg*.top` written by MINT as PDF, and example grids. |
| `docs/JC_Notes.pdf` | Handwritten notes of a journal club talk on the topic. |
| `Makefile` | Build instructions. |

## Build and run

Requires `gfortran` and `make`.

```sh
make            # optimised build
make DEBUG=1    # debug build with floating point traps and run time checks
./integrate
```

The number of dimensions `ndim` is set in `src/integrate_functions.f` (at most `ndimmax`, set in `src/ndimmax.inc`, default 20).
With `flg_sphere = .true.` the volume of an `ndim`-dimensional ball is integrated and compared with the exact value;
otherwise the preset in `func_wrap` is used.

> **Note:** the sphere only works reliably up to `ndim = 12`. MINT can only adapt its grid to points that hit the ball,
> and their fraction of the cube shrinks quickly with `ndim`. For `ndim >= 15` the grid setup does not find the ball
> (even with 3e8 calls per iteration at `ndim = 20`) and the result is wrong or NaN.

By default only the grid setup is run. Enable the integration, 2D integration or event generation by setting the
`flg_*` flags in `src/integrate_functions.f` and recompiling. The run writes one grid file `xg<i>.top` per grid iteration;
generated events go to `events.dat`.

## Grid visualisation

Requires Python 3 with `matplotlib`.

```sh
python3 plotgrid/plot_topdrawer_grid.py xg*.top              # PDFs next to the input files, one page per dimension
python3 plotgrid/plot_topdrawer_grid.py xg5.top -f 2 5 10    # also mark where the folding blocks start
python3 plotgrid/plot_topdrawer_grid.py --classic xg5.top    # plot exactly what the topdrawer file describes
```

Each page shows the cumulative distribution, the old and the new grid, and a reconstruction of the integrand
(`~ PDF`, normalised to its maximum). MINT integrates the absolute value of the function, so this is |f|, not f.
The folding factors (divisors of 50) are given in reverse dimension order, as in POWHEG; see `--help`.

Example: `plotgrid/example_grid/` contains six grid files of a POWHEG-BOX run (7 dimensions each) and the resulting PDFs,
created with

```sh
python3 plotgrid/plot_topdrawer_grid.py plotgrid/example_grid/*.top
```

Adding `-f 2 5 10` additionally marks the folding blocks of the last three dimensions (here dimensions 5 to 7).

## Third-party code and license

`third_party/powheg/` contains `mint-integrator.f`, `random.f` and `cernroutines.f` from the
[POWHEG BOX](https://powhegbox.mib.infn.it/) / [MINT page](https://virgilio.mib.infn.it/~nason/POWHEG/FNOpaper/) by
Paolo Nason and collaborators (`cernroutines.f` consists of routines from CERNLIB). They keep their original copyright
and are distributed under the GPL version 2 and the [MCNET guidelines](https://powhegbox.mib.infn.it/MCNET-GUIDELINES),
see [the license note](https://powhegbox.mib.infn.it/AAAREADME-LICENSE). If you use this code, please cite
[arXiv:0709.2085](https://arxiv.org/abs/0709.2085) and [arXiv:0709.2092](https://arxiv.org/abs/0709.2092).
**`mint-integrator.f` is a modified version**; details are in `third_party/powheg/README.md` and in the file header.

The remaining files are licensed under the GNU General Public License version 2 (see `LICENSE`).
