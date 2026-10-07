# MINT Integrator example

A small, self-contained example of how to use the Fortran implementation of the MINT (Monte Carlo INTegrator) algorithm
for adaptive multidimensional integration and unweighted event generation, together with a script to visualise the
adapted integration grids.

MINT is described in [arXiv:0709.2085](https://arxiv.org/abs/0709.2085) and is a variant of the Vegas algorithm
([Lepage 1978](https://doi.org/10.1016/0021-9991(78)90004-9), [manual](https://inspirehep.net/literature/153221)).
For the hit-and-miss event generation see [arXiv:0709.2092](https://arxiv.org/abs/0709.2092); for a gentler
introduction to it see [hep-ph/0606275](https://arxiv.org/abs/hep-ph/0606275).
`JC_Notes.pdf` are handwritten notes of a journal club talk on the topic.

## Files

| File | Content |
| --- | --- |
| `integrate_functions.f` | Main program: grid setup, integration, 2D integration, event generation (select with the flags at the top). |
| `functions.f` | The integrand and the variable mappings; edit the preset block in `func_wrap` to change the integrand. |
| `ndimmax.inc` | Maximum number of dimensions, shared by `functions.f` and `mint-integrator.f`. |
| `Makefile` | Build instructions. |
| `plotgrid/plot_topdrawer_grid.py` | Plots the grid files `xg*.top` written by MINT as PDF. |
| `mint-integrator.f`, `random.f`, `cernroutines.f` | Third-party code (`mint-integrator.f` modified), see below. |

## Build and run

Requires `gfortran` and `make`.

```sh
make            # optimised build
make DEBUG=1    # debug build with floating point traps and run time checks
./integrate
```

The number of dimensions `ndim` is set in `integrate_functions.f` (at most `ndimmax`, set in `ndimmax.inc`, default 20).
With `flg_sphere = .true.` the volume of an `ndim`-dimensional ball is integrated and compared with the exact value;
otherwise the preset in `func_wrap` is used.

> **Note:** the sphere only works reliably up to `ndim = 12`. MINT can only adapt its grid to points that hit the ball,
> and their fraction of the cube shrinks quickly with `ndim`. For `ndim >= 15` the grid setup does not find the ball
> (even with 3e8 calls per iteration at `ndim = 20`) and the result is wrong or NaN.

By default only the grid setup is run. Enable the integration, 2D integration or event generation by setting the
`flg_*` flags in `integrate_functions.f` and recompiling. The run writes one grid file `xg<i>.top` per grid iteration;
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

`mint-integrator.f`, `random.f` and `cernroutines.f` are taken from the
[POWHEG BOX](https://powhegbox.mib.infn.it/) / [MINT page](https://virgilio.mib.infn.it/~nason/POWHEG/FNOpaper/) by
Paolo Nason and collaborators; `cernroutines.f` consists of routines from CERNLIB. They keep their original copyright.
The POWHEG software is distributed under the GPL version 2 and subject to the
[MCNET guidelines](https://powhegbox.mib.infn.it/MCNET-GUIDELINES), see
[the license note](https://powhegbox.mib.infn.it/AAAREADME-LICENSE). If you use this code, please cite
[arXiv:0709.2085](https://arxiv.org/abs/0709.2085) and [arXiv:0709.2092](https://arxiv.org/abs/0709.2092).

- `random.f` and `cernroutines.f` are included as obtained, not modified in this repository.
- **`mint-integrator.f` is a modified version.** The changes (extended comments, one grid file `xg<i>.top` per
  iteration for the plots, and `ndimmax` raised from 6 to 20, set in `ndimmax.inc`) are listed in the header of the file.

The remaining files are licensed under the GNU General Public License version 2 (see `LICENSE`).
