# Third-party code: POWHEG BOX / MINT

These files were not written by the author of this repository. They are redistributed under the GNU General Public
License version 2, as stated in the [POWHEG license note](https://powhegbox.mib.infn.it/AAAREADME-LICENSE), and are
subject to the [MCNET guidelines](https://powhegbox.mib.infn.it/MCNET-GUIDELINES). Copyright notices and references in
the files must not be removed.

| File | Origin | Status |
| --- | --- | --- |
| `mint-integrator.f` | MINT integrator by Paolo Nason, [arXiv:0709.2085](https://arxiv.org/abs/0709.2085), [POWHEG MINT page](https://virgilio.mib.infn.it/~nason/POWHEG/FNOpaper/) | **modified**, see below |
| `random.f` | [POWHEG BOX](https://powhegbox.mib.infn.it/) (random number generator based on RM48 from CERNLIB) | not modified here |
| `cernroutines.f` | [POWHEG BOX](https://powhegbox.mib.infn.it/), routines from CERNLIB (DDILOG, DZERO, RM48, ...) and a few POWHEG utilities | not modified here |

The exact POWHEG BOX release the files were taken from is not recorded; `random.f` and `cernroutines.f` are included as
they were obtained.

## Modifications of `mint-integrator.f`

These are also listed at the top of the file, as required by the MCNET guidelines (reason, author, date):

- Jakob Linder, 2025-06-06: comments converted to free form and extended, typo fixes, `x` is initialised, debugging
  output of `vtot`, `etot` removed.
- Jakob Linder, 2025-06-11: the grid of every iteration is written to its own topdrawer file `xg<iteration>.top`, to
  visualise the grid adaption with `plotgrid/plot_topdrawer_grid.py`.
- Jakob Linder, 2026-10-07: maximum number of dimensions `ndimmax` raised from 6 to 20 and defined once in
  `src/ndimmax.inc`, to allow integrals in up to 20 dimensions in the example.

## References

- P. Nason, *MINT: a computer program for adaptive Monte Carlo integration and generation of unweighted
  distributions*, [arXiv:0709.2085](https://arxiv.org/abs/0709.2085)
- [arXiv:0709.2092](https://arxiv.org/abs/0709.2092)
