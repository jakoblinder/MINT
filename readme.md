# MINT Integrator

This project provides an example of how to use the original Fortran implementation of the MINT (Monte Carlo INTegrator) algorithm, along with utilities for grid setup, integration, event generation, and visualization of integration grids.
The MINT Fortran implementation is taken from [https://virgilio.mib.infn.it/~nason/POWHEG/FNOpaper/](https://virgilio.mib.infn.it/~nason/POWHEG/FNOpaper/).

## Features

- Adaptive Monte Carlo integration using the MINT algorithm.
- Support for multidimensional integration and event generation.
- Grid visualization using Python and matplotlib.
- Example functions and routines for physics applications.

## Structure

- `integrate_functions.f`: Main Fortran program for grid setup, integration, and event generation.
- `functions.f`: Module were the integrated functions and used mappings are defined.
- `mint-integrator.f`: MINT integrator itself.
- `random.f`: Random number generator.
- `cernroutines.f`, `random.f`: Supporting Fortran source files.
- `Makefile`: Build instructions.
- `plotgrid/plot_topdrawer_grid_classic.py`: Python script to visualize topdrawer grid files as PDF - exactly as it is written in the topdrawer file.
- `plotgrid/plot_topdrawer_grid_useful.py`:  Python script to visualize topdrawer grid files as PDF, using the information from the topdrawer file but showing more useful stuff.
- `LICENSE`: GNU GPL v3 license.

## Build Instructions

1. Ensure you have `gfortran` and `make` installed.
2. In this directory, run:
   ```sh
   make
   ```
   This will build the `integrate` executable.

## Usage

### Fortran Integration

Run the main program:
```sh
./integrate
```
This will perform grid setup and integration as configured in `integrate_functions.f`. You can modify the flags in the source to enable event generation or 2D integration.

### Grid Visualization

To plot the integration grid(s) from a topdrawer file:
```sh
python3 plotgrid/plot_topdrawer_grid.py xg*.top
```
This will create a multi-page PDF with each page showing the grid and function for one dimension.

## License

This project is licensed under the GNU General Public License v3. See the `LICENSE` file for details.
FIXME: This depends on the License under which the MINT integrator and the other not by myself written files are published.

## Contact

For questions or contributions, please contact the project maintainer.
