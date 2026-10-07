# Copyright (C) 2025 Jakob Linder
# SPDX-License-Identifier: GPL-2.0-only
# Usage:
#   make            optimised build of ./integrate
#   make DEBUG=1    debug build (floating point traps, run time checks, backtraces)
#   make clean
#
# -fno-automatic is required by the legacy Fortran 77 code (variables are static).
# Note: -fbounds-check causes a weird error due to the non-lazy evaluation of booleans in gfortran.
FC     = gfortran
FFLAGS = -fno-automatic -ffixed-line-length-none

ifdef DEBUG
FFLAGS += -O -g -ggdb -ffpe-trap=invalid,zero,overflow -fcheck=all -finit-real=nan -fbacktrace -Wall
else
FFLAGS += -O2 -w  # -w: silence warnings of the legacy third-party files
endif

# Sources: own code in src/, third-party POWHEG code in third_party/powheg/
vpath %.f src third_party/powheg

OBJS = integrate_functions.o functions.o mint-integrator.o cernroutines.o random.o

.PHONY: all clean

all: integrate

integrate: $(OBJS)
	$(FC) $(FFLAGS) $^ -o $@

%.o: %.f
	$(FC) $(FFLAGS) -c $<

# integrate_functions uses the module defined in functions.f
integrate_functions.o: functions.o

# the maximum number of dimensions is set in src/ndimmax.inc
functions.o mint-integrator.o: src/ndimmax.inc

clean:
	rm -f *.o *.mod integrate
