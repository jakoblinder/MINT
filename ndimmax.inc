! Copyright (C) 2025 Jakob Linder
! SPDX-License-Identifier: GPL-2.0-only
!
! Maximum number of dimensions of the integral, the only place where it is set.
! Included by mint-integrator.f (mint and gen) and functions.f; integrate_functions.f gets it from there.
      integer ndimmax
      parameter (ndimmax=20)
