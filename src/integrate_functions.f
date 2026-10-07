! Copyright (C) 2025 Jakob Linder
! SPDX-License-Identifier: GPL-2.0-only
      program integrate_functions
         ! Example driver for the MINT integrator: integrates the function defined in functions.f
         ! (func_wrap) in ndim dimensions. Choose the steps to run with the flags below.
         use functions
         implicit none
         integer :: ndim
         common/cndim/ndim
         integer :: ncall1, itmx1, ncall2, itmx2, nevents, ievent, i
         real*8 :: xgrid(0:50,ndimmax), xint, ymax(50,ndimmax)
         real*8 :: intabs_val, intabs_err, estimn, errorn, estimp, errorp
         real*8 :: xgenerated(ndimmax), xtransformed(ndimmax), efficiency(ndimmax)
         real*8 :: xmin, xmax
         common/bounds/xmin, xmax

         integer ifold(ndimmax)
         common/cifold/ifold
         logical negflag ! If true, the function is returning only non zero if it is negative and vice versa.
         common/cnegflag/negflag

         real*8 :: x12(2)
         common/xset/x12

         character(len=30) :: map_type
         common /maptypeblock/ map_type

         character(len=*), parameter :: file_events = 'events.dat'
         integer :: unit_events

         ! Timer variables
         real*8 :: t_start, t_end

         logical :: flg_gridsetup, flg_integration, flg_2dintegration, flg_generation

         ! Number of dimensions of the integral, at most ndimmax (set in ndimmax.inc).
         ndim = 4

         ! Integrate the indicator function of the ndim dimensional ball (sphere) with radius sphere_radius
         ! instead of the function selected in func_wrap. The exact result is the volume of the ball.
         ! Note: this only works reliably up to ndim = 12. The fraction of the cube [-r, r]^ndim covered by the
         ! ball falls quickly with ndim (about 2.5e-8 for ndim = 20), and MINT can only adapt its grid to points
         ! that hit the ball. For larger ndim the result is wrong or NaN.
         flg_sphere    = .false.
         sphere_radius = 1d0

         if (ndim < 1 .or. ndim > ndimmax) stop 'ndim has to be between 1 and ndimmax'

         ! Steps to run (each step needs the grid, so the grid setup is always run):
         flg_gridsetup     = .true.  ! Adapt the grid to |f| and compute the integral of |f|.
         flg_integration   = .false. ! Integrate f itself (positive and negative part separately) on the frozen grid.
         flg_2dintegration = .false. ! Fix the first two variables and integrate over the other two.
         flg_generation    = .false. ! Generate unweighted events (implies the integration step).

         if (flg_gridsetup .or. flg_integration .or. flg_2dintegration .or. flg_generation) then
            ! Number of points to improve the grid:
            ncall1 = 1d5
            ! Number of grid improvement iterations (MaXimum number of ITerations to improve the grid):
            itmx1  = 5

            ! Set up the grid
            call cpu_time(t_start)

            call mint(func_wrap, ndim, ncall1, itmx1, 0, xgrid, xint, ymax, intabs_val, intabs_err)

            call cpu_time(t_end)
            write(*,*) 'Integral over the absolute value of the function:'
            write(*,'(A,G14.6,A,G12.6)') 'Int[ |f| ]: ', intabs_val, ' +- ', intabs_err
            ! Note: xint == intabs_val, since the function is integrated over the absolute value and xint is used as
            !       the initial value for computing the upper bound of the function in the next, imode = 1, step,
            !       where xint will be an input, not an output.
            write(*,'(A,F8.3,A)') 'Grid setup time: ', t_end - t_start, ' s'
            if (flg_sphere) write(*,'(A,I0,A,G14.6)') 'Exact volume of the ', ndim, '-ball: ', ball_volume(ndim, sphere_radius)
         end if


         if (flg_integration .or. flg_generation) then
            ! Setup folding:
            ! Notice that you can only fold by divisors of 50, since this is the number of points in the grid.
            ! Moreover, folding one dimension by e.g. 5 means the integration will take 5 times longer, since it will
            ! call the integrated function 5 times more often in that dimension. Folding by 5 in 2 dimensions will
            ! take 25 times longer, etc..
            ifold(:) = 1

            ! Number of points used for the integration:
            ncall2 = 1d6
            ! Number of integration iterations and upper bound improvements, all done with a number of calls ncall2 to the
            ! integrated function. The different integrand results are combined and only the final result is returned.
            ! Note that this basically corresponds to increasing the number of calls ncall2 by a factor of itmx2.
            itmx2 = 5

            call cpu_time(t_start)

            ! Compute the negative contribution to the integral:
            negflag = .true.
            call mint(func_wrap, ndim, ncall2, itmx2, 1, xgrid, xint, ymax, estimn, errorn)
            ! Compute the positive contribution to the integral. This is done last on purpose: every call of mint
            ! with imode=1 overwrites the upper bounds ymax, and the event generation below needs those of the
            ! positive part.
            negflag = .false.
            call mint(func_wrap, ndim, ncall2, itmx2, 1, xgrid, xint, ymax, estimp, errorp)

            call cpu_time(t_end)

            write(*,'(A,G14.6,A,G12.6)') 'Positive contribution: ', estimp, ' +- ', errorp
            write(*,'(A,G14.6,A,G12.6)') 'Negative contribution: ', estimn, ' +- ', errorn
            write(*,'(A,G14.6,A,G12.6)') 'Total integral:        ', estimp + estimn, ' +- ', sqrt(errorp**2 + errorn**2)

            write(*,'(A,F8.3,A)') 'Integration time: ', t_end - t_start, ' s'
         end if


         if (flg_generation .and. estimp <= 0d0) then
            ! gen would loop forever, since the upper bounds ymax of the positive part are zero.
            write(*,*) 'No positive contribution found, skipping the event generation.'
            flg_generation = .false.
         end if

         if (flg_generation) then
            ! Generate unweighted points distributed according to the positive part of the function (hit and miss).
            ! The routine gen needs the grid and the upper bounds ymax, both set up above.
            ! Arguments of gen(fun, ndim, xgrid, ymax, imode, x):
            !   imode = 0: initialise, imode = 1: generate one point x(1:ndim) in [0,1]^ndim,
            !   imode = 3: return the generation efficiency in x(1).
            nevents = 1000
            write(*,*) 'Generating ', nevents, ' events and writing them to file: ', file_events
            open(newunit=unit_events, file=file_events, status='replace')

            call gen(func_wrap, ndim, xgrid, ymax, 0, xgenerated)
            do ievent = 1, nevents
               call gen(func_wrap, ndim, xgrid, ymax, 1, xgenerated)
               ! Map the generated numbers from [0,1] to the integration region:
               do i = 1, ndim
                  xtransformed(i) = mapping(xgenerated(i), map_type)
               end do
               write(unit_events,*) xtransformed
            end do
            close(unit_events)

            call gen(func_wrap, ndim, xgrid, ymax, 3, efficiency)
            write(*,'(A,F8.4)') 'Generation efficiency: ', efficiency(1)
         end if


         if (flg_2dintegration) then
            ! Since we are able to integrate all dimensions, what about setting the first two to a specific
            ! value and integrating over the remaining ndim-2?
            write(*,*) "Fixing the first two variables, integrating over the remaining ", ndim-2
            if (ndim < 3) stop 'the 2D integration needs ndim >= 3'
            ! Set value of first two x:
            x12 = [0.45d0, 0.55d0]
            ! Note that we remapped all numbers with
            ! xmin + (xmax - xmin) * x(i),
            ! so that even if x(i) lies between 0 and 1, we can still set arbitrary integral boundaries
            ! - they become very fast very unstable though.
            ! This mapping should be inverted here if we want to set the first two values to a specific value:
            x12 = [inverseMapping(x12(1), map_type), inverseMapping(x12(2), map_type)]

            ! The integration grid can be kept! It should be optimal also for the remaining dimensions.
            ! Only redo the integration for the remaining dimensions.
            ifold(:) = 1
            ncall2 = 1000000
            itmx2 = 5

            call cpu_time(t_start)
            negflag = .true.
            call mint(func2d, ndim-2, ncall2, itmx2, 1, xgrid, xint, ymax, estimn, errorn)
            negflag = .false.
            call mint(func2d, ndim-2, ncall2, itmx2, 1, xgrid, xint, ymax, estimp, errorp)
            call cpu_time(t_end)

            write(*,'(A,G14.6,A,G12.6)') 'Positive contribution: ', estimp, ' +- ', errorp
            write(*,'(A,G14.6,A,G12.6)') 'Negative contribution: ', estimn, ' +- ', errorn
            write(*,'(A,G14.6,A,G12.6)') 'Total integral:        ', estimp + estimn, ' +- ', sqrt(errorp**2 + errorn**2)
            write(*,'(A,F8.3,A)') 'Integration time: ', t_end - t_start, ' s'
         end if

         contains
            function func2d(x,weight,ifl)
               ! Version of func_wrap in the style needed by mint with the first two variables fixed to x12 (given
               ! in a common block, since mint cannot pass additional arguments to the function); x are the
               ! remaining ndim-2 variables.
               implicit none
               real*8, intent(in) :: x(*), weight
               integer, intent(in) :: ifl
               real*8 :: func2d, jac
               integer :: i
               integer :: ndim
               common/cndim/ndim

               real*8 :: x12(2)
               common/xset/x12

               character(len=30) :: map_type
               common /maptypeblock/ map_type

               func2d = func_wrap([x12(1), x12(2), x(1:ndim-2)], weight, ifl)
               ! func_wrap multiplies the jacobian of all variables, but the first two are not integrated over.
               jac = 1d0
               do i=1, size(x12)
                  jac = jac * jacobian(x12(i), map_type)
               end do
               func2d = func2d / jac
            end function func2d

      end program integrate_functions
