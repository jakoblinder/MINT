      program integrate_functions
         use functions
         implicit none
         integer :: i, j, pdim, ndim, nevents, ncall1, itmx1, ncall2, itmx2
         parameter (ndim=4, pdim=1)
         real*8 :: xgrid(0:50,ndim), xint, ymax(50,ndim), intabs_val, intabs_err, estimn, errorn, estimp, errorp
         real*8 :: xgenerated(pdim), xtransformed(pdim), xmin, xmax
         common/bounds/xmin, xmax

         integer ifold(ndim)
         common/cifold/ifold
         logical negflag ! If true, the function is returning only non zero if it is negative and vice versa.
         common/cnegflag/negflag

         real*8 :: x12(2)
         common/xset/x12

         character(len=30) :: map_type
         common /maptypeblock/ map_type

         character(len=30) :: file_events
         integer :: unit_events


         logical :: flg_integration, flg_2dintegration, flg_generation

         ! Timer variables
         real*8 :: t_start, t_end

         flg_integration   = .true.
         flg_2dintegration = .false.
         flg_generation    = .false.


         ! Number of points to improve the grid:
         ncall1 = 1000000
         ! Number of grid improvement iterations (MaXimum number of ITerations to improve the grid):
         itmx1  = 5

         ! Set up the grid
         call cpu_time(t_start)
         call mint(func_wrap, ndim, ncall1, itmx1, 0, xgrid, xint, ymax, intabs_val, intabs_err)
         call cpu_time(t_end)
         write(*,*) 'Integral over the absolute value of the function:'
         write(*,*) 'Int[ |f| ]: ', intabs_val, ' +- ', intabs_err
         ! Note: xint == intabs_val, since the function is integrated over the absolute value and xint is used as
         !       the initial value for computing the upper bound of the function in the next, imode = 1, step,
         !       where xint will be an input, not an output.
         write(*,*) 'Grid setup time (s): ', t_end - t_start
         ! write(*,*) xgrid(:,1)

         if (flg_integration) then
            ! Setup folding:
            ! Notice that you can only fold by divisors of 50, since this is the number of points in the grid.
            ! Moreover, folding one dimension by e.g. 5 means the integration will take 5 times longer, since it will
            ! call the integrated function 5 times more often in that dimension. Folding by 5 in 2 dimensions will
            ! take 25 times longer, etc..
            ifold(1)=1
            ifold(2)=1
            ifold(3)=1
            ifold(4)=1

            ! Number of points used for the integration:
            ncall2=100000
            ! Number of integration iterations and upper bound improvements, all done with a number of calls ncall2 to the
            ! integrated function. The different integrand results are combined and only the final result is returned.
            ! Note that this basically corresponds to increasing the number of calls ncall2 by a factor of itmx2.
            itmx2=5

            call cpu_time(t_start)

            ! Compute the positive contribution to the integral:
            negflag = .false.
            call mint(func_wrap, ndim, ncall2, itmx2, 1, xgrid, xint, ymax, estimp, errorp)
            ! Compute the negative contribution to the integral:
            negflag = .true.
            call mint(func_wrap, ndim, ncall2, itmx2, 1, xgrid, xint, ymax, estimn, errorn)
            negflag = .false.

            call cpu_time(t_end)
            write(*,*) estimp,' +- ', errorp
            write(*,*) estimn,' +- ', errorn

            write(*,*) (estimp + estimn),' +- ', sqrt(errorp**2 + errorn**2)
            write(*,*) 'Integration time (s): ', t_end - t_start
         end if





         if (flg_generation) then
            ! Generate random numbers according to the grid.
            ! The function gen is used to generate random numbers according to the function specified in func_wrap.
            ! Inputs:
            ! 1. func_wrap - contains the function according to which the random numbers are generated.
            ! 2. pdim      - dimension of the random numbers to be generated (1 in this case). TODO: Specify which dimension is used.
            ! 3. xgrid     - the integration grid.
            ! 4. ymax      - the maximum value of the function on the grid
            ! 5. ifl - if 0, it generates random numbers, if 1, it generates events, if 2, it generates events with weights, if 3, it generates events with weights and writes them to a file.
            ! Output:
            ! 1. xgenerated - generated random number of dimension pdim.

            write(*,*) "Initialise generation of random numbers:"
            call gen(func_wrap, pdim, xgrid, ymax, 0, xgenerated)

            file_events = 'events.lhe'
            ! open(newunit=unit_events, file=trim(file_events), status='unknown')

            nevents = 1
            write(*,*) "Generating ", nevents, " events and writing them to file: ", trim(file_events)
            call gen(func_wrap, 1, xgrid, ymax, 1, xgenerated)
            write(*,*) xgenerated
            ! do i = 1, nevents
            !    call gen(func_wrap, pdim, xgrid, ymax, 1, xgenerated)
            !    ! do j=1, pdim
            !    !    xtransformed(j) = mapping(xgenerated(j), map_type)
            !    ! end do
            !    xtransformed = xmin + (xmax - xmin) * xgenerated
            !    ! write(unit_events,*) xtransformed
            !    write(*,*) xtransformed
            !    ! if (mod(i, 1000) == 0) then
            !    !    write(*,*) 'Events written: ', i
            !    !    close(unit_events)
            !    !    call sleep(1)
            !    !    open(newunit=unit_events, file=trim(file_events), status='unknown', position='append')
            !    ! end if
            ! end do
            ! ! close(unit_events)

            write(*,*) "Finished generating events."
            call gen(func_wrap, pdim, xgrid, ymax, 3, xgenerated)
         end if


         if (flg_2dintegration) then
            ! Since we are able to integrate all 4 dimensions, what about setting the first two to a specific
            ! value and integrating over the remaining 2?
            write(*,*) "Consider 2D now:"
            ! Set value of first two x:
            x12 = [0.45d0, 0.55d0]
            ! Note that we remapped all numbers with
            ! xmin + (xmax - xmin) * x(i),
            ! so that even if x(i) lies between 0 and 1, we can still set arbitrary integral boundaries
            ! - they become very fast very unstable though.
            ! This mapping should be inverted here if we want to set the first two values to a specific value:
            x12 = [inverseMapping(x12(1), map_type), inverseMapping(x12(2), map_type)]


            ! The integration grid can be kept! It should be optimal also for the remaining two dimensions.

   !          ncall1=1000000
   !          itmx1=5
   !          ! set up the grid
   !          call mint(func2d,2,ncall1,itmx1,0,xgrid,xint,ymax,sigtot,error)
   !          write(*,*) sigtot, error

            ! Only redo the integration for the remaining two dimensions.
            ifold(1)=1
            ifold(2)=1
            ifold(3)=1
            ifold(4)=1
            ncall2=1000000
            itmx2=5

            call cpu_time(t_start)
            negflag=.true.
            call mint(func2d, 2, ncall2, itmx2, 1, xgrid, xint, ymax, estimn, errorn)
            negflag=.false.
            call mint(func2d, 2, ncall2, itmx2, 1, xgrid, xint, ymax, estimp, errorp)
            call cpu_time(t_end)
            write(*,*) estimp,' +- ', errorp
            write(*,*) estimn,' +- ', errorn

            write(*,*) (estimp+estimn),' +- ', sqrt(errorp**2+errorn**2)
            write(*,*) '2D Integration time (s): ', t_end - t_start
         end if

         contains
            function func2d(x,weight,ifl)
               ! This function is a 2D version of the gauss function in the style needed for the mint integration routine.
               ! x12 are given in a common block and are the first two values of x.
               ! I guess they need to be given like that, because the mint routine is not able to pass additional arguments to the function?
               implicit none
               real*8, intent(in) :: x(2), weight
               integer, intent(in) :: ifl
               real*8 :: func2d, jac
               integer :: i

               real*8 :: x12(2)
               common/xset/x12

               character(len=30) :: map_type
               common /maptypeblock/ map_type

               func2d = func_wrap([x12(1), x12(2), x(1), x(2)], weight, ifl)
               ! Since the gauss function is 4 dimensional, also the jacobian is multiplied 4 times.
               ! This is not correct, because we are only integrating over 2 dimensions.
               jac = 1d0
               do i=1, size(x12)
                  jac = jac * jacobian(x(i), map_type)
               end do
               func2d = func2d / jac
            end function func2d

      end program integrate_functions
