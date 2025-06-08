      program integrate_gauss
         use functions
         implicit none
         integer :: i, pdim, ndim, ncall1, itmx1, ncall2, itmx2
         parameter (ndim=4, pdim=1)
         real * 8 :: xgrid(0:50,ndim), xint, ymax(50,ndim), sigtot, error, estimn, errorn, estimp, errorp
         real*8 :: xgenerated(pdim), xtransformed(pdim), xmin, xmax
         common/bounds/xmin, xmax

         integer ifold(ndim)
         common/cifold/ifold
         logical negflag ! If true, the function is returning only non zero if it is negative and vice versa.
         common/cnegflag/negflag

         real*8 :: x12(2)
         common/xset/x12



         abstract interface
            function func (x,weight,ifl)
               implicit none
               real*8 :: func
               integer, parameter :: ndim=4
               real*8, intent(in) :: x(ndim), weight
               integer, intent(in) :: ifl
            end function func
         end interface

         procedure (func), pointer :: f_ptr => null ()

         xmin = -1d10
         xmax =  0d0

         ! Decide which function to integrate:
         f_ptr => gauss

         ! xx = [0.5d0, 0.3d0, 0.9d0, 0.5d0]

         ! print *, f_ptr(xx, 1d0, 1)

         ! Number of points to improve the grid:
         ncall1 = 1000000
         ! #grid iterations:
         itmx1  = 5

         ! Set up the grid
         call mint(f_ptr, ndim, ncall1, itmx1, 0, xgrid, xint, ymax, sigtot, error)
         write(*,*) sigtot, error
         ! write(*,*) xgrid(:,1)

         ! Setup folding:
         ! Notice that with no folding (all ifold=1)
         ! we get around 2 per mill negative contribution to
         ! the integral, for  fixed alphas=0.5; with the
         ! choice below it disappears.
         ifold(1)=1
         ifold(2)=1
         ifold(3)=1
         ifold(4)=1
         ! # of points used for the integration:
         ncall2=100000
         ! TODO: What was this again?
         itmx2=5

         negflag = .true.
         call mint(f_ptr, ndim, ncall2, itmx2, 1, xgrid, xint, ymax, estimn, errorn)
         negflag = .false.
         call mint(f_ptr, ndim, ncall2, itmx2, 1, xgrid, xint, ymax, estimp, errorp)
         write(*,*) estimp,' +- ', errorp
         write(*,*) estimn,' +- ', errorn

         write(*,*) (estimp+estimn),' +- ', sqrt(errorp**2+errorn**2)
         negflag = .false.





         ! call gen(f_ptr, pdim, xgrid, ymax, 0, xgenerated)

         ! open(unit=10, file='events.lhe', status='unknown')
         ! do i = 1, 100000
         !    call gen(f_ptr, pdim, xgrid, ymax, 1, xgenerated)
         !    xtransformed = xmin + (xmax - xmin) * xgenerated
         !    write(10,*) xtransformed
         !    if (mod(i, 1000) == 0) then
         !       write(*,*) 'Events written: ', i
         !       close(10)
         !       open(unit=10, file='events.lhe', status='unknown', position='append')
         !    end if
         ! end do
         ! close(10)

         ! call gen(f_ptr, pdim, xgrid, ymax, 3, xgenerated)


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
         x12 = [inverseMapping(x12(1), "linear"), inverseMapping(x12(2), "linear")]


         ! The integration grid can be kept! It should be optimal also for the remaining two dimensions.

!          ncall1=1000000
!          itmx1=5
!          ! set up the grid
!          call mint(gauss2,2,ncall1,itmx1,0,xgrid,xint,ymax,sigtot,error)
!          write(*,*) sigtot, error

         ! Only redo the integration for the remaining two dimensions.
         ifold(1)=1
         ifold(2)=1
         ifold(3)=1
         ifold(4)=1
         ncall2=1000000
         itmx2=5

         negflag=.true.
         call mint(gauss2, 2, ncall2, itmx2, 1, xgrid, xint, ymax, estimn, errorn)
         negflag=.false.
         call mint(gauss2, 2, ncall2, itmx2, 1, xgrid, xint, ymax, estimp, errorp)
         write(*,*) estimp,' +- ', errorp
         write(*,*) estimn,' +- ', errorn

         write(*,*) (estimp+estimn),' +- ', sqrt(errorp**2+errorn**2)

         contains
            function gauss2(x,weight,ifl)
               ! This function is a 2D version of the gauss function in the style needed for the mint integration routine.
               ! x12 are given in a common block and are the first two values of x.
               ! I guess they need to be given like that, because the mint routine is not able to pass additional arguments to the function?
               implicit none
               real*8, intent(in) :: x(2), weight
               integer, intent(in) :: ifl
               real*8 :: gauss2

               logical negflag
               common/cnegflag/negflag

               real*8 :: xmin, xmax, jac
               common/bounds/xmin, xmax

               real*8 :: x12(2)
               common/xset/x12

               gauss2 = gauss([x12(1), x12(2), x(1), x(2)], weight, ifl)
               ! Since the gauss function is 4 dimensional, also the jacobian is multiplied 4 times.
               ! This is not correct, because we are only integrating over 2 dimensions.
               jac = (xmax - xmin)
               gauss2 = gauss2 / jac**size(x12)
            end function gauss2

      end program integrate_gauss
