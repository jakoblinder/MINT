      module functions
         implicit none
         public

         contains
            function gauss(x,weight,ifl)
               implicit none
               real*8, intent(in) :: x(4), weight
               integer, intent(in) :: ifl
               real*8 :: gauss

               logical negflag
               common/cnegflag/negflag

               real*8 :: accum
               save accum

               real*8, parameter :: pi=3.141592653589793d0
               real*8 :: xmapped(4), xmin, xmax, jac
               common/bounds/xmin, xmax

               real*8 :: mean1, stddev1, mean2, stddev2
               integer :: ndim

               mean1 = 0.0d0
               stddev1 = 0.5d0

               mean2 = 1.0d0
               stddev2 = 0.2d0

               ndim = size(x)

               if(ifl.eq.2) then
                  gauss=accum
                  accum=0

                  if(negflag) then
                     if(gauss.gt.0) gauss=0
                  else
                     if(gauss.lt.0) gauss=0
                  endif
                  return
               endif

               if(ifl.eq.0) then
                  accum=0
               endif

               ! Map x, which is going from 0 to 1, to xmin to xmax.
               xmapped = xmin + (xmax - xmin) * x
               jac     = (xmax - xmin)

               gauss =  (jac * (1.0d0 / (stddev1 * sqrt(2.0d0 * pi))))**ndim * exp(-0.5d0 * sum(((xmapped - mean1) / stddev1)**2))

               gauss = gauss * weight

               accum = accum + gauss
            end function gauss

            function wrapgauss(x,weight,ifl)
               implicit none
               real*8, intent(in) :: x(4), weight
               integer, intent(in) :: ifl
               real*8 :: wrapgauss

               call subgauss(x,weight,ifl,wrapgauss)
            end function wrapgauss

            subroutine subgauss(x,weight,ifl,gauss)
               implicit none
               real*8, intent(in) :: x(4), weight
               integer, intent(in) :: ifl
               real*8, intent(out) :: gauss

               logical negflag
               common/cnegflag/negflag

               real*8 :: accum
               save accum

               ! real*8 :: exp
               real*8, parameter :: pi=3.141592653589793d0
               real*8 :: xmapped, xmin, xmax, jac
               real*8 :: mean1, stddev1, mean2, stddev2
               integer :: ndim, i

               ! write(*,*) x

               mean1 = 0.5d0
               stddev1 = 0.1d0

               mean2 = 1.0d0
               stddev2 = 0.2d0

               ndim = size(x)

               if(ifl.eq.2) then
                  gauss=accum
                  accum=0

                  if(negflag) then
                     if(gauss.gt.0) gauss=0
                  else
                     if(gauss.lt.0) gauss=0
                  endif
                  return
               endif

               if(ifl.eq.0) then
                  accum=0
               endif

               ! xmin = -1d2
               ! xmax = 1d2

               gauss = 1.0d0
               do i = 1, ndim
                  ! Map x, which is going from 0 to 1, to -infinity to +infinity.
                  xmapped = xmin + (xmax - xmin) * x(i)
                  jac     = (xmax - xmin)
                  ! xmapped = log(x(i))
                  ! jac     = 1d0 / x(i)
                  gauss = gauss * jac
     .                    * ((1.0d0 / (stddev1 * sqrt(2.0d0 * pi))) * exp(-0.5d0 * ((xmapped - mean1) / stddev1)**2))
   !   .                          - (1.0d0 / (stddev2 * sqrt(2.0d0 * pi))) * exp(-0.5d0 * ((xmapped - mean2) / stddev2)**2))
               end do
               gauss = gauss * weight

               accum = accum + gauss
            end subroutine subgauss

            function inverseMapping(xmapped)
               implicit none
               real*8, intent(in) :: xmapped
               real*8 :: xmin, xmax
               common/bounds/xmin, xmax
               real*8 :: inverseMapping

               inverseMapping = (xmapped - xmin) / (xmax - xmin)
            end function inverseMapping
      end module functions

