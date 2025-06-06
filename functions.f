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

            function inverseMapping(xmapped)
               implicit none
               real*8, intent(in) :: xmapped
               real*8 :: xmin, xmax
               common/bounds/xmin, xmax
               real*8 :: inverseMapping

               inverseMapping = (xmapped - xmin) / (xmax - xmin)
            end function inverseMapping
      end module functions

