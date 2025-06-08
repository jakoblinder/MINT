      module functions
         implicit none
         public

         contains
            function gauss(x, weight, ifl)
               implicit none
               real*8, intent(in) :: x(4), weight
               integer, intent(in) :: ifl
               real*8 :: gauss

               logical negflag
               common/cnegflag/negflag

               real*8 :: accum
               save accum

               real*8, parameter :: pi=3.141592653589793d0
               real*8 :: xmapped(4), xmin, xmax
               common/bounds/xmin, xmax

               real*8 :: mean1, stddev1, mean2, stddev2
               integer :: i, ndim
               character(len=30) :: map_type

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
               ! xmapped = xmin + (xmax - xmin) * x
               ! jac     = (xmax - xmin)

               map_type = "linear" ! Change this to "linear", "exponential" or "logarithmic" for different mappings.
               do i=1, ndim
                  xmapped(i) = mapping(x(i), map_type)
               end do

               ! gauss =  ((1.0d0 / (stddev1 * sqrt(2.0d0 * pi))))**ndim * exp(-0.5d0 * sum(((xmapped - mean1) / stddev1)**2))
               gauss =  exp(+ sum(xmapped))

               do i=1, ndim
                  gauss = gauss * jacobian(x(i), map_type)
               end do

               gauss = gauss * weight

               accum = accum + gauss
            end function gauss

            function mapping(x, map_type)
               implicit none
               real*8, intent(in) :: x
               character(len=*), intent(in) :: map_type
               real*8 :: xmin, xmax
               common/bounds/xmin, xmax
               real*8 :: mapping

               select case (trim(map_type))
               case ("linear")
                  ! Linear map, mapping x going from 0 to 1, to xmin to xmax.
                  mapping = xmin + (xmax - xmin) * x
               case ("exponential")
                  ! Exponential map, mapping x going from 0 to 1, to xmin to xmax.
                  mapping = exp( log(xmin) + (log(xmax) - log(xmin)) * x )
               case ("logarithmic")
                  ! Logarithmic map, mapping x going from 0 to 1, to xmin to xmax.
                  mapping = log( exp(xmin) + (exp(xmax) - exp(xmin)) * x )
               case default
                  ! Defaults to linear if unknown type.
                  mapping = xmin + (xmax - xmin) * x
               end select
            end function mapping

            function jacobian(x, map_type)
               implicit none
               real*8, intent(in) :: x
               character(len=*), intent(in) :: map_type
               real*8 :: xmin, xmax
               common/bounds/xmin, xmax
               real*8 :: jacobian

               select case (trim(map_type))
               case ("linear")
                  ! Jacobian for linear mapping.
                  jacobian = abs(xmax - xmin)
               case ("exponential")
                  ! Jacobian for exponential mapping.
                  ! jacobian = abs(log(xmax) - log(xmin) * mapping(x, map_type))
                  jacobian = abs(log(xmax) - log(xmin) * exp( log(xmin) + (log(xmax) - log(xmin)) * x ))
               case ("logarithmic")
                  ! Jacobian for logarithmical mapping.
                  jacobian = abs( (exp(xmax) - exp(xmin)) / (exp(xmin) + (exp(xmax) - exp(xmin)) * x))
               case default
                  ! Defaults to linear if unknown type.
                  jacobian = abs(xmax - xmin)
               end select
            end function jacobian

            function inverseMapping(xmapped, map_type)
               implicit none
               real*8, intent(in) :: xmapped
               character(len=*), intent(in) :: map_type
               real*8 :: xmin, xmax
               common/bounds/xmin, xmax
               real*8 :: inverseMapping

               select case (trim(map_type))
               case ("linear")
                  ! Inverse linear map, mapping x going from 0 to 1, to xmin to xmax.
                  inverseMapping = (xmapped - xmin) / (xmax - xmin)
               case ("exponential")
                  ! Inverse exponential map, mapping x going from 0 to 1, to xmin to xmax.
                  inverseMapping = (log(xmapped) - log(xmin)) / (log(xmax) - log(xmin))
               case ("logarithmic")
                  ! Inverse logarithmic map, mapping x going from 0 to 1, to xmin to xmax.
                  inverseMapping = (exp(xmapped) - exp(xmin)) / (exp(xmax) - exp(xmin))
               case default
                  ! Defaults to the inverse linear if unknown type.
                  inverseMapping = (xmapped - xmin) / (xmax - xmin)
               end select
            end function inverseMapping
      end module functions

