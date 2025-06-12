      module functions
         implicit none
         public

         abstract interface
            function func (x)
               implicit none
               real*8 :: func
               real*8, intent(in) :: x
            end function func
         end interface

         contains
            function func_wrap(x, weight, ifl)
               implicit none
               real*8, intent(in) :: x(4), weight
               integer, intent(in) :: ifl
               real*8 :: func_wrap, func_tmp

               procedure (func), pointer :: f_ptr1 => null(), f_ptr2 => null()

               logical negflag
               common/cnegflag/negflag

               real*8 :: accum
               save accum

               real*8, parameter :: pi=3.141592653589793d0
               real*8 :: xmapped(4), xmin, xmax
               common/bounds/xmin, xmax

               integer :: i, ndim
               character(len=30) :: map_type
               common /maptypeblock/ map_type
               logical :: flg_sub2function
               flg_sub2function = .false.

               ndim = size(x)

               ! Set the function pointer to the desired function:
               ! f_ptr1 => gauss1d
               ! xmin = -1d1
               ! xmax =  1d1
               ! [ -Inf, +Inf] => 1

               f_ptr1 => gauss1d
               f_ptr2 => gauss1d2
               flg_sub2function = .true.
               xmin = -5d0
               xmax =  5d0
               ! [ -Inf, +Inf] => 0

               ! f_ptr1 => exponential
               ! xmin = -1d2
               ! xmax =  0d0
               ! [ -Inf, 0] => 1

               ! f_ptr1 => sinus
               ! xmin =  0d0
               ! xmax =  pi
               ! [  0, pi]   => 16
               ! [-pi, pi] =>  0


               if(ifl.eq.2) then
                  func_wrap=accum
                  accum=0

                  if(negflag) then
                     if(func_wrap.gt.0) func_wrap=0
                  else
                     if(func_wrap.lt.0) func_wrap=0
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

               ! Return the 1D functin value at the mapped point exponentiated by the required dimension
               ! and multiply by the jacobian.
               func_wrap = 1d0
               do i=1, ndim
                  func_wrap = func_wrap * f_ptr1(xmapped(i)) * jacobian(x(i), map_type)
               end do

               ! If the second function is used, subtract it from the first one.
               if (flg_sub2function) then
                  func_tmp = 1d0
                  do i=1, ndim
                     func_tmp  = func_tmp  * f_ptr2(xmapped(i)) * jacobian(x(i), map_type)
                  end do
                  func_wrap = func_wrap - func_tmp
               end if

               ! Multipy the function value by the weight.
               func_wrap = func_wrap * weight

               ! Accumulate the function value. Needed for the folding functionality of MINT.
               accum = accum + func_wrap
            end function func_wrap

            function gauss1d(x) result(res)
               ! Gaussian 1D function.
               implicit none
               real*8, intent(in) :: x
               real*8 :: res

               real*8, parameter :: pi=3.141592653589793d0
               real*8 :: mean1, stddev1

               mean1   = -1.5d0
               stddev1 = 0.5d0

               res = ((1.0d0 / (stddev1 * sqrt(2.0d0 * pi)))) * exp(-0.5d0 * ((x - mean1) / stddev1)**2)
            end function gauss1d

            function gauss1d2(x) result(res)
               ! Gaussian 1D function.
               implicit none
               real*8, intent(in) :: x
               real*8 :: res

               real*8, parameter :: pi=3.141592653589793d0
               real*8 :: mean2, stddev2

               mean2   = 1.5d0
               stddev2 = 0.5d0

               res = ((1.0d0 / (stddev2 * sqrt(2.0d0 * pi)))) * exp(-0.5d0 * ((x - mean2) / stddev2)**2)
            end function gauss1d2

            function exponential(x) result(res)
               ! Gaussian 1D function.
               implicit none
               real*8, intent(in) :: x
               real*8 :: res

               res = exp( +x )
            end function exponential

            function sinus(x) result(res)
               ! Gaussian 1D function.
               implicit none
               real*8, intent(in) :: x
               real*8 :: res

               res = sin( +x )
            end function sinus

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

