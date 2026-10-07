! Copyright (C) 2025 Jakob Linder
! SPDX-License-Identifier: GPL-2.0-only
      module functions
         implicit none
         public

         ! Maximum number of dimensions, defined in ndimmax.inc (shared with mint-integrator.f).
         include 'ndimmax.inc'

         ! Set by the main program: integrate the indicator function of a ball (d-dimensional sphere)
         ! of radius sphere_radius instead of the presets in func_wrap.
         logical :: flg_sphere = .false.
         real*8  :: sphere_radius = 1d0

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
               real*8, intent(in) :: x(*), weight
               integer, intent(in) :: ifl
               real*8 :: func_wrap, func_tmp

               procedure (func), pointer :: f_ptr1 => null(), f_ptr2 => null()

               logical negflag
               common/cnegflag/negflag

               real*8 :: accum
               save accum

               real*8, parameter :: pi=3.141592653589793d0
               real*8 :: xmapped(ndimmax), xmin, xmax
               common/bounds/xmin, xmax

               integer :: i, ndim
               common/cndim/ndim ! number of dimensions, set by the main program
               character(len=30) :: map_type
               common /maptypeblock/ map_type
               logical :: flg_sub2function
               flg_sub2function = .false.

               ! Choose the integrand by setting the function pointer(s) and the integration range
               ! [xmin, xmax] (the same for every dimension), and the mapping below. Presets:
               !
               !    f_ptr1 => gauss1d,  f_ptr2 => gauss1d2, flg_sub2function = .true.
               !    xmin = -5, xmax = 5, "linear"
               !       -> product of 1D gaussians at -1.5 minus the same at +1.5 (integral 0,
               !          so positive and negative parts cancel; shows the negflag machinery)
               !    f_ptr1 => gauss1d,  xmin = -10, xmax = 10, "linear"
               !       -> a single gaussian product (integral 1)
               !    f_ptr1 => exponential, xmin = -100, xmax = 0, "linear"
               !       -> integral 1
               !    f_ptr1 => sinus,       xmin = 0, xmax = pi, "linear"
               !       -> integral 2^ndim (with xmin = -pi the integral is 0)
               !
               ! Alternatively the main program can set flg_sphere to integrate the indicator function
               ! of the ndim dimensional ball with radius sphere_radius, see below.
               !
               ! Without the second function set flg_sub2function = .false.. The "exponential" mapping
               ! (see map_type below) needs xmin > 0.

               f_ptr1 => gauss1d
               f_ptr2 => gauss1d2
               flg_sub2function = .true.
               xmin = -5d0
               xmax =  5d0

               if (flg_sphere) then
                  xmin = -sphere_radius
                  xmax =  sphere_radius
               end if

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

               if (flg_sphere) then
                  ! Indicator function of the ball: 1 inside, 0 outside, times the jacobian.
                  ! The exact integral is the volume of the ball, see ball_volume.
                  func_wrap = 0d0
                  if (sum(xmapped(1:ndim)**2) <= sphere_radius**2) func_wrap = 1d0
                  do i=1, ndim
                     func_wrap = func_wrap * jacobian(x(i), map_type)
                  end do
               else
                  ! Return the 1D function value at the mapped point exponentiated by the required dimension
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
               end if

               ! Multiply the function value by the weight.
               func_wrap = func_wrap * weight

               ! Accumulate the function value. Needed for the folding functionality of MINT.
               accum = accum + func_wrap
            end function func_wrap

            function ball_volume(ndim, radius)
               ! Volume of the ndim dimensional ball: pi^(d/2) / Gamma(d/2 + 1) * r^d.
               implicit none
               integer, intent(in) :: ndim
               real*8, intent(in) :: radius
               real*8 :: ball_volume
               real*8, parameter :: pi=3.141592653589793d0

               ball_volume = exp( 0.5d0*ndim*log(pi) - log_gamma(0.5d0*ndim + 1d0) + ndim*log(radius) )
            end function ball_volume

            function gauss1d(x) result(res)
               ! Gaussian 1D function (centred at -1.5).
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
               ! Gaussian 1D function, second one (centred at +1.5).
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
               ! Exponential 1D function.
               implicit none
               real*8, intent(in) :: x
               real*8 :: res

               res = exp( +x )
            end function exponential

            function sinus(x) result(res)
               ! Sine 1D function.
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
                  ! d mapping / d x = (log(xmax) - log(xmin)) * mapping(x, map_type)
                  jacobian = abs(log(xmax) - log(xmin)) * exp( log(xmin) + (log(xmax) - log(xmin)) * x )
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

