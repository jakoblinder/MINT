! MODIFIED VERSION. Original: MINT integrator of the POWHEG BOX (Paolo Nason et al.), see
! arXiv:0709.2085 and https://virgilio.mib.infn.it/~nason/POWHEG/FNOpaper/mint-integrator.f
! Modifications w.r.t. that version, all by Jakob Linder:
!  2025-06-06  (state when added to this repository) Comments converted to free-form style and extended; typo fixes; x is initialised;
!              removed the debugging output of vtot, etot.
!  2025-06-11  The grid of every iteration is written to its own topdrawer file xg<iteration>.top
!              (subroutines regridplotopen, regridplotclose; title line with the dimension added
!              in regrid), to visualise the grid adaption with plotgrid/plot_topdrawer_grid.py.
!  2026-10-07  The maximum number of dimensions ndimmax is raised from 6 to 20 and is now defined once in
!              the include file src/ndimmax.inc (used in mint and gen), to allow the integration of a
!              d-dimensional sphere with d up to 20 in the example.
!
! Integrator Package for POWHEG
! subroutine mint(fun,ndim,ncalls0,nitmax,imode,xgrid,xint,ymax,ans,err)

! ndim = number of dimensions
!
! ncalls0 = # of calls per iteration
!
! nitmax  = # of iterations
!
! fun(xx,www,ifirst): returns the function to be integrated multiplied by www;
!                     xx(1:ndim) are the variables of integration
!                     ifirst=0: normal behaviour
!
! imode: integer flag
!
! imode=0:
!    When called with imode=0 the routine integrates the absolute value of the function
!    and sets up a grid xgrid(0:50,ndim) such that in each ndim-1 dimensional slice
!    (i.e. xgrid(m-1,n)<xx(n)<xgrid(m,n)) the contribution of the integral is the same
!    the array xgrid is setup at this stage; ans and err are the integral and its error
!
! imode=1 (in fact #0)
!    When called with imode=1, the routine performs the integral of the function fun
!    using the grid xgrid. If some number in the array ifold, (say, ifold(n))
!    is different from 1, it must be a divisor of 50, and the 50 intervals xgrid(0:50,n)
!    are grouped into ifold(n) groups, each group containing 50/ifold(n) nearby
!    intervals. For example, if ifold(1)=5, the 50 intervals for the first dimension
!    are divided in 5 groups of 10. The integral is then performed by folding on top
!    of each other these 5 groups. Suppose, for example, that we choose a random point
!    in xx(1) = xgrid(2,1)+x*(xgrid(3,1)-xgrid(2,1)), in the group of the first 5 interval.
!    we sum the contribution of this point to the contributions of points
!    xgrid(2+m*10,1)+x*(xgrid(3+m*10,1)-xgrid(2+m*10,1)), with m=1,...,4.
!    In the sequence of calls to the
!    function fun, the call for the first point is performed with ifirst=0, and that for
!    all subsequent points with ifirst=1, so that the function can avoid to compute
!    quantities that only depend upon dimensions that have ifold=1, and do not change
!    in each group of folded call. The values returned by fun in a sequence of folded
!    calls with ifirst=0 and ifirst=1 are not used. The function itself must accumulate
!    the values, and must return them when called with ifirst=2.
!
! ifold(ndim): in common block cifold
!     If some number in the array ifold, (say, ifold(n))
!     is different from 1, it must be a divisor of 50, and the 50 intervals xgrid(0:50,n)
!     are grouped into ifold(n) groups, each group containing 50/ifold(n) nearby
!     intervals. For example, if ifold(1)=5, the 50 intervals for the first dimension
!     are divided in 5 groups of 10. The integral is then performed by folding on top
!     of each other these 5 groups. Suppose, for example, that we choose a random point
!     in xx(1) = xgrid(2,1)+x*(xgrid(3,1)-xgrid(2,1)), in the group of the first 5 interval.
!     we sum the contribution of this point to the contributions of points
!     xgrid(2+m*10,1)+x*(xgrid(3+m*10,1)-xgrid(2+m*10,1)), with m=1,...,4.
!     ifirst=0,1,2
!     In the folded sequence of calls to the
!     function fun, the call for the first point is performed with ifirst=0, and that for
!     all subsequent points with ifirst=1, so that the function can avoid to compute
!     quantities that only depend upon dimensions that have ifold=1, and do not change
!     in each group of folded call. The values returned by fun in a sequence of folded
!     calls with ifirst=0 and ifirst=1 are not used. The function itself must accumulate
!     the values, and must return them when called with ifirst=2.
!
! xgrid(0:nintervals,ndim)
!     integration grid; initialized and updated with the call to mint with imode=0
!
! xint: real
!     Output value of the integral when called with imode=0,
!     input value of the integral when called with imode=1 (cannot be zero here!)
!     This is the initial value which is used, after exponentiating it with (1/ndim), to compute the upper bounds, i.e. ymax.
!     Since it is determined in at the imode=0 stage, it corresponds to the integral of the absolute value of the function.
!
! xacc(0:nintervals,ndim):
!     distribution of the accumulated value for each dimension; it is used to compute the optimal grid.
!     So, the sum of the array at fixed ndim is the total accumulated value
!
! nhits(0:nintervals,ndim):
!     distribution of number of hits for each dimension; it is used to compute the optimal grid
!
! ymax(nintervals,ndim):
!     integrand upper bound factors, set up by mint when called with imode=1, to be used
!     by the subroutine gen for the generation of unweighted events
!
! ans: real
!     Output value of the integral (both imode=0 and imode=1)
!
! err: real
!     Output value of the error on the integral
!

      subroutine mint(fun,ndim,ncalls0,nitmax,imode,
     #     xgrid,xint,ymax,ans,err)
! imode=0: integrate and adapt the grid
! imode=1: frozen grid, compute the integral and the upper bounds
! others: same as 1 (for now)
      implicit none
      integer nintervals
      parameter (nintervals=50)
      include '../../src/ndimmax.inc'
      integer ncalls0,ndim,nitmax,imode
      real * 8 fun,xgrid(0:nintervals,ndim),xint,ymax(nintervals,ndim),
     #  ans,err
      real * 8 x(ndimmax),vol
      real * 8 xacc(0:nintervals,ndimmax)
      integer icell(ndimmax),ncell(ndimmax)
      integer ifold(ndimmax),kfold(ndimmax)
      common/cifold/ifold
      integer nhits(1:nintervals,ndimmax)
      real * 8 rand(ndimmax)
      real * 8 dx(ndimmax),f,vtot,etot,prod
      integer kdim,kint,kpoint,nit,ncalls,ibin,iret,nintcurr,ifirst
      real * 8 random
      external random,fun
      character(len=30) filename
      if(imode.eq.0) then
         do kdim=1,ndim
            ifold(kdim)=1
            do kint=0,nintervals
               xgrid(kint,kdim)=dble(kint)/nintervals
            enddo
         enddo
      elseif(imode.eq.1) then
         do kdim=1,ndim
            nintcurr=nintervals/ifold(kdim)
            if(nintcurr*ifold(kdim).ne.nintervals) then
               write(*,*)
     # 'mint: the values in the ifold array should be divisors of',
     #  nintervals
               stop
            endif
            do kint=1,nintcurr
               ymax(kint,kdim)=
     #              xint**(1d0/ndim)
            enddo
         enddo
      endif
      ncalls=ncalls0
      nit=0
      ans=0
      err=0
 10   continue
      nit=nit+1
      if(nit.gt.nitmax) then
         if(imode.eq.0) xint=ans
         return
      endif
      if(imode.eq.0) then
         do kdim=1,ndim
            do kint=0,nintervals
               xacc(kint,kdim)=0
               if(kint.gt.0) then
                  nhits(kint,kdim)=0
               endif
            enddo
         enddo
      endif
      vtot=0
      etot=0
      do kpoint=1,ncalls
! find random x, and its random cell
         do kdim=1,ndim
            kfold(kdim)=1
            ncell(kdim)=nintervals/ifold(kdim)*random()+1
            rand(kdim)=random()
         enddo
         f=0
         ifirst=0
 1       continue
         vol=1
         x=0d0
         do kdim=1,ndim
            nintcurr=nintervals/ifold(kdim)
            icell(kdim)=ncell(kdim)+(kfold(kdim)-1)*nintcurr
            ibin=icell(kdim)
            dx(kdim)=xgrid(icell(kdim),kdim)-xgrid(icell(kdim)-1,kdim)
            vol=vol*dx(kdim)*nintcurr
            x(kdim)=xgrid(icell(kdim)-1,kdim)+rand(kdim)*dx(kdim)
            if(imode.eq.0) nhits(ibin,kdim)=nhits(ibin,kdim)+1
         enddo
! contribution to integral
         if(imode.eq.0) then
            f=abs(fun(x,vol,ifirst))+f
         else
! this accumulated value will not be used
            f=fun(x,vol,ifirst)+f
            ifirst=1
            call nextlexi(ndim,ifold,kfold,iret)
            if(iret.eq.0) goto 1
! closing call: accumulated value with correct sign
            f=fun(x,vol,2)
         endif
!
         if(imode.eq.0) then
! accumulate the function in xacc(icell(kdim),kdim) to adjust the grid later
            do kdim=1,ndim
               xacc(icell(kdim),kdim)=xacc(icell(kdim),kdim)+f
            enddo
         else
! update the upper bounding envelope
            prod=1
            do kdim=1,ndim
               prod=prod*ymax(ncell(kdim),kdim)
            enddo
            prod=(f/prod)
            if(prod.gt.1) then
! This guarantees a 10% increase of the upper bound in this cell
               prod=1+0.1d0/ndim
               do kdim=1,ndim
                  ymax(ncell(kdim),kdim)=ymax(ncell(kdim),kdim)
     #          * prod
               enddo
            endif
         endif
         vtot=vtot+f/ncalls
         etot=etot+f**2/ncalls
      enddo
      if(imode.eq.0) then
! iteration is finished; now rearrange the grid
         write(filename, '(A,I0,A)') 'xg', nit, '.top'
         call regridplotopen(filename)
         do kdim=1,ndim
   !          call regrid(xacc(0,kdim),xgrid(0,kdim),
   !   #           nhits(1,kdim),kdim,nintervals,nit)
              call regrid(xacc(0:nintervals,kdim),xgrid(0:nintervals,kdim),nhits(1:nintervals,kdim),kdim,nintervals,nit)
         enddo
         call regridplotclose()
      endif
! the abs is to avoid tiny negative values
      etot=sqrt(abs(etot-vtot**2)/ncalls)
      ! write(*,*) vtot,etot
      if(nit.eq.1) then
         ans=vtot
         err=etot
      else
! prevent annoying division by zero for nearly zero
! integrands
         if(etot.eq.0.and.err.eq.0) then
            if(ans.eq.vtot) then
               goto 10
            else
               err=abs(vtot-ans)
               etot=abs(vtot-ans)
            endif
         elseif(etot.eq.0) then
            etot=err
         elseif(err.eq.0) then
            err=etot
         endif
         ans=(ans/err+vtot/etot)/(1/err+1/etot)
         err=1/sqrt(1/err**2+1/etot**2)
      endif
      goto 10
      end

      subroutine regridplotopen(filename)
      implicit none
      character *(*) filename
      integer iun
      logical iunopen
      common/cregrid/iun,iunopen
      data iunopen/.false./
      open(newunit=iun,file=trim(filename),status='unknown')
      iunopen=.true.
      end

      subroutine regridplotclose
      implicit none
      integer iun
      logical iunopen
      common/cregrid/iun,iunopen
      close(iun)
      iunopen=.false.
      end


      subroutine regrid(xacc,xgrid,nhits,kdim,nint,nit)
      implicit none
      integer  nint,nhits(nint),kdim,nit,iun
      real * 8 xacc(0:nint),xgrid(0:nint)
      real * 8 xn(100),r
      integer kint,jint
      logical iunopen
      common/cregrid/iun,iunopen
      do kint=1,nint
! xacc (xerr) already contain a factor equal to the interval size
! Thus the integral of rho is performed by summing up
         if(nhits(kint).ne.0) then
            xacc(kint) = xacc(kint-1) + abs(xacc(kint))/ nhits(kint)
         else
            xacc(kint) = xacc(kint-1)
         endif
      enddo
      ! xacc(kint) correspond now to the function I_l defined in the
      ! paper, i.e. the sum of the contributions of the integral over
      ! the absolut value of the function divided by the number of
      ! calls in each bin: I_l = \sum_{j=1}^{n} R_j / N_j.
      ! xacc(kint) = \int_{0}^{xgrid(kint)} |f(x)| dx.
      do kint=1,nint
         xacc(kint)=xacc(kint)/xacc(nint)
      enddo
      ! Deviding each xacc(kint) by the total integral:
      ! xacc(kint) = \int_{0}^{xgrid(kint)} |f(x)| dx / \int_{0}^{1} |f(x)| dx
      ! This is the cumulative distribution function (CDF) of the
      ! absolute value of the function.
      write(iun,*) 'set limits x 0 1 y 0 1'
      write(iun,*) ' title top "dim=',kdim,'"'
      write(iun,*) 0, 0
      do kint=1,nint
         write(iun,*) xgrid(kint),xacc(kint)
      enddo
      write(iun,*) 'join 0'

      do kint=1,nint
         ! r = l / m
         !   = 'index of cell' / 'total number of cells'
         r = dble(kint)/ nint

         write(iun,*) 0, r
         write(iun,*) 1, r
         write(iun,*) ' join'

         do jint=1,nint
            if(r.lt.xacc(jint)) then
               xn(kint)=xgrid(jint-1)+(r-xacc(jint-1))
     #        /(xacc(jint)-xacc(jint-1))*(xgrid(jint)-xgrid(jint-1))
               goto 11
            endif
         enddo
         if(jint.ne.nint+1.and.kint.ne.nint) then
            write(*,*) ' error',jint,nint
            stop
         endif
         xn(nint)=1
 11      continue
      enddo
      do kint=1,nint
         xgrid(kint)=xn(kint)
!         xgrid(kint)=(xn(kint)+2*xgrid(kint))/3
!         xgrid(kint)=(xn(kint)+xgrid(kint)*log(dble(nit)))
!     #        /(log(dble(nit))+1)
         write(iun,*) xgrid(kint), 0
         write(iun,*) xgrid(kint), 1
         write(iun,*) ' join'
      enddo
      write(iun,*) ' newplot'
      end

      subroutine nextlexi(ndim,iii,kkk,iret)
! kkk: array of integers 1 <= kkk(j) <= iii(j), j=1,ndim
! at each call iii is increased lexicographycally.
! for example, starting from ndim=3, kkk=(1,1,1), iii=(2,3,2)
! subsequent calls to nextlexi return
!         kkk(1)      kkk(2)      kkk(3)    iret
! 0 calls   1           1           1       0
! 1         1           1           2       0
! 2         1           2           1       0
! 3         1           2           2       0
! 4         1           3           1       0
! 5         1           3           2       0
! 6         2           1           1       0
! 7         2           1           2       0
! 8         2           2           1       0
! 9         2           2           2       0
! 10        2           3           1       0
! 11        2           3           2       0
! 12        2           3           2       1
      implicit none
      integer ndim,iret,kkk(ndim),iii(ndim)
      integer k
      k=ndim
 1    continue
      if(kkk(k).lt.iii(k)) then
         kkk(k)=kkk(k)+1
         iret=0
         return
      else
         kkk(k)=1
         k=k-1
         if(k.eq.0) then
            iret=1
            return
         endif
         goto 1
      endif
      end


      subroutine gen(fun,ndim,xgrid,ymax,imode,x)
! imode=0 to initialize
! imode=1 to generate
! imode=3 store generation efficiency in x(1)
      implicit none
      integer ndim,imode
      integer nintervals
      parameter (nintervals=50)
      include '../../src/ndimmax.inc'
      real * 8 fun,xgrid(0:nintervals,ndim),
     #         ymax(nintervals,ndim),x(ndim)
      real * 8 dx(ndimmax)
      integer icell(ndimmax),ncell(ndimmax)
      integer ifold(ndimmax),kfold(ndimmax)
      common/cifold/ifold
      real * 8 r,f,ubound,vol,random,xmmm(nintervals,ndimmax)
      real * 8 rand(ndimmax)
      external fun,random
      integer icalls,mcalls,kdim,kint,nintcurr,iret,ifirst
      save icalls,mcalls,xmmm
      if(imode.eq.0) then
         do kdim=1,ndim
            nintcurr=nintervals/ifold(kdim)
            xmmm(1,kdim)=ymax(1,kdim)
            do kint=2,nintcurr
               xmmm(kint,kdim)=xmmm(kint-1,kdim)+
     #              ymax(kint,kdim)
            enddo
            do kint=1,nintcurr
               xmmm(kint,kdim)=xmmm(kint,kdim)/xmmm(nintcurr,kdim)
            enddo
         enddo
         icalls=0
         mcalls=0
         return
      elseif(imode.eq.3) then
         if(icalls.gt.0) then
            x(1)=dble(mcalls)/icalls
         else
            x(1)=-1
         endif
         return
      endif
      mcalls=mcalls+1
 10   continue
      do kdim=1,ndim
         nintcurr=nintervals/ifold(kdim)
         r=random()
         do kint=1,nintcurr
            if(r.lt.xmmm(kint,kdim)) then
               ncell(kdim)=kint
               goto 1
            endif
         enddo
 1       continue
         rand(kdim)=random()
      enddo
      ubound=1
      do kdim=1,ndim
         ubound=ubound*ymax(ncell(kdim),kdim)
      enddo
      do kdim=1,ndim
         kfold(kdim)=1
      enddo
      f=0
      ifirst=0
 5    continue
      vol=1
      do kdim=1,ndim
         nintcurr=nintervals/ifold(kdim)
         icell(kdim)=ncell(kdim)+(kfold(kdim)-1)*nintcurr
         dx(kdim)=xgrid(icell(kdim),kdim)-xgrid(icell(kdim)-1,kdim)
         vol=vol*dx(kdim)*nintervals/ifold(kdim)
         x(kdim)=xgrid(icell(kdim)-1,kdim)+rand(kdim)*dx(kdim)
      enddo
      f=f+fun(x,vol,ifirst)
      ifirst=1
      call nextlexi(ndim,ifold,kfold,iret)
      if(iret.eq.0) goto 5
! get final value (x and vol not used in this call)
      f=fun(x,vol,2)
      if(f.lt.0) then
         write(*,*) 'gen: non positive function'
         stop
      endif
      if(f.gt.ubound) then
         call increasecnt
     #       ('upper bound failure in inclusive cross section')
      endif
      ubound=ubound*random()
      icalls=icalls+1
      if(ubound.gt.f) then
         call increasecnt
     #       ('vetoed calls in inclusive cross section')
         goto 10
      endif
      end
