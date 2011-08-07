c+
c NAME:
c     IDL_TEMDEN()
c
c PURPOSE:
c     Iteratively compute the electron temperature or density.
c
c CALLING SEQUENCE:
c     call idl_temden(ionum,den,tem,robs,temden_error)
c
c INPUTS:
c     ionum - integer ion number
c     den         - input/output electron density [cm-3]
c     tem         - input/output electron temperature [K]
c     robs  - observed emission line ratio
c
c OUTPUTS:
c     temden_error - no errors = 0; unphysical conditions 
c
c INTERNAL VARIABLES:
c     den_initial - initial electron density
c     tem_initial - initial electron temperature
c     err         - percentage tolerance in X
c     iter        - iterative index used to find TEM or DEN
c     rcal        - "characteristic" ratio calculated in IDL_SOLVE()
c     x1, x2      - initial (absolute) limits on TEM or DEN
c     xlo, xhi    - current range on iterated TEM or DEN
c     rlo, rhi    - ratios calculated at XLO, XHI
c     x           - current 'best guess' for TEM or DEN
c
c COMMON BLOCKS:
c     None.
c
c PROCEDURES USED:
c     IDL_PARI(), IDL_PARII(), IDL_SOLVE()
c
c INPUT FILES:
c
c COMMENTS:
c
c MODIFICATION HISTORY:
c     J. Moustakas, 2004 July 4, U of A - adopted from de Robertis,
c        Dufour, & Hunt's FIVEL.F
c
c-

      subroutine idl_temden(ionum,den,tem,robs,temden_error)

      implicit none
      
      integer*4 ionum             ! input
      real*8 den, tem, robs       ! input/output

      real*8 err, x1, x2, xlo, xhi, x, rcal, rlo, rhi, ! internal variables
     &     a(5,5), t(5), c(5,5), jl(5,5), pop(5), ncrit(5), 
     &     e(5,5), q(5,5), den_initial, tem_initial
      integer*4 iter, weight(5)
      
      integer*4 temden_error ! output

      temden_error = 0
      den_initial = den
      tem_initial = tem

      x1 = 0.0
      x2 = 0.0
      xlo = 0.0
      xhi = 0.0
      x = 0.0
      rcal = 0.0
      rlo = 0.0
      rhi = 0.0
      
      err = 1.0e-6

c     Note: setting ERR to too high a value can lead to poor results!

c     Temperature limits

c     Note: watch out for machine precision limitations on calculated
c     emissivities at very low temperatures!!

      if (mod(ionum,10) .eq. 2) then
         x1=1.0e3
         x2=4.0e4

c     Density limits

      elseif (mod(ionum,10) .eq. 1) then
         x1=1.0
         x2=1.0e8
      endif

c     initialize atomic parameters      
      
      call idl_pari(ionum,a,t,weight)

c     Permit up to 100 iterations

      do 100 iter=1,100
         if (iter .eq. 1) then
            xlo=x1
            xhi=x2
            x=xlo
         else if (iter .eq. 2) then
            xhi=x2
            x=xhi
         else
            x=sqrt(xlo*xhi)
         endif
c     
         if (mod(ionum,10) .eq. 2) then
            tem=x
         else if (mod(ionum,10) .eq. 1) then
            den=x
         endif

         call idl_parii(ionum,tem,c)
         call idl_solve(den,tem,t,weight,a,c,pop,jl,ncrit,e,q)

c     For new ions, enter appropriate ratio in IF block below

         if (ionum .eq. 7001) then
            rcal=jl(2,1)/jl(3,1)
         else if(ionum .eq. 7012) then
            rcal=(jl(4,3)+jl(4,2))/jl(5,4)
         else if(ionum .eq. 8011) then
            rcal=jl(3,1)/jl(2,1)
         else if(ionum .eq. 8012) then
            rcal=(jl(2,1)+jl(3,1))/(jl(5,3)+jl(5,2)+jl(4,3)+jl(4,2))
         else if(ionum .eq. 8022) then
            rcal=(jl(4,3)+jl(4,2))/jl(5,4)
         else if(ionum .eq. 10022) then
            rcal=(jl(4,1)+jl(4,2))/jl(5,4)
         else if(ionum .eq. 10031) then
            rcal=jl(3,1)/jl(2,1)
         else if(ionum .eq. 10042) then
            rcal=(jl(4,3)+jl(4,2))/jl(5,4)
         else if(ionum .eq. 16011) then
            rcal=jl(3,1)/jl(2,1)
         else if(ionum .eq. 16012) then
            rcal=(jl(2,1)+jl(3,1))/(jl(4,1)+jl(5,1))
         else if(ionum .eq. 16022) then
            rcal=(jl(4,3)+jl(4,2))/jl(5,4)
         else if(ionum .eq. 17021) then
            rcal=jl(3,1)/jl(2,1)
         else if(ionum .eq. 17032) then
            rcal=(jl(4,2)+jl(4,3))/jl(5,4)
         else if(ionum .eq. 18022) then
            rcal=(jl(4,1)+jl(4,2))/jl(5,4)
         else if(ionum .eq. 18031) then
            rcal=jl(3,1)/jl(2,1)
         else if(ionum .eq. 18032) then
            rcal=(jl(2,1)+jl(3,1))/(jl(4,1)+jl(5,1))
         else if(ionum .eq. 18042) then
            rcal=(jl(4,2)+jl(4,3))/jl(5,4)
         endif

c     Test for convergence

         if(iter .eq. 1) then
            rlo=rcal
         else if((iter .eq. 2) .and. ((rcal-robs)/(rlo-robs) .lt. 
     .           0.0)) then
            rhi=rcal
         else if((iter .eq. 2) .and. ((rcal-robs)/(rlo-robs) .gt.
     .           0.0)) then
c              write(6,*) 'Ratio predicts unreasonable conditions - ',
c    &              'restoring initial values.'
            temden_error = 1
            den = den_initial
            tem = tem_initial
            return
         else if(abs(0.5*(xhi-xlo)/(xhi+xlo)) .lt. err) then
            return
         else if((rcal-robs)/(rhi-robs) .lt. 0.0) then
            rlo=rcal
            xlo=x
         else if((rcal-robs)/(rhi-robs) .gt. 0.0) then
            rhi=rcal
            xhi=x
         endif
 100  continue

      return
      end
