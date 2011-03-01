c+
c NAME:
c     IM_TEMDEN()
c
c PURPOSE:
c     Iteratively compute the electron temperature or density.
c
c CALLING SEQUENCE:
c     call im_temden(ionum,den,tem,robs,temden_error)
c
c INPUTS:
c     ionum - integer ion number
c     den   - input/output electron density [cm-3]
c     tem   - input/output electron temperature [K]
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
c     rcal        - "characteristic" ratio calculated in IM_SOLVE()
c     x1, x2      - initial (absolute) limits on TEM or DEN
c     xlo, xhi    - current range on iterated TEM or DEN
c     rlo, rhi    - ratios calculated at XLO, XHI
c     x           - current 'best guess' for TEM or DEN
c
c COMMON BLOCKS:
c     None.
c
c PROCEDURES USED:
c     IM_PARI(), IM_PARII(), IM_SOLVE()
c
c INPUT FILES:
c     tty.dec
c
c COMMENTS:
c
c MODIFICATION HISTORY:
c     J. Moustakas, 2004 July 4, U of A - adopted from de Robertis,
c        Dufour, & Hunt's FIVEL.F
c
c-

      subroutine im_temden(ionum,den,tem,robs,temden_error)

      implicit none
      
      integer*4 ionum             ! input
      real*8 den, tem, robs       ! input/output

      real*8 err, x1, x2, xlo, xhi, x, rcal, rlo, rhi, ! internal variables
     &     a(5,5), t(5), c(5,5), jl(5,5), pop(5), ncrit(5), 
     &     e(5,5), q(5,5), den_initial, tem_initial
      integer*4 iter, weight(5)
      
      integer*4 temden_error ! output

      include 'tty.dec'

      temden_error = 0
      den_initial = den
      tem_initial = tem
      
      err = 0.001

c     Note: setting ERR to too high a value can lead to poor results!

c     Temperature limits

c     Note: watch out for machine precision limitations on calculated
c     emissivities at very low temperatures!!

      if (mod(ionum,10) .eq. 2) then
         x1=1.0e3
         x2=4.0e4
         write(filout,20)
c     
c     Density limits
c     
      elseif (mod(ionum,10) .eq. 1) then
         x1=1.0
         x2=1.0e8
         write(filout,10)
      endif

c     initialize atomic parameters      
      
      call im_pari(ionum,a,t,weight)

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

         call im_parii(ionum,tem,c)
         call im_solve(den,tem,t,weight,a,c,pop,jl,ncrit,e,q)

c     For new ions, enter appropriate ratio in IF block below

         if(ionum .eq. 7001) then
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
c     
c     Output intermediate iterations
c     
         if (mod(ionum,10) .eq. 1) then
            write(filout,30) iter,den,xlo,tem,xhi,(rcal-robs)/robs
         else if(mod(ionum,10) .eq. 2) then
            write(filout,30) iter,tem,xlo,den,xhi,(rcal-robs)/robs
         endif
c     
c     Test for convergence
c     
         if(iter .eq. 1) then
            rlo=rcal
         else if((iter .eq. 2) .and. ((rcal-robs)/(rlo-robs) .lt. 
     .           0.0)) then
            rhi=rcal
         else if((iter .eq. 2) .and. ((rcal-robs)/(rlo-robs) .gt.
     .           0.0)) then
               write(6,*) 'Ratio predicts unreasonable conditions - ',
     &              'restoring initial values.'
            temden_error = 1
            den = den_initial
            tem = tem_initial
            write(filout,40) dmin1(rlo,rcal), dmax1(rlo,rcal)
            write(filout,50) rcal,robs,rlo
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
c     
      write(ttyout,*) 'Did not converge after 100 iterations; stopping'
 10   format(1x,'Iter#',6x,'Density',8x,'XLO',2x,'Temperature',7x,
     1     'XHI',3x,' (RCAL-ROBS)/ROBS')
 20   format(1x,'Iter#',2x,'Temperature',8x,'XLO',6x,'Density',7x,
     1     'XHI',3x,' (RCAL-ROBS)/ROBS')
 30   format(3x,i3,2(3x,1pe10.3,1x,1pe10.3),9x,1pe9.2)
 40   format(' (RCAL-ROBS)/(RLO-ROBS) is > 0, so the entered RATIO',/,
     1     ' predicts unreasonable conditions. Please try again.',/,
     2     ' The allowable RATIO range is:',1pe9.2,'--',1pe9.2)
 50   format(' RCAL = ',f8.3,' ROBS = ',f8.3,' RLO = ',f8.3)
c     
      return
      end
