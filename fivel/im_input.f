c+
c NAME:
c     IM_INPUT()
c
c PURPOSE:
c     Prompt the user for an ion number and one or all of
c     the temperature, the density, or the line ratio.
c
c CALLING SEQUENCE:
c     call im_input(ionum,den,tem,robs)
c
c INPUTS:
c     ionum - integer ion number
c
c OUTPUTS:
c     den   - electron density [cm-3]
c     tem   - electron temperature [K]
c     robs  - observed emission line ratio
c
c INTERNAL VARIABLES:
c     phys  - logical to determine to prompt for ROBS 
c
c COMMON BLOCKS:
c     None.
c
c PROCEDURES USED:
c     None.
c
c INPUT FILES:
c     tty.dec
c
c COMMENTS:
c     Major changes from the original program: (1) Common blocks
c     removed.  Double precision adopted everywhere.
c
c MODIFICATION HISTORY:
c     J. Moustakas, 2004 July 4, U of A - adopted from de Robertis,
c        Dufour, & Hunt's FIVEL.F
c-

      subroutine im_input(ionum,den,tem,robs)

      implicit none
      
      integer*4 ionum           ! input
      logical phys              ! internal
      real*8 den, tem, robs     ! output

      include 'tty.dec'

      phys = .true.

c     Enter relevant physical parameters as required

      if (mod(ionum,10) .eq. 0) then
 100     write(ttyout,10)
         read(kbin,*,err=100,end=100) den,tem
         write(filout,50) int(den),int(tem)
         if(tem .lt. 1.e3 .or. tem .gt. 4.0e4) then
            write(ttyout,*)"Indicated temperature too low/high."
            stop 
         endif
      else if(mod(ionum,10) .eq. 2) then
 110     write(ttyout,20)
         read(kbin,*,err=110,end=110) den
      else
 120     write(ttyout,30) 
         read(kbin,*,err=120,end=120) tem
         if(tem .lt. 1.e3 .or. tem .gt. 3.6e4) then
            write(ttyout,*)"Indicated temperature too low/high."
            stop
         endif
      endif
c     
c     For new ion, line ratio must be included in IF block below
c     
 130  if(ionum .eq. 6021) then
         write(ttyout,*) 'C III] N calculation'
         write(ttyout,*) 'I(1909)/I(1907) = '
      else if(ionum .eq. 7001) then
         write(ttyout,*) '[N I] N calculation'
         write(ttyout,*) 'I(5200)/I(5198) = '
      else if(ionum .eq. 7012) then
         write(ttyout,*) '[N II] T calculation'
         write(ttyout,*) 'I(6548+6583)/I(5755) = '
      else if(ionum .eq. 8011) then
         write(ttyout,*) '[O II] N calculation'
         write(ttyout,*) 'I(3726)/I(3729) = '
      else if(ionum .eq. 8012) then
         write(ttyout,*) '[O II] T calculation'
         write(ttyout,*) 'I(3727)/I(7324) = '
      else if(ionum .eq. 8022) then
         write(ttyout,*) '[O III] T calculation'
         write(ttyout,*) 'I(4959+5007)/I(4363) = '
      else if(ionum .eq. 10022) then
         write(ttyout,*) '[Ne III] T calculation'
         write(ttyout,*) 'I(3869+3968)/I(3343) = '
      else if(ionum .eq. 10031) then
         write(ttyout,*) '[Ne IV] N calculation'
         write(ttyout,*) 'I(2423)/I(2425) = '
      else if(ionum .eq. 10042) then
         write(ttyout,*) '[Ne V] T calculation'
         write(ttyout,*) 'I(3426+3346)/I(2975) = '
      else if(ionum .eq. 14011) then
         write(ttyout,*) '[Si III] N calculation'
         write(ttyout,*) 'I(1883)/I(1893) = '
      else if(ionum .eq. 16011) then
         write(ttyout,*) '[S II] N calculation'
         write(ttyout,*) 'I(6717)/I(6731) = '
      else if(ionum .eq. 16012) then
         write(ttyout,*) '[S II] T calculation'
         write(ttyout,*) 'I(6717+6731)/I(4069+4076) = '
      else if(ionum .eq. 16022) then
         write(ttyout,*) '[S III] T calculation'
         write(ttyout,*) 'I(9069+9532)/I(6312) = '
      else if(ionum .eq. 17021) then
         write(ttyout,*) '[Cl III] N calculation'
         write(ttyout,*) 'I(5518)/I(5538) = '
      else if(ionum .eq. 17032) then
         write(ttyout,*) '[Cl IV] T calculation'
         write(ttyout,*) 'I(7530+8045)/I(5323) = '
      else if(ionum .eq. 18022) then
         write(ttyout,*) '[Ar III] T calculation'
         write(ttyout,*) 'I(7136+7751)/I(5192) = '
      else if(ionum .eq. 18031) then
         write(ttyout,*) '[Ar IV] N calculation'
         write(ttyout,*) 'I(4711)/I(4740) = '
      else if(ionum .eq. 18032) then
         write(ttyout,*) '[Ar IV] T calculation'
         write(ttyout,*) 'I(4711+4740)/I(2854+2868) = '
      else if(ionum .eq. 18042) then
         write(ttyout,*) '[Ar V] T calculation'
         write(ttyout,*) 'I(6435+7006)/I(4626) = '
      else
         phys = .false.
      endif
      if(phys) then
         read(kbin,*,err=130,end=130) robs
         write(filout,40) robs
      endif
 10   format(' Enter density and temperature to find ratio: ')
 20   format(' Enter fixed density to find temperature: ')
 30   format(' Enter fixed temperature to find density: ')
 40   format(' Observed line ratio = ',f8.2)
 50   format(' Electron density = ', i7,' Electron temperature = ', 
     .     i6)

      return
      end
