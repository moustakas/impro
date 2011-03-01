c+
c NAME:
c     IM_OUTPUT()
c
c PURPOSE:
c     Compute emission-line wavelengths, final emissivities, 
c     and the H-beta volume-averaged emissivity.
c
c CALLING SEQUENCE:
c     call im_output(ionum,den,tem,pop,jl,ncrit,t,jhb)
c
c INPUTS:
c     ionum - integer ion number
c     den   - input/output electron density [cm-3]
c     tem   - input/output electron temperature [K]
c     pop   - level population array [5]
c     jl    - emission-line emissivities [5,5]
c     ncrit - critical densities for each transition [5]
c     t     - energy-level separation from ground [1/Angstrom] [5]
c
c OUTPUTS:
c     jhb   - H-beta emissivity
c
c INTERNAL VARIABLES:
c     xl      - wavelength (N.B.- wavelengths > 1e5 
c               Angstroms in microns) [Angstrom]
c     xi      - intensity
c     i, j, k - do loop integers
c
c COMMON BLOCKS:
c     None.
c
c PROCEDURES USED:
c
c INPUT FILES:
c     tty.dec
c
c COMMENTS:
c
c MODIFICATION HISTORY:
c     J. Moustakas, 2004 July 4, U of A - adopted from de Robertis,
c        Dufour, & Hunt's FIVEL.F
c-

      subroutine im_output(ionum,den,tem,pop,jl,ncrit,t,jhb)

      implicit none
      
      integer*4 ionum                                  ! input
      real*8 den, tem, pop(5), jl(5,5), ncrit(5), t(5) ! input
      
      integer*4 i, j, k, dummy           ! internal
      real*8 xl(5), xi(5), t4, e42, jhb2 ! internal

      real*8 jhb                ! output

      include 'tty.dec'

c  Write level populations and critical densities

      write(filout,*) ' '
      write(filout,10) pop(1),(pop(j),ncrit(j), j=2,5)
      write(filout,*) ' '
c
c  Pause statement for screen convenience
c
   50 if (mod(ionum,10) .eq. 0) then
      write(ttyout,*) 'Hit 0 to continue'
      read(kbin,*,err=50,end=50) dummy
c
c  Transform wavelength into microns if > 100000 Angstroms
c  If emissivity is < 1e-38 erg/s, set it to zero.
c
      do 100 i=2,5
         do 200 j=1,4
            if(j .lt. i) then
              xl(j)=1./(t(i)-t(j))
              if(xl(j) .ge. 1.e5) xl(j)=xl(j)*1.e-4
              xi(j)=jl(i,j)/den
              if(xi(j) .lt. 1.e-38) xi(j)=0.
            endif
  200    continue
c
         if(i .eq. 2) then
           write(filout,11) (xl(k),k=1,i-1),(i,k,k=1,i-1),
     1                      (xi(k),k=1,i-1)
         elseif(i .eq. 3) then
           write(filout,12) (xl(k),k=1,i-1),(i,k,k=1,i-1),
     1                      (xi(k),k=1,i-1)
         elseif(i .eq. 4) then
           write(filout,13) (xl(k),k=1,i-1),(i,k,k=1,i-1),
     1                      (xi(k),k=1,i-1)
         elseif(i .eq. 5) then
           write(filout,14) (xl(k),k=1,i-1),(i,k,k=1,i-1),
     1                      (xi(k),k=1,i-1)
         endif
  100 continue

c  Brocklehurst's interpolation formula for H-beta emissivity
c  Accurate for DEN <= 1.e4/cm**3

      t4=tem/1.e4
      e42=1.387/t4**0.983/10.**(0.0424/t4)
      jhb=e42*10.**-25
      jhb2=1.24e-25*(t4**-0.87)
      write(filout,15) jhb
      write(filout,16) jhb2
c
c  Output final results for TEM/DEN ions.
c
      elseif(ionum .eq. 6021) then
        write(filout,*) 'DENS(C++) = ',den
      elseif(ionum .eq. 7001) then
        write(filout,*) 'DENS(N0) = ',den
      elseif(ionum .eq. 7012) then
        write(filout,*) 'TEMP(N+) = ',tem
      elseif(ionum .eq. 8011) then
        write(filout,*) 'DENS(O+) = ',den
      elseif(ionum .eq. 8012) then
        write(filout,*) 'TEMP(O+) = ',tem
      elseif(ionum .eq. 8022) then
        write(filout,*) 'TEMP(O++) = ',tem
      elseif(ionum .eq. 10022) then
        write(filout,*) 'TEMP(Ne++) = ',tem
      elseif(ionum .eq. 10031) then
        write(filout,*) 'DENS(Ne+3) = ',den
      elseif(ionum .eq. 10042) then
        write(filout,*) 'TEMP(Ne+4) = ',tem
      elseif(ionum .eq. 16011) then
        write(filout,*) 'DENS(S+) = ',den
      elseif(ionum .eq. 16012) then
        write(filout,*) 'TEMP(S+) = ',tem
      elseif(ionum .eq. 16022) then
        write(filout,*) 'TEMP(S++) = ',tem
      elseif(ionum .eq. 17021) then
        write(filout,*) 'DENS(Cl++) = ',den
      elseif(ionum .eq. 17032) then
        write(filout,*) 'TEMP(Cl+3) = ',tem
      elseif(ionum .eq. 18022) then
        write(filout,*) 'TEMP(Ar++) = ',tem
      elseif(ionum .eq. 18031) then
        write(filout,*) 'DENS(Ar+3) = ',den
      elseif(ionum .eq. 18042) then
        write(filout,*) 'TEMP(Ar+4) = ',tem
      endif

      write(filout,*) 'Log10 x = ',dlog10(den/sqrt(tem))

   10 format(10x,'Level Populations - Critical Densities (cm**-3):',/,
     110x,'Level 1: ',1pe8.2,/,
     210x,'Level 2: ',1pe8.2,13x,1pe8.2,/,
     310x,'Level 3: ',1pe8.2,13x,1pe8.2,/,
     410x,'Level 4: ',1pe8.2,13x,1pe8.2,/,
     510x,'Level 5: ',1pe8.2,13x,1pe8.2)
   11 format(5x,'Wavelength',5x,f11.2,/,
     1       1x,'Upper--Lower Levels',5x,'(',i1,'--',i1,') ',/,
     2       2x,'Volume Emissivity ',1pe11.3,/)
   12 format(20x,2f11.2,/,20x,2(5x,'(',i1,'--',i1,')'),/,
     1       20x,2(1pe11.3),/)
   13 format(20x,3f11.2,/,20x,3(5x,'(',i1,'--',i1,')'),/,
     1       20x,3(1pe11.3),/)
   14 format(20x,4f11.2,/,20x,4(5x,'(',i1,'--',i1,')'),/,
     1       20x,4(1pe11.3),/)
   15 format(1x,'H-beta volume emissivity (B72 interp.): ',1pe9.2,' 
     1       N(H+)Ne', ' erg/s',/)
   16 format(1x,'H-beta volume emissivity (power law): ',1pe9.2,' 
     1       N(H+)Ne', ' erg/s',/)

      return
      end
