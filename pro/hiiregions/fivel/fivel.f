      program fivel
c
c     /////////////////////////////////////////////////////////////////
c     ////     The LICK 5-level atom program - FIVEL               ////
c     ////                                                         ////
c     ////     Language: FORTRAN 77                                ////
c     ////                                                         ////
c     ////     Certain astrophysically common ions have emission-  ////
c     ////     line intensity ratios which are sensitive to either ////
c     ////     temperature (T) or number density (N). From an      ////
c     ////     observed line ratio (ROBS), then, given a tem-      ////
c     ////     perature, the density can be calculated; or, given  ////
c     ////     the density, the temperature can be calculated.     ////
c     ////     Similarly, given the temperature and density,       ////
c     ////     the emission-line emissivities can be calculated.   ////
c     ////                                                         ////
c     ////     Program Structure:                                  ////
c     ////     1. "Physical conditions"                            ////
c     ////         (a) Specify 'T' option to use observed line     ////
c     ////              ratio and density to find temperature      ////
c     ////         (b) Specify 'N' option to use observed line     ////
c     ////             ratio and temperature to find density       ////
c     ////                                                         ////
c     ////     2. "Line emissivities"                              ////
c     ////         Supply temperature and density to determine     ////
c     ////         the level populations, critical densities,      ////
c     ////         and line emissivities                           ////
c     ////                                                         ////
c     ////     MAIN Program:                                       ////
c     ////          Select option 1 or 2 and ion to investigate.   ////
c     ////                                                         ////
c     ////     Subroutines:                                        ////
c     ////        INPUT:  Supplies data for calculation            ////
c     ////        OUTPUT: Displays final results                   ////
c     ////        PARI:   Contains atomic data independent of T    ////
c     ////        PARII:  Contains atomic data dependent on T      ////
c     ////        SOLVE:  Performs matrix algebra. Finds level     ////
c     ////                populations (POP) and line emissivities  ////
c     ////                (JL) as a function of N and T, and       ////
c     ////                critical densities (NCRIT).              ////
c     ////        TEMDEN: Employs an iterative technique (geo-     ////
c     ////                metric bisection) to find 'T' or 'N'.    ////
c     ////                Used only when T or N is unknown.        ////
c     ////                                                         ////
c     ////     Written May 1987/ M. M. De Robertis, R. J. Dufour,  ////
c     ////     and R. W. Hunt.                                     ////
c     ////---------------------------------------------------------////
c     //// ** The default logical units may have to be changed de- ////
c     //// pending on the computer/operating system. Currently,    ////
c     //// unit '*' is used for input and "interactive output",    ////
c     //// while '6' is used for writing out the calculations.     ////
c     ////                                                         ////
c     /////////////////////////////////////////////////////////////////
c
c     _________________________________________________________________
c     |                                                               |
c     |                        MAIN PROGRAM                           |
c     |                                                               |
c     |    IONUM:   integer code for ion of interest                  |
c     |    ICTG:    category (population or physical conditions)      |
c     |    PHYG:    logical for available "physical condition" ions   |
c     |    POPG:    logical for the available population ions         |
c     |    A:       transition probability matrix                     |
c     |    C:       collision strength matrix                         |
c     |    Q:       collision transition probability matrix           |
c     |    JL:      volume emissivity matrix                          |
c     |    JHB:     H-beta volume emissivity                          |
c     |    NCRIT:   critical density of upper level                   |
c     |    POP:     relative level populations                        |
c     |    WEIGHT:  statistical weights of levels                     |
c     |    T:       energy above ground state matrix (1/Angstroms)    |
c     -----------------------------------------------------------------
c
      real a(5,5),c(5,5),q(5,5),jl(5,5),t(5),pop(5)
      integer kbin,ttyout,filout
      logical phyg,popg

      common /tty/ kbin,ttyout,filout
      data kbin/5/,ttyout/6/,filout/1/

      open(filout,file='fiveout.dat')
      rewind filout

      write(ttyout,10)
  100 write(ttyout,15)
      read(kbin,*,end=100,err=100) ictg
      write(ttyout,*) ' '

c
c Physical conditions
c
      if(ictg .eq. 1) then
c
c Ion code: Atomic number//ionization stage//N(=1),T(=2)
c
        write(ttyout,*)'LICK 5-LEVEL ATOM PHYSICAL CONDITIONS CATEGORY'
        write(ttyout,*)'       *** AVAILABLE IONS ***'
        write(ttyout,*)'   6021: C III]   N;                    '
        write(ttyout,*)'   7001: [N I]    N;   7012: [N II]   T;'
        write(ttyout,*)'   8002: [O I]    T;   8011: [O II]   N;'
        write(ttyout,*)'   8012: [O II]   T;   8022: [O III]  T;'
        write(ttyout,*)'  10022: [Ne III] T;  10031: [Ne IV]  N;'
        write(ttyout,*)'  10042: [Ne V]   T;  16011: [S II]   N;'
        write(ttyout,*)'  16012: [S II]   T;  16022: [S III]  T;'
        write(ttyout,*)'  17021: [Cl III] N;  17032: [Cl IV]  T;'
        write(ttyout,*)'  18022: [Ar III] T;  18031: [Ar IV]  N;'
        write(ttyout,*)'  18032: [Ar IV]  T;  18042: [Ar V]   T;'
  110   write(ttyout,20)
        read(kbin,*,end=110,err=110) ionum
c
c For new ions add appropriate parameters in IONUM and PHYG block
c
        phyg=(ionum .eq. 6021).or.(ionum .eq. 7001).or.
     1       (ionum .eq. 7012).or.(ionum .eq. 8011).or.
     2       (ionum .eq. 8012).or.(ionum .eq. 8022).or.
     3       (ionum .eq.10022).or.(ionum .eq.10031).or.
     4       (ionum .eq.10042).or.(ionum .eq.16011).or.
     5       (ionum .eq.16012).or.(ionum .eq.16022).or.
     6       (ionum .eq.17021).or.(ionum .eq.17032).or.
     7       (ionum .eq.18022).or.(ionum .eq.18031).or.
     8       (ionum .eq.18032).or.(ionum .eq.18042)

        if(phyg) then
          write(filout,*)' '
          write(filout,*)'Physical conditions calculations 
     1                    for ion number: ',ionum
          write(filout,*)' '
          call input(ionum,den,tem,robs)
          call pari(ionum)
          call temden(ionum,den,tem,robs,pop,jl)
          if(robs .ne. 0.) call output(ionum,den,tem,pop,jl)
        else
          write(filout,*)'No parameters for that ion; enter another'
        endif

c
c Level populations, line ratios
c
      else if(ictg .eq. 2) then
        popg=.true.
c
c Ion code: Atomic number//ionization stage//0(=emissivities)
c
        write(ttyout,*)'LICK 5-LEVEL ATOM-POPULATIONS CATEGORY'
        write(ttyout,*)'       *** AVAILABLE IONS ***'
        write(ttyout,*)'   6010: [C II];                  '
        write(ttyout,*)'   6020: C III];      7000: [N I];'
        write(ttyout,*)'   7010: [N II];      7020: [N III];'
        write(ttyout,*)'   8000: [O I];       8010: [O II];'
        write(ttyout,*)'   8020: [O III];    10020: [Ne III];'
        write(ttyout,*)'  10030: [Ne IV];    10040: [Ne V];  '
        write(ttyout,*)'  14020; Si III];    16010: [S II]; '
	write(ttyout,*)'  16020: [S III];    16030: [S IV];  '
        write(ttyout,*)'  17010: [Cl II];    17020: [Cl III];'
        write(ttyout,*)'  17030: [Cl IV];    18020: [Ar III];'
        write(ttyout,*)'  18030: [Ar IV];    18040: [Ar V];  '
  120   write(ttyout,20)
        read(kbin,*,err=120,end=120) ionum
c
c For new ions, add appropriate parameters in IONUM and POPG block.
c
        phyg=(ionum .eq. 6010).or.
     .       (ionum .eq. 6020).or.(ionum .eq. 7000).or.
     .       (ionum .eq. 7010).or.(ionum .eq. 7020).or. 
     1       (ionum .eq. 8000).or.(ionum .eq. 8010).or.
     2       (ionum .eq. 8020).or.(ionum .eq.10020).or.
     3       (ionum .eq.10030).or.(ionum .eq.10040).or.
     3       (ionum .eq.14020).or.
     4       (ionum .eq.16010).or.(ionum .eq.16020).or.
     5       (ionum .eq.16030).or.
     6       (ionum .eq.17010).or.(ionum .eq.17020).or.
     7       (ionum .eq.17030).or.(ionum .eq.18020).or.
     8       (ionum .eq.18030).or.(ionum .eq.18040)

        if(popg) then
          write(filout,*) 'Line emissivity calculations',
     1                    'for ion number ',ionum
          write(filout,*) ' '
          call input(ionum,den,tem,robs)
          call pari(ionum)
          call parii(ionum,tem)
          call solve(den,tem,pop,jl)
          call output(ionum,den,tem,pop,jl)

        else
          write(filout,*) 'No parameters for that ion; try another'
        endif
c
c Stop
c
      else if(ictg .eq. 3) then
        stop
      else
c
c Try again
c 
        write(ttyout,*) 'Enter 1, 2, or 3 only'
      endif
      go to 100
  200 continue

   10 format(//,20x,'<<< LICK 5-LEVEL ATOM PROGRAM >>>',//)
   15 format(/,' Calculate PHYSICAL CONDITIONS(1), LINE RATIOS(2),',
     1' or STOP(3): ')
   20 format(' Select one of the above by typing number: ')
      end
      subroutine input(ionum,den,tem,robs)
c
c     _________________________________________________________________
c     |                                                               |
c     |                        INPUT ROUTINE                          |
c     |                                                               |
c     |   ROBS: the observed line ratio                               |
c     |   DEN:  the electron number density of the gas                |
c     |   TEM:  the electron temperature of the gas                   |
c     |   PHYS: logical to decide to read ROBS                        |
c     |                                                               |
c     -----------------------------------------------------------------
c
      integer kbin,ttyout,filout
      logical phys
      common /tty/ kbin,ttyout,filout

      phys = .true.
c
c  Enter relevant physical parameters as required
c
      if(mod(ionum,10) .eq. 0) then
  100   write(ttyout,10)
        read(kbin,*,err=100,end=100) den,tem
        write(filout,50) int(den),int(tem)
        if(tem .lt. 1.e3 .or. tem .gt. 4.0e4) then
          write(ttyout,*)'Indicated temperature outside permitted range'
          stop
        endif
      else if(mod(ionum,10) .eq. 2) then
  110   write(ttyout,20)
        read(kbin,*,err=110,end=110) den
      else
  120   write(ttyout,30) 
        read(kbin,*,err=120,end=120) tem
        if(tem .lt. 1.e3 .or. tem .gt. 3.6e4) then
          write(ttyout,*)'Indicated temperature outside permitted range'
          stop
        endif
      endif
c
c  For new ion, line ratio must be included in IF block below
c
  130 if(ionum .eq. 6021) then
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
   10 format(' Enter density and temperature to find ratio: ')
   20 format(' Enter fixed density to find temperature: ')
   30 format(' Enter fixed temperature to find density: ')
   40 format(' Observed line ratio = ',f8.2)
   50 format(' Electron density = ', i7,' Electron temperature = ', 
     .       i6)
      return
      end
      subroutine output(ionum,den,tem,pop,jl)
c
c     ________________________________________________________________
c     |                                                              |
c     |                           OUTPUT                             |
c     |                                                              |
c     |     POP:      fraction of electrons in each level            |
c     |     JL:       volume emission coefficient                    |
c     |     JHB:      volume emissivity of H-beta (case B)           |
c     |               (using Brocklehurst's interpolation formula    |
c     |     NCRIT:    critical density                               |
c     |     XL,XI:    wavelength (Angstroms), intensity              |
c     |               (N.B.- wavelengths > 1e5 Angstroms in microns) |
c     |                                                              |
c     ----------------------------------------------------------------
c
      real jl(5,5), pop(5), ncrit, xl(5), xi(5), jhb, jhb2
      integer weight(5)
      integer kbin,ttyout,filout
      common/atom/c(5,5),a(5,5),ncrit(5),t(5),weight
      common/tty/kbin,ttyout,filout
c
c  Write level populations and critical densities
c
      write(filout,*) ' '
      write(filout,10) pop(1),(pop(j),ncrit(j), j=2,5)
      write(filout,*) ' '
c
c  Pause statement for screen convenience
c
   50 if(mod(ionum,10) .eq. 0) then
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
c
c  Brocklehurst's interpolation formula for H-beta emissivity
c  Accurate for DEN <= 1.e4/cm**3
c
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
c
      write(filout,*) 'Log10 x = ',alog10(den/sqrt(tem))
c
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
      subroutine pari(ionum)
c
c     _________________________________________________________________
c     |                                                               |     
c     |                           PAR I                               |
c     |                                                               |
c     |  T:      energy-level separation from ground in 1/Angstroms   |
c     |  WEIGHT: statistical weights for each level                   |
c     |  A:      radiative transition probabilities                   |
c     |                                                               |
c     -----------------------------------------------------------------

      real ncrit
      integer weight(5)
      integer kbin,ttyout,filout
      common /atom/ c(5,5),a(5,5),ncrit(5),t(5),weight
      common /tty/ kbin,ttyout,filout

c
c  Zero arrays
c
      do 100 i=1,5
         do 200 j=1,5
            a(i,j)=0.0
            c(i,j)=0.0
  200    continue
  100 continue
c
      t(1)=0.0
      if((ionum .eq. 6010) .or. (ionum .eq. 6011)
     1                     .or. (ionum .eq. 6012)) then
c
c  C II ion parameters
c
        a(2,1)=2.29e-6
        a(3,1)=55.3
        a(4,1)=1.71
        a(5,1)=0.0
        a(3,2)=65.5
        a(4,2)=5.24
        a(5,2)=43.2
        a(4,3)=2.39e-7
        a(5,3)=3.49e-14
        t(2)=6.342e-7
        t(3)=4.30033e-4
        t(4)=4.30253e-4
        t(5)=4.30536e-4
        weight(1)=2
        weight(2)=4
        weight(3)=2
        weight(4)=4
        weight(5)=6
      elseif((ionum .eq. 6020).or. (ionum .eq. 6021)
     1                     .or. (ionum .eq. 6022)) then
c
c  C III ion parameters
c
        a(2,1)=0.0
        a(3,1)=95.9
        a(4,1)=5.19e-3
        a(5,1)=1.79e9
        a(3,2)=2.39e-7
        a(4,2)=0.0
        a(5,2)=0.0
        a(4,3)=2.41e-6
        a(5,3)=0.0
        a(5,4)=0.0
        t(2)=5.23671e-4
        t(3)=5.23908e-4
        t(4)=5.24471e-4
        t(5)=1.02352e-3
        weight(1)=1
        weight(2)=1
        weight(3)=3
        weight(4)=5
        weight(5)=3
c
      elseif((ionum .eq. 7000) .or. (ionum .eq. 7001)
     1                     .or. (ionum .eq. 7002)) then
c
c  N I ion parameters
c
        a(2,1)=7.27e-6
        a(3,1)=2.02e-5
        a(4,1)=2.71e-3
        a(5,1)=6.58e-3
        a(3,2)=1.27e-8
        a(4,2)=3.45e-2
        a(5,2)=6.14e-2
        a(4,3)=5.29e-2
        a(5,3)=2.76e-2
        t(2)=1.92245e-4
        t(3)=1.92332e-4
        t(4)=2.88389e-4
        t(5)=2.88393e-4
        weight(1)=4
        weight(2)=6
        weight(3)=4
        weight(4)=2
        weight(5)=4
      else if((ionum .eq. 7010) .or. (ionum .eq. 7011)
     1                          .or. (ionum .eq. 7012)) then
c
c  N II ion parameters
c
        a(2,1)=2.08e-6
        a(3,1)=1.16e-12
        a(4,1)=5.35e-7
        a(3,2)=7.46e-6
        a(4,2)=1.01e-3
        a(5,2)=3.38e-2
        a(4,3)=2.99e-3
        a(5,3)=1.51e-4
        a(5,4)=1.12   
        t(2)=4.87e-7
        t(3)=1.308e-6
        t(4)=1.53162e-4
        t(5)=3.26888e-4
        weight(1)=1
        weight(2)=3
        weight(3)=5
        weight(4)=5
        weight(5)=1
      else if((ionum .eq. 7020) .or. (ionum .eq. 7021)
     1                          .or. (ionum .eq. 7022)) then
c
c  N III ion parameters
c
        a(2,1)=4.77e-5
        a(3,1)=339.
        a(4,1)=8.95
        a(3,2)=364.
        a(4,2)=59.0
        a(5,2)=251.
        a(4,3)=0.0
        a(5,3)=0.0
        a(5,4)=0.0
        t(2)=1.744e-6
        t(3)=5.71871e-4
        t(4)=5.72468e-4
        t(5)=5.73279e-4
        weight(1)=2
        weight(2)=4
        weight(3)=2
        weight(4)=4
        weight(5)=6
      else if((ionum .eq. 8000) .or. (ionum .eq. 8001)
     1                          .or. (ionum .eq. 8002)) then
c
c  O I ion parameters
c
        a(2,1)=8.92e-5
        a(4,1)=6.34e-3
        a(5,1)=2.88e-4
        a(3,2)=1.74e-5
        a(4,2)=2.11e-3
        a(5,2)=7.32e-2
        a(4,3)=7.23e-7
        a(5,4)=1.22   
        t(2)=1.583e-6
        t(3)=2.27e-6
        t(4)=1.58679e-4
        t(5)=3.37926e-4
        weight(1)=5
        weight(2)=3
        weight(3)=1
        weight(4)=5
        weight(5)=1

      else if((ionum .eq. 8010) .or. (ionum .eq. 8011)
     1                          .or. (ionum .eq. 8012)) then
c
c  O II ion parameters
c
        a(2,1)=3.82e-5
        a(3,1)=1.65e-4
        a(4,1)=5.64e-2
        a(5,1)=2.32e-2
        a(3,2)=1.2e-7
        a(4,2)=1.17e-1
        a(5,2)=6.15e-2
        a(4,3)=6.14e-2
        a(5,3)=1.02e-1
        a(5,4)=2.08e-11
        t(2)=2.68107e-4
        t(3)=2.68302e-4
        t(4)=4.04675e-4
        t(5)=4.04686e-4
        weight(1)=4
        weight(2)=6
        weight(3)=4
        weight(4)=4
        weight(5)=2
c
      else if((ionum .eq. 8020) .or. (ionum .eq. 8021)
     1                          .or. (ionum .eq. 8022)) then
c
c  O III ion parameters
c
        a(2,1)=2.62e-5
        a(3,1)=3.02e-11
        a(4,1)=2.74e-6
        a(3,2)=9.76e-5
        a(4,2)=6.74e-3
        a(5,2)=2.23e-1
        a(4,3)=1.96e-2
        a(5,3)=7.85e-4
        a(5,4)=1.78
        t(2)=1.132e-6
        t(3)=3.062e-6
        t(4)=2.02733e-4
        t(5)=4.318577e-4
        weight(1)=1
        weight(2)=3
        weight(3)=5
        weight(4)=5
        weight(5)=1
c
      else if((ionum .eq. 10020) .or. (ionum .eq. 10021)
     1                           .or. (ionum .eq. 10022)) then
c
c  Ne III ion parameters
c
        a(2,1)=5.97e-3
        a(3,1)=2.18e-8
        a(4,1)=1.71e-1
        a(5,1)=3.94e-3
        a(3,2)=1.15e-3
        a(4,2)=5.42e-2
        a(5,2)=2.0
        a(4,3)=8.51e-6
        a(5,4)=2.71
        t(2)=6.429e-6
        t(3)=9.205e-6
        t(4)=2.58408e-4
        t(5)=5.57506e-4
        weight(1)=5
        weight(2)=3
        weight(3)=1
        weight(4)=5
        weight(5)=1
c
      else if((ionum .eq. 10030) .or. (ionum .eq. 10031)
     1                           .or. (ionum .eq. 10032)) then
c
c  Ne IV ion parameters
c
        a(2,1)=4.84e-4
        a(3,1)=5.54e-3
        a(4,1)=5.21e-1
        a(5,1)=1.27
        a(3,2)=1.48e-6
        a(4,2)=1.15e-1
        a(5,2)=4.0e-1
        a(4,3)=3.93e-1
        a(5,3)=4.37e-1
        a(5,4)=2.68e-9
        t(2)=4.12346e-4
        t(3)=4.12795e-4
        t(4)=6.24346e-4
        t(5)=6.24413e-4
        weight(1)=4
        weight(2)=6
        weight(3)=4
        weight(4)=2
        weight(5)=4
c
      else if((ionum .eq. 10040) .or. (ionum .eq. 10041)
     1                           .or. (ionum .eq. 10042)) then
c
c  Ne V ion parameters
c
        a(2,1)=1.28e-3
        a(3,1)=5.08e-9
        a(4,1)=2.37e-5
        a(3,2)=4.59e-3
        a(4,2)=1.31e-1
        a(5,2)=4.21
        a(4,3)=3.65e-1
        a(5,3)=6.69e-3
        a(5,4)=2.85
        t(2)=4.124e-6
        t(3)=1.1101e-5
        t(4)=3.02915e-4
        t(5)=6.39136e-4
        weight(1)=1
        weight(2)=3
        weight(3)=5
        weight(4)=5
        weight(5)=1
c
      else if((ionum .eq. 14020) .or. (ionum .eq. 14021)
     1                           .or. (ionum .eq. 14022)) then
c
c  Si III ion parameters
c
        a(2,1)=0.0
        a(3,1)=1.67e4
        a(4,1)=1.2e-2
        a(5,1)=2.60e9
        a(3,2)=3.82e-5
        a(4,2)=3.20e-9
        a(5,2)=1.82e-2
        a(4,3)=2.42e-4
        a(5,3)=2.22
        a(5,4)=2.19e-2
        t(2)=5.27247e-4
        t(3)=5.28533e-4
        t(4)=5.31150e-4
        t(5)=8.28844e-4
        weight(1)=1
        weight(2)=1
        weight(3)=3
        weight(4)=5
        weight(5)=3
c
      else if((ionum .eq. 16010) .or. (ionum .eq. 16011)
     1                           .or. (ionum .eq. 16012)) then
c
c  S II ion parameters
c
        a(2,1)=8.82e-4
        a(3,1)=2.60e-4
        a(4,1)=9.06e-2
        a(5,1)=2.25e-1
        a(3,2)=3.35e-7
        a(4,2)=1.63e-1
        a(5,2)=1.33e-1
        a(4,3)=7.79e-2
        a(5,3)=1.79e-1
        a(5,4)=1.03e-6
        t(2)=1.48530e-4
        t(3)=1.48848e-4
        t(4)=2.45249e-4
        t(5)=2.45718e-4
        weight(1)=4
        weight(2)=4
        weight(3)=6
        weight(4)=2
        weight(5)=4
c
      else if((ionum .eq. 16020) .or. (ionum .eq. 16021)
     1                           .or. (ionum .eq. 16022)) then
c
c  S III ion parameters
c
        a(2,1)=4.72e-4
        a(3,1)=4.61e-8
        a(4,1)=5.82e-6
        a(3,2)=2.07e-3
        a(4,2)=2.21e-2
        a(5,2)=7.96e-1
        a(4,3)=5.76e-2
        a(5,3)=1.05e-2
        a(5,4)=2.22
        t(2)=2.972e-6
        t(3)=8.325e-6
        t(4)=1.1320e-4
        t(5)=2.7163e-4
        weight(1)=1
        weight(2)=3
        weight(3)=5
        weight(4)=5
        weight(5)=1
c
      else if((ionum .eq. 16030) .or. (ionum .eq. 16031)
     1                           .or. (ionum .eq. 16032)) then
c
c  S IV ion parameters
c
        a(2,1)=7.75e-3
        a(3,1)=5.50e+4
        a(4,1)=1.40e+2
        a(3,2)=3.39e+4
        a(4,2)=1.95e+4
        a(5,2)=3.95e+4
        a(4,3)=0.00e-2
        a(5,3)=0.00e-2
        a(5,4)=0.00
        t(2)=9.515e-6
        t(3)=7.118e-4
        t(4)=7.153e-4
        t(5)=7.207e-4
        weight(1)=2
        weight(2)=4
        weight(3)=2
        weight(4)=4
        weight(5)=6
c
      else if((ionum .eq. 17010) .or. (ionum .eq. 17011)
     1                           .or. (ionum .eq. 17012)) then
c
c  Cl II ion parameters
c
        a(2,1)=7.57e-3
        a(3,1)=4.57e-7
        a(4,1)=1.04e-1
        a(5,1)=1.97e-2
        a(3,2)=1.46e-3
        a(4,2)=2.92e-2
        a(5,2)=1.31
        a(4,3)=9.82e-6
        a(5,4)=2.06
        t(2)=6.96e-6
        t(3)=9.965e-6
        t(4)=1.16536e-4
        t(5)=2.78780e-4
        weight(1)=5
        weight(2)=3
        weight(3)=1
        weight(4)=5
        weight(5)=1
c
      else if((ionum .eq. 17020) .or. (ionum .eq. 17021)
     1                           .or. (ionum .eq. 17022)) then
c
c  Cl III ion parameters
c 
        a(2,1)=4.83e-3
        a(3,1)=7.04e-4
        a(4,1)=3.05e-1
        a(5,1)=7.54e-1
        a(3,2)=3.22e-6
        a(4,2)=3.03e-1
        a(5,2)=3.23e-1
        a(4,3)=1.0e-1
        a(5,3)=3.16e-1
        a(5,4)=7.65e-6
        t(2)=1.8053e-4
        t(3)=1.81186e-4
        t(4)=2.9812e-4
        t(5)=2.9907e-4
        weight(1)=5
        weight(2)=3
        weight(3)=1
        weight(4)=5
        weight(5)=1
c
      else if((ionum .eq. 17030) .or. (ionum .eq. 17031)
     1                           .or. (ionum .eq. 17032)) then
c
c  Cl IV ion parameters
c
        a(2,1)=2.14e-3
        a(3,1)=2.70e-7
        a(4,1)=1.54e-5
        a(3,2)=8.25e-3
        a(4,2)=7.23e-2
        a(5,2)=2.47
        a(4,3)=1.79e-1
        a(5,3)=2.62e-2
        a(5,4)=2.80
        t(2)=4.92e-6
        t(3)=1.3419e-5
        t(4)=1.37676e-4
        t(5)=3.2547e-4
        weight(1)=1
        weight(2)=3
        weight(3)=5
        weight(4)=5
        weight(5)=1
c
      else if((ionum .eq. 18020) .or. (ionum .eq. 18021)
     1                           .or. (ionum .eq. 18022)) then
c
c  Ar III ion parameters
c
        a(2,1)=3.08e-2
        a(3,1)=2.37e-6
        a(4,1)=3.14e-1
        a(5,1)=4.17e-2
        a(3,2)=5.17e-3
        a(4,2)=8.32e-2
        a(5,2)=3.91
        a(4,3)=2.21e-5
        a(5,3)=0.0
        a(5,4)=2.59
        t(2)=1.1121e-5
        t(3)=1.5702e-5
        t(4)=1.40100e-4
        t(5)=3.32657e-4
        weight(1)=5
        weight(2)=3
        weight(3)=1
        weight(4)=5
        weight(5)=1
c
      else if((ionum .eq. 18030) .or. (ionum .eq. 18031)
     1                           .or. (ionum .eq. 18032)) then
c
c  Ar IV ion parameters
c Transition probabilities from Kaufman & Sugar 1986, JPCRD 15, 321
c
c       a(2,1)=2.23e-2
c       a(3,1)=1.77e-3
c       a(4,1)=8.62e-1
c       a(5,1)=2.11
c       a(3,2)=2.3e-5
c       a(4,2)=6.03e-1
c       a(5,2)=7.89e-1
c       a(4,3)=0.119
c       a(5,3)=0.598
c       a(5,4)=4.94e-5
        a(2,1)=1.72e-2
        a(3,1)=2.07e-3
        a(4,1)=7.62e-1
        a(5,1)=1.88
        a(3,2)=2.3e-5
        a(4,2)=6.96e-1
        a(5,2)=7.08e-1
        a(4,3)=0.119
        a(5,3)=0.840
        a(5,4)=4.94e-5
        t(2)=2.10904e-4
        t(3)=2.12193e-4
        t(4)=3.48555e-4
        t(5)=3.50326e-4
        weight(1)=4
        weight(2)=4
        weight(3)=6
        weight(4)=2
        weight(5)=4
c
      else if((ionum .eq. 18040) .or. (ionum .eq. 18041)
     1                           .or. (ionum .eq. 18042)) then
c
c  Ar V ion parameters
c
        a(2,1)=7.99e-3
        a(3,1)=1.24e-6
        a(4,1)=3.50e-5
        a(3,2)=2.72e-2
        a(4,2)=0.204
        a(5,2)=6.55
        a(4,3)=0.476
        a(5,3)=5.69e-2
        a(5,4)=3.29
        t(2)=7.639e-6
        t(3)=2.0292e-5
        t(4)=1.62994e-4
        t(5)=3.79125e-4
        weight(1)=1
        weight(2)=3
        weight(3)=5
        weight(4)=5
        weight(5)=1
c
      endif
c
      return
      end
      subroutine parii(ionum,tem)
c
c     ________________________________________________________________
c     |                                                              |
c     |                           PAR II                             |
c     |                                                              |
c     |     Calculate temperature-dependent atomic parameters        |
c     |                                                              |
c     |     C:  collision strengths                                  |
c     |     CJ: coefficients for quadratic fit in T4(T/10000)        |
c     |                                                              |
c     ----------------------------------------------------------------
c
      real cj(3),ncrit,T4
      integer weight(5)
      integer kbin,ttyout,filout
      common /atom/ c(5,5),a(5,5),ncrit(5),t(5),weight
      common /tty/ kbin,ttyout,filout
c
      t4=1.0e-4*tem
      cj(1)=0.0
      cj(2)=0.0
      cj(3)=0.0
c
c  For new ions, appropriate parameters must be included in IF block
c
      if((ionum .eq. 6010) .or. (ionum .eq. 6011)
     1                     .or. (ionum .eq. 6012)) then
c
c  C II ion parameters   (Blum and Pradhan 1992, ApJS 80, 425)
c
        c(2,1)=2.1519
        c(3,1)=0.2425
        c(4,1)=0.3618
        c(5,1)=0.2349
        c(3,2)=0.1771
        c(4,2)=0.4774
        c(5,2)=1.0238
        c(4,3)=0.8237
        c(5,3)=0.8533
        c(5,4)=1.9818
c
c  C II ion parameters (old)
c
c       c(2,1)=2.90
c       c(3,1)=0.276
c       c(4,1)=0.409
c       c(5,1)=0.260
c       c(3,2)=0.197
c       c(4,2)=0.536
c       c(5,2)=1.16
c       c(4,3)=0.911
c       c(5,3)=0.920
c       c(5,4)=2.16
c
      elseif((ionum .eq. 6020) .or. (ionum .eq. 6021)
     1                         .or. (ionum .eq. 6022)) then
c
c  C III ion parameters
c
        cj(1)=1.299-0.426*t4+0.137*t4**2
        c(2,1)=cj(1)/9
        c(3,1)=cj(1)/3
        c(4,1)=cj(1)*5/9
        c(5,1)=3.15+1.61*t4-0.42*t4**2
        c(3,2)=0.783+0.133*t4-0.005*t4**2
        c(4,2)=0.479+0.202*t4-0.004*t4**2
        c(5,2)=0.00
        c(4,3)=2.05+0.63*t4-0.02*t4**2
        c(5,3)=0.0
        c(5,4)=0.0
c
      elseif((ionum .eq. 7000) .or. (ionum .eq. 7001)
     1                     .or. (ionum .eq. 7002)) then
c
c  N I ion parameters
c
        c(2,1)=-1.2124e-2+t4*(0.3616-t4*5.88e-2)
        c(3,1)=-8.386e-3+t4*(0.2419-t4*3.94e-2)
        c(4,1)=-2.601e-3+t4*(0.07003-t4*1.07e-2)
        c(5,1)=-5.074e-3+t4*(0.1395-t4*0.02126)
        c(3,2)=-4.166e-2+t4*(0.368-t4*0.05733)
        c(4,2)=1.227e-2+t4*(0.1046-t4*7.868e-3)
        c(5,2)=4.600e-2+t4*(0.244-t4*0.0240)
        c(4,3)=1.8601e-2+t4*(8.76e-2-t4*0.0092)
        c(5,3)=1.8265e-2+t4*(0.1406-t4*1.187e-2)
        c(5,4)=-3.265e-3+t4*(0.0704-t4*0.00387)
c
      else if((ionum .eq. 7010) .or. (ionum .eq. 7011)
     1                          .or. (ionum .eq. 7012)) then
c
c  N II ion parameters
c
        cj(1)=2.577+t4*(0.137-t4*0.03)
        cj(2)=0.353+t4*(-0.005807+t4*0.006003)
        c(2,1)=0.401
        c(3,1)=0.279
        c(4,1)=cj(1)/9
        c(5,1)=cj(2)/9
        c(3,2)=1.128
        c(4,2)=cj(1)*3/9
        c(5,2)=cj(2)*3/9
        c(4,3)=cj(1)*5/9
        c(5,3)=cj(2)*5/9
        c(5,4)=0.3993+t4*(0.0109+t4*0.001)
c
      else if((ionum .eq. 7020) .or. (ionum .eq. 7021)
     1                          .or. (ionum .eq. 7022)) then
c
c  N III ion parameters    (Blum and Pradhan 1992, ApJS 80, 425, T=15,000 K)
c
        c(2,1)=1.5541
        c(3,1)=0.2041
        c(4,1)=0.3098
        c(5,1)=0.2185
        c(3,2)=0.1621
        c(4,2)=0.4226
        c(5,2)=0.8800
        c(4,3)=1.1400
        c(5,3)=0.6983
        c(5,4)=2.1301
c
c  N III ion parameters   (old)
c
c       c(2,1)=0.701
c       c(3,1)=0.0952
c       c(4,1)=0.139
c       c(5,1)=0.080
c       c(3,2)=0.0616
c       c(4,2)=0.175
c       c(5,2)=0.390
c       c(4,3)=0.695
c       c(5,3)=0.397
c       c(5,4)=1.26
c
      else if((ionum .eq. 8000) .or. (ionum .eq. 8001)
     1                          .or. (ionum .eq. 8002)) then
c
c  O I ion parameters
c
        cj(1)=0.016656+t4*(0.3004-t4*0.02067)
        cj(2)=0.00201+t4*(0.03686-t4*0.002742)
        c(2,1)=-0.00214+t4*(0.09715+t4*0.003707)
        c(3,1)=-0.00109+t4*(0.0332-t4*0.00295)
        c(4,1)=cj(1)*5/9
        c(5,1)=cj(2)*5/9
        c(3,2)=-0.000128+t4*(0.0186+t4*0.00807)
        c(4,2)=cj(1)*3/9
        c(5,2)=cj(2)*3/9
        c(4,3)=cj(1)/9
        c(5,3)=cj(2)/9
        c(5,4)=0.0218+t4*(0.1076-t4*0.02234)
c
      else if((ionum .eq. 8010) .or. (ionum .eq. 8011)
     1                          .or. (ionum .eq. 8012)) then
c
c  O II ion parameters
c
        c(2,1)=0.7894+t4*(0.0098+t4*0.002282)
        c(3,1)=0.5269+t4*(5.204e-3+t4*1.923e-3)
        c(4,1)=0.26+t4*(1.001e-2-t4*2.915e-6)
        c(5,1)=0.1311+t4*(3.445e-3+t4*4.99e-4)
        c(3,2)=1.283+t4*(-0.1389+t4*0.0262)
        c(4,2)=0.7061+t4*(0.0235+t4*4.9e-4)
        c(5,2)=0.285+t4*0.01
        c(4,3)=0.3948+t4*(0.0123+t4*6.262e-4)
        c(5,3)=0.2644+t4*(0.0117-t4*9.16e-4)
        c(5,4)=0.2731+t4*(0.014-t4*2.919e-4)
c
      else if((ionum .eq. 8020) .or. (ionum .eq. 8021)
     1                          .or. (ionum .eq. 8022)) then
c
c  O III ion parameters
c
        cj(1)=1.835+t4*(0.3981-t4*0.06)
        cj(2)=0.2127+t4*(0.0767-0.013*t4)
        c(2,1)=0.4825+t4*(0.0806-t4*0.022)
        c(3,1)=0.2397+t4*(0.0381-t4*0.007)
        c(4,1)=cj(1)/9
        c(5,1)=cj(2)/9
        c(3,2)=1.1325+t4*(0.203-t4*0.05)
        c(4,2)=cj(1)/3
        c(5,2)=cj(2)/3
        c(4,3)=cj(1)*5/9
        c(5,3)=cj(2)*5/9
        c(5,4)=0.3763+t4*(0.3375-t4*0.105)
c
      else if((ionum .eq. 10020) .or. (ionum .eq. 10021)
     1                           .or. (ionum .eq. 10022)) then
c
c  Ne III ion parameters
c  Butler & Zeippen 1994, A&AS 108, 1
        cj(1)=1.357
c	+t4*(0.07967-t4*0.03282)
        cj(2)=0.151
c	+t4*(0.05331-t4*0.01487)
        c(2,1)=0.244
c	+t4*(0.1544-t4*0.05638)
        c(3,1)=0.306
c	+t4*(0.02765-t4*0.012639)
        c(4,1)=cj(1)*5/9
        c(5,1)=cj(2)*5/9
        c(3,2)=0.348
c        +t4*(0.06120-t4*0.02096)
        c(4,2)=cj(1)/3
        c(5,2)=cj(2)/3
        c(4,3)=cj(1)/9
        c(5,3)=cj(2)/9
        c(5,4)=0.269
c	+t4*(0.05962-t4*0.00862)
c old values
c       cj(1)=1.6028+t4*(0.07967-t4*0.03282)
c       cj(2)=0.12956+t4*(0.05331-t4*0.01487)
c       c(2,1)=1.0314+t4*(0.1544-t4*0.05638)
c       c(3,1)=0.2910+t4*(0.02765-t4*0.012639)
c       c(4,1)=cj(1)*5/9
c       c(5,1)=cj(2)*5/9
c       c(3,2)=0.3080+t4*(0.06120-t4*0.02096)
c       c(4,2)=cj(1)/3
c       c(5,2)=cj(2)/3
c       c(4,3)=cj(1)/9
c       c(5,3)=cj(2)/9
c       c(5,4)=0.1739+t4*(0.05962-t4*0.00862)
c
      else if((ionum .eq. 10030) .or. (ionum .eq. 10031)
     1                           .or. (ionum .eq. 10032)) then
c
c  Ne IV ion parameters
c
        c(2,1)=0.8473+t4*(-0.005832-t4*0.002886)
        c(3,1)=0.5652+t4*(-0.00452-t4*0.001535)
        c(4,1)=0.1496+t4*(9.539e-3-t4*3.437e-3)
        c(5,1)=0.2976+t4*(0.0232-t4*0.0088)
        c(3,2)=1.362+t4*(0.0217-t4*0.0186)
        c(4,2)=0.3057+t4*(0.0859-t4*0.026)
        c(5,2)=0.8014+t4*(0.1364-t4*0.0415)
        c(4,3)=0.308+t4*(0.0391-t4*0.0119)
        c(5,3)=0.4291+t4*(0.1103-t4*0.0336)
        c(5,4)=0.2883+t4*(0.0662-t4*0.0127)
c
      else if((ionum .eq. 10040) .or. (ionum .eq. 10041)
     1                           .or. (ionum .eq. 10042)) then
c
c  Ne V ion parameters
c
        cj(1)=1.6175+t4*(0.171-t4*0.01)
        cj(2)=0.3315+t4*(-0.1142+t4*0.034)
        c(2,1)=0.244
        c(3,1)=0.122
        c(4,1)=cj(1)/9
        c(5,1)=cj(2)/9
        c(3,2)=0.578
        c(4,2)=cj(1)/3
        c(5,2)=cj(2)/3
        c(4,3)=cj(1)*5/9
        c(5,3)=cj(2)*5/9
        c(5,4)=1.27-0.02*t4
c
      elseif((ionum .eq. 14020) .or. (ionum .eq. 14021)
     1                         .or. (ionum .eq. 14022)) then
c
c  Si III ion parameters
c
       T4=alog10(t4*10000.)
       cj(1)=539.070-457.121*T4+148.148*T4**2-21.4661*T4**3
     . +1.16536*T4**4
       cj(2)=575.809-522.592*T4+178.762*T4**2-26.8185*T4**3
     . +1.48065*T4**4
       c(2,1)=cj(1)/9
       c(3,1)=cj(1)/3
       c(4,1)=cj(1)*5/9
       c(5,1)=116.465-94.7015*T4+30.0337*T4**2-4.28699*T4**3
     . +0.241284*T4**4
       c(3,2)=21.1105-14.8546*T4+3.71890*T4**2-0.287946*T4**3
     . -3.72024e-3*T4**4
       c(4,2)=-43.3029+38.2889*T4-11.6750*T4**2+1.61607*T4**3
     . -8.92857e-2*T4**4
       c(5,2)=cj(2)/9
       c(4,3)=-41.9423+38.0232*T4-10.5783*T4**2+1.47321*T4**3
     . -9.67262e-2*T4**4
       c(5,3)=cj(2)/3
       c(5,4)=cj(2)*5/9
c
      else if((ionum .eq. 16010) .or. (ionum .eq. 16011)
     1                           .or. (ionum .eq. 16012)) then
c
c  S II ion parameters
c
        c(2,1)=3.06+t4*(-0.29+t4*0.02)
        c(3,1)=4.592+t4*(-0.437+t4*0.03)
        c(4,1)=0.78+t4*(-0.0242+t4*0.01)
        c(5,1)=1.33+t4*(0.264-t4*0.08)
        c(3,2)=8.79+t4*(-1.36+t4*0.16)
        c(4,2)=1.4575+t4*(0.111-t4*0.05)
        c(5,2)=3.282+t4*(0.223-t4*0.13)
        c(4,3)=2.502+t4*(0.1431-t4*0.09)
        c(5,3)=4.613+t4*(0.343-t4*0.17)
        c(5,4)=2.383+t4*(0.043-t4*0.05)
c
      else if((ionum .eq. 16020) .or. (ionum .eq. 16021)
     1                           .or. (ionum .eq. 16022)) then
c
c  S III ion parameters  
c  Tayal SS & Gupta GP 1999 ApJ 526, 544
c
        cj(1)=7.410+t4*(-0.738+t4*0.285)
        cj(2)=1.141+t4*(0.030+t4*0.0075)
        c(2,1)=5.065+t4*(-1.26+t4*0.175)
        c(3,1)=1.591+t4*(-0.407+t4*0.125)
        c(4,1)=cj(1)/9
        c(5,1)=cj(2)/9
        c(3,2)=9.846+t4*(-2.495+t4*0.512)
        c(4,2)=cj(1)/3
        c(5,2)=cj(2)/3
        c(4,3)=cj(1)*5/9
        c(5,3)=cj(2)*5/9
        c(5,4)=1.127+t4*(0.300-t4*0.033)
c
      else if((ionum .eq. 16030) .or. (ionum .eq. 16031)
     1                           .or. (ionum .eq. 16032)) then
c
c  S IV ion parameters
c
c       cj(1)=9.903+t4*(-2.017+t4*0.59)
c       cj(2)=1.135+t4*0.0522
        c(2,1)=8.54
        c(3,1)=0.60
        c(4,1)=1.05 
        c(5,1)=1.17 
        c(3,2)=1.03
        c(4,2)=1.91
        c(5,2)=2.73
        c(4,3)=3.37
        c(5,3)=2.92
        c(5,4)=7.23
c
      else if((ionum .eq. 17010) .or. (ionum .eq. 17011)
     1                           .or. (ionum .eq. 17012)) then
c
c  Cl II ion parameters
c
        c(2,1)=2.17
        c(3,1)=0.443
        c(4,1)=3.86*5/9
        c(5,1)=0.456*5/9
        c(3,2)=0.933
        c(4,2)=3.86/3
        c(5,2)=0.456/3
        c(4,3)=3.86/9
        c(5,3)=0.456/9
        c(5,4)=1.15
c
      else if((ionum .eq. 17020) .or. (ionum .eq. 17021)
     1                           .or. (ionum .eq. 17022)) then
c
c  Cl III ion parameters
c
        c(2,1)=1.1120+t4*(0.3837-t4*0.13750)
        c(3,1)=1.6660+t4*(0.5886-t4*0.21062)
        c(4,1)=0.3912+t4*(0.0085+t4*0.01505)
        c(5,1)=0.7810+t4*(0.0213+t4*0.02784)
        c(3,2)=4.0507+t4*(0.7741-t4*0.29264)
        c(4,2)=1.2051+t4*(0.6197-t4*0.18408)
        c(5,2)=1.8324+t4*(0.4803-t4*0.13531)
        c(4,3)=1.3373+t4*(0.2975-t4*0.08182)
        c(5,3)=3.2157+t4*(1.3672-t4*0.40935)
        c(5,4)=1.7478-t4*(0.0450-t4*0.05217)
c
      else if((ionum .eq. 17030) .or. (ionum .eq. 17031)
     1                           .or. (ionum .eq. 17032)) then
c
c  Cl IV ion parameters
c
        cj(1)=4.702+t4*(0.771-t4*0.01)
        cj(2)=1.712+t4*(0.791-t4*0.25)
        c(2,1)=0.475
        c(3,1)=0.4
        c(4,1)=cj(1)/9
        c(5,1)=cj(2)/9
        c(3,2)=1.5
        c(4,2)=cj(1)/3
        c(5,2)=cj(2)/3
        c(4,3)=cj(1)*5/9
        c(5,3)=cj(2)*5/9
        c(5,4)=0.3388+t4*(1.3214-t4*0.265)
c
      else if((ionum .eq. 18020) .or. (ionum .eq. 18021)
     1                           .or. (ionum .eq. 18022)) then
c
c  Ar III ion parameters
c
        c(2,1)=2.24
        c(3,1)=0.531
        c(4,1)=4.74*5/9
        c(5,1)=0.68*5/9
        c(3,2)=1.18
        c(4,2)=4.74/3
        c(5,2)=0.68/3
        c(4,3)=4.74/9
        c(5,3)=0.68/9
        c(5,4)=0.823
c
      else if((ionum .eq. 18030) .or. (ionum .eq. 18031)
     1                           .or. (ionum .eq. 18032)) then
c
c  Ar IV ion parameters; ZBL87 A&A 188, 251 
c
        c(2,1)=2.2661-t4*(1.2805-t4*0.32167)
        c(3,1)=3.3993-t4*(1.9217-t4*0.48293)
        c(4,1)=0.1637-t4*(0.0351-t4*0.01790)
        c(5,1)=0.3356-t4*(0.0817-t4*0.03930)
        c(3,2)=6.4696-t4*(0.3631-t4*0.04479)
        c(4,2)=1.4523+t4*(0.4098-t4*0.15965)
        c(5,2)=2.3424+t4*(0.2396-t4*0.11250)
        c(4,3)=1.7193+t4*(0.1391-t4*0.07235)
        c(5,3)=3.9465+t4*(0.7872-t4*0.30578)
        c(5,4)=2.1475+t4*(0.1047+t4*0.09440)
c
      else if((ionum .eq. 18040) .or. (ionum .eq. 18041)
     1                           .or. (ionum .eq. 18042)) then
c
c  Ar V ion parameters
c
        cj(1)=5.2075+t4*(-1.985+t4*0.55)
        cj(2)=1.133+t4*(0.127-t4*0.09)
        c(2,1)=0.257
        c(3,1)=0.32
        c(4,1)=cj(1)/9
        c(5,1)=cj(2)/9
        c(3,2)=1.04
        c(4,2)=cj(1)/3
        c(5,2)=cj(2)/3
        c(4,3)=cj(1)*5/9
        c(5,3)=cj(2)*5/9
        c(5,4)=1.27+t4*(-0.02+t4*2.4e-5)
c
      endif
      return
      end
      subroutine temden(ionum,den,tem,robs,pop,jl)
c
c     _________________________________________________________________
c     |                                                               |
c     |                           TEMDEN                              |
c     |                                                               |
c     |     ITER:     iterative index used to find TEM or DEN         |
c     |     RCAL:     'characteristic' ratio calculated in SOLVE      |
c     |     X1,X2:    initial (absolute) limits on TEM or DEN         |
c     |     XLO,XHI:  current range on iterated TEM or DEN            |
c     |     RLO,RHI:  ratios calculated at XLO, XHI                   |
c     |     X:        current 'best guess' for TEM or DEN             |
c     |     ERR:      percentage tolerance in X                       |
c     |                                                               |
c     -----------------------------------------------------------------
c
c     real a(5,5), c(5,5), jl(5,5), cj(3), pop(5), t(5)
      real jl(5,5), pop(5)
      integer kbin,ttyout,filout

      common /tty/ kbin,ttyout,filout

      err = 0.001
c
c  Note: setting ERR to too high a value can lead to poor results!
c
c  Temperature limits
c
c  Note: watch out for machine precision limitations on calculated
c        emissivities at very low temperatures!!
c
      if(mod(ionum,10) .eq. 2) then
        x1=1.0e3
        x2=4.0e4
        write(filout,20)
c
c  Density limits
c
      elseif(mod(ionum,10) .eq. 1) then
        x1=1.0
        x2=1.0e8
        write(filout,10)
      endif
c
c  Permit up to 100 iterations
c
      do 100 iter=1,100
         if(iter .eq. 1) then
           xlo=x1
           xhi=x2
           x=xlo
         else if(iter .eq. 2) then
           xhi=x2
           x=xhi
         else
           x=sqrt(xlo*xhi)
         endif
c
         if(mod(ionum,10) .eq. 2) then
           tem=x
         else if(mod(ionum,10) .eq. 1) then
           den=x
         endif
c
         call parii(ionum,tem)
         call solve(den,tem,pop,jl)
c
c  For new ions, enter appropriate ratio in IF block below
c
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
c  Output intermediate iterations
c
         if(mod(ionum,10) .eq. 1) then
           write(filout,30) iter,den,xlo,tem,xhi,(rcal-robs)/robs
         else if(mod(ionum,10) .eq. 2) then
           write(filout,30) iter,tem,xlo,den,xhi,(rcal-robs)/robs
         endif
c
c  Test for convergence
c
         if(iter .eq. 1) then
           rlo=rcal
         else if((iter .eq. 2) .and. ((rcal-robs)/(rlo-robs) .lt. 
     . 0.0)) then
           rhi=rcal
         else if((iter .eq. 2) .and. ((rcal-robs)/(rlo-robs) .gt.
     . 0.0)) then
           write(filout,40) amin1(rlo,rcal),amax1(rlo,rcal)
           write(filout,50) rcal,robs,rlo
           robs=0.0
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
  100 continue
c
      write(ttyout,*) 'Did not converge after 100 iterations; stopping'
   10 format(1x,'Iter#',6x,'Density',8x,'XLO',2x,'Temperature',7x,
     1'XHI',3x,' (RCAL-ROBS)/ROBS')
   20 format(1x,'Iter#',2x,'Temperature',8x,'XLO',6x,'Density',7x,
     1'XHI',3x,' (RCAL-ROBS)/ROBS')
   30 format(3x,i3,2(3x,1pe10.3,1x,1pe10.3),9x,1pe9.2)
   40 format(' (RCAL-ROBS)/(RLO-ROBS) is > 0, so the entered RATIO',/,
     1    ' predicts unreasonable conditions. Please try again.',/,
     2    ' The allowable RATIO range is:',1pe9.2,'--',1pe9.2)
   50 format(' RCAL = ',f8.3,' ROBS = ',f8.3,' RLO = ',f8.3)
c
      return
      end
      subroutine solve(den,tem,pop,jl)
c
c     _________________________________________________________________
c     |                                                               |
c     |                            SOLVE                              |
c     |                                                               |
c     |     TS:       square root of temperature                      |
c     |     POP1:     storage array for level populations             |
c     |     POP:      level population array                          |
c     |     A:        radiative transition probabilities              |
c     |     C:        collision strengths                             |
c     |     Q:        collision transition probabilities              |
c     |     E:        relative energy separation                      |
c     |     KT,HC,CQ: "constants"                                     |
c     |                                                               |
c     -----------------------------------------------------------------
c
      real e(5,5), q(5,5), a1(5,5), a2(5,5), jl(5,5)
      real pop(5), b1(5), b2(5), dtl(5), pop1(5), ncrit, kt, hc
      integer weight(5)
      integer kbin,ttyout,filout
c
      common /atom/c(5,5),a(5,5),ncrit(5),t(5),weight
      common /tty/kbin,ttyout,filout
c
      kt=(1.3807e-16)*tem
      hc=1.9865e-08
      cq=8.629e-06
      ts=sqrt(tem)
c
c  Transition energy differences
c
      do 100 i=1,5
         do 110 j=i+1,5
            e(j,i)=hc*(t(j)-t(i))
  110    continue
  100 continue
c
c  Collision transition probabilities (if i=j, Q=0)
c
      do 200 i=1,5
         do 210 j=1,5
            if(j .gt. i) then
              q(i,j)=cq*c(j,i)*exp(-e(j,i)/kt)/(weight(i)*ts)
            else if(j .eq. i) then
              q(i,j)=0.0
            else
              q(i,j)=cq*c(i,j)/(weight(i)*ts)
            endif
  210    continue
  200 continue
c
c  Critical density calculation
c
      do 300 i=2,5
         asum=0.0
         qsum=0.0
         do 400 j=1,i-1
            asum=asum+a(i,j)
            qsum=qsum+q(i,j)
  400    continue
         ncrit(i)=asum/qsum
  300 continue
c
c  Matrix manipulation
c
      do 500 i=1,5
         do 600 j=1,5
            if(j .eq. i) then
              a1(i,j)=0.0
            else if(j .gt. i) then
              a1(i,j)=den*q(i,j)
            else
              a1(i,j)=den*q(i,j)+a(i,j)
            endif
  600    continue
  500 continue
c
      do 700 i=1,5
         b1(i)=0.0
         do 800 j=1,5 
            b1(i)=b1(i)+a1(i,j)
  800    continue
  700 continue
c
      do 900 i=1,5
         b2(i)=-a1(1,i)
         a1(i,i)=-b1(i)
  900 continue
c
      do 1000 i=1,5
         do 1100 j=1,5
            do 1200 k=1,5
               a2(j,k)=a1(j,k)
 1200       continue
 1100    continue
         do 1300 j=1,5
            a2(i,j)=b2(j)
 1300    continue
         dtl(i)=1
         do 1400 j=2,4
            t1=a2(j,j)
            if(t1 .eq. 0.0) then
              write(filout,*)'Divide by 0 during T1 normalization
     1 in SOLVE'
              stop
            endif
            dtl(i)=dtl(i)*t1
            do 1500 k=j,5
               a2(j,k)=a2(j,k)/t1
 1500       continue
            do 1600 k=j+1,5
               t2=a2(k,j)
               do 1700 l=j,5 
                  a2(k,l)=a2(k,l)-t2*a2(j,l)
 1700          continue
 1600       continue
 1400    continue
c
         dtl(i)=dtl(i)*a2(5,5)
 1000 continue
c
      pop1(1)=1
      do 1800 i=2,5
         pop1(i)=dtl(i)/dtl(1)
 1800 continue
c
      t1=0
      do 1900 i=1,5
         t1=t1+pop1(i)
 1900 continue
c
c  Finally, level populations, emissivities
c  If jl(i,j) < 1.0e-38, set equal to 1.0e-38 to avoid underflow 
c
      do 2000 i=1,5
         pop(i)=pop1(i)/t1
         do 2100 j=1,5
            jl(i,j)=pop(i)*a(i,j)*e(i,j)
            if(jl(i,j) .lt. 1.0e-38) jl(i,j)=1.0e-38
 2100    continue
 2000 continue
c
      return
      end
