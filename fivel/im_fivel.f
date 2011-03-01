c+
c NAME:
c     IM_FIVEL [The LICK 5-level Atom Program]
c
c PURPOSE:
c     Certain astrophysically common ions have emission- line intensity
c     ratios which are sensitive to either temperature (T) or number
c     density (N). From an observed line ratio (ROBS), then, given a
c     tem- perature, the density can be calculated; or, given the
c     density, the temperature can be calculated.  Similarly, given the
c     temperature and density, the emission-line emissivities can be
c     calculated.
c                                                             
c     Program Structure:                                  
c        1. "Physical conditions"                            
c            (a) Specify 'T' option to use observed line     
c                ratio and density to find temperature      
c            (b) Specify 'N' option to use observed line     
c                ratio and temperature to find density       
c                                                            
c        2. "Line emissivities"                              
c             Supply temperature and density to determine     
c             the level populations, critical densities,      
c             and line emissivities
c
c CALLING SEQUENCE:
c     call im_fivel
c
c INPUTS:
c     See IM_INPUT().
c
c OUTPUTS:
c     See IM_OUTPUT().
c
c INTERNAL VARIABLES:
c
c COMMON BLOCKS:
c     None.
c
c PROCEDURES USED:
c     IM_INPUT()  - Supplies data for calculation
C     IM_OUTPUT() - Displays final resultS
C     IM_PARI()   - Contains atomic data independent of T
C     IM_PARII()  - Contains atomic data dependent on T
C     IM_SOLVE()  - Performs matrix algebra. Finds level populations
C                   (POP) and line emissivities (JL) as a function 
c                   of N and T, and critical densities (NCRIT)
c     IM_TEMDEN() - Employs an iterative technique (geometric 
c                   bisection) to find 'T' or 'N'.  Used only 
c                   when T or N is unknown.
c
c INPUT FILES:
c     tty.dec
c
c COMMENTS:
c     The default logical units may have to be changed de- pending on
c     the computer/operating system. Currently, unit '*' is used for
c     input and "interactive output", while '6' is used for writing out
c     the calculations.
c
c MODIFICATION HISTORY:
c     J. Moustakas, 2004 July 4, U of A - adopted from de Robertis,
c        Dufour, & Hunt's FIVEL.F
c-

      program im_fivel

      implicit none
      
      integer*4 ionum, ictg
      real*8 den, tem, robs

      real*8 a(5,5), c(5,5), q(5,5), jl(5,5), t(5), pop(5), 
     &     ncrit(5), e(5,5), jhb
      integer*4 weight(5), temden_error

      logical phyg, popg

      include 'tty.dec'

      open(filout,file='fiveout.dat')
      rewind filout

      write(ttyout,10)
 100  write(ttyout,15)
      read(kbin,*,end=100,err=100) ictg
      write(ttyout,*) ' '

c     
c     Physical conditions
c     
      if (ictg .eq. 1) then
c     
c     Ion code: Atomic number//ionization stage//N(=1),T(=2)
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
 110     write(ttyout,20)
         read(kbin,*,end=110,err=110) ionum
c     
c     For new ions add appropriate parameters in IONUM and PHYG block
c     
         phyg=(ionum .eq. 6021).or.(ionum .eq. 7001).or.
     1        (ionum .eq. 7012).or.(ionum .eq. 8011).or.
     2        (ionum .eq. 8012).or.(ionum .eq. 8022).or.
     3        (ionum .eq.10022).or.(ionum .eq.10031).or.
     4        (ionum .eq.10042).or.(ionum .eq.16011).or.
     5        (ionum .eq.16012).or.(ionum .eq.16022).or.
     6        (ionum .eq.17021).or.(ionum .eq.17032).or.
     7        (ionum .eq.18022).or.(ionum .eq.18031).or.
     8        (ionum .eq.18032).or.(ionum .eq.18042)

         if (phyg) then

            write(filout,*)' '
            write(filout,*)'Physical conditions calculations 
     1           for ion number: ',ionum
            write(filout,*)' '

            call im_input(ionum,den,tem,robs)
            call im_temden(ionum,den,tem,robs,temden_error)
            call im_pari(ionum,a,t,weight)
            call im_parii(ionum,tem,c)
            call im_solve(den,tem,t,weight,a,c,pop,jl,ncrit,e,q)
            call im_output(ionum,den,tem,pop,jl,ncrit,t,jhb)
               
         else

            write(filout,*)'No parameters for that ion; enter another' 

         endif 

c     
c     Level populations, line ratios
c     
      else if(ictg .eq. 2) then
         popg=.true.
c     
c     Ion code: Atomic number//ionization stage//0(=emissivities)
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
 120     write(ttyout,20)
         read(kbin,*,err=120,end=120) ionum
c     
c     For new ions, add appropriate parameters in IONUM and POPG block.
c     
         phyg=(ionum .eq. 6010).or.
     .        (ionum .eq. 6020).or.(ionum .eq. 7000).or.
     .        (ionum .eq. 7010).or.(ionum .eq. 7020).or. 
     1        (ionum .eq. 8000).or.(ionum .eq. 8010).or.
     2        (ionum .eq. 8020).or.(ionum .eq.10020).or.
     3        (ionum .eq.10030).or.(ionum .eq.10040).or.
     3        (ionum .eq.14020).or.
     4        (ionum .eq.16010).or.(ionum .eq.16020).or.
     5        (ionum .eq.16030).or.
     6        (ionum .eq.17010).or.(ionum .eq.17020).or.
     7        (ionum .eq.17030).or.(ionum .eq.18020).or.
     8        (ionum .eq.18030).or.(ionum .eq.18040)

         if (popg) then

            write(filout,*) 'Line emissivity calculations',
     1           'for ion number ',ionum
            write(filout,*) ' '

            call im_input(ionum,den,tem,robs)
            call im_pari(ionum,a,t,weight)
            call im_parii(ionum,tem,c)
            call im_solve(den,tem,t,weight,a,c,pop,jl,ncrit,e,q)
            call im_output(ionum,den,tem,pop,jl,ncrit,t,jhb)

         else

            write(filout,*) 'No parameters for that ion; try another'

         endif

c     Stop

      else if (ictg .eq. 3) then
         stop
      else
c     
c     Try again
c     
         write(ttyout,*) 'Enter 1, 2, or 3 only'

      endif
      go to 100
 200  continue

 10   format(//,20x,'<<< LICK 5-LEVEL ATOM PROGRAM >>>',//)
 15   format(/,' Calculate PHYSICAL CONDITIONS(1), LINE RATIOS(2),',
     1     ' or STOP(3): ')
 20   format(' Select one of the above by typing number: ')

      end
