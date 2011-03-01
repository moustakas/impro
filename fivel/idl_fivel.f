c+
c NAME:
c     IDL_FIVEL [The LICK 5-level Atom Program]
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
c     call idl_fivel
c
c INPUTS:
c
c OUTPUTS:
c
c INTERNAL VARIABLES:
c
c COMMON BLOCKS:
c     None.
c
c PROCEDURES USED:
C     IDL_PARI()   - Contains atomic data independent of T
C     IDL_PARII()  - Contains atomic data dependent on T
C     IDL_SOLVE()  - Performs matrix algebra. Finds level populations
C                    (POP) and line emissivities (JL) as a function 
c                    of N and T, and critical densities (NCRIT)
c     IDL_TEMDEN() - Employs an iterative technique (geometric 
c                    bisection) to find 'T' or 'N'.  Used only 
c                    when T or N is unknown.
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

      subroutine idl_fivel(ictg,ionum,tem,den,robs,a,c,q,jl,e,
     &     t,pop,ncrit,weight,temden_error)

      implicit none
      
      integer*4 ionum, ictg     ! input
      real*8 den, tem, robs     ! input

      real*8 a(5,5), c(5,5), q(5,5), jl(5,5), e(5,5), t(5), 
     &     pop(5), ncrit(5)
      integer*4 weight(5), temden_error
      
c     Physical conditions

      if (ictg .eq. 1) then

         call idl_temden(ionum,den,tem,robs,temden_error)
         call idl_pari(ionum,a,t,weight)
         call idl_parii(ionum,tem,c)
         call idl_solve(den,tem,t,weight,a,c,pop,jl,ncrit,e,q)

      endif 
         
c     Level populations, line ratios

      if (ictg .eq. 2) then

         call idl_pari(ionum,a,t,weight)
         call idl_parii(ionum,tem,c)
         call idl_solve(den,tem,t,weight,a,c,pop,jl,ncrit,e,q)

      endif

      return 
      end
