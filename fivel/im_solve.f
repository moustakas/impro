c+
c NAME:
c     IM_SOLVE()
c
c PURPOSE:
c     Given the physical conditions of the nebula and the relevant 
c     atomic parameters compute the emission-line emissivities.
c
c CALLING SEQUENCE:
c     call im_solve(den,tem,t,weight,a,c,pop,jl,ncrit,e,q)
c
c INPUTS:
c     den    - electron density [cm-3]
c     tem    - electron temperature [K]
c     t      - energy-level separation [5] (IM_PARI)
c     weight - statistical weights for each level [5] (IM_PARI)
c     a      - radiative transition probabilities [5,5] (IM_PARI)
c     c      - collision strengths [5,5] (IM_PARII)
c
c OUTPUTS:
c     pop   - level population array [5]
c     jl    - emission-line emissivities [5,5]
c     ncrit - critical densities for each transition [5]
c     e     - relative energy separation [5,5]
c     q     - collision transition probabilities [5,5]
c
c INTERNAL VARIABLES:
c     kt, hc, cq - atomic constants
c     i, j, k, l - do loop integers
c     asum, qsum, a1, a2, b1, b2, dtl, t1, t2, pop1, ts
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
c
c-

      subroutine im_solve(den,tem,t,weight,a,c,pop,jl,ncrit,e,q)

      implicit none

      real*8 den, tem, t(5), a(5,5), c(5,5) ! input
      integer*4 weight(5)       ! input

      real*8 asum, qsum, a1(5,5), a2(5,5), b1(5), b2(5), ! internal
     &     dtl(5), pop1(5), t1, t2, kt, hc, cq, ts
      integer*4 i, j, k, l        ! internal
      
      real*8 pop(5), jl(5,5), ncrit(5), e(5,5), q(5,5) ! output
      
      include 'tty.dec'

      kt=(1.3807e-16)*tem
      hc=1.9865e-08
      cq=8.629e-06
      ts=sqrt(tem)

c     Transition energy differences

      do 100 i=1,5
         do 110 j=i+1,5
            e(j,i)=hc*(t(j)-t(i))
 110     continue
 100  continue

c     Collision transition probabilities (if i=j, Q=0)

      do 200 i=1,5
         do 210 j=1,5
            if(j .gt. i) then
               q(i,j)=cq*c(j,i)*exp(-e(j,i)/kt)/(weight(i)*ts)
            else if(j .eq. i) then
               q(i,j)=0.0
            else
               q(i,j)=cq*c(i,j)/(weight(i)*ts)
            endif
 210     continue
 200  continue

c     Critical density calculation

      do 300 i=2,5
         asum=0.0
         qsum=0.0
         do 400 j=1,i-1
            asum=asum+a(i,j)
            qsum=qsum+q(i,j)
 400     continue
         ncrit(i)=asum/qsum
 300  continue

c     Matrix manipulation

      do 500 i=1,5
         do 600 j=1,5
            if(j .eq. i) then
               a1(i,j)=0.0
            else if(j .gt. i) then
               a1(i,j)=den*q(i,j)
            else
               a1(i,j)=den*q(i,j)+a(i,j)
            endif
 600     continue
 500  continue
c     
      do 700 i=1,5
         b1(i)=0.0
         do 800 j=1,5 
            b1(i)=b1(i)+a1(i,j)
 800     continue
 700  continue
c     
      do 900 i=1,5
         b2(i)=-a1(1,i)
         a1(i,i)=-b1(i)
 900  continue
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
     1              in SOLVE'
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

c     Finally, level populations, emissivities
c     If jl(i,j) < 1.0e-38, set equal to 1.0e-38 to avoid underflow 

      do 2000 i=1,5
         pop(i)=pop1(i)/t1
         do 2100 j=1,5
            jl(i,j)=pop(i)*a(i,j)*e(i,j)
            if(jl(i,j) .lt. 1.0e-38) jl(i,j)=1.0e-38
 2100    continue
 2000 continue

      return
      end
