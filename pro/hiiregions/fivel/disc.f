      program disc
      character id*4
      conv=0.017453292

      open(6,file='coords.txt',status='old')
      open(7,file='metal.out',status='unknown')
      write(7,*)'#ID  Gal. dist. (")   Z/Zsun 12+log(O/H) rho/rho0'

      do 100 while (.true.)

      read(6,10)id,ra1,ra2,ra3,d1,d2,d3
 10   format(a4,f3.0,f3.0,f6.3,2x,f4.0,f3.0,f5.2)


c     galactic center coordinates: NCG 300

      dec0=-37.68305556*conv
      ra0=0.914888888*15.*conv*cos(dec0)
      theta=111.*conv
      xincl=43.65*conv

c     galactic metallicity gradient (Urbaneja et al. 2004)
c      rho0 in arcsec
c     oh0: central 12+log(O/H)
c     zgrad: dex/rho0 metallicity gradient

      rho0=585
      oh0=8.52
      zgrad=-0.28

      if(d1.gt.0) then
         dec=(d1+d2/60.+d3/3600.)*conv
      elseif(d1.lt.0) then
         dec=(d1-d2/60.-d3/3600.)*conv
       endif

      ra=(ra1+ra2/60.+ra3/3600.)*15.*conv*cos(dec0)

c     calculate position angle with respect to galaxy center: phi
c     and projected galactocentric distance: P
      x1=(ra-ra0)
      y1=(dec-dec0)

      P=sqrt(x1**2+y1**2)
      angle=abs(atan(y1/x1))

      if (x1.ge.0..and.y1.ge.0.) then
         phi=90.*conv-angle
      elseif (x1.ge.0..and.y1.lt.0.) then
         phi=90.*conv+angle
      elseif (x1.lt.0..and.y1.lt.0.) then
         phi=270.*conv-angle
      elseif (x1.lt.0..and.y1.ge.0.) then
         phi=270.*conv+angle
      endif


      x=sqrt((P*cos(phi-theta))**2)
      y=sqrt(P**2-x**2)/cos(xincl)

      Gdist=sqrt(x**2+y**2)
      Gdist=Gdist*3600./conv
      P=P*3600./conv
      xxx=x*3600./conv
      yyy=y*3600./conv


c     calculate abundance from galactocentric distance and gradient parameters
      f=Gdist/rho0
      deltaz=f*zgrad
      z=10**((oh0+deltaz)-8.69)

      write(7,20)id,Gdist,z,(oh0+deltaz),f


 100  continue

 20   format(a4,f7.1,10x,f7.2,2f9.2)

      end

