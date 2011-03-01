CDOC ...................................................................
C    Name           hubble
C    
C    Synopsis       call hubble(l,b,vh,vcor,triple)
C    
C    Application    
C       This routine computes the distance given the radial velocity (VH)
C       and the galactic coordinates (L) and (B).
C       It uses the three components velocity field model from 
C       Burstein et al. 1989; Faber+Burstein 1988, Large-scale Motions:
C               Great Attractor: l=309, b=18, d=4200km/s
C                                V(LG)=535, n=1.7, c=1350/4200
C               Virgo infall:    l(vir)=284. b(vir)=75, d=1360km/sec
C                                V(Virgo)=100 n=3, c=0.5
C               Local Group anomaly: bulk motion 360 km/s
C                                toward l=199 b=0
C               Reduction to CMB  (Lubin et al. 1985)
C                                l(cmb)=267, b(cmb)=50, V(cmb)=377
C
C       The triple valued region extend 20 degrees around each attractor.     
C       When a galaxy is found in a triple valued region:
C             1- (TRIPLE) take a non-zero value (see below)
C             2- the `nearby' solution is computed
C    Note
C       In average the corrected velocity (VCOR) is greater than the 
C       heliocentric velocity. For distance larger than an attractor distance,
C       the mean effect of this attractor in a galactocentric infall (obvious).
C       Only for distance smaller than the distance of the attractor, and
C       directions close to the direction of the attractor the effect is
C       opposite.
C
C    Auteur         PhP
C
C    Methode
C
C    Arguments
C       L      R4 In  [deg] Galactic longitude
C       B      R4 In  [deg] Galactic latitude
C       VH     R4 In  [km/s] Heliocentric cz
C       VCOR   R4 Out [km/s] Corrected cz
C       TRIPLE R4 Out 0= Not in a triple valued zone
C                     1= In GA Triple valued Zone
C                     2= In Virgo Triple Valued Zone
C                    -1= In the core of GA
C                    -2= In the core of Virgo
C
C    Sous_Programmes none
C
C    References (J. Moustakas, 2003 July 29, U of A):
C       Kolatt et al 1995, MNRAS, 275, 797
C       Burstein 2000, ASP Conf. Ser. 201: Cosmic Flows Workshop, p. 178
C
CDOC ...................................................................
      subroutine hubble(l,b,vh,vcor,triple)
c      include 'plparam.f'    !Modified by C.T.
c
      real*4 l,b,vh,vcor
      real*4 lr,br
      real*4 lcmb,bcmb,vcmb
      real*4 lvir,bvir,rv,nv
      real*4 lga,bga,rga,nga
      real*4 llg,blg,vlg
c     
      triple=0.
      us=0.
      usv=0.
      PLPrad = 180.0 / 3.1415926  !Added by C.T.

c     write(*,*) 'INPUT:  ', l, b, vh, vcor, triple
c
c  Virgo Infall center

c the following quantities were updated by J. Moustakas (03may21)
      
      rv=1330.0                 ! distance to Virgo [Huchra 1996]
      vv=220.0                  ! local infall to Virgo [Burstein 1999]
      lvir=283.8/PLPrad         ! [Burstein 1999]
      bvir=74.5/PLPrad          ! [Burstein 1999]

c     rv=1360
c     lvir=284/PLPrad
C     bvir=75/PLPrad
c     vv=100.

      cv=0.25*rv
      nv=1.5
      xv=rv*cos(bvir)*cos(lvir)	!x y z position in helioc coord
      yv=rv*cos(bvir)*sin(lvir)
      zv=rv*sin(bvir)
c
c  Great Attractor center

      rga=4000.0                ! Distance to the GA [Kolatt 1995]
      lga=320.0/PLPrad          ! [Kolatt 1995]
      bga=0.0/PLPrad            ! [Kolatt 1995]
      vga=535.0                 ! local infall to the GA

c     rga=4200.
c     lga=309./PLPrad
c     bga=18./PLPrad
c     vga=535. 

      nga=1.7
      cga=0.34*rga
      xga=rga*cos(bga)*cos(lga)	!x y z in helioc coord
      yga=rga*cos(bga)*sin(lga)
      zga=rga*sin(bga)
c
c  Local Group anomaly

      vlg=-273.0                ! [Burstein 1999]
      llg=222.9/PLPrad          ! [Burstein 1999]
      blg=-4.7/PLPrad           ! [Burstein 1999]

c     vlg=-360.
c     llg=199./PLPrad
c     blg=-4./PLPrad
c
c
c  heliocentric velocity 

      br=b/PLPrad
      lr=l/PLPrad
      x0=vh*cos(br)*cos(lr)
      y0=vh*cos(br)*sin(lr)
      z0=vh*sin(br)
c
c
c  reduction to CMB

      vcmb=368.6                ! [Burstein 1999]
      lcmb=264.7/PLPrad         ! [Burstein 1999]
      bcmb=48.2/PLPrad          ! [Burstein 1999]

c     vcmb=377.                 !Lubin et al. 1985
c     lcmb=267./PLPrad
c     bcmb=50./PLPrad

      x1=vcmb*cos(bcmb)*cos(lcmb)
      y1=vcmb*cos(bcmb)*sin(lcmb)
      z1=vcmb*sin(bcmb)
      vcorcmb=(x1*x0+y1*y0+z1*z0)/vh

c     
c  Local Group bulk motion
      xlg=0.
      ylg=0.
      zlg=0.
      if (vh.lt.700) then
         vvlg=vlg*(1.-abs(vh)/700.)**3 !smooth transition to 0
      else
         vvlg=0.
      endif
      xlg=vvlg*cos(blg)*cos(llg)
      ylg=vvlg*cos(blg)*sin(llg)
      zlg=vvlg*sin(blg)
      vcorlg=(xlg*x0+ylg*y0+zlg*z0)/vh
c     
c First guess for GA effect:
      r=(vh+vcorcmb+vcorlg)/vh
      x4=x0*r-xga
      y4=y0*r-yga
      z4=z0*r-zga
      r4=sqrt(x4**2+y4**2+z4**2)
      if(r4.ge.500.)then
         du=-vga*(r4/rga)*
     s        ((rga**2+cga**2)/(r4**2+cga**2))**((nga+1.)/2.)
      else
         du=0.
         triple=-1
         vcor=rga
         return
      endif
      dv4=du*(x4*x0+y4*y0+z4*z0)/vh/r4
c     
c First guess for Virgo effect:
      r=(vh+vcorcmb+vcorlg-dv4)/vh
      x3=x0*r-xv
      y3=y0*r-yv
      z3=z0*r-zv
      r3=sqrt(x3**2+y3**2+z3**2)
      if(r3.ge.cv)then
         du=-vv*((rv**2+cv**2)/(r3**2+cv**2))**(nv/2.)
      else
         du=0.                  ! Core of Virgo
         triple=-2
         vcor=rv
         return
      endif
      dv3=du*(x3*x0+y3*y0+z3*z0)/vh/r3
c     
c  Diagnotize the existence of a Triple valued zone for the GA:
      vcor=vh+vcorcmb+vcorlg-dv3
      rs=(x0*xga+y0*yga+z0*zga)/vh ! rs is the distance for 0-effect
      drs=0.
      if(rs.gt.0.and.triple.eq.0) then
         xs=x0*rs/vh
         ys=y0*rs/vh
         zs=z0*rs/vh
         x4=xs-xga
         y4=ys-yga
         z4=zs-zga
         r4=sqrt(x4**2+y4**2+z4**2)
         r4=amax1(r4,cga/sqrt(1.*nga))
         us=vga*(r4/rga)*
     s        ((rga**2+cga**2)/(r4**2+cga**2))**((nga+1.)/2.)
         if(abs(vcor-rs).le.us) then
            rsd=rs+1.
            xs=x0*rsd/vh
            ys=y0*rsd/vh
            zs=z0*rsd/vh
            x4=xs-xga
            y4=ys-yga
            z4=zs-zga
            r4=sqrt(x4**2+y4**2+z4**2)
            usd=-vga*(r4/rga)*
     s           ((rga**2+cga**2)/(r4**2+cga**2))**((nga+1.)/2.)
            drs=usd*(x4*x0+y4*y0+z4*z0)/vh/r4
            if(drs.lt.-1.) then
               triple=1.
c         print*,' in the GA triple valued zone ',drs
            endif
         endif
      endif
c
c Diagnotize the existence of a Triple valued zone for Virgo:
      vcor=vh+vcorcmb+vcorlg-dv4
      rsv=(x0*xv+y0*yv+z0*zv)/vh
      drv=0.
      if(rsv.gt.0.and.triple.eq.0) then
         xs=x0*rsv/vh
         ys=y0*rsv/vh
         zs=z0*rsv/vh
         x3=xs-xv
         y3=ys-yv
         z3=zs-zv
         r3=sqrt(x3**2+y3**2+z3**2)
         usv=vv*
     s        ((rv**2+cv**2)/(r3**2+cv**2))**(nv/2.)
         if(abs(vcor-rsv).le.usv) then
            rsdv=rsv+1.
            xs=x0*rsdv/vh
            ys=y0*rsdv/vh
            zs=z0*rsdv/vh
            x3=xs-xv
            y3=ys-yv
            z3=zs-zv
            r3=sqrt(x3**2+y3**2+z3**2)
            usdv=-vv*
     s           ((rv**2+cv**2)/(r3**2+cv**2))**(nv/2.) 
            drsv=usdv*(x3*x0+y3*y0+z3*z0)/vh/r3
            if(drsv.lt.-1.) then
c     print*,' in VIRGO triple valued zone '
               triple=2.
            endif
         endif
      endif
c     
c  first guess for distance
c  (given this guess the model predicts the velocity perturbation
c   due to Virgo infall and Great Attractor, this peculiar
c   velocity is then added to the observed velocity to produce a
c   new guess of the distance)
      r=(vcorcmb+vcorlg-dv3-dv4)/vh
c	r=amin1(abs(vh)/7000.,1.)*vcorcmb/vh
      if(triple.eq.1)then
         r=(rs-us-10.)/vh-1.
         dv4=0.
         do while((vcorcmb+vcorlg-dv3-dv4)/vh.gt.r)
c              vest=(r+1)*vh-vcorcmb-vcorlg+dv3
c              print*,vest,vest+dv4,vh
            r=r+10./vh
            x4=x0*(1.+r)-xga
            y4=y0*(1.+r)-yga
            z4=z0*(1.+r)-zga
            r4=sqrt(x4**2+y4**2+z4**2)
            if(r4.ge.500.)then
               du=-vga*(r4/rga)*
     s              ((rga**2+cga**2)/(r4**2+cga**2))**((nga+1.)/2.)
            else
               du=0.
            endif
            dv4=du*(x4*x0+y4*y0+z4*z0)/vh/r4
         enddo
      elseif(triple.eq.2)then
         r=(rsv-usv)/vh-1.
      endif
c
c  iterations for estimation of peculiar streaming
      do k=1,10
         x2=x0*(1.+r)
         y2=y0*(1.+r)
         z2=z0*(1.+r)
         x3=x2-xv
         y3=y2-yv
         z3=z2-zv
         r3=sqrt(x3**2+y3**2+z3**2)
         if(r3.ge.cv)then
            du=-vv*((rv**2+cv**2)/(r3**2+cv**2))**(nv/2.)
         else
            du=0.
            vcor=rv
            triple=-2
            return
         endif
         dv3=du*(x3*x0+y3*y0+z3*z0)/vh/r3
         x4=x2-xga
         y4=y2-yga
         z4=z2-zga
         r4=sqrt(x4**2+y4**2+z4**2)
         if(r4.ge.500.)then
            du=-vga*(r4/rga)*
     s           ((rga**2+cga**2)/(r4**2+cga**2))**((nga+1.)/2.)
         else
            du=0.
            vcor=rga
            triple=-1
            return
         endif
         dv4=du*(x4*x0+y4*y0+z4*z0)/vh/r4
         r=(vcorcmb+vcorlg-dv3-dv4)/vh     
      enddo
c     
      vcor=vh+vcorcmb+vcorlg-dv3-dv4
c     
c       if(triple.ne.0)then
c       print'(a,4f9.2)','HUBBLE test: cmb lg virgo GA',
c    s       vcorcmb,vcorlg,-dv3,-dv4
c       endif
c       print*, 'OUTPUT: ', l, b, vh, vcor, triple
c     
      return
      end



