ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c PROGRAMMA CHE CALCOLA LA DISTANZA DELLE GALASSIE ISOLATE DEL GPBG.
c LE DISTANZE, CON I DATI DELLE GALASSIE, SONO RIPORTATE IN ''daticampo''
c IN CUI LE GALASSIE A TV SONO CARATTERIZZATE DA UN *.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      program prova
      implicit none
      integer nt,tab(8600),i,nnnn
      character*13 nome,nomeg(8600),t*2,tg(8600)*2,nomefile*30
      real ll,bb,czg,dout(3),lg(8600),bg(8600),czgg(8600),appmag,vrot,A(4)
      real pi,l(4),b(4),r(4),rc(4),vmb,lmb,bmb,aa,appmagg(8600),vrotg(8600)
      real appmagcg(8600),appmagc,O
      common/parametri/pi,l,b,r,rc,vmb,lmb,bmb
      common/vista/ll,bb,czg,aa
      common/out/dout
      common/normal/A,O


      pi=2.*acos(0.)



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Model Parameters: EDIT THIS SECTION OF THE CODE UDSING THE  
c                               BEST-FITTING PARAMETERS FOR THE MULTIATTRACTOR MODEL 
c                               QUOTED IN TABLE 4 of Marinoni et al ApJ, 1998, 505, 484 

c  As an example I use here the parameters obtained fitting peculiar velocity field models
c  to the data Mark III* (i.e column 9 of table 4) 

c Coordinate della Vergine ( Not fit parameters in our model)

      r(1) = 1350.
      l(1) = 284.*pi/180.
      b(1) = 74.*pi/180.

c Coordinate del Grande Attrattore 

      r(2) = 4200.
      l(2) = 300.*pi/180.
      b(2) = 25.*pi/180.

c Coordinate di Perseus-Pisces ( Not fit parameters in our model)

      r(3)=5000.
      l(3)=120.*pi/180.
      b(3)=-30.*pi/180.

c Coordinate del Superammmasso di  Shapley 
c        (r(4) is  not a fit parameter in our model)
      r(4)=14000.      
      l(4)=308.*pi/180.
      b(4)=2.*pi/180.

c Parametri del dipolo del CMB ( Not fit parameters in our model)

      vmb=627.
      lmb=277.*pi/180.
      bmb=30.*pi/180.

c NORMALIZATION

c     Virgo
      A(1)=4.0 

c     Great Attractor
      A(2)=1.5

c     Perseus Pisces
      A(3)=3.4

c     Shapley
      A(4)=4.6


c Omega matter
      O=0.6


c CORE RADIUS

c     Virgo
      rc(1)=480

c    Great Attractor
      rc(2)=1820.

c     Perseus Pisces
      rc(3)=1410.

c     Shapley
      rc(4)=2930.
 

c Omega matter
      O=0.6


       nomefile='inputdata'
       nnnn=2810


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c      print*,'Inserisc il file contenente i seguenti dati:'
c      print*,' '
c      print*,'nomegal, coord l, coord b, app_mag, app_magcorretta, vrot,'
c      print*,'tipomorf, czlg'
c      print*,' '
c      print*,'con il seguente formato'
c      print*,'1x,A13,2x,f6.2,1x,f6.2,3x,f5.2,1x,f5.2,3x,f5.3,3x,a2,3x,f5.0'
c      print*,'************************************************************'

c      print*,'input  file'
c      read*,nomefile
c      print*,'# of galaxies in inputfile?'
c      read*, nnnn

		

c     LEGGE COORDINATE E REDSHIFT DELLE GALASSIE
    
      open(10,file=nomefile,status='old')
      open(20,file='outputdata',status='unknown')
      open(30,file='tv_output',status='unknown')
      nt=0
      aa=0.

c      do i=1,3381
c         read(15,100) b1(i),b2(i),vr(i)
c      enddo

      write(20,*),'       nome            l       b      mag   magc  vrot   czlg  tm   d'
      do i=1,nnnn

c         read(10,100,end=233) nome,ll,bb,appmag,appmagc,vrot,t,czg
          read(10,100,end=233) ll,bb,czg
	
         lg(i)=ll
         bg(i)=bb
         czgg(i)=czg
         appmagg(i)=0.
         appmagcg(i)=0.
         vrotg(i)=0.
         nomeg(i)=' '
         tg(i)='1'

c     trasforma l e b in radianti
        
         ll=ll*pi/180.
         bb=bb*pi/180.

c    chiama la subroutine che ,invertendo la curva di C. trova le dist.     
       
         call dist


c    calcola il n di galassie nella regione a triplo valore

         if (aa.eq.1.) then
          write(20,400) i,nomeg(i),lg(i),bg(i),appmagg(i),appmagcg(i),
     $           vrotg(i),czgg(i),tg(i),dout(1),dout(2),dout(3)
          nt=nt+1
          aa=0.
          tab(nt)=i
            elseif (aa.eq.0.) then
               write(20,500)i,nomeg(i),lg(i),bg(i),appmagg(i),appmagcg(i),
     $           vrotg(i),czgg(i),tg(i),dout(1)
         endif


      enddo
 233  close(10)

c    stampa i gruppi a triplo valore

c      print*,'numero galassie nella zona a triplo valore',nt
      do i=1,nt
            write(30,160)tab(i),lg(tab(i)),bg(tab(i)),czgg(tab(i))
      enddo


      close(20)

c 100  format(1x,A13,2x,f6.2,1x,f6.2,3x,f5.2,1x,f5.2,3x,f5.3,3x,a2,3x,f5.0)
 100  format(16x,f6.2,1x,f6.2,30x,f5.0)
 160  format(1x,'gal. n',i4,1x,'cord.',f6.2,1x,f6.2,1x,'czlg',f5.0)
 170  format(1x,f6.2,1x,f6.2)
 400  format('*',i4,1x,a13,3x,f6.2,1x,f6.2,3x,f5.2,1x,f5.2,1x,f5.3,3x,
     $     f5.0,1x,a2,3x,3(f5.0,1x))
 500  format(1x,i4,1x,a13,3x,f6.2,1x,f6.2,3x,f5.2,1x,f5.2,1x,f5.3,3x,
     $     f5.0,1x,a2,3x,f5.0)
      end
     


      subroutine dist

c inverte date l b la curva di C. e a seconda del valore del r-s del gruppo
c da 1 o 3 valori della distanza

      integer nn
      parameter(nn=201)
      real deltad
      parameter (deltad=7000./float(nn-1))
      real d(nn),z(nn),z1(nn),z2(nn),z3(nn),cz,aa
      real d1(nn),d2(nn),d3(nn),ind(2),czg,ll,bb
      real y21(nn),y22(nn),y23(nn),de(nn),dout(3)
      integer n(3)
      common/vista/ll,bb,czg,aa
      common/out/dout

      do i=1,nn
         d(i)=300.+(i-1)*deltad
         z(i)=cz(d(i))
c         print*,d(i),z(i)
         if (i.ne.1.and.i.ne.2) then
            de(i-1)=.5*(z(i)-z(i-2))
         endif
      enddo
      de(1)=z(2)-z(1)
      de(nn)=z(nn)-z(nn-1)
      ind(1)=0.
      ind(2)=0.
      n(1)=0
      n(2)=0
      n(3)=0
      do i=1,nn-1
         if ((de(i).ge.0).and.(de(i+1).le.0)) then
            ind(1)=max(z(i+1),z(i))
            n(1)=(i+1)
         elseif ((de(i).le.0).and.(de(i+1).ge.0)) then
            ind(2)=min(z(i+1),z(i))
            n(2)=(i+1)-n(1)
            goto 200
         endif
      enddo
 200  n(3)=nn-n(1)-n(2)


      
      if (n(3).eq.nn) then
         call spline(z,d,nn,1.E30,1.E30,y21)
         call splint(z,d,y21,nn,czg,dout(1))
c         write(*,500) dout(1)
      elseif (n(3).lt.nn) then
         
         do i=1,n(1)
            z1(i)=z(i)
            d1(i)=d(i)
         enddo
         if ((z1(n(1))-z1(n(1)-1)).le..5) then
            n(1)=n(1)-1
         endif

         
         do i=1,n(2)
            z2(i)=z(i+n(1))
            d2(i)=d(i+n(1))
         enddo
         call inv(n(2),z2)
         call inv(n(2),d2)
         

         do i=1,n(3)
            z3(i)=z(i+n(1)+n(2))
            d3(i)=d(i+n(1)+n(2))
         enddo
         
         call spline(z1,d1,n(1),1.E30,1.E30,y21)
         call spline(z2,d2,n(2),1.E30,1.E30,y22)
         call spline(z3,d3,n(3),1.E30,1.E30,y23)
      


         if ((czg.gt.0).and.(czg.lt.ind(2))) then
            call splint (z1,d1,y21,n(1),czg,dout(1))
c            write(*,400) dout(1)
         elseif ((czg.ge.ind(2)).and.(czg.le.ind(1))) then
            call splint(z1,d1,y21,n(1),czg,dout(1))
            call splint(z2,d2,y22,n(2),czg,dout(2))
            call splint(z3,d3,y23,n(3),czg,dout(3))
            aa=1.
c            print*,'ECCOLA'
            do j=1,nn
c               print*,d(j),z(j)
            enddo
            
c            write(*,300) dout(1),dout(2),dout(3)
         elseif (czg.gt.ind(1)) then
            call splint(z3,d3,y23,n(3),czg,dout(1))
c            write(*,400) dout(1)
         endif
      endif




 300  format(2x,'zona a triplo valore',2x,'d=',3(f6.1,2x))
 400  format(2x,'zona a singolo valore',2x,'d=',f6.1)
 500  format('la curva non ha t.v d=',2x,f6.1)
 550  format(1x,f6.2,1x,f6.2)
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine inv(n,x)
      real x(n),temp
      integer n,m
      if ((n/2)*2.eq.n)  then
         m=n/2
      else
         m=(n+1)/2
      endif
      do i=1,m
         temp=x(i)
         x(i)=x((n+1)-i)
         x((n+1)-i)=temp
      enddo
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      
      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      INTEGER n,NMAX
      REAL yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=500)
      INTEGER i,k
      REAL p,qn,sig,un,u(NMAX)
      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+
     *1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*
     *u(i-1))/p
11    continue
      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      SUBROUTINE splint(xa,ya,y2a,n,x,y)
      INTEGER n
      REAL x,y,xa(n),y2a(n),ya(n)
      INTEGER k,khi,klo
      REAL a,b,h
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) pause 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**
     *2)/6.
      return
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc




      function cz(d)

      real coseno(4),xr(4),rc(4),d,aa,A(4),pi
      real DD(4),V(4)
      real vpec,l(4),b(4),r(4),O,lmb,bmb,vmb,tmb
      real ll,bb,czgr,cz
      common/parametri/pi,l,b,r,rc,vmb,lmb,bmb
      common/vista/ll,bb,czgr,aa
      common/normal/A,O
 
      tmb=cos(bmb)*cos(lmb)*cos(bb)*cos(ll)
     $        +cos(bmb)*sin(lmb)*cos(bb)*             
     $                  sin(ll)+sin(bmb)*sin(bb)
c Parametri A

c      A(1)=4.0
c      A(2)=1.5
c      A(3)=3.4
c      A(4)=4.6
c Omega
c      O=0.6

c      print*,'PARAMETRI',pi,l,b,r,rc,vmb,lmb,bmb,A,O
      vpec=0.
c      print*,A
      do i=1,4
         coseno(i)=cos(b(i))*cos(l(i))*cos(bb)*cos(ll)+cos(b(i))*sin(l(i))  
     $        *cos(bb)*sin(ll)+sin(b(i))*sin(bb)           
         
c         print*,'coseno',coseno(i)
         xr(i)=(r(i)*r(i)+d*d-2.*d*r(i)*coseno(i))
c     
c     
         DD(i)=3.*A(i)*rc(i)**3./(xr(i)**(1.5))
     $        *(log((xr(i)**(0.5)+
     $        sqrt(rc(i)**2.+xr(i)))/rc(i))-
     $        xr(i)**(0.5)/sqrt(rc(i)**2.+xr(i)))

         V(i)=1./3.*(O**0.6)*xr(i)**(0.5)*
     $        (DD(i)-0.189*O**(0.08)*DD(i)**(2.)+0.081*O**(0.11)*
     $        DD(i)**(3.)-0.045*O**(0.12)*DD(i)**(4.)+0.03*O**(0.15)*
     $        DD(i)**(5.))*(r(i)*coseno(i)-d)/(xr(i)**(0.5))   

c         print*,'velocita',V(i)
         vpec=vpec+V(i)
c         print*,'vpec',vpec
      enddo
      
c###############################################################
c     
c     ECCO IL MODELLO: redshift per una singola galassia
c     generato da 4 attrattori Virgo,GA,PP e Shapley
      
      cz=d+vpec-vmb*tmb
c      print*,'funzione',cz
c     
c###############################################################
c     

      return 
      end
