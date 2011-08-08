c     program catalog.for
c     this program will take as input a data file which contains the id#, 
c     ra, dec, apparent mag, and redshift. 

      double precision ras, decs, mv, z 
      double precision hourra, degra, degdec, radra, raddec
      character id*9
      integer lines, n, rah, ram, decd, decm

      open(4,file='hinz0402.cat.in')
      open(2,file='hinz0402.cat')
   
c     the number of lines in the input file
      lines = 12

      rah = 0
      ram = 0
      ras = 0
      decd = 0
      decm = 0
      decs = 0
      mv = 0
      z = 0
      PI = 3.14159265359

      do n=1, lines
         
c     the input file must have all spaces filled in the id designation and 
c     all colons must br removed from the ra and dec.  The magnitudes must 
c     be positive.
         
c     3       read (4,50,err=3,end=30)id,rah,ram,ras,decd,decm,decs
         
c     50      format(a9,x,i2,x,i2,x,f5.2,x,i3,x,i2,x,f4.1)
         
 3       read (4,50,err=3,end = 30)id,rah,ram,ras,decd,decm,decs,mv
         
 50      format(a9,x,i2,x,i2,x,f4.1,x,i3,x,i2,x,i2,x,f5.2)
         
c     print*, 'goto',n 
c     This converts the ra and dec into radians
         
         
         if (decd .lt. 0) then
            decm = -decm
            decs = -decs
         endif
         
         hourra = rah+(ram+ras/60.)/60.
         degra = hourra * 360/24
         degdec = decd+(decm+decs/60.)/60.
         radra = degra * (PI)/180
         raddec = degdec * (PI)/180
         mv=mv*100.
c     print*,mv
         
         if (raddec .lt. 0) then
            write(2,20)n,int(mv),radra,raddec,id
         elseif (raddec .ge. 0) then
            write(2,25)n,int(mv),radra,raddec,id
         endif

c     if (raddec .lt. 0) then
c     write(2,20)n,radra,raddec,id
c     elseif (raddec .ge. 0) then
c     write(2,25)n,radra,raddec,id
c     endif

      enddo
      continue
      
 20   format(x,i3,x,i4,x,f12.10,x,f13.10,6x,'0',6X, 
     *     '0',x,A9,56X,'2000.00   ')
      
 25   format(x,i3,x,i4,x,f12.10,x,,'+',f12.10,6x,'0',6X, 
     *     '0',x,A9,56X,'2000.00   ')
 
 30   close(unit = 4)
      close(unit = 2)


      end
