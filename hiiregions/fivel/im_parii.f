c+
c NAME:
c     IM_PARII()
c
c PURPOSE:
c     Calculate temperature-dependent atomic parameters. 
c
c CALLING SEQUENCE:
c     call im_parii(ionum,tem,c)
c
c INPUTS:
c     ionum - integer ion number
c     tem   - electron temperature [K]
c
c OUTPUTS:
c     c - collision strengths [5,5]
c
c INTERNAL VARIABLES:
c     cj   - coefficients for quadratic fit in T4 = T/10000 [3]
c     i, j - do loop integers
c     t4   - tem/1E4
c     logt - alog10(tem)
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

      subroutine im_parii(ionum,tem,c)

      implicit none

      integer*4 ionum             ! input
      real*8 tem                  ! input
      
      integer*4 i, j              ! internal
      real*8 t4, logt, cj(3)      ! internal

      real*8 c(5,5)               ! output

      include 'tty.dec'

      t4=1.0e-4*tem
      cj(1)=0.0
      cj(2)=0.0
      cj(3)=0.0

c     Zero arrays

      do 100 i=1,5
         do 200 j=1,5
            c(i,j)=0.0
 200     continue
 100  continue


c     For new ions, appropriate parameters must be included in IF block

      if ((ionum .eq. 6010) .or. (ionum .eq. 6011)
     1     .or. (ionum .eq. 6012)) then

c     C II ion parameters   (Blum and Pradhan 1992, ApJS 80, 425)

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
c     C II ion parameters (old)
c     
c     c(2,1)=2.90
c     c(3,1)=0.276
c     c(4,1)=0.409
c     c(5,1)=0.260
c     c(3,2)=0.197
c     c(4,2)=0.536
c     c(5,2)=1.16
c     c(4,3)=0.911
c     c(5,3)=0.920
c     c(5,4)=2.16
c     
      elseif ((ionum .eq. 6020) .or. (ionum .eq. 6021)
     1        .or. (ionum .eq. 6022)) then
c     
c     C III ion parameters
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
     1        .or. (ionum .eq. 7002)) then
c     
c     N I ion parameters
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
     1        .or. (ionum .eq. 7012)) then
c     
c     N II ion parameters
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
     1        .or. (ionum .eq. 7022)) then
c     
c     N III ion parameters    (Blum and Pradhan 1992, ApJS 80, 425, T=15,000 K)
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
c     N III ion parameters   (old)
c     
c     c(2,1)=0.701
c     c(3,1)=0.0952
c     c(4,1)=0.139
c     c(5,1)=0.080
c     c(3,2)=0.0616
c     c(4,2)=0.175
c     c(5,2)=0.390
c     c(4,3)=0.695
c     c(5,3)=0.397
c     c(5,4)=1.26
c     
      else if((ionum .eq. 8000) .or. (ionum .eq. 8001)
     1        .or. (ionum .eq. 8002)) then
c     
c     O I ion parameters
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
     1        .or. (ionum .eq. 8012)) then
c     
c     O II ion parameters
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
     1        .or. (ionum .eq. 8022)) then
c     
c     O III ion parameters
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
     1        .or. (ionum .eq. 10022)) then
c     
c     Ne III ion parameters
c     Butler & Zeippen 1994, A&AS 108, 1
         cj(1)=1.357
c     +t4*(0.07967-t4*0.03282)
         cj(2)=0.151
c     +t4*(0.05331-t4*0.01487)
         c(2,1)=0.244
c     +t4*(0.1544-t4*0.05638)
         c(3,1)=0.306
c     +t4*(0.02765-t4*0.012639)
         c(4,1)=cj(1)*5/9
         c(5,1)=cj(2)*5/9
         c(3,2)=0.348
c     +t4*(0.06120-t4*0.02096)
         c(4,2)=cj(1)/3
         c(5,2)=cj(2)/3
         c(4,3)=cj(1)/9
         c(5,3)=cj(2)/9
         c(5,4)=0.269
c     +t4*(0.05962-t4*0.00862)
c     old values
c     cj(1)=1.6028+t4*(0.07967-t4*0.03282)
c     cj(2)=0.12956+t4*(0.05331-t4*0.01487)
c     c(2,1)=1.0314+t4*(0.1544-t4*0.05638)
c     c(3,1)=0.2910+t4*(0.02765-t4*0.012639)
c     c(4,1)=cj(1)*5/9
c     c(5,1)=cj(2)*5/9
c     c(3,2)=0.3080+t4*(0.06120-t4*0.02096)
c     c(4,2)=cj(1)/3
c     c(5,2)=cj(2)/3
c     c(4,3)=cj(1)/9
c     c(5,3)=cj(2)/9
c     c(5,4)=0.1739+t4*(0.05962-t4*0.00862)
c     
      else if((ionum .eq. 10030) .or. (ionum .eq. 10031)
     1        .or. (ionum .eq. 10032)) then
c     
c     Ne IV ion parameters
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
     1        .or. (ionum .eq. 10042)) then
c     
c     Ne V ion parameters
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
     1        .or. (ionum .eq. 14022)) then
c     
c     Si III ion parameters
c     
         logt=dlog10(tem)
         cj(1)=539.070-457.121*logt+148.148*logt**2-21.4661*logt**3
     .        +1.16536*logt**4
         cj(2)=575.809-522.592*logt+178.762*logt**2-26.8185*logt**3
     .        +1.48065*logt**4
         c(2,1)=cj(1)/9
         c(3,1)=cj(1)/3
         c(4,1)=cj(1)*5/9
         c(5,1)=116.465-94.7015*logt+30.0337*logt**2-4.28699*logt**3
     .        +0.241284*logt**4
         c(3,2)=21.1105-14.8546*logt+3.71890*logt**2-0.287946*logt**3
     .        -3.72024e-3*logt**4
         c(4,2)=-43.3029+38.2889*logt-11.6750*logt**2+1.61607*logt**3
     .        -8.92857e-2*logt**4
         c(5,2)=cj(2)/9
         c(4,3)=-41.9423+38.0232*logt-10.5783*logt**2+1.47321*logt**3
     .        -9.67262e-2*logt**4
         c(5,3)=cj(2)/3
         c(5,4)=cj(2)*5/9
c     
      else if ((ionum .eq. 16010) .or. (ionum .eq. 16011)
     1        .or. (ionum .eq. 16012)) then
c     
c     S II ion parameters
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
     1        .or. (ionum .eq. 16022)) then
c     
c     S III ion parameters  
c     Tayal SS & Gupta GP 1999 ApJ 526, 544
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
      else if ((ionum .eq. 16030) .or. (ionum .eq. 16031)
     1        .or. (ionum .eq. 16032)) then
c     
c     S IV ion parameters
c     
c     cj(1)=9.903+t4*(-2.017+t4*0.59)
c     cj(2)=1.135+t4*0.0522
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
     1        .or. (ionum .eq. 17012)) then
c     
c     Cl II ion parameters
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
     1        .or. (ionum .eq. 17022)) then
c     
c     Cl III ion parameters
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
     1        .or. (ionum .eq. 17032)) then
c     
c     Cl IV ion parameters
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
     1        .or. (ionum .eq. 18022)) then
c     
c     Ar III ion parameters
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
     1        .or. (ionum .eq. 18032)) then
c     
c     Ar IV ion parameters; ZBL87 A&A 188, 251 
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
     1        .or. (ionum .eq. 18042)) then
c     
c     Ar V ion parameters
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
