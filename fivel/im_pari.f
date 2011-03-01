c+
c NAME:
c     IM_PARI()
c
c PURPOSE:
c     Initialize atomic parameters.
c
c CALLING SEQUENCE:
c     call im_pari(ionum,a,t,weight)
c
c INPUTS:
c     ionum - integer ion number
c
c OUTPUTS:
c     a      - radiative transition probabilities [5,5]
c     t      - energy-level separation from ground [1/Angstrom] [5]
c     weight - statistical weights for each level [5]
c
c INTERNAL VARIABLES:
c     i, j - do loop integers
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

      subroutine im_pari(ionum,a,t,weight)

      implicit none

      integer*4 ionum             ! input
      integer*4 i, j              ! internal

      real*8 a(5,5), t(5)         ! output
      integer*4 weight(5)         ! output

      include 'tty.dec'

c     Zero arrays

      do 100 i=1,5
         do 200 j=1,5
            a(i,j)=0.0
 200     continue
 100  continue
c     
      t(1)=0.0
      if((ionum .eq. 6010) .or. (ionum .eq. 6011)
     1     .or. (ionum .eq. 6012)) then

c     C II ion parameters

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
     1        .or. (ionum .eq. 6022)) then
c     
c     C III ion parameters
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
     1        .or. (ionum .eq. 7002)) then
c     
c     N I ion parameters
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
     1        .or. (ionum .eq. 7012)) then
c     
c     N II ion parameters
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
     1        .or. (ionum .eq. 7022)) then
c     
c     N III ion parameters
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
     1        .or. (ionum .eq. 8002)) then
c     
c     O I ion parameters
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
     1        .or. (ionum .eq. 8012)) then
c     
c     O II ion parameters
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
     1        .or. (ionum .eq. 8022)) then
c     
c     O III ion parameters
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
     1        .or. (ionum .eq. 10022)) then
c     
c     Ne III ion parameters
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
     1        .or. (ionum .eq. 10032)) then
c     
c     Ne IV ion parameters
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
     1        .or. (ionum .eq. 10042)) then
c     
c     Ne V ion parameters
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
     1        .or. (ionum .eq. 14022)) then
c     
c     Si III ion parameters
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
     1        .or. (ionum .eq. 16012)) then
c     
c     S II ion parameters
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
     1        .or. (ionum .eq. 16022)) then
c     
c     S III ion parameters
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
     1        .or. (ionum .eq. 16032)) then
c     
c     S IV ion parameters
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
     1        .or. (ionum .eq. 17012)) then
c     
c     Cl II ion parameters
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
     1        .or. (ionum .eq. 17022)) then
c     
c     Cl III ion parameters
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
     1        .or. (ionum .eq. 17032)) then
c     
c     Cl IV ion parameters
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
     1        .or. (ionum .eq. 18022)) then
c     
c     Ar III ion parameters
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
     1        .or. (ionum .eq. 18032)) then
c     
c     Ar IV ion parameters
c     Transition probabilities from Kaufman & Sugar 1986, JPCRD 15, 321
c     
c     a(2,1)=2.23e-2
c     a(3,1)=1.77e-3
c     a(4,1)=8.62e-1
c     a(5,1)=2.11
c     a(3,2)=2.3e-5
c     a(4,2)=6.03e-1
c     a(5,2)=7.89e-1
c     a(4,3)=0.119
c     a(5,3)=0.598
c     a(5,4)=4.94e-5
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
     1        .or. (ionum .eq. 18042)) then
c     
c     Ar V ion parameters
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
