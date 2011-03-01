      subroutine csp_convol(sfr,flux,cspf,npix,nage)

      integer npix, nage
      integer i, j, k
      real*4 sfr(npix,nage)
      real*4 cspf(npix,nage) 
      real*4 sfr(npix,nage)
      real*4 sarr(npix,nage)

      do j = 2, nage
         do k = 1, j
            do l = 1, j
               sarr(l,k) = sfr(j-k)*flux(l,k)
            enddo
         enddo
         do i = 1, npix
            cspf(i,j) = trapez(age(1:j),sarr(i,
         enddo

      enddo
      
      return
      end

