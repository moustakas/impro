function psfex_psf_recon,xc,yc,psffile


psf = mrdfits(psffile,1,hdr1,/silent)
psf = psf.psf_mask

xzero = sxpar(hdr1,'POLZERO1')
yzero = sxpar(hdr1,'POLZERO2')
xscale = sxpar(hdr1,'POLSCAL1')
yscale = sxpar(hdr1,'POLSCAL2')
sample_factor = sxpar(hdr1,'PSF_SAMP')
order = sxpar(hdr1,'POLDEG1')

x = (xc-xzero)/xscale
y = (yc-yzero)/yscale

psf_image = psf[*,*,0]
for n = 1,order do begin
   for ny = 0,n do begin
      nx = n-ny
      k = nx+ny*(order+1)-(ny*(ny-1))/2.
      psf_image += x^nx * y^ny * psf[*,*,k]
   endfor
endfor
   
return,psf_image
end
