function galex_read_psf, fuv=fuv
; jm10apr25ucsd - read the GALEX PSF; by default, read the NUV PSF
    psfpath = getenv('CATALOGS_DIR')+'/galex/'
    if keyword_set(fuv) then band = 'fuv' else band = 'nuv'

    psffile = psfpath+'TedPSF_'+band+'.fits'
    splog, 'Reading '+psffile
    psf = mrdfits(psffile,/silent)
    psf = psf/total(psf)
return, psf
end
    
