pro build_igmtau_grid
; jm09feb11nyu - build a lookup table of IGM attenuation as a function
;   of wavelength and redshift
; jm09aug18ucsd - use a finer redshift and wavelength grid spacing 

    dz = 0.05
    zmin = dz
    zmax = 6.0
    zgrid = findgen((zmax-zmin)/dz+1)*dz+zmin
    nz = n_elements(zgrid)

    wavemin = 1.0
    wavemax = 1220.0*(1.0+zmax)
    dwave = 1.0
    wave = findgen((wavemax-wavemin)/dwave+1)*dwave+wavemin
    nw = n_elements(wave)

    grid = {wave: wave, zgrid: zgrid, igm: fltarr(nw,nz)}
    for iz = 0L, nz-1L do grid.igm[*,iz] = exp(-lm_igmtau(wave,zgrid[iz]))

    outfile = getenv('IMPRO_DIR')+'/dust/igmtau_grid.fits'
    im_mwrfits, grid, outfile

; example of usage:
;   mywave = findgen(6000)+2000.0 ; [3000-8000 A]
;   myz = [1.5,2.5,3.5,4.5]
;   for jj = 0L, n_elements(myz)-1L do begin
;      plot, mywave, interpolate(grid.igm,findex(grid.wave,mywave),$
;        findex(grid.zgrid,myz[jj]),/grid), xsty=3, ysty=3
;      cc = get_kbrd(1)
;   endfor
    
return
end
    
