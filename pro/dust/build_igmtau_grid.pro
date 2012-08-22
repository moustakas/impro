;+
; NAME:
;   BUILD_IGMTAU_GRID
;
; PURPOSE:
;   Build a lookup table of IGM attenuation as a function of
;   wavelength and redshift. 
;
; INPUTS: 
;
; OUTPUTS: 
;   Lookup table written to ${IMPRO_DIR}/etc
;
; EXAMPLES:
;   Example of how to use the output:
;
;      grid = mrdfits($IMPRO_DIR/etc/igmtau_grid.fits.gz',1)
;      mywave = findgen(6000)+2000.0 ; [3000-8000 A]
;      myz = [1.5,2.5,3.5,4.5]
;      for jj = 0L, n_elements(myz)-1L do begin
;         plot, mywave, interpolate(grid.igm,findex(grid.wave,mywave),$
;           findex(grid.zgrid,myz[jj]),/grid), xsty=3, ysty=3
;         cc = get_kbrd(1)
;      endfor
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2009 Feb 11, NYU
;   jm09aug18ucsd - use a finer redshift and wavelength grid spacing 
;
; Copyright (C) 2009, John Moustakas
; 
; This program is free software; you can redistribute it and/or modify 
; it under the terms of the GNU General Public License as published by 
; the Free Software Foundation; either version 2 of the License, or
; (at your option) any later version. 
; 
; This program is distributed in the hope that it will be useful, but 
; WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
; General Public License for more details. 
;-

pro build_igmtau_grid, clobber=clobber

    dz = 0.02 ; 0.05
    zmin = dz
    zmax = 15.0 ; 6.0
    zgrid = findgen((zmax-zmin)/dz+1)*dz+zmin
    nz = n_elements(zgrid)

    wavemin = 1.0
    wavemax = 1220.0*(1.0+zmax)
    dwave = 1.0
    wave = findgen((wavemax-wavemin)/dwave+1)*dwave+wavemin
    nw = n_elements(wave)

    grid = {wave: wave, zgrid: zgrid, igm: fltarr(nw,nz)}
    for iz = 0L, nz-1L do grid.igm[*,iz] = exp(-lm_igmtau(wave,zgrid[iz]))

    outfile = getenv('IMPRO_DIR')+'/etc/igmtau_grid.fits'
    im_mwrfits, grid, outfile, clobber=clobber

return
end
    
