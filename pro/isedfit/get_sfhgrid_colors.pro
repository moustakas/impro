;+
; NAME:
;   GET_SFHGRID_COLORS()
;
; PURPOSE:
;   Synthesize photometry in arbitrary bandpasses given an iSEDfit
;   grid number generated using BUILD_ISEDFIT_SFHGRID.
;
; INPUTS: 
;   sfhgrid - SFH grid number
;   filterlist - filter curves to use 
;
; OPTIONAL INPUTS - see BUILD_ISEDFIT_SFHGRID for details:
;   imf - initial mass function (default 'chab')
;   synthmodels - population synthesis models (default 'bc03')
;   redcurve - reddening curve (default 1)
;   sfhgrid_paramfile - SFH grid parameter file
;   isedfit_sfhgrid_dir - pathname indicating the location of the
;     grids 
;   band_shift - band-shifting factor to apply to the input filters
;     (default 0.0) 
;
; KEYWORD PARAMETERS: 
;
; OUTPUTS: 
;   colors - output data structure that contains the details of every
;     model and the rest-frame AB magnitudes in each of the bandpasses
;     in FILTERLIST 
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Jan 28, UCSD
;
; Copyright (C) 2010, John Moustakas
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

function get_sfhgrid_colors, sfhgrid, filterlist=filterlist, imf=imf, $
  synthmodels=synthmodels, redcurve=redcurve, sfhgrid_paramfile=sfhgrid_paramfile, $
  isedfit_sfhgrid_dir=isedfit_sfhgrid_dir, band_shift=band_shift

    if (n_elements(sfhgrid) eq 0) then begin
       doc_library, 'get_sfhgrid_colors'
       return, -1
    endif

    nfilt = n_elements(filterlist)
    if (nfilt eq 0) then begin
       splog, 'FILTERLIST input required'
       return, -1
    endif

    if (n_elements(synthmodels) eq 0) then synthmodels = 'bc03'
    if (n_elements(imf) eq 0) then imf = 'chab' ; 'salp'
    if (n_elements(redcurve) eq 0) then redcurve = 1

; read the parameter file describing each of the grids and get the
; reddening curve
;   params = read_sfhgrid_paramfile(sfhgrid,sfhgrid_paramfile=sfhgrid_paramfile)
    redcurvestring = redcurve2string(redcurve)

    if (n_elements(isedfit_sfhgrid_dir) eq 0) then isedfit_sfhgrid_dir = $
      '${ISEDFIT_SFHGRID_DIR}/'
    sfhgridstring = 'sfhgrid'+string(sfhgrid,format='(I2.2)')
    sfhgridpath = isedfit_sfhgrid_dir+sfhgridstring+$
      '/'+synthmodels+'/'+redcurvestring+'/'

    chunkinfofile = sfhgridpath+imf+'_chunkinfo.fits.gz'
    if (file_test(chunkinfofile) eq 0) then begin
       splog, 'No SFH grid information file found!'
       return, -1
    endif

; read the chunks and synthesize the photometry    
    chunkinfo = mrdfits(chunkinfofile,1,/silent)
    chunkfiles = sfhgridpath+strtrim(chunkinfo.chunkfiles,2)
    nchunk = n_elements(chunkfiles)

    for ichunk = 0L, nchunk-1L do begin
       splog, 'Reading '+chunkfiles[ichunk]
       info1 = mrdfits(chunkfiles[ichunk],1)
       nmodel = n_elements(info1)
       colors1 = struct_trimtags(info1,except=['wave','flux'])
       colors1 = struct_addtags(temporary(colors1),$
         replicate({abmag: fltarr(nfilt,n_elements(info1[0].age))},nmodel))
       
       kwave = k_lambda_to_edges(info1[0].wave)
       for ii = 0L, nmodel-1 do begin
          k_projection_table, rmatrix, info1[ii].flux, kwave, $
            /silent, zvals, filterlist, zmin=0.0, zmax=0.0, nz=1, $
            band_shift=band_shift
          colors1[ii].abmag = -2.5*alog10(transpose(reform(rmatrix)))
       endfor
       if (ichunk eq 0L) then colors = colors1 else $
         colors = [temporary(colors),colors1]
    endfor 
    
return, colors
end
