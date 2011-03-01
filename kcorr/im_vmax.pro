;+
; NAME:
;   IM_VMAX()
;
; PURPOSE:
;   Compute Vmax given the magnitude and redshift limits of a sample.
;
; INPUTS: 
;   redshift - galaxy redshift [NGAL]
;   mag - apparent magnitude in the selection band [NGAL]
;   coeffs - coefficients calculated by KCORRECT [5,NGAL]
;   bright - bright magnitude cut-off of the survey
;   faint - faint magnitude cut-off of the survey
;   filter - filter function name corresponding to MAG [NGAL]
;
; OPTIONAL INPUTS: 
;   vname - basis set used to calculate K-corrections
;
; KEYWORD PARAMETERS: 
;
; OUTPUTS: 
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Feb 07, UCSD
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

pro im_vmax, coeffs

    code incomplete!
    
    
    select_filter = 'ndwfs_I.par'
    aa = read_ages(/ancillary)

    cut = where((aa.i_obs gt 15.0) and (aa.i_obs lt 20.0) and $
      (aa.z gt 0.05) and (aa.z lt 0.75) and (aa.mass gt 0.0) and $
      (vv.kcorr_flag eq 1))
    aa = aa[cut]

    vname = 'default.nolines'
    coeffs = aa.coeffs

    kcorrect, dummaggies, dummaggies_ivar, 0.0, dumk, band_shift=band_shift, $
      filterlist=select_filter, rmatrix=rmatrix, zvals=zvals, coeffs=coeffs[*,0],$
      vname=vname, /silent

    for ii = 0, n_elements(aa)-1 do begin
       k_reconstruct_maggies, coeffs[*,ii], zvals, recm, $
         band_shift=band_shift, zvals=zvals, rmatrix=rmatrix, $
         filterlist=select_filter

    endfor

    
    
return
end
    
