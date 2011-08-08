;+
; NAME:
;       HB_SFR()
;
; PURPOSE:
;       Estimate the star-formation rate (SFR) from the B-band
;       luminosity and the observed H-beta luminosity according to
;       Moustakas, Kennicutt, & Tremonti (2006).
;
; INPUTS:
;       loglb  - B-band luminosity [L(B)_sun]
;       loglhb - observed (stellar-absorption corrected, but not
;                reddening-corrected) H-beta luminosity [erg/s]
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;       log - return the log-base-10 SFR and error
;
; OUTPUTS:
;       sfr - star-formation rate [M_sun/yr]
;
; OPTIONAL OUTPUTS:
;       sfr_err - uncertainty in SFR based on the standard-deviation
;                 of the empirical calibration [M_sun/yr]
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2005 Aug 26, U of A - written
;       jm06apr10uofa - also return SFR_ERR; documented
;
; Copyright (C) 2005-2006, John Moustakas
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

function hb_sfr, loglb, loglhb, sfr_err=sfr_err, log=log

    nloglb = n_elements(loglb)
    nhb = n_elements(loglhb)
    
    if (nloglb eq 0L) or (nhb eq 0L) then begin
       doc_library, 'hb_sfr'
       return, -1L
    endif

    if (nloglb ne nhb) then begin
       splog, 'LOGLB and LOGLHB must have the same number of elements!'
       return, -1L
    endif

; Table 1 from Moustakas, Kennicutt, & Tremonti (2006)
    
    hb_sfr = {loglb: 0.0, mb: 0.0, p25: 0.0, p50: 0.0, p75: 0.0, mean: 0.0, stddev: 0.0}
    hb_sfr = replicate(hb_sfr,8)

    hb_sfr.loglb  = [7.25,7.75,8.25,8.75,9.25,9.75,10.25,10.75]
    hb_sfr.mb     = [-12.6781,-13.9281,-15.1781,-16.4281,-17.6781,-18.9281,-20.1781,-21.4281]
    hb_sfr.p25    = [0.400381,0.371458,0.401395,0.349113,0.426179,0.489636,0.639466,0.819381]
    hb_sfr.p50    = [0.448798,0.397380,0.418290,0.399117,0.494880,0.547120,0.766920,0.892183]
    hb_sfr.p75    = [0.538937,0.434633,0.497325,0.450828,0.649410,0.669939,0.890293,1.031510]
    hb_sfr.mean   = [0.435313,0.403330,0.431358,0.413791,0.525946,0.583498,0.774705,0.926661]
    hb_sfr.stddev = [0.0793577,0.0418002,0.0568947,0.0800698,0.147184,0.156696,0.188477,0.189738]

;   hb_sfrpath = atlas_path(/projects)+'sfrs/'
;   hb_sfrfile = 'hb_sfr.fits'
;   hb_sfr = mrdfits(hb_sfrpath+hb_sfrfile,1,/silent)

;   sfr = interpol(hb_sfr.mean,hb_sfr.loglb,loglb) + loglhb - 41.0
    sfr     = interpol(hb_sfr.p50,hb_sfr.loglb,loglb) + loglhb - 41.0
    sfr_err = interpol(hb_sfr.stddev,hb_sfr.loglb,loglb)
;   ploterror, loglb, sfr, sfr_err, ps=3, errstyle=1

    if (not keyword_set(log)) then begin
       sfr = 10.0^sfr
       sfr_err = sfr*sfr_err*alog(10.0)
    endif

return, sfr
end
