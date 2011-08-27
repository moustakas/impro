;+
; NAME:
;       OII_SFR()
;
; PURPOSE:
;       Estimate the star-formation rate (SFR) from the B-band
;       luminosity and the observed [O II] luminosity according to 
;       Moustakas, Kennicutt, & Tremonti (2006).
;
; INPUTS:
;       loglb   - B-band luminosity [L(B)_sun]
;       logloii - observed (not reddening-corrected) [O II] 3727
;                 luminosity [erg/s]
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

function oii_sfr, loglb, logloii, sfr_err=sfr_err, log=log

    nloglb = n_elements(loglb)
    noii = n_elements(logloii)
    
    if (nloglb eq 0L) or (noii eq 0L) then begin
       doc_library, 'oii_sfr'
       return, -1L
    endif

    if (nloglb ne noii) then begin
       splog, 'LOGLB and LOGLOII must have the same number of elements!'
       return, -1L
    endif

; Table 2 from Moustakas, Kennicutt, & Tremonti (2006)
    
    oii_sfr = {loglb: 0.0, mb: 0.0, p25: 0.0, p50: 0.0, p75: 0.0, mean: 0.0, stddev: 0.0}
    oii_sfr = replicate(oii_sfr,7)

    oii_sfr.loglb  = [   7.75000D,  8.25000D,   8.75000D,   9.25000D,  9.75000D, 10.2500D, 10.7500D]
    oii_sfr.mb     = [  -13.9281D, -15.1781D,  -16.4281D,  -17.6781D, -18.9281D,-20.1781D,-21.4281D]
    oii_sfr.p25    = [-0.0280972D,-0.136150D, -0.139559D,-0.0455551D,0.0563614D,0.296572D,0.516810D]
    oii_sfr.p50    = [ 0.0978434D,0.0385868D,-0.0735958D, 0.0920122D, 0.201360D,0.529993D,0.762256D]
    oii_sfr.p75    = [  0.413915D, 0.291297D,  0.110720D,  0.281650D, 0.395194D,0.781435D,0.972045D]
    oii_sfr.mean   = [  0.156174D, 0.108897D,-0.0173798D,  0.136472D, 0.261247D,0.545016D,0.749037D]
    oii_sfr.stddev = [  0.239896D, 0.313019D,  0.185411D,  0.292403D, 0.347343D,0.331816D,0.285470D]

;   oii_sfrpath = atlas_path(/projects)+'sfrs/'
;   oii_sfrfile = 'oii_sfr.fits'
;   oii_sfr = mrdfits(oii_sfrpath+oii_sfrfile,1,/silent)

;   sfr = interpol(oii_sfr.mean,oii_sfr.loglb,loglb) + logloii - 41.0
    sfr     = interpol(oii_sfr.p50,oii_sfr.loglb,loglb) + logloii - 41.0
    sfr_err = interpol(oii_sfr.stddev,oii_sfr.loglb,loglb)
;   ploterror, loglb, sfr, sfr_err, ps=3, errstyle=1

    if (not keyword_set(log)) then begin
       sfr = 10.0^sfr
       sfr_err = sfr*sfr_err*alog(10.0)
    endif

return, sfr
end
