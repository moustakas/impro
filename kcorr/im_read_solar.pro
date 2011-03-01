;+
; NAME:
;       IM_READ_SOLAR()
;
; PURPOSE:
;       Read the solar spectrum provided by the CALSPEC database. 
;
; CALLING SEQUENCE:
;       solar = im_read_solar()
;
; INPUTS:
;       None
;
; OPTIONAL INPUTS:
;       None
;
; KEYWORD PARAMETERS:
;       None
;
; OUTPUTS:
;       solar
;          wave - wavelength [A]
;          flux - flux density [erg/s/cm2/A]
;
; OPTIONAL OUTPUTS:
;       None
;
; COMMON BLOCKS:
;       None
;
; PROCEDURES USED:
;       READCOL
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2004 Mar 21, U of A
;-

function im_read_solar

    sedpath = filepath('',root_dir=getenv('LMFILTERZ_DIR'))
    sedfile = 'calspec_solar_spectrum.fits'

    if file_test(sedpath+sedfile) then begin
       sun = mrdfits(sedpath+sedfile,1,h,/silent)
       solar = {wave: sun.wavelength, flux: sun.flux}
    endif else begin
       print, 'Solar spectrum '+sedpath+sedfile+' not found.'
       return, -1L
    endelse

; normalize

    norm = 1.0
    solar.flux = solar.flux / norm[0]

;   plot, solar.wave, solar.flux, ps=10, xsty=3, ysty=3, $
;     xrange=[3000,10000], charsize=2.0, charthick=2.0, $
;     xthick=2.0, ythick=2.0, xtitle='Wavelength', ytitle='Flux Density'
    
return, solar
end    
