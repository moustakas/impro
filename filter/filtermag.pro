;+
; NAME:
;   filtermag
;
; PURPOSE:
;   Compute the AB mag of a spectrum in a given filter
;
; CALLING SEQUENCE:
;   mag = filtermag(wave, flux_flam, filter_wave, filter_trans)
;
; INPUTS:
;   wave         - wavelength array of spectrum in units of angstroms
;   flux_flam    - flux array of spectrum in units of ergs/s/cm^2/A
;   filter_wave  - wavelength array of filter throughput
;   filter_trans - transmission of filter at each wavelength 
;
; OUTPUTS:
;   AB magnitude of the input spectrum in the filter
;
; PROCEDURES CALLED:
;   linterp
;
; COMMENTS:
;   Spectrum can have an irregular wavelength grid.
;
; BUGS:
;   Not sure what happens if the wavelength range of the spectrum is smaller
;   than that of the filter -- should check.
;-
;---------------------------------------------------------------------------------

function filtermag, wave, flux_flam, filter_wave, filter_trans

   ; Determine width of each pixel in log-lambda units
   nx = n_elements(wave)
   pix = findgen(nx)
   logwave = alog10(wave)
   logdiff = abs(logwave[1:*] - logwave)  ; diff from 1 pixel center to the next
   linterp, pix[0:nx-2] + 0.5, logdiff, pix, pixwidth 

   ; Interpolate filter transmission to match the wavelength array of 
   ; the input spectrum 
   linterp, filter_wave, filter_trans, wave, filter_img
   
   ; set to zero outside of measured bounds
   bad = where(wave lt min(filter_wave) or wave gt max(filter_wave))
   if bad[0] ne -1 then filter_img[bad] = 0

   ; Convert flux from units of erg/s/cm^2/A to erg/s/cm^2/Hz 
   flux_fnu = flux_flam * wave * wave / 3.0e18 

   ; Add up the flux transmitted by the filter
   filter_norm = total(filter_img * pixwidth)
   filter_flux = total(flux_fnu * filter_img * pixwidth) / filter_norm

   ;**********************
   ; NOTE!  We are really summing f_nu * dnu / nu  * filter_trans
   ; The trick is that dlog-lambda = dnu / nu / ln(10)
   ; The reason that each pixel is divided by nu is that we want to integrate
   ; in units of photons rather than units of energy

   filter_abmag = -2.5 * alog10(filter_flux) - 48.6 

   return, filter_abmag
end

