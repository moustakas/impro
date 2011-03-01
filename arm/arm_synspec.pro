;+
; NAME: ARM_SYNSPEC.pro
;       
; CATEGORY: astronomy
;
; PURPOSE: create a synthetic spectrum
;
; CALLING SEQUENCE: 
;   spectrum = ARM_SYNSPEC(wave, line, [cont=, ew=, fwhm=])
;
; INPUT: 
;   wave - wavelength array of spectrum
;   line - wavelengths of spectral lines
;
; OPTIONAL INPUTS:
;   cont - continuum values for spectrum (default=[1,1...])
;   ew   - equivalent widths of spectral lines (default=[1,1...]),
;          negative values indicate emission lines
;   fwhm - full width at half max of spectral lines 
;          (default=[2,2...])
;
; OUTPUTS: 
;   spec - continuum plus/minus spectral lines
;   wave - wavelength array
;
; MODIFICATION HISTORY:
;    written by A.R.Marble, Steward Obs., July 2003
;    error checking improved, ARM, 2004 June 15
;-

function ARM_SYNSPEC, wave, line, cont=cont, ew=ew, fwhm=fwhm, dispersion=dispersion

; default values & error-checking

    if N_ELEMENTS(dispersion) eq 0L then dispersion = 1.0 else $
      dispersion = FLOAT(dispersion)
    if N_ELEMENTS(line) eq 0L then begin
       MESSAGE, 'WAVE and LINE must be defined.'
       return, -1
    endif

    nl = N_ELEMENTS(line)       ; number of spectral lines

    if N_ELEMENTS(cont) eq 0L then cont = wave * 0d0 + 1d0 $
    else if N_ELEMENTS(cont) ne N_ELEMENTS(wave) then $
      MESSAGE, 'wavelength and continuum arrays differ in size'

    if N_ELEMENTS(ew) eq 0 then ew = line * 0d0 + 1d1 $
    else if N_ELEMENTS(ew) ne nl then $
      MESSAGE, 'line wavelength and ew arrays differ in size'

    if N_ELEMENTS(fwhm) eq 0 then fwhm = line * 0d0 + 2d1 $
      else if N_ELEMENTS(fwhm) ne nl then $
      MESSAGE, 'line wavelength and fwhm arrays differ in size'
    
; create synthetic spectrum

    var = (0.42466 * fwhm)^2.   ; variance
    
    spec = cont * 0d + 1d       ; unity spectrum
    
;  subtract opacity array for each absorption line

    if line[0] ne -1 then for i=0,nl-1 do spec = $
      spec - (ew[i] / sqrt(2.*!pi*var[i])) * $
      exp(-1. * ((wave-line[i])^2.) / (2.*var[i]))

    spec = spec * cont          ; multiply continuum by opacity map 

    return, spec

 end


