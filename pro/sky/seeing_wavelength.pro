;+
; NAME:
;   SEEING_WAVELENGTH
; PURPOSE:
;   Plot the variation of the seeing disk with wavelength.  
; KEYWORD PARAMETERS: 
;   postscript - write out a postscript file
; COMMENTS:
;   Based on notes on Adaptive Optics from L. Close.
; MODIFICATION HISTORY:
;   J. Moustakas, 2002 Feb 15, U of A
;-

pro seeing_wavelength, postscript=postscript

    minwave = 1000.0 ; [Angstrom]
    maxwave = 5E4    ; 5 micron [Angstrom]
    dwave = 10.0     ; Angstrom per pixel

    wave = findgen((maxwave-minwave)/dwave)*dwave+minwave

    r_0 = 10.0 * (wave / 0.5E4)^(6./5) ; Fried parameter [cm]

    seeing = 206265.0 * 1E-8 * wave / r_0 ; proportional to wave^(1./5) [arcsec]

    if keyword_set(postscript) then begin
       ps_open, 'seeing_wavelength', /ps_fonts
       device, /inches, /times
    endif else window, 0, xs=550, ys=550
    plot, wave, seeing, xsty=1, ysty=1, /xlog, xrange=[minwave,maxwave], $
      xtitle='Wavelength (\AA)', ytitle='Seeing (arcsec)', $
      charsize=2.0, charthick=2.0, xthick=2.0, ythick=2.0, thick=2.0, $
      title='Seeing Disk Wavelength Dependence'
    oplot, [3500,3500], !y.crange, line=2, thick=1.5
    oplot, [7000,7000], !y.crange, line=2, thick=1.5
    legend, ['Optical Window'], linestyle=2, /right, /top, box=0, $
      charsize=1.3, charthick=2.0
    if keyword_set(postscript) then ps_close

return
end
