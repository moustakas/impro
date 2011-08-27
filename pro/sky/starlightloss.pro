;+
; NAME:
;   starlightloss
; 
; PURPOSE: 
;   Compute light losses from slit spectra as a function of seeing, 
;     slit width, wavelength, paralactic angle, airmass.
;
; CALLING SEQUENCE: 
;   starlightloss, seeing_fwhm, slit_width, theta=, airmass=, 
;       wave0=, max_arcsec=, pixels_per_arcsec=, noplot=
;
; INPUTS:
;   seeing_fwhm  - FWHM of seeing in arcseconds at the wavelength of interest
;   slit_width   - Width of slit in arcseconds
;
; OPTIONAL KEYWORDS:
;
;   theta        - Paralactic angle  (angle of slit relative to the horizon)
;                  Defaults to 90 degrees  (perpendicual to the horizon)
;   airmass      - airmass of observation, defaults to 1.2 
;   wave0        - central wavelength in Angstroms, defaults to 5500 A.
;                  Atmospheric dispersion is computed relative to 5500 A.
;   max_arcsec   - maximum distance in arcseconds over which the calculation
;                  is computed.  Defaults to 10 arcsec.
;   pixels_per_arcsec - pixels per arcsecond, defaults to 20
;   noplot       - toggles plotting, defalt is to plot
;
;  OUTPUTS:
;    Returns fraction of light lost from the slit
;  
;  EXAMPLES:  
;    lost = starlightloss(2.0, 4.5, theta = 70, airmass=1.3, wave0 = 3800)
;
;  BUGS: 
;    Probably a slight dependence on the resolution (arcsec/pix)
;
;  PROCEDURES CALLED:
;    compute_refraction
;
;  REVISION HISTORY:
;    20-April-2003  Written by C. Tremonti, Steward Observatory
;-

function starlightloss, seeing_fwhm, slit_width, theta = theta, $
         airmass = airmass, wave0 = wave0, max_arcsec = max_arcsec, $
         pixels_per_arcsec = pixels_per_arcsec, noplot = noplot, $
         silent = silent

;-------------
; Set defaults
;-------------

if not keyword_set(theta) then theta = 0
if not keyword_set(airmass) then airmass = 1.2
if not keyword_set(wave0) then wave0 = 5500
if not keyword_set(max_arcsec) then max_arcsec = 10
if not keyword_set(pixels_per_arcsec) then pixels_per_arcsec = 20

theta = double(theta)

;---------------------
; Set up in 2d arrays
;---------------------

r_arcsec = findgen(max_arcsec * pixels_per_arcsec) / pixels_per_arcsec

npix = long(n_elements(r_arcsec) * 2 - 1)
xaxis = [-1 * reverse(r_arcsec), r_arcsec[1:*]]
yaxis = [-1 * reverse(r_arcsec), r_arcsec[1:*]]

x = fltarr(npix, npix) + rebin(xaxis, npix, npix)
y = fltarr(npix, npix) + transpose(rebin(yaxis, npix, npix))

;-----------------------------
; Compute displacement due to atmospheric dispersion
;-----------------------------

displacement = compute_refraction(airmass, refwave = 5500, lambda = lambda) 
linterp, lambda, displacement, wave0, disp_wave0

xshift = x + disp_wave0
r = sqrt(xshift^2 + y^2)

;------------------
; Define star profile
;-----------------

sigma = seeing_fwhm / (2.0 * sqrt(2.0 * alog(2.0)))
star_img = 1.0 * exp (-r^2/(2.0 * sigma^2))

;--------------
; Define slit aperture
;--------------

vertical = (theta mod 90 eq 0) and (theta mod 180 ne 0)

if vertical then begin
   slit_aperture = where(x gt -slit_width/2 and x lt slit_width/2)
endif else begin
  ymax = tan(theta * !DTOR)*x + slit_width/2 / abs(cos(theta * !DTOR))
  ymin = tan(theta * !DTOR)*x - slit_width/2 / abs(cos(theta * !DTOR))
  slit_aperture = where(y gt ymin and y lt ymax) 
endelse

inslit = total(star_img[slit_aperture]) / total(star_img)

;-------------------------------------------------------------------------------
; Plot star image + slit
;-------------------------------------------------------------------------------

if not keyword_set(noplot) then begin

  sz = size(star_img)
  if !D.NAME eq 'PS' then begin
     smdim = !D.Y_SIZE < !D.X_VSIZE 
     sz = [2, 0.8 * smdim, 0.8 * smdim]
  endif else begin
     window, 0, xsize = sz[1] * 1.2, ysize = sz[2] * 1.2
  endelse
  xi = 0.10 * !D.X_VSIZE
  yi = 0.10 * !D.Y_VSIZE

  plot_img = fltarr(npix, npix) + 0.45
  plot_img[slit_aperture] = star_img[slit_aperture] 
  if !D.NAME eq 'PS' then $
    tv, bytscl(plot_img, top = 220), xi, yi, xsize = sz[1], ysize = sz[2] $
  else tv, bytscl(plot_img, top = 220), xi, yi

  plot, xaxis, yaxis, /xstyle, /ystyle, /nodata, /noerase, $
        position = [xi, yi, xi + sz[1] - 1, yi + sz[2] - 1], /device, $
        xtitle = 'Arcseconds', ytitle = 'Arcseconds' 
  oplot, [-100, 100], [0, 0]
  oplot, [0, 0], [-100, 100]
endif

;----------------
; Print results
;----------------

if not keyword_set(silent) then $
print, 'Percent of light lost = ' + $
        string((1 - inslit) * 100, format = '(F6.2)') + '%'

return, (1 - inslit)

end
