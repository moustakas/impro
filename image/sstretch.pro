function sstretch, im, npix=npix, silent=silent
;+
; NAME:
;	SSTRETCH (that is, "smart stretch")
;
; PURPOSE:
;	To return the median, standard deviation, minimum and maximum
;	of an image (but see COMMENTS below) for computing a general
;	image stretch based on a subset of the pixels.  
;
; INPUTS:
;	im - two-dimensional image or an array
;
; OPTIONAL INPUTS:
;	npix - number of image pixels to use for computing the stretch
;              (default = 500)
;	
; KEYWORD PARAMETERS:
;	silent - suppress printing image statistics to the screen
;
; OUTPUTS:
;	sstats - structure containing the calculated stretch values
;
; COMMENTS:
;	This routine can also return statistics of a one- or
;       multi-dimensional array.  This routine was inspired by a
;       section of Marc Buie's CCDPHOT
;       (http://www.lowell.edu/users/buie/idl/idl.html#categ12). 
;
; EXAMPLE:
;	Given an image, im, try the following:
;
;	IDL> iminfo = sstretch(im)
;	IDL> tv, bytscl(im,min=iminfo.min,max=iminfo.max)
;
; PROCEDURES USED:
;	STRN(), STDEV()
;
; MODIFICATION HISTORY:
;	John Moustakas, 2001 April 28, U of A
;-

    npixtot = n_elements(im)

    if not keyword_set(npix) then npix = (500L < npixtot)
    
; choose npix random pixels from the image and sort the pixel values

    pix = floor(randomu(seed,npix)*(npixtot-1L))
    sim = (im[pix])[sort(im[pix])]

; compute the image statistics.  define a "one-sigma" interval as +/-20%
    
    x1 = 0.2*npix
    x2 = 0.8*npix

    med = sim[npix/2L]
    stdv = stdev(sim[x1:x2])

    minim = med - 3.0*stdv ; image minimum (3-sigma)
    maxim = med + 5.0*stdv ; image maximum (5-sigma)

    if (minim eq maxim) then maxim = minim + 1.0

    if not keyword_set(silent) then begin
       
       print
       print, 'Median : '+strn(med,format='(G0.0)')
       print, 'Sigma  : '+strn(stdv,format='(G0.0)')
       print, 'Minimum: '+strn(minim,format='(G0.0)')
       print, 'Maximum: '+strn(maxim,format='(G0.0)')
       print
       
    endif

; pack everything into a structure and return
    
    sstats = {median: float(med), sigma: float(stdv), $
              min: float(minim), max: float(maxim)}
    
return, sstats
end
