;+
; NAME: arm_plotthicken
;       
; CATEGORY: plotting
;
; PURPOSE: thicken all lines on a given plot
;
; CALLING SEQUENCE: arm_plotthicken, infile, 
;                      [outfile=, datapath=, /writeover] 
;
; INPUTS:
;   infile - string name of FITS file to thicken
;       
; OPTIONAL INPUTS:
;   outfile  - string name of outputted FITS file
;              (default='arm_plotthicken.fits') 
;   datapath - string name of path to data 
;   npix     - number of neighboring pixels to darken 
;
; KEYWORDS: 
;   writeover - do not prompt before writing over existing output file
;   negative  - interpret high/low values as dark/light
;   bw        - force all values to be black or white
;
; OUTPUTS:
; 
; OPTIONAL OUTPUTS:
;
; EXAMPLE:
;
; PROCEDURES USED:
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;    written by A.R.Marble & J.Murphy, Steward Obs., 2004 May 12
;-

pro arm_plotthicken, infile, outfile=outfile, datapath=datapath, $
       writeover=writeover, negative=negative, bw=bw, npix=npix

; defaults and error checking

    if N_ELEMENTS(datapath) eq 0L then datapath = '' else $
      if FILE_SEARCH(datapath) eq '' then begin
       PRINT & PRINT, 'Invalid pathname. Aborting...'
       return
    endif

    if FILE_SEARCH(datapath+infile) eq '' then begin
       PRINT & PRINT, 'Invalid input file name.  Aborting...'
       return
    endif
    
    if N_ELEMENTS(outfile) eq 0L then outfile = 'arm_plotthicken.fits'
    
    if not KEYWORD_SET(writeover) then begin
       
       if FILE_SEARCH(datapath+outfile) ne '' then begin
          
          PRINT & PRINT, 'The file '+datapath+outfile+' already exists.  Write over? (Y/N)'
          
          key = ''
          while key ne 'N' and key ne 'Y' do key = STRUPCASE(GET_KBRD(1))
          if key eq 'N' then begin
             PRINT, 'arm_plotthicken ABORTED...'
             return
          endif
          
       endif
       
    endif
    
    if N_ELEMENTS(npix) eq 0L then npix = 1L

; read in and process the image

    im = READFITS(datapath+infile, header)
    im = im[*,*]                ; elliminate higher dimensionality
    
    nx = n_elements(im[*,0])    ; number of columns
    ny = n_elements(im[0,*])    ; number of rows

    lo = MIN(im)                ; original lowest value
    hi = MAX(im)                ; original highest values

    im = im - lo
    im = im * 1d0 / MAX(im)

    if KEYWORD_SET(bw) then im = ROUND(im)

; identify dark pixels and darken neighbors

    indices = WHERE(im eq 0L, count)

    if count gt 0L then begin
    
       yindex = indices / FIX(nx)
       xindex = indices - yindex * nx

       for i = 0L, count - 1L do $
         im[(xindex[i] - npix) > 0 : (xindex[i] + npix) < (nx - 1), $
            (yindex[i] - npix) > 0 : (yindex[i] + npix) < (ny - 1)] = $
         im[xindex[i], yindex[i]]
       
    endif

    if KEYWORD_SET(negative) then im = 1 - im

    im = im * (hi - lo)
    im = im + lo

    WRITEFITS, datapath+outfile, im, header

end
