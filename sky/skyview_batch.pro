;+
; NAME:
;	SKYVIEW_BATCH
;
; PURPOSE:
;	Retrieve images and data from SKYVIEW.
;
; CALLING SEQUENCE:
;
; INPUTS:
;	filename - output FITS file name for the retrieved image
;	vcoord   - [name/value] coordinate position or name of the
;                  object (see below for format)
;
; OPTIONAL INPUTS:
;	survey   - [name] survey name [default: Digitized Sky Survey]
;	perlpath - full path name to the perl scripts SKVBATCH and
;                  WEBQUERY [default: ${SKYVIEW_DIR}]
;	datapath - full path name to output the data [default: PWD]
;	scoord   - coordinate system [Equatorial,Galactic,Ecliptic]
;	equinx   - equinox
;	maproj   - map projection
;	sfactr   - image size [degrees]
;	pixelx   - x image size [pixels]
;	pixely   - y image size [pixels]
;	nameres  - name resolver method
;	
; KEYWORD PARAMETERS:
;	silent   - do not print web query message
;
; OUTPUTS:
;	
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;	See http://skyview.gsfc.nasa.gov/docs/batchpage.html for much
;	more detailed documentation.  This IDL script does not have
;	the full flexibility of the PERL script.
;
;	The PERL scripts called by this routine should be located in a
;	directory defined in your .idlenv or .cshrc file (e.g., setenv
;	SKYVIEW_PATH ${HOME}/idl/impro/skyview).  
;
; EXAMPLE:
;
; PROCEDURES USED:
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2002 Feb 18, U of A
;	jm02jun14uofa - vectorized and added more keywords
;-

pro skyview_batch, filename, vcoord, survey=survey, perlpath=perlpath, datapath=datapath, $
  imtype=imtype, scoord=scoord, equinx=equinx, maproj=maproj, sfactr=sfactr, $
  pixelx=pixelx, pixely=pixely, nameres=nameres

    nfiles = n_elements(filename)
    nvcoord = n_elements(vcoord)

    if n_params() ne 2L then begin
       print, 'Syntax - skyview_batch, filename, vcoord, [survey=, perlpath=], $'
       print, '   [datapath=, imtype=, scoord=, equinx=, maproj=, sfactr=], $'
       print, '   [pixelx=, pixely=]'
       return
    endif

    if nfiles ne nvcoord then begin
       print, 'FILENAME and VCOORD do not have the same number of elements.'
       return
    endif
    
    if not keyword_set(perlpath) then perlpath = filepath('',root_dir=getenv('SKYVIEW_DIR'))
    if not keyword_set(datapath) then begin
       spawn, ['pwd'], datapath & datapath = datapath[0]+'/'
    endif
    if n_elements(survey) eq 0L then survey = 'Digitized Sky Survey'
    if n_elements(nameres) eq 0L then nameres = 'NED/SIMBAD'

    nsfactr = n_elements(sfactr)
    if (nsfactr eq 0L) then sfactr = replicate(5./60,nfiles) else begin
       if nsfactr eq 1L then sfactr = replicate(sfactr,nfiles) else $
         if nsfactr ne nfiles then begin
          print, 'SFACTR and VCOORD do not have the same number of elements.'
          return
       endif
    endelse

    if nfiles gt 1L then begin

       for k = 0L, nfiles-1L do begin
       
          skyview_batch, filename[k], vcoord[k], survey=survey, perlpath=perlpath, $
            datapath=datapath, imtype=imtype, scoord=scoord, sfactr=sfactr[k], $
            pixelx=pixelx, pixely=pixely, nameres=nameres

       endfor

    endif
       
;   command = "./skvbatch file='"+datapath+filename+"' VCOORD='"+vcoord+"' SURVEY='"+survey+"'"
    command = "./skvbatch file='"+datapath+filename+"' VCOORD='"+vcoord+"' SURVEY='"+survey+"' PXLCNT='Yes'"

    if keyword_set(imtype) then command = command+" RETURN='"+strcompress(string(imtype),/remove)+"'"
    if keyword_set(scoord) then command = command+" SCOORD='"+strcompress(string(scoord),/remove)+"'"
    if keyword_set(sfactr) then command = command+" SFACTR='"+strcompress(string(sfactr),/remove)+"'"
    if keyword_set(pixelx) then command = command+" PIXELX='"+strcompress(string(pixelx),/remove)+"'"
    if keyword_set(pixely) then command = command+" PIXELY='"+strcompress(string(pixely),/remove)+"'"
    if keyword_set(nameres) then command = command+" NAMERES='"+strcompress(string(nameres),/remove)+"'"
    
    pushd, perlpath
    if not keyword_set(silent) then print, command
    spawn, [command], /sh
    popd

return
end
