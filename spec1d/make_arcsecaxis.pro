;+
; NAME:
;	MAKE_ARCSECAXIS()
;
; PURPOSE:
;	Generate the spatial axis in arcsec for two-dimensional
;	spectra. 
;
; CALLING SEQUENCE:
;       arcsecaxis = make_arcsecaxis(header,npix=,axis_0=,refpix=,pscale=)
;
; INPUTS:
;	header - spectrum header
;
; OUTPUTS:
;       arcsecaxis - spatial axis vector
;
; OPTIONAL OUTPUTS:
;       npix    - number of pixels (NAXIS2)
;       axis_0  - spatial position at REFPIX [arcsec] (CRVAL2)
;       pscale  - plate scale [arcsec/pixel] (CD2_2)
;       refpix  - reference pixel number [pixel] (CRPIX2)
;
; COMMENTS:
;
; PROCEDURES USED:
;	SXPAR(), SPLOG
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2002 Mar 1, U of A
;       jm03apr14uofa, added error checking
;       jm03may01uofa, added more keyword outputs
;-

function make_arcsecaxis, header, npix=npix, axis_0=axis_0, refpix=refpix, pscale=pscale

    nhead = n_elements(header)
    if nhead eq 0L then begin
       print, 'Syntax - arcsecaxis = make_arcsecaxis(header,[npix=,axis_0=,refpix=,pscale=])'
       return, -1L
    endif
    
    npix = sxpar(header,'NAXIS2',count=count1)
    axis_0 = sxpar(header,'CRVAL2',count=count2)
    refpix = sxpar(header,'CRPIX2',count=count3)-1L ; IDL is zero-indexed
    pscale = sxpar(header,'CD2_2',count=count4)

    if (float(count1+count2+count3+count4))[0] ne 4.0 then begin
       splog, 'Insufficient spatial header information.'
       return, -1L
    endif
    
    indx = lindgen(npix)-refpix
    arcsecaxis = axis_0 + indx*pscale
    
return, arcsecaxis
end
