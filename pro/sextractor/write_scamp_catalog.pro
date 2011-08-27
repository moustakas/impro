;+
; NAME:
;   WRITE_SCAMP_CATALOG
;
; PURPOSE:
;   Write a SCAMP-compatible (in FITS_LDAC format) astrometric
;   reference catalog. 
;
; INPUTS: 
;   ra - right ascension (double-precision, decimal degrees)
;   dec - declination (double-precision, decimal degrees)
;   mag - photometric magnitude (arbitrary bandpass)
;
; OPTIONAL INPUTS: 
;   ra_err - error in RA (default: 50 mas)
;   dec_err - error in DEC (default: 50 mas)
;   mag_err - error in MAG (default: 0.1 mag); note that MAG_ERR not
;     seem to be used by the current version (1.4.6) of SCAMP
;   outfile - output file name
;
; KEYWORD PARAMETERS: 
;
; OUTPUTS: 
;   Writes a FITS_LDAC file to OUTFILE, which can then be passed to
;   SCAMP using the following SCAMP parameters:
;     ASTREF_CATALOG = 'FILE'
;     ASTREFCAT_NAME = 'OUTFILE' [where OUTFILE is defined above]
;     ASTREFCENT_KEYS = 'RA,DEC'
;     ASTREFERR_KEYS = 'ERR_A,ERR_B'
;     ASTREFMAG_KEY = 'MAG'
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2009 May 04, NYU - written, based in large part on a
;     piece of C code (and documentation) written by D. Lang (Toronto)
;
; Copyright (C) 2009, John Moustakas
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

pro write_scamp_catalog, ra, dec, mag, flags=flags, ra_err=ra_err, $
  dec_err=dec_err, mag_err=mag_err, outcat=outcat, $
  outfile=outfile

    nobj = n_elements(ra)
    if (n_params() ne 3) or (nobj eq 0) then begin
       doc_library, 'write_scamp_catalog'
       return
    endif

    if (n_elements(ra_err) eq 0) then ra_err = ra*0.0+0.05/3600.0    ; 50 mas
    if (n_elements(dec_err) eq 0) then dec_err = dec*0.0+0.05/3600.0 ; 50 mas
    if (n_elements(mag_err) eq 0) then mag_err = mag*0.0+0.1   
    if (n_elements(flags) eq 0) then flags = long(mag*0.0)
    if (n_elements(outfile) eq 0L) then outfile = 'scamp.cat'
    
    if (n_elements(dec) ne nobj) then begin
       splog, 'Dimensions of RA, DEC, and MAG do not match'
       return
    endif
    if (n_elements(mag) ne nobj) then begin
       splog, 'Dimensions of RA, DEC, and MAG do not match'
       return
    endif
    if (n_elements(ra_err) ne nobj) then begin
       splog, 'Dimensions of RA and RA_ERR do not match'
       return
    endif
    if (n_elements(dec_err) ne nobj) then begin
       splog, 'Dimensions of DEC and DEC_ERR do not match'
       return
    endif
    if (n_elements(mag_err) ne nobj) then begin
       splog, 'Dimensions of MAG and MAG_ERR do not match'
       return
    endif
    if (n_elements(flags) ne nobj) then begin
       splog, 'Dimensions of RA, DEC, and FLAGS do not match'
       return
    endif

; generate the generic field header card; note that every header entry
; must be 80 characters long, but this is taken care of by MKHDR
    mkhdr, hdr, 0, extend=0
    sxdelpar, hdr, 'DATE'
    sxdelpar, hdr, 'COMMENT'
    sxdelpar, hdr, 'EXTEND'
    sxaddpar, hdr, 'BITPIX', 0
    sxaddpar, hdr, 'NAXIS', 2
    sxaddpar, hdr, 'NAXIS1', 0
    sxaddpar, hdr, 'NAXIS2', 0
    nhdr = n_elements(hdr)
    out1 = {field_header_card: strarr(nhdr)}
    out1.field_header_card = hdr

; now make the source list FITS table; at minimum, you need RA, DEC,
; RA_ERR, DEC_ERR, and MAG; scamp does not seem to use MAG_ERR; use
; the default scamp keywords
    outcat = {x_world: 0.0D, y_world: 0.0D, erra_world: 0.0D, $
      errb_world: 0.0D, errtheta_world: 0.0D, mag: 0.0, magerr: 0.0, $
      flags: 0}
    outcat = replicate(outcat,nobj)
    outcat.x_world = ra
    outcat.y_world = dec
    outcat.erra_world = ra_err
    outcat.errb_world = dec_err
    outcat.mag = mag
    outcat.magerr = mag_err
    outcat.flags = flags

; write out; we have to perform some FITS table juggling to get the
; EXTNAME header keywords (required by SCAMP) written out
    splog, 'Writing '+outfile
    mwrfits, 0, outfile, /create
    mwrfits, out1, outfile
    mwrfits, outcat, outfile

    j1 = mrdfits(outfile,1,h1,/silent)
    j2 = mrdfits(outfile,2,h2,/silent)
    sxaddpar, h1, 'EXTNAME', 'LDAC_IMHEAD'
    sxaddpar, h2, 'EXTNAME', 'LDAC_OBJECTS'

    mwrfits, 0, outfile, /create
    mwrfits, out1, outfile, h1
    mwrfits, outcat, outfile, h2

return
end
    
