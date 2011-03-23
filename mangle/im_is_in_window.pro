;+
; NAME:
;   IM_IS_IN_WINDOW()
;
; PURPOSE:
;   Much faster alternative to IS_IN_WINDOW() that relies on the
;   mangle 'polyid' command.  
;
; INPUTS: 
;   polyfile - polygon file name
;   ra,dec - celestial coordinates (decimal degrees, preferably
;     double precision) [NGALS]
;
; OPTIONAL INPUTS: 
;   tolerance - intersection tolerance (default 1D-8 arcsec)
;   scheme - pixelization scheme; default is the native SDSS scheme,
;     -Pd0,6 
;   tmpdir - output directory for temporary, junk files (default /tmp) 
;
; KEYWORD PARAMETERS: 
;   pixelize - pixelize the polygon file before running polyid
;     (default is to assume that 'pixelize' has already been run) 
;
; OUTPUTS: 
;   outflag - Boolean flag: 1=in window, 0=not in window [NGALS] 
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Given a polygon file name and a list of coordinates (RA,DEC), this
;   routine returns a Boolean flag indicating whether or not each
;   object in the input list falls within any of the polygons in
;   POLYFILE.  Unlike IS_IN_WINDOW(), however, this routine uses
;   polyid, which is much fast for large numbers of polygons and/or
;   galaxies. 
; 
;   Loosely based on J. Aird's IS_IN_WIN_PIX. 
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Sep 09, UCSD
;   jm11mar15ucsd - optionally return the polyid
;
; Copyright (C) 2010-2011, John Moustakas
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

function im_is_in_window, polyfile, ra=ra, dec=dec, pixelize=pixelize, $
  scheme=scheme, tmpdir=tmpdir, tolerance=tolerance, polyid=polyid, $
  silent=silent

    if (n_elements(polyfile) eq 0) then begin
       doc_library,'im_is_in_window'
       return, 0
    endif
    if (file_test(polyfile) eq 0) then begin
       splog, 'Polygon file '+polyfile+' not found!'
       return, 0
    endif

    ngals = n_elements(ra)
    if (ngals eq 0L) then begin
       splog, 'RA and DEC inputs required!'
       return, 0
    endif
    if (ngals ne n_elements(dec)) then begin
       splog, 'RA and DEC must have the same number of elements!'
       return, 0
    endif

    if (n_elements(tmpdir) eq 0) then tmpdir = '/tmp/'
    if (n_elements(scheme) eq 0) then scheme = '-Pd0,6'
    if (n_elements(tolerance) eq 0) then tolerance = 1D-8 ; [arcsec]
    tol = strtrim(string(tolerance,format='(E16.2)'),2)+'s '

    t0 = systime(1)
; create the coordinate file list
    coordfile = tmpdir+'tmp_coords.dat'
    if (keyword_set(silent) eq 0) then splog, $
      'Writing temporary coordinate file list '+coordfile
    forprint, ra, dec, textout=coordfile, /nocomment, $
      format='E16.8,2x,E16.8', /silent

; pixelize and then snap, if necessary
    infile = polyfile
    if keyword_set(pixelize) then begin
       tmp_outfile = tmpdir+'tmp_pixelize.ply'
       cmd = 'pixelize -m'+tol+' '+scheme+' '+infile+' '+tmp_outfile
       if (keyword_set(silent) eq 0) then splog, cmd
       spawn, cmd, /sh
       infile = tmp_outfile

       tmp_outfile = tmpdir+'tmp_snapfile.ply'
       cmd = 'snap -m'+tol+' '+infile+' '+tmp_outfile
       if (keyword_set(silent) eq 0) then splog, cmd
       spawn, cmd, /sh
       infile = tmp_outfile
    endif
    
;; run snap and then polyid
;    snapfile = tmpdir+'tmp_snapfile.ply'
;    polyidfile = tmpdir+'tmp_polyidfile.ply'
;    cmd = 'snap -m'+tol+' '+tmp_polyfile+' '+snapfile
;    splog, cmd & spawn, cmd, /sh
;    spawn, 'polyid '+snapfile+' '+coordfile+' '+polyidfile, /sh

; run polyid
    polyidfile = tmpdir+'tmp_polyidfile.ply'
    cmd = 'polyid '+infile+' '+coordfile+' '+polyidfile
    if (keyword_set(silent) eq 0) then splog, cmd
    spawn, cmd, /sh

; parse the results using some junk files; a non-blank third column
; means the galaxy overlaps a particular polygon
    tmpfile1 = tmpdir+'tmp_out1.ply'
    spawn, "awk < "+polyidfile+" '{print $3}' > "+tmpfile1, /sh

    junk = ''
    openr, lun, tmpfile1, /get_lun
    readf, lun, junk ; skip the header
    strflag = strarr(ngals)
    readf, lun, strflag
    free_lun, lun

    polyid = strcompress(strflag,/remove)
    none = where(polyid eq '')
    if (none[0] ne -1L) then polyid[none] = -1
    polyid = long(polyid)

    outflag = polyid ne -1 ; 1=in, 0=out
;   outflag = fix(strcompress(strflag,/remove) ne '') 

    if (keyword_set(silent) eq 0) then splog, 'Total time ', $
      (systime(1)-t0)/60.0

return, outflag
end
