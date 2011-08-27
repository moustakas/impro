;+
; NAME:
;   IM_CALL_MANGLE()
;
; PURPOSE:
;   Simple wrapper on various mangle scripts.
;
; INPUTS: 
;   polyfiles - one or more polygon files to process using mangle
;   outfile - output file name (required)
;
; OPTIONAL INPUTS: 
;   tolerance - intersection tolerance (default 1D-8 arcsec)
;   scheme - pixelization scheme (note the native SDSS scheme is
;     -Pd0,6)  
;   minarea - discard polygons with area<MINAREA using poly2poly
;     (default 1D-13 sr)
;   tmpdir - output directory for temporary, junk files (default /tmp) 
;
; KEYWORD PARAMETERS: 
;   Optional mangle commands (see COMMENTS):
;     selfsnap
;     snap
;     pizelize
;     balkanize
;     unify
;     poly2poly
;   silent - suppress messages to STDOUT
;
; OUTPUTS: 
;   Mangled polygon file is written to OUTFILE.
;
; COMMENTS:
;   See http://space.mit.edu/~molly/mangle for lots of handy info.  
; 
;   This routine is a simple wrapper on mangle, (very!) loosely based
;   on K. Wong's CALL_MANGLE.  It basically just does 'snap -S', 
;   'pixelize', 'snap', 'balkanize', 'unify' (optional), and finally
;   'poly2poly'.  
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Sep 09, UCSD
;
; Copyright (C) 2010, John Moustakas
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

pro do_mangle, cmd, silent=silent
    if (keyword_set(silent) eq 0) then splog, cmd
    t0 = systime(1)
    spawn, cmd, /sh
    if (keyword_set(silent) eq 0) then splog, 'Total time = ', $
      (systime(1)-t0)/60.0
return
end

pro im_call_mangle, polyfiles, outfile=outfile, tolerance=tolerance, $
  scheme=scheme, minarea=minarea, tmpdir=tmpdir, selfsnap=selfsnap, $
  pixelize=pixelize, snap=snap, balkanize=balkanize, unify=unify, $
  poly2poly=poly2poly, silent=silent

    npoly = n_elements(polyfiles)
    if (npoly eq 0L) then begin
       doc_library, 'im_call_mangle'
       return
    endif
    join_polyfiles = strjoin(polyfiles,' ') ; join the file names

    if (n_elements(outfile) eq 0) then begin
       splog, 'Output file name OUTFILE required'
       return
    endif

    if (n_elements(tmpdir) eq 0) then tmpdir = '/tmp/'
    if (n_elements(scheme) eq 0) then scheme = '' ; '-Pd0,6'
    if (n_elements(minarea) eq 0) then minarea = 1D-13 ; [sr]
    if (n_elements(tolerance) eq 0) then tolerance = 1D-8 ; 1E-5 ; [arcsec]
    tol = strtrim(string(tolerance,format='(E16.2)'),2)+'s '
    
; do it!
    infile = join_polyfiles
    if keyword_set(selfsnap) then begin
       tmp_outfile = tmpdir+'junk_snapfile1.ply'
       do_mangle, 'snap -S -m'+tol+' '+infile+' '+tmp_outfile, silent=silent
       infile = tmp_outfile
    endif

    if keyword_set(pixelize) then begin
       tmp_outfile = tmpdir+'junk_pixelfile1.ply'
       do_mangle, 'pixelize -m'+tol+' '+scheme+' '+infile+' '+tmp_outfile, silent=silent
       infile = tmp_outfile
    endif
    
    if keyword_set(snap) then begin
       tmp_outfile = tmpdir+'junk_snapfile2.ply'
       do_mangle, 'snap -m'+tol+' '+infile+' '+tmp_outfile, silent=silent
       infile = tmp_outfile
    endif

    if keyword_set(balkanize) then begin
       tmp_outfile = tmpdir+'junk_balkanfile1.ply'
       do_mangle, 'balkanize -m'+tol+' '+infile+' '+tmp_outfile, silent=silent
       infile = tmp_outfile
    endif

    if keyword_set(unify) then begin
       tmp_outfile = tmpdir+'junk_unifyfile1.ply'
       do_mangle, 'unify -m'+tol+' '+infile+' '+tmp_outfile, silent=silent
       infile = tmp_outfile
    endif

    if keyword_set(poly2poly) then begin
       tmp_outfile = tmpdir+'junk_polyfile1.ply'
       do_mangle, 'poly2poly -m'+tol+' -k'+strtrim(string(minarea,$
         format='(E16.2)'),2)+' '+infile+' '+tmp_outfile, silent=silent
    endif

    if (n_elements(tmp_outfile) eq 0) then begin
       splog, 'Nothing done!!'
    endif else begin
       splog, 'Writing '+outfile
       spawn, '/bin/cp -f '+tmp_outfile+' '+outfile, /sh
    endelse

return
end
