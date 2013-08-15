;+
; NAME:
;   WRITE_PEGASE_IMF
;
; PURPOSE:
;   Write a Pegase-compatible IMF file.
;
; INPUTS: 
;   mass - stellar mass boundaries of the desired IMF [NPIECE+1] 
;   slope - logarithmic slope between each mass interval [NPIECE] 
;
; OPTIONAL INPUTS: 
;   outfile - output file name (default 'IMF_custom.dat')
;
; KEYWORD PARAMETERS: 
;
; OUTPUTS: 
;   The appropriate output file
;
; OPTIONAL OUTPUTS:
;
; EXAMPLES:
;   First try the Salpeter IMF:
;      IDL> write_pegase_imf, [0.1,100.0], 1.0-2.35, $
;      IDL>   outfile='IMF_mySalpeter.dat'

;   Alternatively, try the "cosmic IMF" of Wilkins+08 (eq. 2):
;      IDL> write_pegase_imf, [0.01,0.08,0.5,150.0], $
;      IDL>   1.0-[0.3,1.3,2.35], outfile='IMF_wilkins.dat'
;
; COMMENTS:
;   See pg. 15 of the Pegase manual for details.  In particular, note
;   that Pegase defines the *logarithmic* slope (i.e., the Salpeter
;   slope is 1.35 in this case, not 2.35 as in the linear definition
;   of the IMF).  Also remember that the IMF filename must be added to
;   the 'list_IMFs.dat' file in order to be recognized! 
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Mar 12, UCSD
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

pro write_pegase_imf, mass, slope, outfile=outfile

    if (n_elements(mass) eq 0) or (n_elements(slope) eq 0) then begin
       doc_library, 'write_pegase_imf'
       return
    endif

    if (n_elements(mass) ne n_elements(slope)+1) then begin
       splog, 'MASS must an NPIECE+1 element array!'
       return
    endif
    
    npiece = n_elements(slope)
    if (n_elements(outfile) eq 0) then outfile = 'IMF_custom.dat'
    splog, 'Writing '+outfile
    openw, lun, outfile, /get_lun
    printf, lun, string(npiece,format='(I0)')
    for ii = 0, npiece-1 do printf, lun, $
      strtrim(string(mass[ii],format='(F12.2)'),2)+' '+$
      strtrim(string(slope[ii],format='(F12.2)'),2)
    printf, lun, strtrim(string(mass[npiece],format='(F12.2)'),2)
    printf, lun, ' '
    free_lun, lun

return
end
    
