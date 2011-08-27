;+
; NAME:
;       IM_READ_VIZIER_TSV()
;
; PURPOSE:
;       Read "|" separated ASCII files downloaded from Vizier and
;       outputs a data structure. 
;
; INPUTS:
;       filename - see READ_FMR()
;
; OPTIONAL INPUTS:
;       prefix - string to prepend to each tag in the output data
;                structure 
;
; OUTPUTS:
;       data - parsed data structure
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2007 Nov 02, NYU, written
;
; Copyright (C) 2007, John Moustakas
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

function im_read_vizier_tsv, filename, prefix=prefix

    if (n_elements(filename) ne 1L) then begin
       doc_library, 'im_read_vizier_tsv'
       return, -1L
    endif

    if (n_elements(prefix) eq 0L) then prefix = '' else prefix = prefix+'_'

    alldata = djs_readlines(filename)
    alldata = alldata[where((strmatch(strcompress(alldata,/remove),'*#*') eq 0B) and $ ; remove comments & blank lines
      (strcompress(alldata,/remove) ne''))]

    head = alldata[0L] ; column names
    data = alldata[3L:n_elements(alldata)-1L] ; first three rows are the column names, units, and a bunch of dashes 
    ngalaxy = n_elements(data)

    head_tags = strsplit(head,'|',/extract)
    nhead_tags = n_elements(head_tags)
    tags = strarr(nhead_tags)
    for ii = 0L, nhead_tags-1L do tags[ii] = idl_validname(head_tags[ii],/convert_all)
;   tags = idl_validname(strsplit(head,'|',/extract),/convert_all)
    ntags = n_elements(tags)

    tempfile = '/tmp/junk.junk'
    
    openw, lun, tempfile, /get_lun
    for it = 0L, ntags-1L do printf, lun, '# '+string(it+1L,format='(I2)')+' '+tags[it]
    for ig = 0L, ngalaxy-1L do begin
       line = strcompress(strsplit(data[ig],'|',/extract),/remove)
       for iline = 0L, n_elements(line)-1L do if (line[iline] eq '') then $ ; replace blank lines with junk
         line[iline] = '-999'
       printf, lun, line, format='('+string(ntags,format='(I0)')+'A25)'
    endfor
    free_lun, lun

    struct = rsex(tempfile)
    rmfile, tempfile
    
return, struct
end
    
