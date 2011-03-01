pro im_ps2pdf, gz=gz
; jm08aug28nyu - convert landscape PS plots to properly rotated PDF

    if keyword_set(gz) then begin
       psgzname = file_search('*.ps.gz',count=nps)
       psname = repstr(psgzname,'.ps.gz','.ps')
    endif else psname = file_search('*.ps',count=nps)
    pdfname = repstr(psname,'.ps','.pdf')

    t0 = systime(1)
    for ii = 0L, nps-1L do begin
       if keyword_set(gz) then begin
          splog, 'Unzipping '+psgzname[ii]
          spawn, 'gunzip -f '+psgzname[ii]
       endif
       splog, 'Converting '+psname[ii]+'-->'+pdfname[ii]
       spawn, 'ps2pdf '+psname[ii]+' '+pdfname[ii]
       spawn, 'pdftk '+pdfname[ii]+' cat 1-endW output /tmp/junk.pdf'
       spawn, '/bin/mv -f /tmp/junk.pdf '+pdfname[ii]
    endfor
    splog, 'Total time = ', (systime(1)-t0)/60.0
       
return
end
    
