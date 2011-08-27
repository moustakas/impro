;+
; NAME:
;   IM_PS2HTML
;
; PURPOSE:
;   Convert all PS and EPS files in a directory to PNG files and
;   render a web page.
;
; INPUTS:
;   htmlbase - base name for the HTML page and the subdirectory
;              (relative to HTML_PATH) where the postscript files
;              have been written
;
; OPTIONAL INPUTS:
;   pslist    - input list of PS files; default it to search for
;               all EPS files in HTMLBASE (or all PS files if
;               PSFIND=1) 
;   npscols   - number of columns in the HTML page (default 3)
;   html_path - postscript files must exist in the path
;               HTML_PATH+HTMLBASE and the web page is written to
;               HTML_PATH 
;
; KEYWORD PARAMETERS:
;   psfind    - search for PS files; default is to only search for
;               EPS files
;   pngfind   - search for PNG files; default is to only search for
;               EPS files (this keyword is stronger than PSFIND); if
;               set, then force: ONLY_HTML=1, ONLY_PNG=0, CLEANPNG=0 
;   only_html - do not generate PNG files, just re-render the html 
;   only_png  - do not render the html, just generate PNG files 
;   cleanpng  - remove all PNG files from HTML_PATH+HTMLBASE 
;   silent    - suppress messages to STDOUT 
;
; OUTPUTS:
;   An html file called HTML_PATH+HTMLBASE+'.html'.  PNG files are
;   also generated from the EPS or PS files found in
;   HTML_PATH+HTMLBASE. 
;
; COMMENTS:
;   Uses spawn to run the UNIX command CONVERT.
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2004 Jan 08, U of A - written
;   jm04nov16uofa - added menu at the top of the web page 
;   jm05jan18uofa - added image captions and better
;                   cross-referencing information
;   jm05jan21uofa - added support for .EPS files 
;   jm05mar31uofa - fixed image centering 
;   jm05nov25uofa - added PSFIND keyword; default is to only look
;                   for EPS files 
;   jm09jan22nyu  - added PNGFIND keyword; default is to only look 
;                   for EPS files 
;
; Copyright (C) 2004-2005, 2009, John Moustakas
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

pro im_ps2html, htmlbase, pslist=pslist, pnglist=pnglist, npscols=npscols, $
  psfind=psfind, pngfind=pngfind, only_html=only_html, only_png=only_png, $
  html_path=html_path, cleanpng=cleanpng, silent=silent
  
    nhtmlbase = n_elements(htmlbase)
    if (nhtmlbase ne 1L) then begin
       doc_library, 'im_ps2html'
       return
    endif
    
    if (n_elements(html_path) eq 0L) then html_path = './'
    if (n_elements(npscols) eq 0L) then npscols = 3L

    fullpath = html_path+htmlbase+'/'
    if file_test(fullpath,/directory) eq 0L then begin
       splog, 'Path '+html_path+htmlbase+' does not exist.'
       return
    endif

    if keyword_set(pngfind) then begin
       only_html = 1
       only_png = 0
       cleanpng = 0
       suff = 'PNG' 
    endif else suff = 'PS'

    if (n_elements(pslist) eq 0L) then begin
       psfiles = file_search(fullpath+'*.eps')
       if keyword_set(psfind) then psfiles = file_search(fullpath+'*.ps')
       if keyword_set(pngfind) then psfiles = file_search(fullpath+'*.png')
       npsfiles = n_elements(psfiles)
       if (npsfiles eq 0L) then begin
          splog, 'No '+suff+' files found in '+fullpath
          return
       endif else psfiles = file_basename(psfiles)
    endif else begin
       psfiles = pslist
       tempfiles = file_search(fullpath+pslist,count=npsfiles)
       if (npsfiles ne n_elements(pslist)) then begin
          splog, 'Some files in PSLIST do not exist in '+fullpath
          return
       endif
    endelse

    if (not keyword_set(silent)) then splog, 'Found '+$
      string(npsfiles,format='(I0)')+' '+suff+' files'

; clean up PNG files from calling this routine previously 

    if keyword_set(cleanpng) then begin
       if not keyword_set(silent) then splog, 'Deleting all PNG files in '+fullpath+'.'
       spawn, 'rm -f '+html_path+htmlbase+'/*.png*', /sh
;      spawn, 'rm -f '+html_path+htmlbase+'/*.ps', /sh
    endif

    pngfiles = repstr(repstr(psfiles,'.ps','.png'),'.eps','.png')
    if not keyword_set(only_html) then for i = 0L, npsfiles-1L do begin
       if not keyword_set(silent) then splog, 'Converting '+psfiles[i]+' to '+pngfiles[i]
       spawn, 'convert '+fullpath+psfiles[i]+' '+fullpath+pngfiles[i], /sh
    endfor

; generate the web page HTML page

    rootpsfiles = repstr(repstr(psfiles,'.ps',''),'.eps','') ; NOTE! does not include EPS files

    nmenucols = npscols

    if (not keyword_set(only_png)) then begin

       xwidth = string(fix(100/float(npscols)),format='(I0)')
;      xwidth = string(fix(100/float(npscols)-npscols*3.0-npscols*15.0),format='(I0)')
;      xwidth = string(fix(100/float(npscols)-npscols*5.0),format='(I0)')
;      xwidth = '100'
       
       if (not keyword_set(silent)) then splog, 'Writing '+html_path+htmlbase+'.html'
       
       openw, lun1, html_path+htmlbase+'.html',/get_lun
       printf, lun1, '<html>'
       printf, lun1, '<head>'
       printf, lun1, '<title>Postscript Output for '+strupcase(htmlbase)+'</title>'
       printf, lun1, '<style type="text/css">'
       printf, lun1, '<!--'
       printf, lun1, 'body {'
       printf, lun1, '   color: white;'
       printf, lun1, '   background-color: black;'
       printf, lun1, '   text-align: center;'
       printf, lun1, '   border: 0;'
       printf, lun1, '   font-size: 100%;'
       printf, lun1, '   font-weight: bold;'
       printf, lun1, '}'
       printf, lun1, 'td {'
       printf, lun1, '   font-size: 100%;'
       printf, lun1, '   font-weight: bold;'
       printf, lun1, '   text-align: center;'
       printf, lun1, '}'
       printf, lun1, 'a {'
       printf, lun1, '   color: cyan;'
       printf, lun1, '   visited: cyan;'
       printf, lun1, '}'
       printf, lun1, '.center {'
       printf, lun1, '   text-align: center;'
       printf, lun1, '}'
       printf, lun1, '.smaller {'
       printf, lun1, '   font-size: 90%;'
       printf, lun1, '}'
       printf, lun1, 'h1 {color: gold;}'
       printf, lun1, '-->'
       printf, lun1, '</style>'
       printf, lun1, '</head>'
       printf, lun1, '<body>'
       printf, lun1, '<h1>Postscript Output for '+strupcase(htmlbase)+'</h1>'
       printf, lun1, '</br>'

; generate the menu       
       
       printf, lun1, '<table border="1" align="center" frame="box" width="100%">'
       printf, lun1, '<tbody>'
       for i = 0L, ceil(npsfiles/float(nmenucols))-1L do begin

          printf, lun1, '<tr width="100%">'
          for j = 0L, nmenucols-1L do begin
             indx = j + i*nmenucols
             if (indx le npsfiles-1L) then begin
                printf, lun1, '<td width="3%">'+string(j+1,format='(I0)')+':'+$
                  string(i+1,format='(I0)')+'</td>'
                printf, lun1, '<td width="10%"><a href="'+htmlbase+'/'+pngfiles[indx]+'">png</a>'
                if file_test(fullpath+rootpsfiles[indx]+'.ps',/regular) then printf, lun1, $
                  '<a href="'+htmlbase+'/'+rootpsfiles[indx]+'.ps"> ps</a>' 
                if file_test(fullpath+rootpsfiles[indx]+'.eps',/regular) then printf, lun1, $
                  '<a href="'+htmlbase+'/'+rootpsfiles[indx]+'.eps"> eps</a>'
                printf, lun1, '</td>'
                printf, lun1, '<td>'+rootpsfiles[indx]+'</td>'
             endif
;            print, i, psfiles[indx]
          endfor
          printf, lun1, '</tr>'

       endfor
       printf, lun1, '</tbody>'
       printf, lun1, '</table>'
       printf, lun1, '<br />'
       printf, lun1, '<table align="center" cellpadding="10" border="2" cellspacing="10" width="100%">'
       printf, lun1, '<tbody>'

       for i = 0L, ceil(npsfiles/float(npscols))-1L do begin

; images          
          
          printf, lun1, '<tr>'
          for j = 0L, npscols-1L do begin
             indx = j + i*npscols
             if (indx le npsfiles-1L) then begin
                printf, lun1, '<td width="'+xwidth+'%"><a href="'+htmlbase+'/'+pngfiles[indx]+'">'+$
                  '<img width="100%" align="center" valign="center" src="'+htmlbase+'/'+pngfiles[indx]+'"></a></td>'
             endif
          endfor
          printf, lun1, '</tr>'

; image captions          
          
          printf, lun1, '<tr>'
          for j = 0L, npscols-1L do begin
             indx = j + i*npscols
             if (indx le npsfiles-1L) then begin
                printf, lun1, '<td width="'+xwidth+'%" class="center">'+string(j+1,format='(I0)')+':'+$
                  string(i+1,format='(I0)')
                printf, lun1, rootpsfiles[indx]+' <a href="'+htmlbase+'/'+pngfiles[indx]+'">png</a>'
                if file_test(fullpath+rootpsfiles[indx]+'.ps',/regular) then printf, lun1, $
                  '<a href="'+htmlbase+'/'+rootpsfiles[indx]+'.ps"> ps</a>' 
                if file_test(fullpath+rootpsfiles[indx]+'.eps',/regular) then printf, lun1, $
                  '<a href="'+htmlbase+'/'+rootpsfiles[indx]+'.eps"> eps</a>'
                printf, lun1, '</td>'
             endif
          endfor
          printf, lun1, '</tr>'
          
       endfor 

       printf, lun1, '</tbody>'
       printf, lun1, '</table>'
       printf, lun1, '</body>'
       printf, lun1, '</html>'
       free_lun, lun1

    endif 

return
end

