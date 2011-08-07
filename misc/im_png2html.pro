;+
; NAME:
;   IM_PNG2HTML
;
; PURPOSE:
;   Render a web page from a list of PNG (or JPG) files.
;
; INPUTS:
;   htmlbase - base name for the HTML page and the subdirectory
;              (relative to HTML_PATH) where the PNG files have been
;              written
;
; OPTIONAL INPUTS:
;   pnglist   - input list of PNG files **including the full path
;               name**; default it to search for all PNG files in
;               HTMLBASE 
;   npngcols  - number of columns in the HTML page (default 3)
;   html_path - files must exist in the path HTML_PATH+HTMLBASE and
;               the web page is written to HTML_PATH 
;   title     - title of the HTML page
;
; KEYWORD PARAMETERS:
;   silent    - suppress messages to STDOUT 
;
; OUTPUTS:
;   An html file called HTML_PATH+HTMLBASE+'.html'.  
;
; COMMENTS:
;   Uses spawn to run the UNIX command CONVERT.
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2009 Jan 26, NYU - based entirely on IM_PS2HTML 
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

pro im_png2html, htmlbase, pnglist=pnglist, npngcols=npngcols, $
  html_path=html_path, title=title, silent=silent
  
    nhtmlbase = n_elements(htmlbase)
    if (nhtmlbase ne 1L) then begin
       doc_library, 'im_png2html'
       return
    endif
    
    if (n_elements(html_path) eq 0L) then html_path = './'
    if (n_elements(npngcols) eq 0L) then npngcols = 3L

    fullpath = html_path+htmlbase+'/'
    if file_test(fullpath,/directory) eq 0L then begin
       splog, 'Path '+html_path+htmlbase+' does not exist.'
       return
    endif

    if (n_elements(pnglist) eq 0L) then begin
       if keyword_set(or_jpeg) then $
         pngfiles = file_search(fullpath+'*.png') else $
           pngfiles = [file_search(fullpath+'*.png'),$
         file_search(fullpath+'*.jpg')]
       npngfiles = n_elements(pngfiles)
       if (npngfiles eq 0L) then begin
          splog, 'No PNG files found in '+fullpath
          return
       endif else pngfiles = file_basename(pngfiles)
    endif else begin
       pngfiles = pnglist
       npngfiles = n_elements(pngfiles)
;      tempfiles = file_search(fullpath+pnglist,count=npngfiles)
;      if (npngfiles ne n_elements(pnglist)) then begin
;         splog, 'Some files in PNGLIST do not exist in '+fullpath
;         return
;      endif
    endelse

    if (not keyword_set(silent)) then splog, 'Found '+$
      string(npngfiles,format='(I0)')+' PNG files'

; generate the web page HTML page

    rootpngfiles = repstr(pngfiles,'.png','')
    nmenucols = npngcols

    if (not keyword_set(only_png)) then begin

       xwidth = string(fix(100/float(npngcols)),format='(I0)')
;      xwidth = string(fix(100/float(npscols)-npscols*3.0-npscols*15.0),format='(I0)')
;      xwidth = string(fix(100/float(npscols)-npscols*5.0),format='(I0)')
;      xwidth = '100'
       
       if (not keyword_set(silent)) then splog, 'Writing '+html_path+htmlbase+'.html'
       
       openw, lun1, html_path+htmlbase+'.html',/get_lun
       printf, lun1, '<html>'
       printf, lun1, '<head>'
       if (n_elements(title) eq 0L) then title = $
         'Postscript Output for '+strupcase(htmlbase)
       printf, lun1, '<title>'+title+'</title>'
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
       printf, lun1, '</br>'
       printf, lun1, '<h1>'+title+'</h1>'
       printf, lun1, '</br>'

; generate the menu       
       
       printf, lun1, '<table border="1" align="center" frame="box" width="100%">'
       printf, lun1, '<tbody>'
       for i = 0L, ceil(npngfiles/float(nmenucols))-1L do begin

          printf, lun1, '<tr width="100%">'
          for j = 0L, nmenucols-1L do begin
             indx = j + i*nmenucols
             if (indx le npngfiles-1L) then begin
                printf, lun1, '<td width="3%">'+string(j+1,format='(I0)')+':'+$
                  string(i+1,format='(I0)')+'</td>'
                printf, lun1, '<td width="10%"><a href="'+htmlbase+'/'+pngfiles[indx]+'">png</a>'
                printf, lun1, '</td>'
                printf, lun1, '<td>'+rootpngfiles[indx]+'</td>'
             endif
;            print, i, pngfiles[indx]
          endfor
          printf, lun1, '</tr>'

       endfor
       printf, lun1, '</tbody>'
       printf, lun1, '</table>'
       printf, lun1, '<br />'
       printf, lun1, '<table align="center" cellpadding="10" border="2" cellspacing="10" width="100%">'
       printf, lun1, '<tbody>'

       for i = 0L, ceil(npngfiles/float(npngcols))-1L do begin

; images          
          
          printf, lun1, '<tr>'
          for j = 0L, npngcols-1L do begin
             indx = j + i*npngcols
             if (indx le npngfiles-1L) then begin
                printf, lun1, '<td width="'+xwidth+'%"><a href="'+htmlbase+'/'+pngfiles[indx]+'">'+$
                  '<img width="100%" align="center" valign="center" src="'+htmlbase+'/'+pngfiles[indx]+'"></a></td>'
             endif
          endfor
          printf, lun1, '</tr>'

; image captions          
          
          printf, lun1, '<tr>'
          for j = 0L, npngcols-1L do begin
             indx = j + i*npngcols
             if (indx le npngfiles-1L) then begin
                printf, lun1, '<td width="'+xwidth+'%" class="center">'+string(j+1,format='(I0)')+':'+$
                  string(i+1,format='(I0)')
                printf, lun1, rootpngfiles[indx]+' <a href="'+htmlbase+'/'+pngfiles[indx]+'">png</a>'
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

