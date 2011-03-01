;+
; NAME:
;       IMAGES2HTML
;
; PURPOSE:
;       Generate web pages for lists of images.
;
; CALLING SEQUENCE:
;       
;
; INPUTS:
;
;
; OPTIONAL INPUTS:
;
;
; KEYWORD PARAMETERS:
;       nothumbs - do not generate thumbnails
;
; OUTPUTS:
;
;
; OPTIONAL OUTPUTS:
;
;
; COMMON BLOCKS:
;
;
; PROCEDURES USED:
;
;
; COMMENTS:
;
;
; EXAMPLES:
;
;
; MODIFICATION HISTORY:
;
;-

pro images2html, htmlfile=htmlfile, htmlpath=htmlpath, profilefile=profilefile, $
  profilepath=profilepath, datfile=datfile, nothumbs=nothumbs, title=title, $
  textcolor=textcolor, bgcolor=bgcolor

    if n_elements(htmlfile) eq 0L then htmlfile = 'index.html'
    if n_elements(htmlpath) eq 0L then htmlpath = './'
    if n_elements(profilefile) eq 0L then profilefile = 'profile.dat'
    if n_elements(profilepath) eq 0L then profilepath = './'
    if n_elements(datfile) eq 0L then datfile = 'images.dat'

    if file_test(htmlpath+htmlfile) then begin
       splog, 'Making backup of '+htmlpath+htmlfile+'.'
       spawn, ['cp -f '+htmlpath+htmlfile+' '+htmlpath+htmlfile+'.old'], /sh
    endif

    if n_elements(title) eq 0L then title = ''
    if n_elements(textcolor) eq 0L then textcolor = 'red'
    if n_elements(bgcolor) eq 0L then bgcolor = 'black'

; read the images data file

    if file_test(htmlpath+datfile) eq 0L then begin
       print, 'Images file '+htmlpath+datfile+' not found.'
       return
    endif
    
;   readcol, htmlpath+datfile, imlist, captions, format='A,A', $
;     /silent, delimiter='|', comment='#'

    dattext = djs_readlines(htmlpath+datfile)
    keep = where(strmatch(dattext,'#*') ne 1B,nkeep)
    if (nkeep ne 0L) then dattext = dattext[keep]
    
    imlist = '' & captions = '' 
    for i = 0L, n_elements(dattext)-1L do begin
       info = strsplit(dattext[i],'|',/extract,count=ninfo)
       imlist = [imlist,info[0]] 
       if ninfo eq 1L then captions = [captions,''] else captions = [captions,info[1]]
    endfor 
    imlist = strcompress(imlist[1L:n_elements(imlist)-1L],/remove)
    captions = strtrim(captions[1L:n_elements(captions)-1L],2)
    nimlist = n_elements(imlist)

; read the profile data file

    profile = {$
      name:       '', $
      email:      '', $
      url:        ''}

    if file_test(profilepath+profilefile) then begin
       protext = djs_readlines(profilepath+profilefile)
       profile.name = protext[0]
       profile.email = protext[1]
       profile.url = protext[2]
    endif
    
; generate a log file

    openw, lun, 'images2html.log', /get_lun
    printf, lun, 'HMTLPATH    - '+htmlpath
    printf, lun, 'HMTLFILE    - '+htmlfile
    printf, lun, ''
    printf, lun, 'PROFILEPATH - '+profilepath
    printf, lun, 'PROFILEFILE - '+profilefile
    printf, lun, '   name     - '+profile.name
    printf, lun, '   email    - '+profile.email
    printf, lun, '   url      - '+profile.url
    printf, lun, ''
    printf, lun, 'DATFILE     - '+datfile
    free_lun, lun
    
; generate smaller images

    thumblist = 'thumb_'+imlist
    if not keyword_set(nothumbs) then begin

       iconsize = 250
       for k = 0L, nimlist-1L do begin
          print, 'Writing '+thumblist[k]
;         spawn, ['mkdir -p thumbnails'], /sh
          spawn, ['convert -geometry x'+strn(iconsize)+' '+imlist[k]+' '+thumblist[k]], /sh
       endfor
    
    endif
    
    ncolumn = 2
    nrow = long(ceil(nimlist/float(ncolumn)))

    openw, lun1, htmlpath+htmlfile, /get_lun
    printf, lun1, '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://cerebus.as.arizona.edu/~ioannis">'
    printf, lun1, '<html>'
    printf, lun1, '<head>'
    printf, lun1, '<meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1"></meta>'
    printf, lun1, '<title>'+title+'</title>'
    printf, lun1, '<style type="text/css">'
    printf, lun1, 'body {'
    printf, lun1, '   background-color: '+bgcolor+';'
    printf, lun1, '   font-family: "Times New Roman" ;'
    printf, lun1, '   font-size: 100%;'
    printf, lun1, '   font-weight: bold;'
    printf, lun1, '   color: '+textcolor+';'
    printf, lun1, '   text-align: center;'
;   printf, lun1, '   margin-left: 5%;'
;   printf, lun1, '   margin-right: 5%;'
;   printf, lun1, '   margin-top: 5%;'
    printf, lun1, '}'
    printf, lun1, 'a:link {'
    printf, lun1, '   color: '+textcolor+';'
    printf, lun1, '   text-decoration: none       ;'
    printf, lun1, '}'
    printf, lun1, 'a:visited {'
    printf, lun1, '   color: '+textcolor+';'
    printf, lun1, '   text-decoration: none       ;'
    printf, lun1, '}'
;    printf, lun1, 'img {'
;    printf, lun1, 'height: 300'
;    printf, lun1, '}'
    printf, lun1, 'div.base {'
    printf, lun1, '   margin: 1ex;'
    printf, lun1, '   padding: 0;'
    printf, lun1, '   width: auto;'
    printf, lun1, '}'
    printf, lun1, 'div.row {'
    printf, lun1, '   border: none;'
    printf, lun1, '   margin-top: 0;'
    printf, lun1, '   margin-right: auto;'
    printf, lun1, '   margin-bottom: 0;'
    printf, lun1, '   margin-left: auto;'
    printf, lun1, '   padding: 0;'
    printf, lun1, '   width: 100%;'
    printf, lun1, '}'
    printf, lun1, 'div.left {'
    printf, lun1, '   border: none;'
    printf, lun1, '   float: left;'
    printf, lun1, '   margin: 0;'
    printf, lun1, '   margin-bottom: 2ex;'
    printf, lun1, '   padding: 0;'
    printf, lun1, '   width: 49%;'
    printf, lun1, '}'
    printf, lun1, 'div.right {'
    printf, lun1, '   border: none;'
    printf, lun1, '   float: right;'
    printf, lun1, '   margin: 0;'
    printf, lun1, '   margin-bottom: 2ex;'
    printf, lun1, '   padding: 0;'
    printf, lun1, '   width: 49%;'
    printf, lun1, '}'
    printf, lun1, '</style>'
    printf, lun1, '</head>'
    printf, lun1, '<body>'
    printf, lun1, '<h1>'+title+'</h1>'

    printf, lun1, '<div class="base">'

    for i = 0L, nrow-1L do begin

       printf, lun1, '<div class="row">' ; open the row

       printf, lun1, '<div class="left">'
       printf, lun1, '<a href="'+imlist[2*i]+'"><img src="'+thumblist[2*i]+'" '+$
         'alt = "'+imlist[2*i]+'" /></a><br /><p>'+captions[2*i]+'</p>'
       printf, lun1, '</div>'

       if ((2.0*i+1) lt nimlist) then begin
          printf, lun1, '<div class="right">'
          printf, lun1, '<a href="'+imlist[2*i+1]+'"><img src="'+thumblist[2*i+1]+'" '+$
            'alt = "'+imlist[2*i+1]+'" /></a><br /><p>'+captions[2*i+1]+'</p>'
          printf, lun1, '</div>'
       endif

       printf, lun1, '</div>' ; close the row

    endfor
    
    printf, lun1, '</div>' ; close the base
    printf, lun1, '<br /><br /><br clear="bottom">'
    
; ---------------------------------------------------------------------------
;    printf, lun1, '<table width="100%" border="0" cellspacing="0" cellpadding="0" bordercolor="#000000">'
;    printf, lun1, '<tbody>'
;
;    for i = 0L, nrow-1L do begin
;
;       printf, lun1, '<tr>'
;       for j = 0L, ncolumn-1L do begin
;          indx = j+i*ncolumn
;          if indx lt nimlist then $
;            printf, lun1, '<td><a href="'+imlist[indx]+'"><img src="'+thumblist[indx]+'" '+$
;              'alt = "'+imlist[indx]+'" /></a><br />'+captions[indx]+'</td>'
;       endfor
;       printf, lun1, '</tr>'
;
;    endfor
;
;    printf, lun1, '</tbody>'
;    printf, lun1, '</table>'
; ---------------------------------------------------------------------------
    
;   printf, lun1, '<p>Web page generated '+im_today()+'.</p>'
    printf, lun1, '</body>'
    printf, lun1, '</html>'
    free_lun, lun1

return
end    
