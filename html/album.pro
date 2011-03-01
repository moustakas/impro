pro album_format, file, col1, col2, col3, col4, backup=backup, path=path

    if not keyword_set(path) then path=cwd()
    pushd,path

    col1 = arm_strfmt(col1) & l1 = strn(strlen(col1[0]))
    col2 = arm_strfmt(col2) & l2 = strn(strlen(col2[0]))
    col3 = arm_strfmt(col3) & l3 = strn(strlen(col3[0]))
    col4 = arm_strfmt(col4) & l4 = strn(strlen(col4[0]))

    rearrange, reverse(sort(col1)), col1, col2, col3, col4

    if keyword_set(backup) then spawn,'cp '+file+' '+file+'.old'

    close,2
    openw,2,file
    for i=0,n_elements(col1)-1 do printf,2, $
      f='(A'+l1+',x,A1,x,A'+l2+',x,A1,x,A'+l3+',x,A1,x,A'+l4+')', $
    col1[i],'|',col2[i],'|',col3[i],'|',col4[i]
    close,2
    
    return

 end

pro album_header, fnum, title=title, folder=folder, bgcolor=bgcolor

    if not keyword_set(title) then title = 'TITLE'
    if not keyword_set(folder) then folder = 'FOLDER'

    printf,fnum,'<html>'
    printf,fnum,'<head><title>'+title+'</title></head>'
    printf,fnum,'<body bgcolor="#'+bgcolor+'" text="#000000" '+ $
      'link="#000000" vlink="#000000" alink="#000000">'
    printf,fnum,'<table border=0 width=100% height=100% border=0 bgcolor="#'+bgcolor+'"><tr><td>'
    printf,fnum,'<table border=0 align=center valign=center bgcolor="#'+bgcolor+'" cellpadding=10>'
    printf,fnum,'<tr><td align=center valign=center>'
    printf,fnum,'<font size=150%><b>'+title+'</b></font><br>'
    printf,fnum,'<i>'+folder+'</i><p>'
    printf,fnum,'</table></td></tr>'
    printf,fnum

    return

 end

pro album_main, fnum, image, caption

    printf, fnum, '<tr><td align=center><table border=0 bgcolor="white" cellpadding=10 height=450>'
    printf, fnum, '<tr><td align=center valign=center>'
    printf, fnum, '<img src="images/'+image+'" height=450 align=center><br>'
    printf, fnum, '</td></tr></table><table border=0><tr><td align=center>'+caption
    printf, fnum, '</td></tr></table></td>'


    return

 end

pro album_body, fnum, images, captions

    for i=0,n_elements(images)-1,2 do begin

       printf, fnum, '<tr><td align=center><table width=880>'
       printf, fnum, '<tr><td align=center><table border=0 bgcolor="white" cellpadding=10>'
       printf, fnum, '<tr><td align=center valign=center>'
       printf, fnum, '<img src="images/'+images[i]+'" height=300 align=center><br>'
       printf, fnum, '</td></tr></table><table border=0><tr><td align=center>'+captions[i]
       printf, fnum, '</td></tr></table></td>'

       if n_elements(images)-1 gt i then begin
          printf,fnum

;          printf,fnum,'<td><table border=0 width=420 height=340 bgcolor="white">'
;          printf,fnum,'<tr><td align=center valign=center>'
;          printf,fnum,'<img src="images/'+images[i+1]+'" height=300 align=center><br>'
;
;          printf,fnum,captions[i+1]
;          printf,fnum,'</td></tr></table></td></tr>'
       printf, fnum, '<td align=center><table border=0 bgcolor="white" cellpadding=10>'
       printf, fnum, '<tr><td align=center valign=center>'
       printf, fnum, '<img src="images/'+images[i+1]+'" height=300 align=center><br>'
       printf, fnum, '</td></tr></table><table border=0><tr><td align=center>'+captions[i+1]
       printf, fnum, '</td></tr></table></td></tr>'
          printf,fnum
          printf,fnum
       endif else printf,fnum,'<td></td></tr>'

       printf, fnum, '</table></td></tr>'

    endfor
    
    return

 end

pro album_folders, fnum, folder2, type2, names=names

;stop

; output html for "listed" folders

;    old = findfile('folder_*.html')

    names= arm_unique(folder2[where(type2 ne 'hide')], indices=indices)
    names = ['MAIN PAGE',names]

    count=n_elements(names)

    printf, fnum, '<tr><td align=center><table border=0 cellpadding=10>'

    if count gt 0 then begin
       if count gt 1 then word='photos' else word='photo'

       num = ceil(count/3.)

       printf, fnum, '<tr><td align=center valign=top>'
       for i=0,num-1 < (count-1) do begin
          if i ne 0 then begin
             total = n_elements(where(folder2 eq names[i]))
             printf,fnum,'<a href="folder_'+strcompress(names[i],/remove_all)+'.html">'
             printf,fnum,names[i]+' ('+strn(total)+' '+word+')</a><br>'
          endif else begin
             printf,fnum,'<a href="index.html">'
             printf,fnum,names[i]+'</a><br>'
          endelse
       endfor
       printf, fnum, '</td><td align=center valign=top>'
       if count-1 ge num then for i=num,2*num-1<(count-1) do begin
          total = n_elements(where(folder2 eq names[i]))
          printf,fnum,'<a href="folder_'+strcompress(names[i],/remove_all)+'.html">'
          printf,fnum,names[i]+' ('+strn(total)+' '+word+')</a><br>'
       endfor
       printf, fnum, '</td><td align=center valign=top>'
       if count-1 ge 2.*num then for i=2.*num,3.*num-1<(count-1) do begin
          total = n_elements(where(folder2 eq names[i]))
          printf, fnum,'<a href="folder_'+strcompress(names[i],/remove_all)+'.html">'
          printf, fnum,names[i]+' ('+strn(total)+' '+word+')</a><br>'
       endfor
    endif
    
    printf,fnum,'</td></tr></table></td></tr>'

    return
    
 end

pro album_footer, fnum

    printf,fnum,'</td></tr></table></body></html>'

    return

 end

pro album

    path = '/home/amarble/public_html/emma/'
    filename = 'images.dat'

; update data file

    pushd,path+'images/'
    all = findfile('*_*.jpg',count=nall)
    readcol,filename,f='A,A,A,A',delimiter='|',image,type,folder,caption,/silent
    image0 = strtrim(image)
    folder0 = strtrim(folder)
    new = 0L
    for i=0,nall-1 do begin
       if (where(strcompress(image,/remove_all) eq $  ; add new entries
                 strcompress(all[i],/remove_all)))[0] eq -1 then begin
          image   = [image,all[i]]
          type    = [type,'']
          folder  = [folder,'']
          caption = [caption,'']
          new = 1L
       endif
    endfor
    album_format, filename, image, type, folder, caption, /backup
    if new then begin
       print & print,'MODIFY images.dat AND HIT ANY KEY WHEN DONE...'
       input='' & read,input
    endif
    readcol,filename,f='A,A,A,A',delimiter='|',image2,type2,folder2,caption2,/silent
    album_format, filename, image2, type2, folder2, caption2
    image2   = strtrim(image2,2)
    type2    = strtrim(type2,2)
    folder2  = strtrim(folder2,2)
    caption2 = strtrim(caption2,2)
    landscape = bytarr(n_elements(image2))
    whland = where(stregex(type2,'\*') ne -1, count)
    if count gt 0 then begin
       landscape[whland] = 1
       for i=0,count-1 do type2[whland[i]] = strsplit(type2[whland[i]],'\*',/regex,/extract)
    endif
    popd

; add date to caption

    dates = strarr(n_elements(image2))
    for i=0,n_elements(image2)-1 do begin
       text = (strsplit(image2[i],'_',/extract))[0]
       dates[i] = num_month(strmid(text,2,2),/reverse)+' '+ $
         strmid(text,4,2)+', 20'+strmid(text,0,2)
       caption2[i] = dates[i]+' - '+caption2[i]
    endfor

; create html file

    pushd, path
    spawn, 'cp index.html index.html.old' ; make backup of webpage
    close,2
    openw,2,path+'index.html'
    
    album_header, 2, 'born May 24, 2003', ncol=1

; output html for "show" images

    wh = where(type2 eq 'show', count)
    if count gt 0 then album_body, 2, image2[wh], caption2[wh]

    album_footer, 2
    close,2

; create folder pages
    
    for i=0,n_elements(names)-1 do begin
       close,3
       openw,3,path+'folder_'+strcompress(names[i],/remove_all)+'.html'
       
       album_header, 3, names[i]
       
; output html for "show" images
       
       wh = where(folder2 eq names[i] and type2 ne 'hide', count)
       if count gt 0 then album_body, 3, image2[wh], caption2[wh], landscape=landscape[wh]

; output html for "listed" folders

       printf,3,'<tr><td colspan=2 align=center valign=center>'
       printf,3,'<a href="index.html">'
       printf,3,'<b>MAIN PAGE</b></a><br>'
       for j=0,n_elements(names)-1 do begin
          if names[j] ne names[i] then begin
             wh = where(folder2 eq names[j] and type2 ne 'hide', count)
             if count gt 1 then word='photos' else word='photo'
             printf,3,'<a href="folder_'+strcompress(names[j],/remove_all)+'.html">'
             printf,3,'<b>'+names[j]+' ('+strn(count)+' '+word+')</b></a><br>'
          endif
       endfor
       printf,3,'</td></tr>'

       album_footer, 3
       close,3
    endfor

; clean out old folders

    if old[0] ne '' then begin
       for i=0,n_elements(old)-1 do begin
          wh = where('folder_'+strcompress(names,/remove_all)+'.html' eq old[i], count)
          if count eq 0 then spawn,'\rm '+old[i]+'*'
       endfor
    endif

    if new then begin

       index = -1
       for i=0,n_elements(image2)-1 do begin
          wh = where(strmatch(image0,image2[i]) eq 1, count)
          if count eq 0 and type2[i] ne 'hide' then index = [index,i]
       endfor
       index = index[1:n_elements(index)-1]
       nfldr = arm_unique(folder2[index])

    endif
    
    return

 end









