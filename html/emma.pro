pro emma

    path = '/home/amarble/public_html/emma/'
    filename = 'images.dat'
    bgcolor='FFFF93'

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
    whland = where(stregex(type2,'\*') ne -1, count)
    if count gt 0 then begin
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
    
    album_header, 2, bgcolor=bgcolor, $
      title='Emmages Unlimited', $
      folder='bringing images of Emma Heather Marble to the world'

; output html for "show" images

    wh = where(type2 eq 'main', count)
    if count gt 0 then album_main, 2, image2[wh[0]], caption2[wh[0]]

; output html for "listed" folders

    album_folders, 2, folder2, type2, names=names

    old = findfile('folder_*.html')
;
;    unique, folder2, unique=names, indices=indices
;
;    printf,2,'<tr><td colspan=2 align=center valign=center>'
;    for i=0,n_elements(names)-1 do begin
;
;       wh = where(folder2 eq names[i] and type2 ne 'hide', count)
;       if count gt 0 then begin
;          if count gt 1 then word='photos' else word='photo'
;          printf,2,'<a href="folder_'+strcompress(names[i],/remove_all)+'.html">'
;          printf,2,'<b>'+names[i]+' ('+strn(count)+' '+word+')</b></a><br>'
;       endif
;    endfor
;    printf,2,'</td></tr>'

    album_footer, 2
    close,2

; create folder pages
    
    for i=0,n_elements(names)-1 do begin
       close,3
       openw,3,path+'folder_'+strcompress(names[i],/remove_all)+'.html'
       
       album_header, 3, title='Emma Heather Marble', folder=names[i], bgcolor=bgcolor
       
; output html for "show" images
       
       wh = where(folder2 eq names[i] and type2 ne 'hide', count)
       if count gt 0 then album_body, 3, image2[wh], caption2[wh]


; output html for "listed" folders

;       printf,3,'<tr><td colspan=2 align=center valign=center>'
;       printf,3,'<a href="index.html">'
;       printf,3,'<b>MAIN PAGE</b></a><br>'
;       for j=0,n_elements(names)-1 do begin
;          if names[j] ne names[i] then begin
;             wh = where(folder2 eq names[j] and type2 ne 'hide', count)
;             if count gt 1 then word='photos' else word='photo'
;             printf,3,'<a href="folder_'+strcompress(names[j],/remove_all)+'.html">'
;             printf,3,'<b>'+names[j]+' ('+strn(count)+' '+word+')</b></a><br>'
;          endif
;       endfor
;       printf,3,'</td></tr>'

       album_folders, 3, folder2, type2

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
       unique, folder2[index], unique=nfldr

       close,4
       openw,4,'/home/amarble/public_html/emma/files/sigfile'
       printf,4,systime()
       printf,4
       whmain = where(type2[index] eq 'show', mainpage)
       if mainpage gt 0 then begin
          if mainpage eq 1 then text = 'PHOTO.' else text = 'PHOTOS.'
          printf,4,'THE MAIN PAGE FEATURES '+digitize(mainpage,2)+' NEW '+text
       endif
       for i=0,n_elements(nfldr)-1 do begin
          n = n_elements(where(folder2[index] eq nfldr[i]))
          if n eq 1 then text = ' PHOTO HAS ' else text = ' PHOTOS HAVE '
          printf,4,digitize(n,2)+text+'BEEN ADDED TO FOLDER: '+nfldr[i]
       endfor
       printf,4
       printf,4,'This message has been automatically mailed to you by:'
       printf,4
       printf,4,'                  Emmages Unlimited'
       printf,4
       printf,4,'"bringing images of Emma Heather Marble to the world"'
       printf,4
       printf,4,'    http://shapley.as.arizona.edu/~amarble/emma/'
       close,4
       
       if n_elements(index) eq 1 then word='photo' else word='photos'
       spawn,"pine -signature-file=/home/amarble/public_html/emma/files/sigfile "+ $
         "-url 'mailto:amarble@u.arizona.edu,marble@umocm.com,sallye4th@yahoo.com,"+ $
         "lemarble@juno.com,j_andy_neill@yahoo.com,eneill@kc.rr.com,mynaptime@hotmail.com"+ $
         "?subject="+strn(n_elements(index),f='(i)')+"%20new%20"+word+"'"
       
    endif
    
    return

 end









