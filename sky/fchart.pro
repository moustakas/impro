pro fchart, objectinfo, date, radius=radius, scoord=scoord, $
            obsname=obsname, datapath=datapath, cleanup=cleanup
; jm02jun14uofa
; generate a finding chart for observations

    nobject = n_elements(objectinfo)
    if nobject eq 0L then begin
       print, 'Syntax - fchart, objectinfo, date, [radius=, scoord=], $'
       print, '           [obsname=], cleanup=cleanup'
       return
    endif

    if n_elements(date) eq 0L then message, 'Variable DATE not defined.'

    if tag_exist(objectinfo,'object') eq 0L then $
      message, 'Structure tag OBJECT does not exist.'
    if tag_exist(objectinfo,'ra') eq 0L then $
      message, 'Structure tag RA does not exist.'
    if tag_exist(objectinfo,'dec') eq 0L then $
      message, 'Structure tag DEC does not exist.'

    if n_elements(radius) eq 0L then radius = 5.0/60.0 ; 10 arcmin [degrees]
    if n_elements(datapath) eq 0L then begin
       spawn, ['pwd'], datapath & datapath = datapath[0]+'/'
    endif

    pushd, datapath

    htmlfile = strlowcase(objectinfo.object)+'.html'
    
; generate the index HTML page

    openw, lun, 'index.html', /get_lun
    printf, lun, '<html>'
    printf, lun, '<head>'
    printf, lun, '<title>Finding Charts</title>'
    printf, lun, '<FRAMESET COLS="20%,*" >'
    printf, lun, '<FRAME SRC="contents.html" NAME="contents">'
    printf, lun, '<FRAME SRC="main.html" NAME="main">'
    printf, lun, '</FRAMESET>'
    printf, lun, '</head>'
    printf, lun, '<body text="wheat" bgcolor="black" link="cyan" vlink="cyan" alink="cyan">'
    printf, lun, '</body>'
    printf, lun, '</html>'
    free_lun, lun
    
; generate the main HTML page

    openw, lun, 'main.html', /get_lun
    printf, lun, '<html>'
    printf, lun, '<body text="wheat" bgcolor="black" link="cyan" vlink="cyan" alink="cyan">'
    printf, lun, '<p>'
    printf, lun, '</body>'
    printf, lun, '</html>'
    free_lun, lun
    
; generate the table of contents HTML page

    openw, lun, 'contents.html', /get_lun
    printf, lun, '<html>'
    printf, lun, '<head>'
    printf, lun, '<title>Finding Chart Index</title>'
    printf, lun, '</head>'
    printf, lun, '<body text="wheat" bgcolor="black" link="cyan" vlink="cyan" alink="cyan">'
    printf, lun, '<h1 align=center>Finding Chart Index</h1>'
    printf, lun, '<table align=center cellpadding=2 border>'
    printf, lun, '<tbody>'
    printf, lun, '<tr align=center><th>Object</th><th>RA</th><th>Dec</th>'

    for k = 0L, nobject-1L do begin

       printf, lun, '<tr>'
       printf, lun, '<td><a href="'+htmlfile[k]+'" target="main">'+objectinfo[k].object+'</td>'+$
         '<td>'+objectinfo[k].ra+'</td><td>'+objectinfo[k].dec+'</td>'
       printf, lun, '</tr>'

    endfor

    printf, lun, '</tbody>'
    printf, lun, '</table>'
    printf, lun, '</body>'
    printf, lun, '</html>'
    free_lun, lun

    for i = 0L, nobject-1L do begin

       info = objectinfo[i]
       
; grab the DSS image for this object, display, and write out

       root = strlowcase(info.object)
       title = strupcase(root)
       imname = root+'.fits'
       radec = info.ra+', '+info.dec
    
       skyview_batch, imname, radec, sfactr=radius, datapath=datapath, imtype='fits'
       
       image = readfits(imname,h,/silent)
       im = scaleimage(image)   ; byte scale the image

       extast, h, astr          ; astrometry structure

       pixscale = abs(astr.cdelt[0])*3600.0 ; pixel scale [arcsec/pixel]

       xsize = 2.0*astr.crpix[0]-1.0 & ysize = 2.0*astr.crpix[1]-1.0
       xaxis = (findgen(xsize)-astr.crpix[0]-1.0)*pixscale/60.0
       yaxis = (findgen(ysize)-astr.crpix[1]-1.0)*pixscale/60.0

       window, 0, xsize=550, ysize=550
       display, im, xaxis, yaxis, charsize=1.2, charthick=2.0, xthick=2.0, ythick=2.0, $
         /xsty, /ysty, xtitle='x (arcmin)', ytitle='y (arcmin)'
    
; write out

       imout = tvrd()
       write_jpeg, root+'.jpg', imout

; generate the airmass plot for this object and write out

       airmass_plots, date, info.ra, info.dec, info.object, $
         obsname=obsname, imjpg=imjpg

       write_jpeg, root+'_airmass.jpg', imjpg

; open up the html file

       openw, lun, root+'.html', /get_lun
       printf, lun, '<html>'
       printf, lun, '<head>'
       printf, lun, '<title>Finding Chart for '+title+'</title>'
       printf, lun, '</head>'
       printf, lun, '<body text="wheat" bgcolor="black" link="cyan" vlink="cyan" alink="cyan">'
       printf, lun, '<center><h1><font size=+3> Finding Chart for <b>'+title+'</b></font></h1></center>'
       printf, lun, '<hr size=+4 rule=+2>'
       printf, lun, '<table border=5 cellspacing=3 cellpadding=2 width=30%>' ; information table
       printf, lun, '<tr><td>RA</td><td align=center>'+info.ra+'</td></tr>'
       printf, lun, '<tr><td>DEC</td><td align=center>'+info.dec+'</td></tr>'
       printf, lun, '</table>'  ; close information table
       printf, lun, '<center>'
       printf, lun, '<table border=5 cellspacing=3 cellpadding=2 width=100%>' ; image table
       printf, lun, '<tr>'
       printf, lun, '<td width=410 align=middle><img src="'+root+'.jpg" width=400 height=400</td>'
       printf, lun, '<td width=410 align=middle><img src="'+root+'_airmass.jpg" width=400 height=400></td>'
       printf, lun, '</tr>'
       printf, lun, '</table>' ; close image table
;      printf, lun, '<p>Page created by <a href="mailto:jmoustakas@as.arizona.edu">John</a>.'
       printf, lun, '</body>'
       printf, lun, '</html>'
       free_lun, lun
       
; cleanup un-necessary images    

    endfor
    popd
    
stop
    
return
end
