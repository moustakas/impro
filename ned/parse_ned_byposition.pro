pro parse_ned_byposition, nedfile, outfile=outfile, outpath=outpath
; jm02may24uofa
; this routine will parse a ned text file returned by searching near a
; sky position
    
    if n_elements(nedfile) eq 0L then begin
       print, 'Syntax - parse_ned_byposition, nedfile, outfile=outfile, outpath=outpath'
       return
    endif
    
    cspeed = 2.99792458E5 ; [km/s]

; read the whole file and crop

    if file_test(nedfile) eq 0L then begin
       print, 'File '+nedfile+' not found.'
       return
    endif

    nedata = djs_readilines(nedfile)
    lstart = (where(strmatch(nedata,'*SEARCH RESULTS*') eq 1B))[0]
    lend = (where(strmatch(nedata,'*all the objects found*') eq 1B))[0]

    nedata = nedata[lstart:lend]
    
; now find out how many objects are in the file by doing a search on BYNAME

    objname = where(strmatch(nedata,'NEARPOSN *') eq 1B,nobj)
    print, 'There are '+strn(nobj)+' objects in '+nedfile+'.'

    ginfo = strarr(8,nobj)
    data = create_struct('galaxy', '', 'ra', '', 'dec', '', $
                         'index', -99L, 'cz', -99.0, 'morph', '', 'mag', -99.0, $
                         'dmaj', -99.0, 'dmin', -99.0)
    data = replicate(data,nobj)

    for i = 0L, nobj-1L do begin

       if strmatch(nedata[objname[i]+1L],'*No object is found*') ne 1B then begin

          line = strjoin(nedata[objname[i]+4:objname[i]+5])
       
          ginfo[0,i] = strn(strmid(line,7,16)+strmid(line,89,12)) ; NED galaxy name
          ginfo[1,i] = strn(strmid(line,24,11)) ; ra
          ginfo[2,i] = strn(strmid(line,36,10)) ; dec

          z = strmid(line,46,9)
          if strtrim(z) ne '' then begin
             if (abs(float(z)) gt 1.0) then cz = z else cz = cspeed*float(z) 
          endif else cz = ''
          ginfo[3,i] = strn(cz) ; cz (km/s)
          
          ginfo[4,i] = strcompress(strmid(line,60,20),/remove) ; morphology
          ginfo[5,i] = strn(strmid(line,140,8)) ; magnitude
          ginfo[6,i] = strn(strmid(line,148,5)) ; major axis diameter [arcmin]
          ginfo[7,i] = strn(strmid(line,155,5)) ; minor axis diameter [arcmin]

; parse the ra and dec
       
          ginfo[1,i] = strjoin(strsplit(ginfo[1,i],'hms',/extr),':')
          ginfo[2,i] = strjoin(strsplit(ginfo[2,i],'dms',/extr),':')

       endif

    endfor

; fill in blanks

    cz = strcompress(reform(ginfo[3,*]),/remove)
    w = where(cz eq '',nw)
    cz[w] = '-99.0'
    
; fill the data structure

    data.galaxy = strcompress(reform(ginfo[0,*]),/remove)
    data.ra     = strcompress(reform(ginfo[1,*]),/remove)
    data.dec    = strcompress(reform(ginfo[2,*]),/remove)
    data.cz     = float(cz)
    data.morph  = strcompress(reform(ginfo[4,*]),/remove)
    indx = where(strcompress(ginfo[5,*],/remove) ne '') & data[indx].mag = reform(float(ginfo[5,indx]))
    indx = where(strcompress(ginfo[6,*],/remove) ne '') & data[indx].dmaj = reform(float(ginfo[6,indx]))
    indx = where(strcompress(ginfo[7,*],/remove) ne '') & data[indx].dmin = reform(float(ginfo[7,indx]))

    good = where(data.galaxy ne '',ngood)
    if ngood ne 0L then data = data[good]

    data.index = good
    
; sort by RA

;   srtra = sort(im_hms2dec_arr(data.ra))
;   data = data[srtra]
    
; write a binary fits table
    
    if not keyword_set(outpath) then outpath = cwd()
    if not keyword_set(outfile) then outfile = 'outfile.fits'
;   if not keyword_set(outfile) then outfile = 'outfile.txt'; nedfile+'.parsed'
    print, 'Writing '+outpath+outfile+'.'

    mwrfits, data, outpath+outfile, /create
    
;   space = '  |  '
;   openw, lun, outfile, /get_lun
;   for j = 0L, nobj-1L do begin
;      outline = ginfo[*,j]
;      printf, lun, outline[0], space, outline[1], space, outline[2], space, outline[3], $
;        space, outline[4], space, outline[5], space, outline[6], space, outline[7], $
;        format='(A28,A5,A10,A5,A10,A5,A7,A5,A20,A5,A6,A5,A4,A5,A4)'
;   endfor
;   free_lun, lun
    
return
end
