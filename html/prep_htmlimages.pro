pro prep_htmlimages, inpath=inpath, outpath=outpath, refsize=refsize
; jm03nov20uofa

    if n_elements(inpath) eq 0L then inpath = './'
    if n_elements(outpath) eq 0L then outpath = './'

    if n_elements(refsize) eq 0L then refsize = 1024L ; reference size
    strefsize = string(refsize,format='(I0)')
    
    pushd, inpath
    imlist = file_search('*',count=imcount)
    outimlist = 'www_'+imlist

    for nim = 0L, imcount-1L do begin

; determine the image dimensions
       
       spawn, ['/usr/bin/identify '+imlist[nim]], iminfo, /sh
       imsplit = strsplit(iminfo,' ',/extract)
       imstrdims = imsplit[2]
       if strmatch(imstrdims,'*x*') eq 0B then begin

          splog, 'Image '+imlist[nim]+' is not a proper image.'
          
       endif else begin

          imdims = long(strsplit(imstrdims,'x',/extract))
          imratio = imdims[0]/float(imdims[1])

; ------------------------------------------------------------------
; this bit of code figures out if the image is portrait or landscape
; and the proper scaling factor
; ------------------------------------------------------------------
          if (imratio gt 1.0) then $
            newdims = fix(refsize*[1.0,1.0/imratio]) else $ ; landscape
            newdims = fix(refsize*[imratio,1.0])            ; portrait
          newstrdims = strjoin(string(newdims,format='(I0)'),'x')
          
          splog, 'Writing '+imlist[nim]+' '+imstrdims+' --> '+newstrdims
          spawn, ['convert '+'-resize '+newstrdims+' '+imlist[nim]+' '+outpath+outimlist[nim]], /sh

; ------------------------------------------------------------------
; we can achieve the same effect using convert
; ------------------------------------------------------------------

;         if (imratio gt 1.0) then prefix = '' else prefix = 'x'
;         splog, 'Writing '+imlist[nim]
;         spawn, ['convert '+'-scale '+prefix+strefsize+' '+imlist[nim]+' '+outpath+outimlist[nim]], /sh
          
       endelse
       
    endfor

    popd
    
return
end    
