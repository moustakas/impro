;+
; NAME:
;       WSEX
;
; PURPOSE:
;       Write out arbitrary SExtractor format catalogs, using native
;       header information & data themselves.  Correctly reads longs,
;       strings, and doubles.  
;
; INPUTS:
;       A SExtractor-format catalog
;
; OUTPUTS:
;       Prints out a catalog file for the given structure
;
; COMMON BLOCKS:
;       None.
;
; RESTRICTIONS:
;       None.
;
; PROCEDURE:
;       Use syntax
;       wset,catalog,outfile='catalog.cat'
;
; COMMENTS:
;       To do -- 
;       . better error checking, e.g. for existing output filename,
;       etc. 
;
; PROCEDURES USED:
;       VALID_NUM
;
; MODIFICATION HISTORY:
;       L. Moustakas '04feb11 - created
;       J. Moustakas, 09jul08, UCSD - check for arrays in the input
;         structure 
;-

pro wsex,cat,outfile=outfile

     nbody = n_elements(cat)

; Check that an argument has been passed     
     if n_params() le 0 then begin 
         print,'wsex,cat,outfile=''catalog.cat'''
         return
     endif 

     if n_elements(outfile) eq 0 then begin
; herein, assign an output catalog name
         outfile='wsexout.cat'
     endif

; check to see if any of the tags are arrays, and expand the tag names
; accordingly
     tags1  = tag_names(cat)
     ntags1 = n_elements(tags1)
     narray = lonarr(ntags1)
     for itag = 0L, ntags1-1L do begin
        narray[itag] = n_elements(cat[0].(itag))
        if (n_elements(tags) eq 0L) then begin
           if (narray[itag] eq 1L) then tags = tags1[itag] else $
             tags = tags1[itag]+string(lindgen(narray[itag]),format='(I0)')
        endif else begin
           if (narray[itag] eq 1L) then tags = [tags,tags1[itag]] else $
             tags = [tags,tags1[itag]+string(lindgen(narray[itag]),format='(I0)')]
        endelse
     endfor
     ntags = n_elements(tags)

; figure out the format
     formarr = strarr(ntags1)
     for i = 0L, ntags1-1L do begin
        case size((cat[0].(i))[0],/tname) of
           'BYTE'   : formarr[i] = '(I10)'
           'INT'    : formarr[i] = '(I10)'
           'LONG'   : formarr[i] = '(I10)'
           'FLOAT'  : formarr[i] = '(G16.8)'
           'DOUBLE' : formarr[i] = '(E16.8)'
           'STRING' : formarr[i] = '(A20)'
           'UINT'   : formarr[i] = '(I10)'
           'ULONG'  : formarr[i] = '(I10)'
           'LONG64' : formarr[i] = '(I10)'
           'ULONG64': formarr[i] = '(I10)'
           else: begin
              print, 'Unsupported data type.'
              return
           end 
        endcase
     endfor 

; write out     
     openw,u0,outfile, /get_lun
     for i=0l,ntags-1 do $
       printf,u0,'# '+string(i+1,format='(I5)')+' '+strupcase(tags[i])

     for kk = 0L, nbody-1L do begin          
        tstr=''
        for ii = 0L, ntags1-1L do for jj = 0L, narray[ii]-1L do begin
           tstr = tstr+string((cat[kk].(ii))[jj],format=formarr[ii])
        endfor
        printf, u0, tstr
     endfor
     free_lun, u0
     
;    tind    = intarr(nhead)
;    tindint = intarr(nhead)
;
;    for i=0l,nhead-1 do tind[i]    = valid_num(cat[0].(i))
;    for i=0l,nhead-1 do tindint[i] = valid_num(cat[0].(i),/int)
;
;    formstr='(a20)'
;    formlon='(i10)'
;    formdbl='(g14.7)' ; jm07aug18nyu
;    formarr=strarr(nhead)
;    for i=0l,nhead-1 do $
;      if tind[i] eq 0 then $
;      formarr[i]=formstr else $
;      if tindint[i] eq 1 then $
;      formarr[i]=formlon else $
;      formarr[i]=formdbl
;
;    for k=0l,nbody-1 do begin 
;        tstr=''
;        for i=0l,nhead-1 do tstr=tstr+string(cat[k].(i),format=formarr[i])
;        printf,u0,tstr
;    endfor
;    close,u0,/force

return
end

