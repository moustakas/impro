FUNCTION rsex,catalog, use_row=use_row
;+
; NAME:
;       RSEX
;
; PURPOSE:
;       Read in arbitrary SExtractor format catalogs, using native
;       header information & data themselves.  Correctly reads longs,
;       strings, and doubles.  
;
; INPUTS:
;       A SExtractor-format catalog
;
; OPTIONAL INPUTS:
;       use_row - use this row of the catalog to determine the format
;                 of the output data structure (zero-indexed; default
;                 0); this is necessary because occasionally a column
;                 in the first row is a different format (e.g., LONG)
;                 than the rest of the rows (e.g., DOUBLE)
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       Returns a structure with all catalog entries, using field
;       names for tagnames. 
;
; COMMON BLOCKS:
;       None.
;
; RESTRICTIONS:
;       None.
;
; PROCEDURE:
;       Use syntax
;       cs = rsex('catalog.cat')
;
; COMMENTS:
;
; PROCEDURES USED:
;       FINDFILE
;       CREATE_STRUCT
;       VALID_NUM
;       VALID_NUM_ARR
;
; MODIFICATION HISTORY:
; LAM = L. Moustakas
; JM  = J. Moustakas
;       LAM '04may02 - fixed problem case of there being an array of
;                      values based on the last header tag position.
;                      will now work with that case after special
;                      check. 
;       LAM '04feb04 - converted from the old lrsex.pro; adapted to
;                      detect and read longs, strings, and doubles.  
;       JM  '04sep22 - skip comment lines in the header designated
;                      with a double-hash (##)
;       JM  '07aug18 - added USE_ROW optional input
;       JM  '08jul25 - use a variable logical unit number
;       LAM+JM '08dec11 - vetted - no changes made
;-

; Check that an argument has been passed     
     IF n_params() LE 0 THEN BEGIN 
         print,'cat=rsex(catalog)'
         return,-1
     ENDIF 

; See if the catalog file exists at all

     if (file_test(catalog,/regular) eq 0L) then begin ; jm04sep20uofa
         print,'error - the specified catalog does not seem to exist'
         return,-1
     endif

;    jnk = findfile(catalog,count=catexist)
;    if catexist eq -1 then begin
;        message,'error - the specified catalog does not seem to exist'
;        return,-1
;    endif

; Count catalog lines
     spawn,'wc -l < '+catalog, nlines, /sh
     nlines = long(nlines[0]);-1

; Open the catalog to read
     openr,lun,catalog,error=err,/get_lun
     if err ne 0 then begin
         message,'error opening catalog - '+!error_state.msg ; jm04sep20uofa
;        message,'error opening catalog - ',!error_state.msg
         free_lun, lun
         return,-1
     endif

     if (n_elements(use_row) eq 0L) then use_row = 0L ; jm07aug18nyu
     
; Read in the header
     head = ''
     numb = 0l
     cstr=''
     tag = 1
     nheadcomment = 0L ; jm04sep20uofa
     while tag do begin         ; while '#' tag is true
         readf,lun,cstr
         if (strmatch(cstr,'##*') eq 0B) then begin ; skip comment lines
            if strpos(cstr,'#') ne -1 then begin
               head = [head, strupcase((str_sep(strcompress(cstr),' '))[2])]
               numb = [numb, (str_sep(strcompress(cstr),' '))[1]]
            endif else tag=0
         endif else nheadcomment = nheadcomment + 1L ; jm04sep20uofa
     endwhile
     nhead = n_elements(head)   ; number of header entries
     head = head[1:(nhead-1)]
     numb = numb[1:(nhead-1)]
     nhead = n_elements(head)
     free_lun, lun

; Check header for implicit array entries, and expand relevant tagnames
     tothead=''
     for i=0l,nhead-2 do begin
         tothead = [tothead,head[i]]
         mlen = numb[i+1]-numb[i]-1
         if mlen ne 0 then $
           for j=1,mlen do $
           tothead=[tothead,head[i]+strcompress(j,/remove_all)]
     endfor

; The header info so far
     ntothead = n_elements(tothead)
     tothead=[tothead[1:(ntothead-1)],head[(nhead-1)]]

; Need to check one more thing -- whether the last header field is
; actually the first one of an array...

; Read the first data line.  If it has more entries than the total
; number of header fields so far, we've missed an array at the end. 
     junk=strarr(nhead+nheadcomment) ; jm04sep20uofa
     openr,lun,catalog,/get_lun
     readf,lun,junk
     readf,lun,cstr
     free_lun, lun
     nstr = n_elements(str_sep(strcompress(strtrim(cstr,2)),' '))

     if nstr gt ntothead then begin
         mlen = nstr - ntothead
         for j=1,mlen do $
           tothead=[tothead,head[(nhead-1)]+strcompress(j,/remove_all)]
     endif
     ntothead=n_elements(tothead)

; That should do it!  Onwards, to read the data...
     nbody=nlines-nhead-nheadcomment ; jm04sep20uofa
     junk=strarr(nhead+nheadcomment) ; jm04sep20uofa
     openr,lun,catalog
     readf,lun,junk
     for ijunk = 0L, use_row-1L do readf,lun,junk2 ; jm07aug18nyu
     readf,lun,cstr
     free_lun, lun
     body=strarr(nbody)
     openr,lun,catalog
     readf,lun,junk
     readf,lun,body
     free_lun, lun
     bodyslam=strarr(ntothead,nbody)
     for i=0l,nbody-1 do $
;      bodyslam[*,i]=(str_sep(strcompress(strtrim(body[i],2)),' '))[0L:nhead-1L]
       bodyslam[*,i]=(str_sep(strcompress(strtrim(body[i],2)),' '))[0L:ntothead-1L] ; jm04nov22uofa
     carr=str_sep(strcompress(strtrim(cstr,1)),' ')

     tind    = valid_num_arr(carr) ; 0=string 1=number
     tindint = valid_num_arr(carr,/int) ; 0=non-integer 1=integer
     if tind[0] eq 0 then $
       type='' else $
       if tindint[0] eq 1 then $
       type=0l else $
       type=0.d0
     cs=create_struct(tothead[0],type)

     for i=1l,ntothead-1 do begin
         if tind[i] eq 0 then $
           type='' else $
           if tindint[i] eq 1 then $
           type=0l else $
           type=0.d0
         cs=create_struct(cs,tothead[i],type)
     endfor
     cs=replicate(cs,nbody)

     for i=0l,ntothead-1 do begin
        for j=0l,nbody-1 do begin
           cs[j].(i)=bodyslam[i,j]
;          print, bodyslam[i,j], i, j
;          print, i, j
        endfor
;stop
     endfor
     
     return,cs
END 
