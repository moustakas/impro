function html_structure_parse, instruct, tags=tags, tagformats=tagformats, $
  tagindices=tagindices, blank=blank, keep_error=keep_error, keep_upperlimits=keep_upperlimits
; jm05apr05uofa - written

; take a data structure and parse the entries of interest so that they
; can be incorporated into an HTML table

; TAGS - tag names of interest
; TFORMATS - string formatting codes for each tag
; TAGINDICES - output list of indices in case not every tag in TAGS is
;              in INSTRUCTU

    if (n_elements(instruct) eq 0L) or (n_elements(tags) eq 0L) or (n_elements(tagformats) eq 0L) then begin
       print, 'Syntax - outstruct = html_structure_parse(instruct,tags=,tagformats=,tagindices=)'
       return, -1L
    endif
    
    if (n_elements(blank) eq 0L) then blank = '&nbsp'
    
    nstruct = n_elements(instruct)

    doit = match_string(tags,tag_names(instruct),/exact,index=match,/silent)
    tagindices = where(match ne -1L)
    tags = tags[tagindices]
    tagformats = tagformats[tagindices]
    
    routstruct = struct_trimtags(instruct,select=tags)
    outstruct = struct_trimtags(instruct,select=tags,format=tagformats)
    ntags = n_tags(outstruct)

    for itag = 0L, ntags-1L do begin
       property = reform(strtrim(reform(outstruct.(itag)),2))
       rproperty = routstruct.(itag)
       dims = size(outstruct.(itag),/dimension)
       if (dims[0] eq 2L) then begin
          no = where((property[1,*] lt 0.0),nno,comp=go,ncomp=ngo) ; flag on the error
          if (not keyword_set(keep_error)) then property[1,*] = blank
          if (nno ne 0L) then if (not keyword_set(keep_upperlimits)) then property[0,no] = blank
          outstruct.(itag) = strtrim(reform(property,2,nstruct),2)
       endif else begin
          if (size(rproperty,/type) eq 7L) then begin
             outstruct.(itag) = strtrim(reform(property,1,nstruct),2)
          endif else begin
             no = where((property lt -900.0),nno,comp=go,ncomp=ngo)
             if (nno ne 0L) then property[no] = blank
             outstruct.(itag) = strtrim(reform(property,1,nstruct),2)
          endelse
       endelse
    endfor

return, outstruct
end

