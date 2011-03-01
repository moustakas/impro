function im_empty_structure, oldstruct, ncopies=ncopies, $
  empty_value=empty_value, empty_string=empty_string, _extra=extra
; jm05apr07uofa - empty a data structure; use Schlegel's trick
; jm06aug28uofa - added _extra 
; jm07aug22nyu - added EMPTY_VALUE and EMPTY_STRING optional inputs 

    newstruct = oldstruct[0]
    struct_assign, {junk: 0}, newstruct

    if (n_elements(empty_value) ne 0L) then $
      for itag = 0L, n_tags(newstruct)-1L do $
        if size(newstruct.(itag),/tname) ne 'STRING' then $
          newstruct.(itag) = empty_value
    if (n_elements(empty_string) ne 0L) then $
      for itag = 0L, n_tags(newstruct)-1L do $
        if size(newstruct.(itag),/tname) eq 'STRING' then $
          newstruct.(itag) = empty_string
    
    newstruct = im_struct_trimtags(newstruct,_extra=extra)
    
    if (n_elements(ncopies) ne 0L) then begin
       if (ncopies gt 1L) then newstruct = replicate(newstruct,ncopies)
    endif

return, newstruct
end    
