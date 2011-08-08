;+
; NAME:
;   IM_EMPTY_STRUCTURE()
;
; PURPOSE:
;   Empty a data structure.
;
; INPUTS: 
;   oldstruct - input structure
;
; OPTIONAL INPUTS: 
;   ncopies - number of times to copy the output structure
;   empty_value - replace numbers with this value (e.g., -999) 
;   empty_string - replace strings with this string (e.g., '...') 
;
; OUTPUTS: 
;   newstruct - emptied structure
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2005 Apr 07, U of A
;   jm06aug28uofa - added _extra 
;   jm07aug22nyu - added EMPTY_VALUE and EMPTY_STRING optional inputs  
;
; Copyright (C) 2005-2007, John Moustakas
; 
; This program is free software; you can redistribute it and/or modify 
; it under the terms of the GNU General Public License as published by 
; the Free Software Foundation; either version 2 of the License, or
; (at your option) any later version. 
; 
; This program is distributed in the hope that it will be useful, but 
; WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
; General Public License for more details. 
;-

function im_empty_structure, oldstruct, ncopies=ncopies, $
  empty_value=empty_value, empty_string=empty_string, _extra=extra

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
