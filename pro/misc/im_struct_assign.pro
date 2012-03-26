;+
; NAME:
;   IM_STRUCT_ASSIGN()
;
; PURPOSE:
;   Simple wrapper on struct_assign to get around the
;   pass-by-reference pain-in-the-neck.
;
; INPUTS: 
;   in - assign data *from* this input structure
;   out - assign data *to* this output structure
;
; KEYWORD PARAMETERS: 
;   nozero - do not zero out data in OUT that does not exist in IN
;     (see STRUCT_ASSIGN)
;
; OUTPUTS: 
;   out - (modified)
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2009 Mar 01, NYU
;
; Copyright (C) 2009, John Moustakas
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

function im_struct_assign, in, out1, nozero=nozero
    out = out1 ; make sure OUT1 is not affected
    struct_assign, in, out, nozero=nozero
return, out
end
