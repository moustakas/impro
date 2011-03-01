function rnd,invalue,round,up=up,down=down,out=out,in=in
;+
; NAME:
;	RND
;
; PURPOSE:
;	This function will round the given value up/down/in or out using the
;	second parameter as the modulo (default is 1)
;
; CATEGORY:
;	utilities
;
; CALLING SEQUENCE:
;	rounded_value = RND(value [,round])
;
; INPUTS:
;	VALUE = value to round
;
;   OPTIONAL INPUT PARAMETERS:
;	ROUND = option parameter to round to nearest of. (Can be non-integer)
;  
;   KEYWORD PARAMETERS:
;	UP   = if set then round VALUE up (towards greater +ve values)
;	DOWN = if set then round VALUE down (towards greater -ve values)
;	OUT  = if set then round VALUE out (away from zero)
;	IN   = if set then round VALUE in (towards zero)
;
; OUTPUTS:
;	ROUNDED_VALUE = value rounded to the nearest ROUND (default is 1)
;
; COMMON BLOCKS:
;	none.
; SIDE EFFECTS:
;	none.
; MODIFICATION HISTORY:
;	Written by: Trevor Harris, Physics Dept., University of Adelaide,
;		July, 1990.
;
;-

    on_error,2
    
    value = double([invalue])
    fnctn = value*0.0
    nonzero = where(value ne 0,count)

    if (count gt 0) then begin

       if (n_params() lt 2) then round = 1.0
       round = abs(round)
       val=value[nonzero]+round*0.5

       if (keyword_set(up)) then $
         val = value[nonzero]+round*0.999999 
       if (keyword_set(down)) then $
         val = float(long((value[nonzero]+round)/float(round))*round)-round*0.999999
       if (keyword_set(out)) then $
         val = value[nonzero]+float(value[nonzero])/abs(value[nonzero])*round*0.999999
       if (keyword_set(in)) then $
         val = value[nonzero]-float(value[nonzero])/abs(value[nonzero])*round*0.000001

       fnctn[nonzero] = float(long(val/float(round))*round)

    endif 

    sz = size(invalue)
    if (sz(0) eq 0) then fnctn = fnctn[0] ;data originally a scalar

return, fnctn
end
