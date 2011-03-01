function im_string_stats, arr, type=type, ndecimal=ndecimal, _extra=extra
; jm09mar07nyu - given an array, compute statistics and build a string
; that looks nice as a legend label; for example, in most cases ARR
; will be residuals of some type

    if (n_elements(arr) eq 0L) then return, ''
    if (n_elements(ndecimal) eq 0L) then ndecimal = 2
    if (n_elements(type) eq 0L) then type = 1

    nd = string(ndecimal,format='(I0)')
    ss = im_stats(arr,_extra=extra)
    case type of
       1: str = strtrim(string(ss.median,format='(F12.'+nd+')'),2)+$
         ' ('+strtrim(string(ss.mean,format='(F12.'+nd+')'),2)+'\pm'+$
         strtrim(string(ss.sigma,format='(F12.'+nd+')'),2)+')'
       2: str = strtrim(string(ss.median,format='(F12.'+nd+')'),2)+$
         '\pm'+strtrim(string(ss.sigma,format='(F12.'+nd+')'),2)
       else: message, 'STRING type not supported'
    endcase
    
return, str
end
