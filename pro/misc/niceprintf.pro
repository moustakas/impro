;+
; NAME:
;   NICEPRINTF
; PURPOSE:
;   Write to an existing LUN arrays of numbers in a nice column
;   format. 
; MODIFICATION HISTORY:
;   Leonidas Moustakas, ???
;-

pro niceprintf, u0, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11

    on_error, 2
    
if (keyword_set(v11)) then begin 
    for i=0,n_elements(v1)-1 do $
     printf,u0,string(v1[i])+string(v2[i])+string(v3[i])+$
     string(v4[i])+string(v5[i])+string(v6[i])+string(v7[i])+$
     string(v8[i])+string(v9[i])+string(v10[i])+string(v11[i])
    return
endif 
if (keyword_set(v10)) then begin 
    for i=0,n_elements(v1)-1 do $
     printf,u0,string(v1[i])+string(v2[i])+string(v3[i])+$
     string(v4[i])+string(v5[i])+string(v6[i])+string(v7[i])+$
     string(v8[i])+string(v9[i])+string(v10[i])
    return
endif 
if (keyword_set(v9)) then begin 
    for i=0,n_elements(v1)-1 do $
     printf,u0,string(v1[i])+string(v2[i])+string(v3[i])+$
     string(v4[i])+string(v5[i])+string(v6[i])+string(v7[i])+$
     string(v8[i])+string(v9[i])
    return
endif 
if (keyword_set(v8)) then begin 
    for i=0,n_elements(v1)-1 do $
     printf,u0,string(v1[i])+string(v2[i])+string(v3[i])+$
     string(v4[i])+string(v5[i])+string(v6[i])+string(v7[i])+$
     string(v8[i])
    return
endif 
if (keyword_set(v7)) then begin 
    for i=0,n_elements(v1)-1 do $
     printf,u0,string(v1[i])+string(v2[i])+string(v3[i])+$
     string(v4[i])+string(v5[i])+string(v6[i])+string(v7[i])
    return
endif
if (keyword_set(v6)) then begin 
    for i=0,n_elements(v1)-1 do $
     printf,u0,string(v1[i])+string(v2[i])+string(v3[i])+$
     string(v4[i])+string(v5[i])+string(v6[i])
    return
endif
if (keyword_set(v5)) then begin 
    for i=0,n_elements(v1)-1 do $
     printf,u0,string(v1[i])+string(v2[i])+string(v3[i])+$
     string(v4[i])+string(v5[i])
    return
endif
if (keyword_set(v4)) then begin 
    for i=0,n_elements(v1)-1 do $
     printf,u0,string(v1[i])+string(v2[i])+string(v3[i])+$
     string(v4[i])
    return
endif
if (keyword_set(v3)) then begin 
    for i=0,n_elements(v1)-1 do $
     printf,u0,string(v1[i])+string(v2[i])+string(v3[i])
    return
endif
if (keyword_set(v2)) then begin 
   for i=0,n_elements(v1)-1 do $
    printf,u0,string(v1[i])+string(v2[i])
   return
endif
if (keyword_set(v1)) then begin 
   for i=0,n_elements(v1)-1 do $
    printf,u0,string(v1[i])
   return
endif

return
end
