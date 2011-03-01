pro free_pointer, pointer
; jm01jan12uofa
; handles arrays of pointers
    
; array input

    np=n_elements(pointer)
    if np gt 1L then begin
       for i=0L, np-1L do free_pointer, pointer[i]
       return
    endif

; if it's a valid pointer, free it, otherwise we don't care

    if ptr_valid(pointer) then ptr_free, pointer
    
return
end
