pro get_element, x, value, position
; jm01jan28uofa
; a generalization of GETELEMENT_VECTOR, this routine will also accept
; an array of positions, and return an array of indices

    position = long(value-value)
    for i = 0L, n_elements(value)-1L do begin
       array_value = min((abs(x-value[i])),temp)
       position[i] = temp
    endfor
      
return
end
