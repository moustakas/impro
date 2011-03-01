function uniq_indx, nindx, xmax, sort=sort
; jm06mar13uofa - return an index array of NINDX unique elements (for
;                 subscripting a unique subset of elements of an
;                 array); XMAX is the number of elements in the array
;                 you're trying to subscript

    if (n_elements(nindx) eq 0L) then nindx = 1L else nindx = nindx>1L
    if (n_elements(xmax) eq 0L) then xmax = 1L else xmax = xmax>1L
    
    indx = long(randomu(seed,nindx)*xmax)
    while n_elements(uniq(indx,sort(indx))) ne long(nindx) do begin
       indx = long(randomu(seed,nindx)*xmax)
    endwhile

    if keyword_set(sort) then indx = indx[sort(indx)]

return, indx
end
