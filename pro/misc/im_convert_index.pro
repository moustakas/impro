;+
; NAME:
;       IM_CONVERT_INDEX
; PURPOSE:
;       Converts a one dimensional array index into a vector of indices, whose
;       length is equal to the number of dimensions of the array. This is
;       useful when wanting to know, for instance, what row and column element
;       10034 corresponds to in a 200x150 2-D array. The routine is general and
;       can handle arrays with any number of array dimensions, up to the IDL
;       maximum of 7.
; CALLING SEQUENCE:
;       new_index = CONVERT_INDEX(index, array)
; INPUTS:
;   INDEX
;       A 1 dimensional array index to be converted. IDL can reference
;       multidimensional arrays using a simple 1 dimensional index.
;       Such an index is obtained, for instance from functions such as
;       MAX, MIN and WHERE.
;   ARRAY
;       The array to which this index applies. This routine only uses this
;       parameter to determine the array dimensions, it does not actually use
;       the data stored in the array.
; OUTPUTS:
;   NEW_INDEX
;       The function returns an array of indices, in increasing array index
;       order. NEW_INDEX has a maximum lenght of 7, since IDL arrays are limited
;       to 7 dimensions.
; EXAMPLE:
;       If ARRAY is a 4x3 array and INDEX=7 then this function will return
;       [3,1], since array element 7 (when ARRAY is viewed as a
;       one-dimensional array) is actually column 3, row 1 when ARRAY is viewed
;       as a 2-dimensional array.
; MODIFICATION HISTORY:
;       Created October 1990 by Mark Rivers
;       J. Moustakas, 2005 Jun 22, U of A - vectorized
;-

function im_convert_index, index, array

    nindx = n_elements(index)
    if (nindx eq 0L) or (n_elements(array) eq 0L) then begin
       print, 'Syntax - result = im_convert_index(index,array)'
       return, -1L
    endif

    if (nindx gt 1L) then begin
       for iindx = 0L, nindx-1L do begin
          result1 = im_convert_index(index[iindx],array)
          if (iindx eq 0L) then result = result1 else $
            result = [ [result], [result1] ]
       endfor
       return, result
    endif
    
    nd = size(array)
    ndims = nd[0]

    denom = 1L

    for i = 1L, ndims do denom = denom * nd[i]
    result = lonarr(7)

    for i = ndims, 0L, -1L do begin
       result[i] = index / denom
       index = index mod denom
       denom = denom / nd[i]
    endfor

return, result[0:ndims-1L]
end
