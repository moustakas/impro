function im_d2dmod, dist, err_dist=err_dist, err_dmod=err_dmod
; jm08jun06nyu - convert a distance and error (in Mpc) to a distance
;                modulus and error
    
    ndist = n_elements(dist)
    if (ndist eq 0L) then begin
       doc_library, 'im_d2dmod'
       return, -1L
    endif
    
    dmod = 5.0*alog10(dist)+25.0
    if arg_present(err_dmod) then begin
       if (n_elements(err_dist) eq 0L) then err_dist = dist*0.0 else begin
          if (n_elements(err_dist) ne ndist) then begin
             print, 'Dimensions of ERR_DIST and DIST must agree'
             return, -1L
          endif
       endelse
       err_dmod = 5.0*err_dist/dist/alog(10.0)
    endif 
    
return, dmod
end
