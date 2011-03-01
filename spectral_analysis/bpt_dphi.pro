function bpt_dphi, niiha, oiiihb, phi=phi
; jm08feb24nyu - based on a routine by C. Tremonti

    if (n_elements(niiha) eq 0L) or (n_elements(oiiihb) eq 0L) then begin
       doc_library, 'bpt_dphi'
       return, -1L
    endif
    
    x = niiha + 0.45
    y = oiiihb + 0.5
    d = sqrt(x^2 + y^2)
    phi = atan(x/y)*!radeg

return, d
end
