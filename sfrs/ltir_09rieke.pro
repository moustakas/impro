function ltir_09rieke, log_l24
; jm09nov23ucsd - convert various input luminosities to L(TIR); deal
; with errors later

    log_ltir = 1.445D + 0.945D*log_l24 ; eq (A6)
    
return, log_ltir
end
    
