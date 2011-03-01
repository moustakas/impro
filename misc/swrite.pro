; save wrapper
; jm00feb20ucb

pro swrite, data, filename

    if n_params() eq 0L then begin
       print, 'Syntax - swrite, data, filename'
       return
    endif

    save, data, filename=filename

return
end

