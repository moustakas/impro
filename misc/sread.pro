; restore wrapper
; jm00feb20ucb

function sread, filename

    if n_params() eq 0L then begin
       print, 'Syntax - data = sread(filename)'
       return, -1
    endif

    restore, filename

return, data
end

