
; jm00june15ucb

; read one column of data very quickly (see readfast, too) (floating point)

pro read1col, filename, array

	on_error, 2

        spawn, ['wc -l '+filename], result
        result = strtrim(result[0], 2)
        lines = long(strmid(result,0,strpos(result,' ')))

        array = fltarr(lines)

        openr, lun1, filename, /get_lun
        readf, lun1, array
        free_lun, lun1

return
end
