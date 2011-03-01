pro look_for_crap_in_gandalf, spec

    for ii = 0, n_elements(spec)-1 do begin
       lines = strtrim(spec[0].linename,2)
       for ll = 0, n_elements(lines)-1 do begin
          indx1 = tag_indx(spec[0],lines[ll])
          indx2 = tag_indx(spec[0],lines[ll]+'_amp')
          indx3 = tag_indx(spec[0],lines[ll]+'_sigma')
          indx4 = tag_indx(spec[0],lines[ll]+'_linez')
; if flux>0 and flux_err=0 then bad
          if ((spec[ii].(indx1))[0] gt 0.0) and $
            ((spec[ii].(indx1))[1] eq 0.0) then begin
             splog, 'FLUX_ERR=0 on index '+strtrim(ii,2)+', line '+lines[ll]
;            print, (spec[ii].(indx2))[0]/(spec[ii].(indx2))[1]
          endif
; if flux=0 and amp>0 then bad
          if ((spec[ii].(indx1))[0] eq 0.0) and $
            ((spec[ii].(indx2))[0] gt 0.0) then begin
             splog, lines[ll], ii & stop
          endif
; if amp>0 and amp_err=0 then bad
          if ((spec[ii].(indx2))[0] gt 0.0) and $
            ((spec[ii].(indx2))[1] eq 0.0) then begin
             splog, 'AMP_ERR=0 on index '+strtrim(ii,2)+', line '+lines[ll]
          endif
; if z>0 and z_err=0 then bad
          if ((spec[ii].(indx4))[0] gt 0.0) and $
            ((spec[ii].(indx4))[1] eq 0.0) then begin
             splog, 'ZERR=0 on index '+strtrim(ii,2)+', line '+lines[ll]
;            print, (spec[ii].(indx1))[0]/(spec[ii].(indx1))[1], $
;              (spec[ii].(indx2))[0]/(spec[ii].(indx2))[1]
          endif
       endfor
    endfor

return
end
