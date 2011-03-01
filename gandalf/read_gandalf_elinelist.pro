function read_gandalf_elinelist, linefile, actionfit=actionfit, $
  fitlines=fitlines, nlines=nlines
; jm09nov12ucsd - read a GANDALF-format emission-line file

    if (file_test(linefile) eq 0) then begin
       splog, 'Emission-line file '+linefile+' not found'
       return, -1
    endif

    readcol, linefile, eml_i, eml_name, eml_lambda, eml_action, $
      eml_kind, eml_a, eml_v, eml_s, eml_fit, eml_fit_iter2, skip=2, $
      format='I,A,F,A,A,F,F,F,A,A', comment='#', /silent
    nlines = n_elements(eml_i)

    names = ['i','name','lambda','action','kind','a','v','s','fit','fit_iter2']
    values = ['0',"''",'0.0',"''","''",'0.0','0.0','0.0',"''","''"]
    linepars = mrd_struct(names,values,nlines)

    linepars.i = eml_i
    linepars.name = eml_name
    linepars.lambda = eml_lambda
    linepars.action = eml_action
    linepars.kind = eml_kind
    linepars.a = eml_a
    linepars.v = eml_v
    linepars.s = eml_s
    linepars.fit = eml_fit
    linepars.fit_iter2 = eml_fit_iter2
    
;   linepars = create_struct('i',eml_i, 'name', eml_name, 'lambda',$
;     eml_lambda, 'action', eml_action, 'kind', eml_kind, 'a', $
;     eml_a,'v', eml_v,'s', eml_s,'fit', eml_fit)

; if /ACTIONFIT then change lines flagged as 'm' (mask) to 'f' (fit)
; and also optionally return just the lines being fitted
    if keyword_set(actionfit) then begin
       these = where((linepars.action eq 'm') and (linepars.name ne 'sky') and $
         (linepars.name ne 'NaI'),nthese)
       if (nthese ne 0) then begin
          fitlines = linepars[these]
          fitlines.action = 'f'
;         linepars.action[these] = 'f'
;         fitlines = create_struct('i', linepars.i[these], 'name', $
;           linepars.name[these], 'lambda', linepars.lambda[these], $
;           'action', linepars.action[these], 'kind', linepars.kind[these], $
;           'a', linepars.a[these], 'v', linepars.v[these],'s', $
;           linepars.s[these], 'fit', linepars.fit[these])
       endif
    endif
    
    junk = where(linepars.kind eq 'l',nlines)

return, linepars
end
    
