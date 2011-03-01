function logarithmic_ages, agemin=agemin, agemax=agemax, $
  dlogage=dlogage, silent=silent
; jm03sep12uofa - return a grid of ages logarithmic in age; all ages
;                 in Myr
; jm08apr19nyu - added AGEMIN, AGEMAX optional inputs    

    if n_elements(agemin) eq 0L then agemin = 1.0D
    if n_elements(agemax) eq 0L then agemax = 13500D
    if n_elements(dlogage) eq 0L then dlogage = 0.20D ; logarithmic age bin

    logage = findgen((alog10(agemax)-alog10(agemin))/dlogage+1.0)*dlogage+alog10(agemin)
    age = 10.0^logage
    age = fix(age[uniq(fix(age))])

    splog, agemin, agemax, dlogage
    
    if (not keyword_set(silent)) then begin
       splog, 'Logarithmic age interval in Myr:'
       niceprint, age
;      niceprint, strjoin(string(age,format='(G0.0)'),', ')
    endif
    
return, age
end    
