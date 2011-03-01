function morph2type, morph1, help=help
; jm02may21uofa - written; convert a morphology (string) to RC3
;                 numerical type 
; jm06oct10nyu - generalized a little

    morph = morph1
    morph = repstr(repstr(repstr(morph,'SAB','S'),'SA','S'),'SB','S')
    morph = repstr(repstr(repstr(morph,'IAB','I'),'IA','I'),'IB','I')
    morph = repstr(repstr(morph,'p',''),'?','')

    e_early = where(strmatch(morph,'*E0*',/fold) or strmatch(morph,'*E1*',/fold) or $
      strmatch(morph,'*E2*',/fold),ne_early)
    if (ne_early ne 0L) then morph = repstr(repstr(repstr(morph,'E0','E0'),'E1','E0'),'E2','E0')
    e_late = where((strmatch(morph,'*E5*',/fold) or strmatch(morph,'*E6*',/fold) or $
      strmatch(morph,'*E7*',/fold)) and (strmatch(morph,'*E0*') eq 0),ne_late)
    if (ne_late ne 0L) then morph = repstr(repstr(repstr(morph,'E5','E+'),'E6','E+'),'E7','E+')
    e_mid = where((strmatch(morph,'*E3*',/fold) or $ ; this must be last
      strmatch(morph,'*E4*',/fold)) and (strmatch(morph,'*E0*') eq 0) and (strmatch(morph,'*E+*') eq 0),ne_mid)
    if (ne_mid ne 0L) then morph = repstr(repstr(morph,'E3','E0-1'),'E4','E0-1')
;   e_mid = where((strmatch(morph,'*E*',/fold) or strmatch(morph,'*E3*',/fold) or $ ; this must be last
;     strmatch(morph,'*E4*',/fold)) and (strmatch(morph,'*E0*') eq 0) and (strmatch(morph,'*E+*') eq 0),ne_mid)
;   if (ne_mid ne 0L) then morph = repstr(repstr(repstr(morph,'E3','E0-1'),'E4','E0-1'),'E','E0-1')
    
    if keyword_set(help) then begin

       print, 'The available morphological types are:'
       
       print, ['cE','E0','E0-1','E+','E','S0-','S0',$ ; added E and I
         'S0+','S0/a','Sa','Sab','Sb','Sbc',$
         'Sc','Scd','Sd','Sdm','Sm','I0',$
         'Im','I','cI','Pec']
       return, -99L

    endif
    
    nmorph = n_elements(morph)
    if (nmorph eq 0L) or (size(morph,/type) ne 7L) then begin
       print, 'Syntax - type = morph2type(morph,/help)'
       return, -99L
    endif

    type = lonarr(nmorph)
    for i = 0L, nmorph-1L do begin
    
       case strlowcase(strcompress(morph[i],/remove)) of
          'ce'  : type[i] = -6
          'e0'  : type[i] = -5
          'e0-1': type[i] = -5
          'e'   : type[i] = -5 ; added - assumed "mid-elliptical"
          'e+'  : type[i] = -4
          's0-' : type[i] = -3
          's0'  : type[i] = -2
          's0+' : type[i] = -1
          's0/a': type[i] = 0
          'sa'  : type[i] = 1
          'sab' : type[i] = 2
          'sb'  : type[i] = 3
          'sbc' : type[i] = 4
          'sc'  : type[i] = 5
          'scd' : type[i] = 6
          'sd'  : type[i] = 7
          'sdm' : type[i] = 8
          'sm'  : type[i] = 9
          'i0'  : type[i] = 90
          'i'   : type[i] = 90 ; added - assumed irregular
          'im'  : type[i] = 10
          'ci'  : type[i] = 11
          'pec' : type[i] = 99
          else: begin
             print, 'Unrecognized type '+morph[i]+'...classifying as Pec.'
             type[i] = 99
          end
       endcase

    endfor

return, type
end    
