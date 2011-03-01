function type2morph, type
; jm02apr21uofa
; convert RC3 numeric morphological type to string type for plot
; labeling 

    ntype = n_elements(type)
    if ntype eq 0L then message, 'Variable TYPE not defined.'

    type = long(type)
    morph = strarr(ntype)
    
    for i = 0L, ntype-1L do begin

       case type[i] of
          -6: morph[i] = 'cE' 
          -5: morph[i] = 'E0' 
          -4: morph[i] = 'E+' 
          -3: morph[i] = 'S0-'
          -2: morph[i] = 'S0' 
          -1: morph[i] = 'S0+'
          0: morph[i] = 'S0/a'
          1: morph[i] = 'Sa'  
          2: morph[i] = 'Sab' 
          3: morph[i] = 'Sb'  
          4: morph[i] = 'Sbc' 
          5: morph[i] = 'Sc'  
          6: morph[i] = 'Scd' 
          7: morph[i] = 'Sd'  
          8: morph[i] = 'Sdm' 
          9: morph[i] = 'Sm'  
          90: morph[i] = 'I0' 
          10: morph[i] = 'Im' 
          11: morph[i] = 'cI' 
          99: morph[i] = 'Pec'
          else: begin
             print, 'Unrecognized type '+type[i]+'...classifying as Pec.'
             morph[i] = 'pec'
          endelse
       endcase 
       
    endfor 

return, morph
end
