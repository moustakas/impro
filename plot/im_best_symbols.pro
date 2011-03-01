function im_best_symbols, fill=fill, psize=psize
; jm04jul13uofa - return the "best" symbols to use in a plot from
;                 IM_SYMBOLS
; jm05jun22uofa - the priority of the symbols was improved    

    psym =  [104,105,106,108,115,122,104,105,106,108,101,115,122,102,109,109,116,119,121,123,123,125]
    fill =  [  0,  1,  0,  1,  1,  1,  1,  0,  1,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  0,  1,  1]
    psize = [1.0,1.7,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0]
    
return, psym
end
    
