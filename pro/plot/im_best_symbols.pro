;+
; NAME:
;   IM_BEST_SYMBOLS()
; PURPOSE:
;   Return a priority-sorted list of the 22 "best" SYMCAT() symbols.
; OPTIONAL INPUTS:
;   fill - fill array (defaults to reasonable values)
;   psize - symsize array (defaults to reasonable values)
; OUTPUT:
;   psym - plot symbol to use in PLOT (i.e., psym=psym)
; MODIFICATION HISTORY:
;   J. Moustakas, 2004 Jul 13, U of A
;   jm05jun22uofa - the priority of the symbols was improved 
;-

function im_best_symbols, fill=fill, psize=psize
    psym =  [104,105,106,108,115,122,104,105,106,108,101,115,122,102,109,109,116,119,121,123,123,125]
    fill =  [  0,  1,  0,  1,  1,  1,  1,  0,  1,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  0,  1,  1]
    psize = [1.0,1.7,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0]
return, psym
end
    
