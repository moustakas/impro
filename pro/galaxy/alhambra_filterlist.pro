;+
; NAME:
;   ALHAMBRA_FILTERLIST()
; PURPOSE:
;   Return the ALHAMBRA filters.
; MODIFICATION HISTORY:
;   J. Moustakas, 2014 Jul 6, Siena
;-

function alhambra_filterlist, shortfilt=shortfilt
    shortfilt = [$
      'f365w',$ 
      'f396w',$ 
      'f427w',$ 
      'f458w',$ 
      'f489w',$ 
      'f520w',$ 
      'f551w',$ 
      'f582w',$ 
      'f613w',$ 
      'f644w',$ 
      'f675w',$ 
      'f706w',$ 
      'f737w',$ 
      'f768w',$ 
      'f799w',$ 
      'f830w',$ 
      'f861w',$ 
      'f892w',$ 
      'f923w',$ 
      'f954w',$ 
      'H',$ 
      'J',$ 
      'KS']      
    return, 'alhambra_'+shortfilt+'.par' 
end
