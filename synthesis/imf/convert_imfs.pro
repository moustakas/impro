;+
; NAME:
;       CONVERT_IMFS()
;
; PURPOSE:
;       Return the conversion between a variety of initial-mass
;       functions.  
;
; CALLING SEQUENCE:
;
;
; INPUTS:
;
;
; OPTIONAL INPUTS:
;
;
; KEYWORD PARAMETERS:
;
;
; OUTPUTS:
;
;
; OPTIONAL OUTPUTS:
;
;
; COMMON BLOCKS:
;
;
; PROCEDURES USED:
;
;
; COMMENTS:
;
;
; EXAMPLES:
;
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2005 Mar 22, U of A
;-

function convert_imfs, minmass=minmass, maxmass=maxmass

    if (n_elements(minmass) eq 0L) then minmass = 0.01
    if (n_elements(maxmass) eq 0L) then maxmass = 150.0

; initialize the mass vector and the fiducial Salpeter IMF     
    
    dlogmass = 0.05
    nmass = round(1.0 + (alog10(maxmass)-alog10(minmass))/dlogmass)
    logmass = alog10(minmass) + dlogmass*findgen(nmass)
    mass = 10^logmass

    x = -1.35
    salpeter = imf(mass,x,[0.09999,100.0])

; Bell & de Jong "diet" Salpeter

    x = -1.35
    bell = imf(mass,x,[0.09999,125.0])
    

stop       
    
    
return, constant
end
