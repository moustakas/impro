;+
; NAME: hubble_distance
; 
; PURPOSE: 
;   Compute distances of nearby galaxies using the velocity field model
;   of Burstein et al. 1989
;
; CALLING SEQUENCE:
;    distance = hubble_distance(name, ra, dec, cz, h0= )
;  
; INPUTS:
;   name - string array containing names of galaxies
;   ra   - string array of RA values in hours (either decimal or sexigesimal)   
;   dec  - string array of DEC values (either decimal or sexigesimal)    
;   cz   - observed galaxy velocity in km/s
; 
; OPTIONAL INPUTS:
;   h0   - Hubble constant in km/s/Mpc, defaults to 70
;
; OUTPUTS:
;   Returns model distances in Mpc
;
; COMMENTS: 
;   Calls fortran code "hubble.so" compiled from hubble.f
;
; BUGS:
;   Fortran code sensitive to I/O formatting
;
; REVISION HISTORY:
;   29-April-2003  Written by C. Tremonti, Steward Observatory
;-

function hubble_distance, name, ra, dec, cz, h0 = h0

    if not keyword_set(h0) then h0 = 70.0 ; km/s/Mpc

    ngal = n_elements(ra)
    distance = fltarr(ngal) - 999.0
    tripple_name = ' '

    srcdir = filepath('',root_dir=getenv('IMPRO_DIR'),subdirectory='src')

; loop through galaxies

    for igal = 0, ngal - 1 do begin

       get_coords, radec, instring = ra[igal] + '  ' + dec[igal]
       glactc, radec[0], radec[1], 2000.0, lgal, bgal, 1
       
; call FORTRAN program 

       vcor = 0.0
       tripple = 0.0
       flag = call_external(srcdir+'hubble.so', 'hubble_wrapper_', $
         float(lgal), float(bgal), float(cz[igal]), vcor, tripple, unload=0)

; tripple valued distances are left as -999.0

       if tripple eq float(0) then distance[igal] = vcor / h0 else $
         tripple_name = tripple_name + name[igal]

    endfor  

; print out the data on the tripple valued distances for good measure
    
;   print, ' ' 
;   print, 'GALAXIES WITH TRIPPLE VALUED DISTANCES: ' + tripple_name
 
return, distance 
end

