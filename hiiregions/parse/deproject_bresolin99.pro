function deproject_bresolin99, diam=diam
; jm05nov08uofa - return the de-projected galactocentric radii for the
;                 HII regions in this paper

    data = mrdfits('raw_bresolin99.fits',1,/silent)
    ned = mrdfits('raw_bresolin99_ned.fits',1,/silent)
    if (n_elements(diam) eq 0L) then diam = mrdfits('raw_bresolin99_diameters.fits',1,/silent)

; compute the inclination angles

    d25_maj = float(diam.rc3_major_axis)
    d25_min = float(diam.rc3_minor_axis)
    pa = diam.rc3_posangle
    
    ratio = d25_min/d25_maj ; "b" divided by "a"
    quantity = sqrt((1.0-ratio^2)/0.96)
    incl = asin(quantity<1)*!radeg

    nhii = n_elements(data)
    radius = fltarr(nhii)-999.0

    result = {$
      id:         0L, $
      galaxy:     '', $
      hii:        '', $
      ra:         '', $
      dec:        '', $
      d25_maj:   0D0, $
      d25_min:   0D0, $
      incl:      0.0, $
      pa:        0.0, $
      radius: -999.0}
    result = replicate(result,nhii)

    result.id = lindgen(nhii)
    result.galaxy = strcompress(data.name,/remove)
    result.hii = strcompress(data.hiiname,/remove)

;   for i = 53, 53 do begin
    for i = 0L, nhii-1L do begin

       match = where((strtrim(data[i].name,2) eq strtrim(ned.nedgalaxy,2)) or $
         (strtrim(data[i].name,2) eq strtrim(ned.galaxy,2)),nmatch)

       result[i].ra = ned[match].ra
       result[i].dec = ned[match].dec
       result[i].d25_maj = d25_maj[match]
       result[i].d25_min = d25_min[match]
       result[i].incl = incl[match]
       result[i].pa = pa[match]

;      if (abs(data[i].raoff) eq 0.0) and (abs(data[i].deoff) eq 0.0) then result[i].radius = 0.0
       
       if (pa[match] gt -900.0) then begin
       
          result[i].radius = im_hiiregion_deproject(ned[match[0]].ra,ned[match[0]].dec,$
            incl[match[0]],pa[match[0]],float(data[i].raoff),float(data[i].deoff))

          result[i].radius = result[i].radius/(result[i].d25_maj/2.0)

       endif

    endfor

    struct_print, result

return, result
end
    
