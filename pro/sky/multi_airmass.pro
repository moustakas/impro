;+
; NAME:
;   MULTI_AIRMASS()
;
; PURPOSE:
;   Generate a multi-page plot of airmass and PA versus time. 
;
; INPUTS: 
;   date - date of the observation in FITS format [YYYY,MM,DD] (e.g.,
;     '2002-02-10') 
;   ra, dec - celestial coordinates of object(s)
;
; OPTIONAL INPUTS: 
;   object - object name
;   obsname - observatory name (default KPNO)
;   psname - postscript file name 
;
; KEYWORD PARAMETERS: 
;
; OUTPUTS: 
;
; MODIFICATION HISTORY:
;   J. Moustakas ???
;
; Copyright (C) 20??, John Moustakas
; 
; This program is free software; you can redistribute it and/or modify 
; it under the terms of the GNU General Public License as published by 
; the Free Software Foundation; either version 2 of the License, or
; (at your option) any later version. 
; 
; This program is distributed in the hope that it will be useful, but 
; WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
; General Public License for more details. 
;-

pro multi_airmass, date, ra, dec, object=object, obsname=obsname, psname=psname

    ndate = n_elements(date)
    nra = n_elements(ra)
    ndec = n_elements(dec)

    if (nra eq 0L) or (ndec eq 0L) then begin
       splog, 'Inputs RA and/or DEC not defined.'
       return
    endif

    for ii = 0, ndec -1 do $
       if strmid(dec[ii], 0, 1) ne '-' and $
          strmid(dec[ii], 0, 1) ne '+' then dec[ii] = '+' + dec[ii]

    if (n_elements(object) eq 0L) then object = ra + dec
    nobj = n_elements(object)

    if (nra ne ndec) and (ndec ne nobj) then begin
       splog, 'Incompatible dimensions of RA, DEC, and OBJECT.'
       return
    endif

    ; Set up plotting

    if keyword_set(psname) then begin
        splog, 'Writing '+psname+'.'
        dfpsplot, psname, /square, /color
    endif

    if (!d.name eq 'PS') then postthick = 6.0 else begin
       window, 0, xs=600, ys=600
       postthick = 2.0
    endelse

    colors = djs_icolor(['red', 'orange', 'green', 'dark blue', $
             'navy', 'dark magenta', 'purple', 'dark cyan', 'light red', $
             'pink', 'light green', 'light yellow', 'blue'])

    lines = [0, indgen(nobj - 1) + 2]


    if (n_elements(obsname) eq 0L) then obsname = 'KPNO'

;   splog, 'Selecting observatory '+obsname+'.'
    observatory, obsname, obsinfo

    mdy = strsplit(date,'-',/extract) ; [Year,Month,Day]
    
; time vector from 18:00 to 06:00 hours in 10 minute intervals
    
    jd_ut = timegen(6,12,start=julday(mdy[1],mdy[2],mdy[0],18+obsinfo.tz),$
      unit='minute',step=10)
    jd_ut = reform(jd_ut,n_elements(jd_ut))

    jd = timegen(6,12,start=julday(mdy[1],mdy[2],mdy[0],18),$
      unit='minute',step=10)
    jd = reform(jd,n_elements(jd))

    ct2lst, lst, (360.0-obsinfo.longitude), junk, jd_ut ; convert to LST

    airmass = fltarr(n_elements(jd), nobj)
    pa = fltarr(n_elements(jd), nobj)
    galaxy = strarr(nobj)
    ;hmlabel = label_date(date_format=['%H:%I'])

    for iobj = 0L, nobj-1L do begin

       galaxy[iobj] = strtrim(object[iobj],2)
       radd = 15.0*hms2dec(ra[iobj]) & radd = radd[0] ; RA [degree]
       ddec = hms2dec(dec[iobj]) & ddec = ddec[0]     ; DEC [degree]

; calculate sunset time and sunrise time

       ha = lst*15.0D - radd    ; hour angle [degrees]

       airmass[*, iobj] = compute_airmass(obsinfo.latitude,ddec,ha,zd=zd)

       ; parallactic angle       
       pa[*, iobj] = compute_parangle(obsinfo.latitude,ddec,ha) 

    endfor
; generate the plots          
          
    plot, jd_ut, airmass[*,0], xsty=11, /ysty, yrange=[0,3], $
          ytitle='Airmass', charsize=1.7, charthick=postthick, $
          xthick=postthick, ythick=postthick, thick=postthick, $
          position=[0.16,0.5,0.95,0.9], xtickname = replicate(' ',20), $
          /nodata, xminor = 8

    for iobj=0L, nobj-1L do oplot, jd_ut, airmass[*,iobj], $
        linestyle=lines[iobj], color=colors[iobj], thick=postthick

    axis, xaxis=1, /xsty, xtickformat='LABEL_DATE', xtickunits='TIME', $
          xrange=[jd[0],jd[n_elements(jd)-1]], charsize=1.7, $
          charthick=postthick, xthick=postthick, xtitle='Local Time', $
          xminor = 8

    oplot, [0, 2e7], [1, 1], thick=postthick, linestyle=1
    oplot, [0, 2e7], [2, 2], thick=postthick, linestyle=1

    xyouts, jd_ut[1], 0.3, obsname + ': ' + date, $
            charthick=postthick, charsize=1.7

    legend, [strupcase(galaxy)], /right, /bot, box=0, charsize=1.2, $
            charthick=postthick, thick=postthick, $
            linestyle = lines, color = colors
            
    plot, jd_ut, pa, xtickformat='LABEL_DATE', xsty=11, /ysty, $
          xtickunits='TIME', ytitle='Parallactic Angle', xtitle='UT Time', $
          charsize=1.7, charthick=postthick, $
          xthick=postthick, ythick=postthick, thick=postthick, $
          position=[0.16,0.1,0.95,0.5], /noerase, /nodata, xminor=8
          ;xgridstyle=1, ygridstyle=1, xticklen=1, yticklen=1, /nodata

    for iobj=0L, nobj-1L do  oplot, jd_ut, pa[*,iobj], $
        linestyle=lines[iobj], color=colors[iobj], thick=postthick

    if keyword_set(psname) then dfpsclose
    
return

end
