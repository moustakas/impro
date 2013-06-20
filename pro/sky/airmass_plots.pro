;+
; NAME:
;   AIRMASS_PLOTS
;
; PURPOSE:
;   Plot of airmass versus time for astronomical observations. 
;
; INPUTS:
;   date - date of the observations in FITS standard string
;          format, [YYYY-MM-DD] (e.g., date = '2002-02-10')
;   ra   - right ascension in H:M:S or decimal degrees [NOBJECT] 
;   dec  - declination in D:M:S or decimal degrees [NOBJECT]
;
; OPTIONAL INPUTS:
;   object  - name of object [NOBJECT]
;   obsname - observatory name [default: KPNO]
;   
; KEYWORD PARAMETERS:
;   postscript    - generate postscript output for each object 
;   bigpostscript - generate a single, large postscript file
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; TODO:
;   [1] Add an epoch keyword and to allow precession of
;       coordinates.
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2001 Nov 2, U of A - written
;   jm02feb04uofa - generalized and documented
;   jm04mar01uofa - cleaned up; OBJECT is now an optional input;
;                   gzip postscript output 
;   jm04nov07uofa - added PSNAME and BIGPOSTSCRIPT keywords;
;                   documentation cleaned up
;   jm05mar01uofa - added a check to see if the PS device has been
;                   opened; this allows an airmass postscript file
;                   to be opened by another routine
;   jm09apr21nyu  - checks whether RA and DEC are already in degrees 
;
; Copyright (C) 2001-2002, 2004-2005, 2009, John Moustakas
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

pro airmass_plots, date, ra, dec, object=object, obsname=obsname, $
  minairmass=minairmass, psname=psname, postscript=postscript, $
  bigpostscript=bigpostscript, pdf=pdf, silent=silent

    ndate = n_elements(date)
    if (ndate eq 0L) then begin
       doc_library, 'airmass_plots'
       return
    endif

    nra = n_elements(ra)
    ndec = n_elements(dec)
    if (nra eq 0L) or (ndec eq 0L) then begin
       splog, 'Inputs RA and/or DEC not defined.'
       return
    endif

    if (n_elements(object) eq 0L) then begin
       object = strarr(nra)
       for iobj = 0L, nra-1L do object[iobj] = 'object_'+$
         string(iobj,format='(I3.3)')
    endif
    nobj = n_elements(ra)

    if (nra ne ndec) and (ndec ne nobj) then begin
       splog, 'Incompatible dimensions of RA, DEC, and OBJECT.'
       return
    endif

    if (n_elements(obsname) eq 0L) then obsname = 'KPNO'
    if (keyword_set(silent) eq 0) then $
      splog, 'Selecting observatory '+obsname
    observatory, obsname, obsinfo

    mdy = strsplit(date,'-',/extract) ; [Year,Month,Day]
    
; time vector from 18:00 to 06:00 hours in 10 minute intervals
    jd = timegen(6,12,start=julday(mdy[1],mdy[2],$
      mdy[0],18),unit='minute',step=10)
    jd_ut = timegen(6,12,start=julday(mdy[1],mdy[2],$
      mdy[0],18+obsinfo.tz),unit='minute',step=10)
    jd_ut = reform(jd_ut,n_elements(jd_ut))
;   caldat, jd_ut, month, day, year, hour_ut, minute_ut ; convert to UT time

    jd = timegen(6,12,start=julday(mdy[1],mdy[2],$
      mdy[0],18),unit='minute',step=10)
    jd = reform(jd,n_elements(jd))
;   caldat, jd, month, day, year, hour, minute ; convert to local time

    ct2lst, lst, (360.0-obsinfo.longitude), junk, jd_ut ; convert to LST

    if (n_elements(psname) eq 0) then psname = 'airmass_plots.ps'
    if keyword_set(postscript) or keyword_set(bigpostscript) then $
      postthick = 4 else postthick = 2
    if keyword_set(bigpostscript) then begin
       postscript = 0
       if (keyword_set(silent) eq 0) then splog, 'Writing '+psname
       dfpsplot, psname, /square
    endif

    if (!d.name eq 'PS') then postthick = 4 else begin
       window, 0, xs=600, ys=600
       postthick = 2.0
    endelse

    minairmass = fltarr(nobj)
    for iobj = 0L, nobj-1L do begin
       galaxy = strtrim(object[iobj],2)
       if (size(ra[iobj],/type) eq 7) then begin ; RA [degree]
          radd = 15.0*hms2dec(ra[iobj])
          rastr = ra[iobj]
       endif else begin
          radd = ra[iobj]
          rastr = dec2hms(ra[iobj]/15.0,/colon)
       endelse
       if (size(dec[iobj],/type) eq 7) then begin ; DEC [degree]
          ddec = hms2dec(dec[iobj])
          decstr = dec[iobj]
       endif else begin
          ddec = dec[iobj]
          decstr = dec2hms(dec[iobj],/colon)
       endelse

; calculate sunset time and sunrise time
       ha = lst*15.0D - radd[0]    ; hour angle [degrees]

       airmass = compute_airmass(obsinfo.latitude,ddec[0],ha,zd=zd)
       
       minair = min(airmass,mindx)
       minairmass[iobj] = minair
       best_jd = jd[mindx]
       caldat, best_jd, month, day, $
         year, hour, minute, second
       best_time = strmid(strjoin(strsplit((dec2hms($
         hour+minute/(60D)+second/(3600D)))[0],' ',/extract),':'),0,5)
       thedate = strjoin(strcompress(date,/remove),'-')

       if (minair gt 3.0) then begin
          if (keyword_set(silent) eq 0) then splog, $
            galaxy+' is not observable on '+thedate+$
            ' (minimum airmass = '+string(airmass[mindx],format='(F5.3)')+')'
       endif else begin
          if keyword_set(postscript) then begin
             psfile = strlowcase(strcompress(galaxy,/remove))+'_airmass.ps'
             if (keyword_set(silent) eq 0) then splog, 'Writing '+psfile
             dfpsplot, psfile, /square
          endif 

          pa = compute_parangle(obsinfo.latitude,ddec[0],ha) ; parallactic angle

; generate the plots          
          hmlabel = label_date(date_format=['%H:%I'])
          plot, jd_ut, airmass, xsty=9, ysty=1, yrange=[0,3], $
            ytitle='Airmass', charsize=1.7, charthick=postthick, $
            xthick=postthick, ythick=postthick, thick=postthick, $
            position=[0.16,0.4,0.95,0.9], xtickname = replicate(' ',10), $
            xgridstyle=1, ygridstyle=1, yticklen=1, xticklen=1
          axis, xaxis=1, xsty=1, xtickformat='LABEL_DATE', xtickunits='TIME', $
            xrange=[jd[0],jd[n_elements(jd)-1]], charsize=1.7, charthick=postthick, $
            xthick=postthick, xtitle='Local Time'
          legend, [galaxy,rastr,decstr], /right, /top, box=0, $
            charsize=1.7, charthick=postthick, margin=0
          legend, ['Minimum airmass '+strn(minair,length=4)+' at '+best_time+' Local Time',$
            thedate+' at '+strupcase(obsname)], $
            /left, /bottom, box=0, charsize=1.7, charthick=postthick, margin=0
          
          plot, jd_ut, pa, xtickformat='LABEL_DATE', xsty=9, /ysty, $
            xtickunits='TIME', ytitle='Parallactic Angle', xtitle='UT Time', $
            charsize=1.7, charthick=postthick, xthick=postthick, ythick=postthick, $
            thick=postthick, position=[0.16,0.1,0.95,0.4], /noerase, $
            xgridstyle=1, ygridstyle=1, xticklen=1, yticklen=1

          if keyword_set(postscript) then begin
             dfpsclose
             spawn, 'gzip -f '+psfile, /sh
          endif else begin
             if (nobj gt 1L) and (not keyword_set(bigpostscript)) then begin
                splog, strupcase(galaxy)+' - press any key to continue.'
                cc = get_kbrd(1)
             endif else if (keyword_set(silent) eq 0) then $
               splog, strupcase(galaxy)
          endelse
       endelse
    endfor

    if keyword_set(bigpostscript) then begin
       dfpsclose
       if keyword_set(pdf) then begin
          spawn, 'ps2pdf '+psname+' '+repstr(psname,'.ps','.pdf'), /sh
          rmfile, psname
       endif else spawn, 'gzip -f '+psname, /sh
    endif
    
return
end
