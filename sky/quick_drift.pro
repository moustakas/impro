;+
; NAME: 
;      QUICK_DRIFT
;
; PURPOSE:
;      Create a PS file showing the DSS image of a galaxy, the drift
;      boxes, and all the necessary driftscan parameters (scan rates, etc)
;
; CALLING SEQUENCE:
;      quick_drift, galaxy, ra, dec, scanlen, exptime, scan_pa=, $
;        offset_along=, offset_perp=, nscan=, date=, psname=, imsize=, $
;        notes=
; INPUTS:
;      galaxy - galaxy name
;      ra -- right ascension [HH:MM:SS.S]
;      dec -- declination [DD:MM:SS.S]
;      scanlength -- length of driftscan in arcseconds
;      exptime -- exposure time in minutes
;
; OPTIONAL INPUTS:
;      scan_pa -- position angle of the slit, default = 90 degrees
;      offset_along -- offset in arcsec along the slit, default = 0
;        (if this parameter is a vector, multiple slits will be used)
;      offset_perp -- offset in arcsec perpendicular to the slit, default = 0
;      nscan -- number of scans
;      date - date of the observations in YYYY-MM-DD format, defaults to today
;      psname -- name of postscript output, defaults to 'quick_drift.ps'
;      imsize -- size in arcmin of DSS image returned, default=9'
;      notes -- comments to appear in lower right of plot 
;
; EXAMPLE:
;      quick_drift, 'UGC07321', '12:17:34.0', '22:32:23', 40, 10, $
;         scan_pa = 81, offset_along = [-90,90], notes='11Mpc high priority'
;
;-
;------------------------------------------------------------------------------


pro quick_drift, galaxy, ra, dec, scanlen, exptime, scan_pa = scan_pa, $
    offset_along = offset_along, offset_perp = offset_perp, nscan = nscan, $
    date=date, psname = psname, imsize=imsize, notes=notes

;------------------
; Set up defaults

if not keyword_set(scan_pa) then scan_pa = 90
if not keyword_set(offset_perp) then offset_perp = 0
if not keyword_set(offset_along) then offset_along = 0
nslits = n_elements(offset_along)
if not keyword_set(nscan) then nscan=10
if not keyword_set(date) then spawn, 'date +%F', date
if not keyword_set(notes) then notes = ''
if not keyword_set(imsize) then imsize = 9.0
if not keyword_set(psname) then psname = 'quick_drift.ps'

;---------------------------------
; Compute parameters of the drift

for ii = 0, nslits -1 do begin
  scan1 = scan_rates(long(scan_pa), dec, long(scanlen), long(exptime*60), $
       offset_along_slit=offset_along[ii], offset_perp_slit=offset_perp, $ 
       nscan=nscan, object=galaxy, verbose=0)
  if ii eq 0 then scan = scan1 else scan = [scan, scan1]
endfor

driftlabel = [ $
  'RA RATE = ' + string(scan[0].ra_rate, format='(F6.3)'), $
  'DEC_RATE = ' + string(scan[0].dec_rate, format='(F6.3)'), $
  'SCANTIME = ' + string(scan[0].scantime, format='(I4)'), $
  'EXPTIME = ' + string(scan[0].exptime, format='(I5)'), $
  'NSCAN = ' + string(scan[0].nscan, format='(I3)'), $
  'RA OFFSET  = [' +strjoin(string(scan.ra_offset, form='(F5.1)'), ',') + ']', $
  'DEC OFFSET = [' +strjoin(string(scan.dec_offset, form='(F5.1)'), ',') + ']']


;---------------------
; Query the DSS for an image

dra = 15.0D*im_hms2dec(ra)
ddec = im_hms2dec(dec)
querydss, [dra,ddec], dssimage, hdss, imsize=imsize, survey='2r', /ned
;writefits, galaxy+'_dss.fits', dssimage, hdss, /create
;image = readfits(galaxy+'_dss.fits',hdss,/silent)

;-----------------------
; Determine the image size in physical and WCS coords

extast, hdss, astr
imsize = size(dssimage,/dimension)
xsize = imsize[0] & xcen = xsize/2.0  
ysize = imsize[1] & ycen = ysize/2.0
xpixscale = (astr.pltscl*1E-3*astr.xsz)/60 ; [arcmin/pixel]
ypixscale = (astr.pltscl*1E-3*astr.ysz)/60 ; [arcmin/pixel]
xaxis = (findgen(xsize)-xcen)*xpixscale ; centered on the image [arcsec]
yaxis = (findgen(ysize)-ycen)*ypixscale ; centered on the image [arcse

;----------------------
; Scale the image for display

topvalue = !d.table_size-2
minvalue = -40L
xbox = fix(xsize*0.20) & ybox = fix(ysize*0.20)
djs_iterstat, dssimage[xcen-xbox:xcen+xbox,ycen-ybox:ycen+ybox], $
  sigma=rms, mean=mean

lo = -2.0 & hi = 3.0
img = imgscl(dssimage,min=(mean+lo*rms)>min(dssimage),$
             max=(mean+hi*rms)<max(dssimage),top=topvalue)
img = bytscl(topvalue-img,min=minvalue,top=topvalue)

;-------------------
; Plot the image 

dfpsplot, psname, /square, /color
postthick = 6.0

plotimage, img, /preserve_aspect, position=pos, /normal, $
  imgxrange=minmax(xaxis), imgyrange=minmax(yaxis), charsize=1.8, $
  charthick=postthick, xthick=postthick, ythick=postthick, $
  xtitle=textoidl('\Delta\alpha [arcmin]'), $
  ytitle=textoidl('\Delta\delta [arcmin]'), title=galaxy+' - DSS'

;--------------------
; Find PA, diameters, etc, in RC3

rc3 = read_rc3()
rc3ra = 15.0D*im_hms2dec(rc3.ra)
rc3dec = im_hms2dec(rc3.dec)
searchrad = 125.0 ; search radius
ntot = im_djs_angle_match(dra,ddec,rc3ra,rc3dec,dtheta=searchrad/3600.0,$
       units='degrees',mcount=mcount,mindx=mindx,mdist=mdist)

if mindx[0] ne -1 then begin
  pa = rc3[mindx[0]].pa
  d25_min = rc3[mindx[0]].d25_min
  d25_maj = rc3[mindx[0]].d25_maj
  bmag = rc3[mindx[0]].bmag
  im_oplot_box, d25_min, d25_maj, pa, line=2, color=djs_icolor('red'), $
    thick=postthick
endif else begin
  pa = 0
  d25_min = 0
  d25_max = 0
  bmag = 0
endelse

;----------------------------
; Add legends

  label1 = [$
    '\alpha = '+ ra, $
    '\delta = '+ dec, $
    'B = '+strtrim(string(bmag,format='(F5.1)'),2),$
    'Dmaj = '+strtrim(string(d25_maj,format='(F5.1)'),2)+"'",$
    'Dmin = '+strtrim(string(d25_min,format='(F5.1)'),2)+"'",$
    '\Delta\theta = '+strtrim(string(pa,format='(I3)'),2)]

  legend, textoidl(label1), /left, /top, box=0, charsize=charsize, $
    charthick=postthick, textcolor=djs_icolor('black'), /clear

  label2 = [$ 
    'Scan Length = '+strtrim(string(scanlen,format='(I3)'),2)+'"', $
    'Slit PA = '+strtrim(string(scan_pa,format='(I3)'),2), $
    'Rotator PA = '+strtrim(string(181.2-scan_pa,format='(I3)'),2), $
    'Offset Perp = '+strtrim(string(offset_perp,format='(I5)'),2)+'"']
  for ii = 0, nslits -1 do label2 = [label2, $
    'Offset'+string(ii+1, format='(I0)')+' Along = '+$
           strtrim(string(offset_along[ii],format='(I5)'),2)+'"']

  legend, label2, /right, /top, box=0, charsize=charsize, $
     charthick=postthick, textcolor=djs_icolor('black'), /clear

  legend, driftlabel, /left, /bottom, box=0, charsize=charsize, $
     charthick=postthick, textcolor=djs_icolor('black'), /clear

  legend, string(notes), /right, /bottom, box=0, charsize=charsize, $
     charthick=postthick, textcolor=djs_icolor('black'), /clear
  
;------------------------------
; Draw drift boxes

colors = ['orange', 'red', 'yellow', 'cyan', 'green', 'blue']

for ii = 0, nslits - 1 do $
  im_oplot_box, scanlen/60.0, 3.5, scan_pa, corners = corners, $
    line = 0, color=djs_icolor(colors[ii]), thick=postthick+2, $
    xoffset=offset_perp/60.0, yoffset=offset_along[ii]/60.0 

;-------------------------------
; Add Airmass Plot

airmass_plots, date, ra, dec, object=galaxy, obsname='KPNO' 


dfpsclose

end
