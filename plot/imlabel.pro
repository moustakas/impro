pro imlabel,h,pflag,key=key,xdis=xdis,ydis=ydis,chsize=chsize
;+
; NAME:
;       IMLABEL
; PURPOSE:
;       To label an image which has been displayed on the workstation
;       with pixel values, or right ascension and declination if
;       astrometry info is present in the header. Will also write to a 
;       postscript file.
;
; CALLING SEQUENCE:
;       IMLABEL, h, pflag, [ KEY , XDIS =, YDIS =, CHSIZE = ]
;
; INPUTS:
;       h       - FITS image header array
;
; OPTIONAL INPUTS:
;       pflag   - if present and nonzero, display pixels only
;
; OPTIONAL KEYWORD INPUTS:
;       key   - optional FITS keyword to be extracted for title, scalar string
;       xdis  - approx. dist. in pixels for x tics (def=20)
;       ydis  - approx. dist. in pixels for y tics (def=20)
;       chsize - scale factor for size of character labels (def = 1.0)
;
; OUTPUTS:
;       None
; SYSTEM VARIABLES: 
;       The system variable !GRID should be set to 1 to display 
;       full grid lines, or set to 0 to see just tic marks.
;
; RESTRICTIONS:
;       The following image present difficulties for IMLABEL
;               (1) Declinations greater than 70 or less than -70 degrees.
;               (2) Rotation angles near 90 or 270 degrees
;       The RA and Dec labels will only be approximate for large
;       declinations.
;
; SIDE EFFECTS:
;       IMLABEL will mask the outside 50 pixels of the screen.  If you need
;       To see part of an image hidden under the mask, suppress the graphics
;       overlay (graphon).  The color and intensity of the coordinate overlay
;       may be modified with the procedures GRCOL and GRINT, respectively.
;       NOTE: IF writing to a postscript file, the border will NOT be masked.
;       Use procedure Border to create the border first.
;
; REVISON HISTORY:
;       written by B. Pfarr, STX, 4/14/87
;       REVISED G. HENNESSY FOR IVAS SEP 1988
;       Modified for use with a workstation.  M. Greason, STX, May 1990
;       Modified to use normal coordinates.  Will work with postscript. 
;                B. PFarr, STX, 1/91
;       Modified to also work on the IVAS and DEANZA, BP, 2/91
;       converted key,xdis,ydis to keywords BP, 4/91
;       added chsize keyword, BP, 4/91
;       Updated to use ASTROMETRY structure.  J.Offenberg, Hughes STX, Jan 1993
;-
 On_error,2

 common tv,chan,zoom,xroam,yroam
 common images,x00,y00,xsize,ysize

; Check for required param

 npar = N_params()
 if npar eq 0 then begin
   print,'Syntax - IMLABEL, h, pflag, [ KEY , XDIS =, YDIS =, CHSIZE = ]'
   print,'   Keyword variables - key,xdis,ydis'
   return
 endif

 if !D.WINDOW EQ -1 then begin
   message,'ERROR - No image windows active'
   return
 endif

 zparcheck,'IMLABEL',h,1,7,1,'FITS header array'

;                       Initialize.

 sv_fancy=!FANCY                        ;save system variable fancy
 !FANCY = 1
 NPIX=50
 NPIX1=NPIX+1
 color = 255
 if (!d.name eq 'PS') then color = 0

 grchan = chan
 if (!d.name eq 'DEANZA') then GRCHAN=4
 if (!d.name eq 'IVAS')   then GRCHAN=3
;
;                                    Make sure IMAGES common block exists
 if N_elements(x00) eq 0 then begin
   print,string(7B),'IMLABEL: WARNING - No previous image loaded with CTV'
   x00 = intarr(11)  &  y00 = x00
   xsize = x00 + !D.X_VSIZE   &  ysize = y00 + !D.Y_VSIZE
 endif

;                                       check for pixels flag
 if npar lt 2 then pflag = 0       

;     get keyword for optional label
 if keyword_set(key) then keyword = key else keyword = '' 

;                                       load default tic spacings
 if keyword_set( XDIS ) then pixx = xdis else pixx = 20
 if keyword_set( YDIS ) then pixy = ydis else pixy = 20
 if not keyword_set( CHSIZE ) then chsize = 1.0

;                       Compute corners and limits

   if (!d.name eq 'PS') then begin
      xdis = xsize(chan)
      ydis = ysize(chan)
      x0 =   NPIX
      y0 =   NPIX
      x1 = xsize(chan)-npix
      y1 = ysize(chan)-npix
      xtv = 0 
      ytv = 0
      x0_im = x0
      y0_im = y0
      x1_im = x1
      y1_im = y1
      xmid=xsize(chan)/2
      ymid=ysize(chan)/2
   endif else begin
      xdis= !d.x_vsize
      ydis = !d.y_vsize
      x0 = x00(chan) > NPIX
      y0 = y00(chan) > NPIX
      x1 = x0+xsize(chan) < (!D.X_VSIZE - NPIX1)
      y1 = y0+ysize(chan) < (!D.Y_VSIZE - NPIX1)
      x0_im = (x0 - x00(chan)) > 0
      y0_im = (y0 - y00(chan)) > 0
      xtv = x0 - x0_im 
      ytv = y0 - y0_im  ;Offset between image & TV coords
      xsiz =  x1-x0     &    ysiz = y1 - y0
      x1_im = x0_im + xsiz & y1_im = y0_im + ysiz
      xmid = x0_im + xsiz/2.  &   ymid = y0_im + ysiz/2.
 endelse

;blank out border if not postscript
 if (!d.name ne 'PS') then begin
   if (x0 le npix) or (y0 le npix) or (x1 ge X1+1) or (y1 ge Y1+1) then begin
     print,'erasing image edges'
     tv,bytarr((!D.X_VSIZE/2)*2,((!D.Y_VSIZE-y1)/2)*2)+1,0,((y1+1)/2)*2,CHAN
     tv,bytarr(!D.X_VSIZE,(y0/2)*2)+1,0,0,CHAN
     tv,bytarr(((!D.X_VSIZE-x1)/2)*2+3,!D.Y_VSIZE)+1,(x1/2)*2,0,CHAN
     tv,bytarr((x0/2)*2,!D.Y_VSIZE) + 1,0,0,CHAN
   endif
 endif

;                       Draw box around image

 xp = [x0-1,x0-1, x1+1, x1+1, x0-1]
 yp = [y0-1, y1+1,y1+1, y0-1, y0-1]
 xp = float(xp)/xdis
 yp = float(yp)/ydis
 plots, xp, yp, /normal, COLOR = color, CHAN = grchan

;                       Try to get astrometry from header

 extast, h, astr, noparams
 proj = strmid( astr.ctype(0),5,3)
 if noparams ge 0 then begin     ; Determine RA and Dec of center
   if proj EQ 'UIT' then $
           uit_xy2ad, xmid, ymid, astr, a, d else $
           xy2ad, xmid, ymid, astr, a, d
   radec, a, d, i1, i2, i3, i4, i5, i6
   acen = string(i1,format='(I2)') + ':' + string(i2,format='(I2)') + ':' + $
          string(i3,format='(F5.1)')
   dcen = string(i4,format='(I3)') + ':' + string(i5,format='(I2)') + ':' + $
          string(fix(i6),format='(I2)')
   cen = 'CENTER PIXEL:  R A: ' + acen + '  DEC: ' + dcen
   xlen = strlen(cen) * 12   
   yc = (y0 - 70) > 3
   xyouts, float(xdis/2)/xdis, float(yc)/ydis, cen, $
    alignment=0.5,size=chsize*1.00,color=color,/normal,chan=grchan
 endif 

 if (noparams lt 0) or (pflag ge 1) then begin

; Make labels for tic marks (every other tic)

   numtica = (x1 - x0) / pixx + 1
   numticd = (y1 - y0) / pixy + 1
   ticlabx = strarr(numtica)
   ticlaby = strarr(numticd)
;
   for i=0,numtica/2 do ticlabx(i)=string((i*2*pixx)+x0,'(i3)')
   for i=0,numticd/2 do ticlaby(i)=string((i*2*pixy)+y0,'(i3)')

;                       First tics at first pixels

   xtic1 = 0
   ytic1 = 0

;                       Set units = pixels

   xunits = 'PIXELS'
   yunits = 'PIXELS'

  rot = 0.0                     ;No rotation angle

 endif else begin                 

;                       Determine rotation angle

    getrot,h,rot
    while rot gt 90 do rot = rot - 180
    while rot lt -90 do rot = rot + 180

; Determine ra and dec max and mins of bottom and left axes

   if proj EQ 'UIT' then begin
        UIT_xy2ad, x0_im, y0_im, astr, botmin, leftmin
        UIT_xy2ad, x1_im, y0_im, astr, botmax, d 
        UIT_xy2ad, x0_im, y1_im, astr, a, leftmax 
   endif else begin 
        xy2ad, x0_im, y0_im, astr, botmin, leftmin
        xy2ad, x1_im, y0_im, astr, botmax, d
        xy2ad, x0_im, y1_im, astr, a, leftmax
   endelse

   if (!debug gt 0) then print,'IMLABEL: leftmin,leftmax=',leftmin,leftmax
   if (!debug gt 0) then print,'IMLABEL: botmin,botmax=',botmin,botmax

; Determine tic size and label units

   tics, botmin, botmax, x1-x0, pixx, raincr,/RA
   tics, leftmin, leftmax, y1-y0, pixy, decincr
   numtica = (x1 - x0) / pixx + 1
   numticd = (y1 - y0) / pixy + 1

; Determine pos and value at 1st tic

   tic_one, botmin, pixx, raincr, botmin2, xtic1, /RA
   tic_one, leftmin, pixy, decincr, leftmin2, ytic1

; Create tic labels

   ticlabels, botmin2,  numtica, raincr, ticlabx, /RA
   ticlabels, leftmin2,  numticd, decincr, ticlaby

; Label tic marks in ra and dec

   xunits = 'RIGHT ASCENSION'
   yunits = 'DECLINATION'

 endelse

; Set up x, y pixel grids

 rot_rad = rot / !radeg
 gridx = [indgen(numtica)*pixx+x0, x1]
 gridy = [indgen(numticd)*pixy+y0, y1]

; If !grid = 1 draw grid lines

 if (!grid eq 1) then begin

; Draw x grid lines for pixels

 if (noparams lt 0) or (pflag ge 1) then begin
        for i = x0+xtic1, x1, pixx do begin
           xest = replicate(i,numticd+1)
           xx= [i,xest]
           yy= [y0,gridy]
       plots, float(xx)/xdis,float(yy)/ydis, /normal, color=color,chan=grchan
        end

; Draw y grid lines for pixels

        for i = y0+ytic1, y1, pixy do begin
           yest = replicate(i,numtica+1)
           xx=[x0,gridx]
           yy=[i,yest]
           plots, float(xx)/xdis, float(yy)/ydis, /normal,  $
             color=color, chan=grchan
         end
    endif else begin

        gridx_im = gridx - xtv
        gridy_im = gridy - ytv

;                       Draw x grid lines for ra display

        if rot gt 0 then xedge = x0 else xedge = x1
        for i = x0+xtic1, x1, pixx do begin
           ii = i - xtv
           if proj EQ 'UIT' then uit_xy2ad, ii, y0, astr, a,d else $
           xy2ad, ii, y0, astr, a, d
           newx = cons_ra(a, gridy_im, astr) + xtv
           if (newx(numticd) gt x1) or (newx(numticd) lt x0) then begin 
                  yeff = interpol(gridy, newx, fltarr(1)+xedge) 
                  good = where ( (newx ge x0) and (newx le x1))
                  newx = [newx(good),xedge] 
                  newy = [gridy(good),yeff]
           endif else newy = gridy
           xx=[i,newx]
           yy=[y0,newy]
           plots, float(xx)/xdis, float(yy)/ydis, /normal, $
             color=color,chan=grchan
        end
        if rot gt 0 then begin
                x_1 = i & x_2 = 2000 & delta = pixx & xedge = x1
        endif else begin
                x_1 = x0+xtic1-pixx & x_2 = -2000 & delta = -pixx & xedge = x0
        endelse
        for i = x_1, x_2, delta do begin
           ii = i - xtv
           if proj EQ 'UIT' then uit_xy2ad, ii, y0_im, astr, a, d else $
           xy2ad, ii, y0_im, astr,a , d
           newx = cons_ra(a,gridy_im,astr) + xtv
           yeff1 = interpol(gridy,newx,fltarr(1) + xedge)
           yeff2 = y1
           if  ( (rot le 0) and (newx(numticd) lt x0) )  or $
               ( (rot ge 0) and (newx(numticd) gt x1) ) then goto,donex 
           if (newx(numticd) gt x1) then begin
               yeff2 = interpol(gridy,newx,fltarr(1) + x1)
               xeff2 = x1           
           endif else if (newx(numticd) lt x0) then begin
               yeff2 = interpol(gridy,newx,fltarr(1) + x0)
               xeff2 = x0
           endif else xeff2 = newx(numticd)
           good = where ( (newx ge x0) and (newx le x1))
           newx = [xedge,newx(good),xeff2] 
           newy = [yeff1,gridy(good),yeff2]
           
           plots, float(newx)/xdis, float(newy)/ydis,$
                /normal, color=color, chan=grchan
        end

; Draw y grid lines 

donex:  if rot gt 0 then yedge = y1 else yedge =y0
        for i = y0+ytic1, y1, pixy do begin
            ii = i - ytv
            if proj EQ 'UIT' then uit_xy2ad, x0_im, ii, astr, a, d else $
            xy2ad, x0_im, ii, astr, a, d
            newy = cons_dec(d,gridx_im,astr) + ytv
            if (newy(numtica) lt y0) or (newy(numtica) gt y1) then begin
                    xeff = interpol(gridx,newy,fltarr(1)+yedge)
                    good = where( (newy ge y0) and (newy le y1) )
                    newx = [gridx(good),xeff]
                    newy = [newy(good),yedge]
            endif else newx = gridx
            xx=[x0,newx]
            yy=[i,newy]
            plots, float(xx)/xdis, float(yy)/ydis, $
                    /normal, color=color,chan=grchan
        end
        if rot lt 0 then begin
            y_1 = i & y_2 = 2000 & yedge = y1 & delta = pixy
        endif else begin
            y_1 = y0 + ytic1 - pixy & y_2 = -2000 & yedge = y0 & delta = -pixy
        endelse
        for i = y_1,y_2,delta do begin
            ii = i - ytv
            if proj EQ 'UIT' then uit_xy2ad, x0_im, ii, astr, a, d else $
            xy2ad, x0_im, ii, astr, a, d
            newy = cons_dec(d,gridx_im,astr) + ytv
            xeff1 = interpol(gridx,newy,fltarr(1)+yedge)
            xeff2 = x1
            if ( (rot le 0) and (newy(numtica) gt y1) ) or $
               ( (rot ge 0) and (newy(numtica) lt y0) ) then goto,doney
            if (newy(numtica) gt y1) then begin
                         yeff2 = y1
                         xeff2 = interpol(gridx,newy,fltarr(1)+y1)
            endif else if (newy(numtica) lt y0) then begin
                         yeff2 = y0
                         xeff2 = interpol(gridx,newy,fltarr(1)+y0)
            endif else yeff2 = newy(numtica)
            good = where( (newy ge y0) and (newy le y1) )
            newx = [xeff1,gridx(good),xeff2] 
            newy = [yedge,newy(good),yeff2] 
            plots, float(newx)/xdis,float(newy)/ydis,$
                        /normal, color=color, chan=grchan
        end
   endelse
doney:
endif else begin

; If !grid = 0 compute tics using rot angle

   xtic=cos(rot_rad)*10.0
   ytic=sin(rot_rad)*10.0

; Draw x tic marks

   for i = x0+xtic1, x1, pixx do begin
       plots, float([i, i-ytic])/xdis, float([y0, y0+xtic])/ydis,$
            /normal, color=color, chan=grchan
       plots, float([i, i+ytic])/xdis, float([y1, y1-xtic])/ydis,$
            /normal, color=color, chan=grchan
   end

; Draw y tic marks

   for i = y0+ytic1, y1, pixy do begin
       plots, float([x0, x0+xtic])/xdis, $
          float( [i, i+ytic])/ydis, /normal, color=color, chan=grchan
       plots, float([x1, x1-xtic])/xdis, $
         float([i, i-ytic])/ydis,/normal, color=color, chan=grchan
   end
endelse

; Extract name and label and write to screen

 name = sxpar(h,'OBJECT')
 if (!err ne -1) then xyouts, float(x0)/xdis,  $
    float(y1+12)/ydis, name,size=chsize*1.6,color=color,/normal,chan=grchan
 if (keyword ne '') then begin
   remchar,keyword,' '
   date = sxpar(h,keyword)
   if (!err ne -1) then xyouts, float(x0+100)/xdis, $ 
                           float( y1+12)/ydis, date, $
                           size=chsize*1.5,color=color,/normal,chan=grchan
 endif

; Label tic marks (every other tic)

 k = 0
 for i = xtic1, x1-x0, 2*pixx do begin
     label = strtrim(ticlabx(k),2)
     xyouts, float(x0+i)/xdis,float( y0-20)/ydis,$
       label, size=chsize*0.9,color=color,alignment=0.5,/normal,chan=grchan
     k = k + 1
 endfor

 k = 0
 for i = ytic1, y1-y0, 2*pixy do begin
    label = strtrim(ticlaby(k),2)
    xyouts, float(x0-5)/xdis, float(y0+i-2)/ydis, label, $
         size=chsize*0.9,color=color,alignment=1.0,/normal,chan=grchan
    k = k + 1
 end

; Label x and y units 

 xlen = strlen(xunits)*12
 ylen = strlen(yunits)*22
 xyouts, float(xdis/2)/xdis,float( y0-32)/ydis, $
    xunits,size=chsize*1.00,color=color,alignment=0.5,/normal,chan=grchan
 xyouts, float((x0-30)>12)/xdis,float( ydis/2)/ydis, $
    yunits, size=chsize*1.00,orientation=90.0,color=color, $
    alignment=0.5,/normal,chan=grchan

 !fancy = sv_fancy

 return
 end
