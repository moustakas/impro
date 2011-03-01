; Doug Finkbeiner, April 5, 1997
; revised 9 October to include display of galactic and equatorial
; coords (see below)

; image display tool that does exactly what I want with no frills. 
; see comments below for dm


PRO dm_helpful_hints

;    dm, image, align=align, min=imin, max=imax, same=same, $
;        exit_flag=exit_flag, hemi=hemi, epoch=epoch, $
;        issa=issa, plate=plate, mark=mark, file=file, $
;        overplot=overplotin, wdelete=dm_wdelete, help=dm_help
  print, ' '
  print, '--------------------------------------------------------------------'
  print, ' Display Manager  20-Jan-1998 version                D. Finkbeiner  '
  print, '--------------------------------------------------------------------'
  print, ' MOUSE controls:  L = LEFT, M = MIDDLE, R = RIGHT                   '
  print, '   L   - hold down and move for brightness and contrast             '
  print, '   L+M - reassign color map (press L first)                         '
  print, '   M   - click for recenter, flick up zoom in,  flick down zoom out '
  print, '   R+L - exit (press R first to preserve stretch)                   '
  print, '   R+M - add image to blink ring (press R first)                    '
  print, ' (flick means press button and move mouse before releasing button)  '
  print, '                                                                    '
  print, ' KEYWORDS:                                                          '
  print, '   min, max  - minimum and maximum values for color stretch         '
  print, '   same      - same color stretch as last display                   '
  print, '   align     - align image with previous display                    '
  print, '   exit_flag - exit immediately after displaying image              '
  print, '   wdelete   - delete dm windows                                    '
  print, '   help      - this page                                            '
  print, '--------------------------------------------------------------------'
  print, ' COORDINATES                                                        '
  print, ' Lambert projections used by Schlegel, Finkbeiner & Davis are       '
  print, ' recognized.  for such projections the following keywords are used  '
  print, '   hemi      - 1 for NGP, 0 for SGP (NGP default)                   '
  print, '   epoch     - set to epoch for RA, dec display                     '
  print, '   issa      - set flag for ISSA plate number and position on plate '
  print, '   mark      - cause mouse R+M to print CR, to keep coord info      '
  print, '   file      - set to filename for list of coordinates from mark    '
  print, '   overplot  - [n, 2] array is (l, b); [n,4] is (l,b,pol,theta)     '
  print, '--------------------------------------------------------------------'
  print, ' '
  return

END


FUNCTION smartbytscl, image, min=min, max=max, top=top

sz = size(image)

IF sz(3) EQ 6 THEN BEGIN 
    print, 'Converting Complex to Real'
    imageb = bytscl(abs(image), min=min, max=max, top=top)
ENDIF ELSE BEGIN 
    imageb = bytscl(image, min=min, max=max, top=top)
ENDELSE

return, imageb
END



PRO dm_init_main_window, mx, my, mainxsize, mainysize, mainimage, image, imageb, imin, imax, ncolors, sz, colorscale, scale, overplot=overplot

  mainxsize = !d.x_size
  mainysize = !d.y_size
  
  imagesz = size(image)

  sx2 = mainxsize/(2*scale)
  sy2 = mainysize/(2*scale)
  cor = [mx-sx2, mx+sx2, my-sy2, my+sy2]
  realcorners = cor
  cor = cor > 0
  cor(0) = cor(0) < (imagesz(1)-1)
  cor(1) = cor(1) < (imagesz(1)-1)
  cor(2) = cor(2) < (imagesz(2)-1)
  cor(3) = cor(3) < (imagesz(2)-1)
  
  dsx = fix(((cor(0:1)-mx)*scale)+mainxsize/2)
  dsy = fix(((cor(2:3)-my)*scale)+mainysize/2)

; CONVENTION:   indices range from cor(0):cor(1)-1 and
; dsx(0):dsx(1)-1 etc.  This makes rebinning much easier. 

;  mainimage = image(mx:mx+mainxsize-1, my:my+mainysize-1)
;  imageb = bytscl(mainimage, min=imin, max=imax, top=ncolors)

  cor = fix(cor*scale)/scale
  imagebpiece = imageb(cor(0):cor(1)-1, cor(2):cor(3)-1)

  szb = size(imagebpiece)
  IF (scale EQ 1) THEN $
    scalimageb = imagebpiece ELSE $
    scalimageb = rebin(imagebpiece, szb(1)*scale, szb(2)*scale, sample=1)
;  mainimage = bytarr(mainxsize, mainysize)

;  print, dsx(1)-dsx(0), dsy(1)-dsy(0), mainxsize,  mainysize
  mainimage = scalimageb


  print, 'Displaying image...'
  IF (dsx(0) NE 0) OR (dsy(0) NE 0) THEN erase
  tv, mainimage, dsx(0), dsy(0)
  IF keyword_set(overplot) THEN BEGIN

; set up coordinates for overplot
      plot, [0, mainxsize-1], [0, mainysize-1], xstyle=5, ystyle=5, $
        /nodata, /noerase, xmargin=[0, 0], ymargin=[0, 0]
      
      overpltx = (overplot[*, 0]-mx)*scale + mainxsize/2
      overplty = (overplot[*, 1]-my)*scale + mainysize/2
      IF (size(overplot))[2] EQ 4 THEN BEGIN ; polarization

;          doverpltx = (overplot[*, 2]-mx)*scale + mainxsize/2
;          doverplty = (overplot[*, 3]-my)*scale + mainysize/2
          doverpltx = overplot[*, 2]
          doverplty = overplot[*, 3]

          FOR oind=0, n_elements(overpltx)-1 DO BEGIN 
              oplot, overpltx(oind)+doverpltx(oind)*[-1, 1], $
                     overplty(oind)+doverplty(oind)*[-1, 1]
          ENDFOR 
;          oplot, overpltx, overplty, ps=4
;          oplot, overpltx+doverpltx, overplty+doverplty, ps=3
      ENDIF ELSE BEGIN 
          oplot, overpltx, overplty, ps=4
      ENDELSE 
  ENDIF

  sz = size(mainimage)
  colorscale = bytscl(lindgen(mainxsize, 10) MOD mainxsize, top=ncolors)
  tv, colorscale
  
  return
END


PRO dmb, clear=clear

COMMON dmwindowinfo, mainwindex, zoomwindex, specwindex, mx, my, scale, $ 
  stretchmin, stretchmax, blinkring

wset, mainwindex
nblink = n_elements(blinkring)-1

IF keyword_set(clear) THEN BEGIN 
    print, 'Clearing blink ring              '
    FOR i=1, nblink DO $
      wdelete, blinkring(i)

    blinkring = -1
    return
ENDIF


nblink = n_elements(blinkring)-1
IF nblink LT 1 THEN BEGIN 
    print, 'NO images in blink ring'
    return
ENDIF

print, nblink, ' images in blink ring - press L+R to quit'

i = 1
mscode = 0

WHILE mscode NE 5 DO BEGIN 
    tvrdc, x, y, 1
    mscode = !err
    IF !err EQ 1 THEN i = i+1
    IF !err EQ 2 THEN i = i-1
    IF i GT nblink THEN i = 1
    IF i LT 1 THEN i = nblink
    device,copy=[0,0,!d.x_size,!d.y_size,0,0,blinkring(i)]
    REPEAT BEGIN 
        tvrdc, x, y, 0
        IF !err EQ 5 THEN return
    ENDREP UNTIL !err EQ 0
ENDWHILE

return
END



PRO overplotsetup, overplot, projsize, projedge, hemi
;+------------------------------------------------------------------------  
; NAME:
;       overplotsetup
;+------------------------------------------------------------------------  
; INPUTS:
;       overplot      - overplot coords
;+------------------------------------------------------------------------  
; PURPOSE:
; conver (l,b) to (x,y) in pixel coords on raw image
;
;+------------------------------------------------------------------------  
;-

  IF (size(overplot))[2] EQ 2 THEN BEGIN ; points
      l = overplot[*, 0]
      b = overplot[*, 1]
      IF hemi EQ 1 THEN BEGIN   ;north
          good = where(b GT 0)
      ENDIF ELSE BEGIN          ;south
          good = where(b LT -0)
      ENDELSE
      
      polar_proj_ns, l[good], b[good], nx, ny, sx, sy
      
      IF hemi EQ 1 THEN BEGIN   ;north
          overx = (nx*projsize/2)+(projsize/2+projedge)
          overy = (ny*projsize/2)+(projsize/2+projedge)
      ENDIF ELSE BEGIN          ;south
          overx = (sx*projsize/2)+(projsize/2+projedge)
          overy = (sy*projsize/2)+(projsize/2+projedge)
      ENDELSE
      overplot = [[overx], [overy]]
      return
  ENDIF 


  IF (size(overplot))[2] EQ 4 THEN BEGIN ; polarization

      b = overplot[*, 1]

      IF hemi EQ 1 THEN BEGIN   ;north
          good = where(b GT 5)
      ENDIF ELSE BEGIN          ;south
          good = where(b LT -5)
      ENDELSE

      l = overplot[good, 0]
      b = overplot[good, 1]
      p = overplot[good, 2]
      th = overplot[good, 3]

      ; this way of doing things will result in a drawn line that
          ; lines up with structure IN THIS projection. 

      dl = 0.01*sin(th/!radeg)/cos(b/!radeg)
      db = 0.01*cos(th/!radeg)
      polar_proj_ns, l, b, nx, ny, sx, sy
      polar_proj_ns, l+dl, b+db, dnx, dny, dsx, dsy
      
      IF hemi EQ 1 THEN BEGIN   ;north
          x  = (nx*projsize/2)+(projsize/2+projedge)
          y  = (ny*projsize/2)+(projsize/2+projedge)
          dx = (dnx*projsize/2)+(projsize/2+projedge)
          dy = (dny*projsize/2)+(projsize/2+projedge)
      ENDIF ELSE BEGIN          ;south
          x  = (sx*projsize/2)+(projsize/2+projedge)
          y  = (sy*projsize/2)+(projsize/2+projedge)
          dx = (dsx*projsize/2)+(projsize/2+projedge)
          dy = (dsy*projsize/2)+(projsize/2+projedge)
      ENDELSE

      dx = (dx-x)*1e2*p
      dy = (dy-y)*1e2*p

      overplot = [[x], [y], [dx], [dy]]

      return
  ENDIF 


END
    


PRO dm, image, align=align, min=imin, max=imax, same=same, $
        exit_flag=exit_flag, hemi=hemi, epoch=epoch, eclip=eclip, $
        issa=issa, plate=plate, mark=mark, file=file, $
        overplot=overplotin, wdelete=dm_wdelete, help=dm_help, $
        hispec=hispec, hirange=hirange
;+------------------------------------------------------------------------  
; NAME:
;       dm
;+------------------------------------------------------------------------  
; INPUTS:
;       image      - image to display
;+------------------------------------------------------------------------  
; KEYWORDS:
;       min,max    - min, max of stretch to display (black=min)
;       hemi       - hemi=1 for north, -1 for south (north is default)
;       epoch      - to display equatorial coords instead of (l,b),
;                      set this to the desired equinox (e.g. 1950)
;  
;+------------------------------------------------------------------------  
; FLAGS:
;       align      - align (same zoom and pan) with previous display
;       same       - display with same stretch as previous display
;       exit_flag  - exit immediately after displaying. 
;+------------------------------------------------------------------------  
; OUTPUTS:
;       just the image on the screen
;+------------------------------------------------------------------------  
; COMMENTS:
;
; This display tool is written to optimize speed with large images. 
;   Two windows are created.  "DM Main" displays the image, "DM Zoom"
;   a zoom box.  Pan and zoom are supported
;
; Square arrays at 1024, 2048, and 4096 are assumed to be Lambert
;   projections in Galactic coordinates.  Northern hemisphere is
;   assumed, unless keyword "hemi" equals -1.  Coordinates are
;   displayed as (l,b) unless the epoch keyword is used, in which 
;   case (RA,dec) for that epoch are displayed. 
;
; 20-Jan-1998 Modified to check whether windows have been deleted,
;  and if so, create windows again.  DPF
; 20-Jan-1998 Added wdelete keyword
; 20-Jan-1998 Added help keyword
;
; 29-Jan-1998 Added hispec, hirange keywords for HI spectrum display.
; 
;+------------------------------------------------------------------------  
;-
  COMMON dmwindowinfo, mainwindex, zoomwindex, specwindex, mx, my, scale, $ 
    stretchmin, stretchmax, blinkring

; first check for /wdelete keyword

  IF keyword_set(dm_wdelete) THEN BEGIN ; remove dm windows
      IF n_elements(mainwindex) EQ 1 THEN BEGIN 
          wdelete, mainwindex, zoomwindex ;delete windows
          mainwindex = -1
          zoomwindex = -1
          return
      ENDIF ELSE BEGIN 
          print, 'There were no DM windows defined yet - cannot delete'
          return
      ENDELSE
  ENDIF


; now check for /help keyword
  IF keyword_set(dm_help) THEN BEGIN 
      dm_helpful_hints
      return
  ENDIF

       
; Now we can begin the program  

  ver = 1-keyword_set(exit_flag) ; verbose keyword

  IF ver THEN print, '---------------------'
  IF ver THEN print, 'BEGIN Display Manager'
  IF ver THEN print, '---------------------'

  imagesize = size(image)
  IF imagesize(0) NE 2 THEN BEGIN
      IF ver THEN print, 'ERROR - image is not two-dimensional!'
      return
  ENDIF
  
; Added 20 Jan 1998 by DPF

  dmwinfault = 0                ; flag eq 1 if a window has been deleted
  IF n_elements(mainwindex) EQ 1 THEN BEGIN ; if mainwindex defined
    IF mainwindex EQ -1 THEN BEGIN ; if windows were previously deleted
      dmwinfault = 1
    ENDIF ELSE BEGIN 
      device, window_state=windowstate   ; get list of valid window numbers
      IF n_elements(windowstate) GT (mainwindex > zoomwindex) THEN BEGIN
          ; if the list is long enough, then check if windows are gone
          IF (windowstate[mainwindex] EQ 0) OR $
             (windowstate[zoomwindex] EQ 0) THEN BEGIN
              print, 'You deleted a window!!!'
              dmwinfault = 1
              wdelete, mainwindex, zoomwindex
          ENDIF 
      ENDIF ELSE BEGIN 
          ; if window list was not long enough, something wrong with IDL
          print, 'Bad window list - will make new windows'
          dmwinfault = 1
      ENDELSE
    ENDELSE 
  ENDIF 

  IF (n_elements(mainwindex) EQ 0) OR dmwinfault THEN BEGIN 
      IF ver THEN print, 'Creating new display windows'
      device, pseudo=8
      window, title='DM Main', /free, xsize=512, ysize=512, xpos=630, ypos=388
      mainwindex = !d.window
      window, title='DM Zoom', /free, xsize=128, ysize=128, xpos=490, $ 
        ypos=772, retain=0      ; retain=0 means no backing store

      zoomwindex = !d.window
      specwindex = -1
      blinkring = -1
      scale = 1.
      mx = imagesize(1)/2       ; coords of middle
      my = imagesize(2)/2
      
      loadct, 0
  ENDIF
  

; This information used for calculating Galactic coordinates for such
; images.  It could also be used for alignment by storing in common block.
projsize = 0
IF imagesize(1) EQ imagesize(2) THEN BEGIN ; square
    CASE imagesize(1) OF 
    1024: BEGIN projsize = 1024 & projedge=0   & ENDCASE
    2048: BEGIN projsize = 2048 & projedge=0   & ENDCASE
    4096: BEGIN projsize = 4096 & projedge=0   & ENDCASE
    1536: BEGIN projsize = 1024 & projedge=256 & ENDCASE
    2248: BEGIN projsize = 2048 & projedge=100 & ENDCASE
    4496: BEGIN projsize = 4096 & projedge=200 & ENDCASE
        ELSE: projsize = 0
    ENDCASE
    IF projsize NE 0 THEN BEGIN
        IF keyword_set(hemi) EQ 0 THEN BEGIN
            hemi = 1
            print, 'Image assumed to be polar projection - NORTH'
            print, '  (Call with hemi=-1 for south)'
        ENDIF ELSE BEGIN 
            IF hemi EQ 1  THEN $
              print, 'Image recognized as polar projection - NORTH'
            IF hemi EQ -1 THEN $
              print, 'Image recognized as polar projection - SOUTH'
        ENDELSE
        print, 'SIZE = ', projsize, '  EDGE = ', projedge
        csys = 1
        IF keyword_set(epoch) THEN BEGIN
            print, 'Coordinates:  RA,dec (', string(epoch,form='(I4)'), ')' 
            csys = 2
        ENDIF
        IF keyword_set(issa) THEN BEGIN 
            print, 'Coordinates:  ISSA plates'
            csys = 3
        ENDIF
        IF keyword_set(plate) THEN BEGIN
            print, 'Coordinates:  l,b ON ISSA plate number ', plate
;            csys = 4
            ; not yet implemented
        ENDIF 
        IF csys EQ 1 THEN BEGIN
            print, 'Coordinates:  Gal l,b'            
        ENDIF

    ENDIF
ENDIF


IF keyword_set(hispec) THEN BEGIN 

    IF specwindex EQ -1 THEN BEGIN 
        window, title='DM Spectrum', /free, xsize=512, ysize=512, $
          xpos=0, ypos=388
        specwindex = !d.window
    ENDIF      

ENDIF


IF NOT keyword_set(align) THEN BEGIN 
    scale = 1
    mx = imagesize(1)/2         ; coords of middle
    my = imagesize(2)/2
ENDIF
    

wshow, zoomwindex
wshow, mainwindex
wset, mainwindex

IF keyword_set(same) EQ 0 THEN BEGIN 
    erase
    IF ver THEN print, 'window cleared'
ENDIF

;imin = min(image)
;imax = max(image)
;print, 'imin,'


IF (keyword_set(same) EQ 1) AND (keyword_set(stretchmin)) THEN BEGIN 
    imin = float(stretchmin)
    imax = float(stretchmax)
ENDIF ELSE BEGIN 

    IF keyword_set(imin) EQ 0 THEN imin = 0. ; these must be real!
    IF keyword_set(imax) EQ 0 THEN imax = 80.
    imin = float(imin)
    imax = float(imax)
    IF ver THEN print, imin, imax
; use min and max for type byte
    IF imagesize(3) EQ 1 THEN imin = min(image, max=imax) 
    stretchmin = imin
    stretchmax = imax
ENDELSE
  
IF keyword_set(overplotin) THEN BEGIN 
    overplot = overplotin
    overplotsetup, overplot, projsize, projedge, hemi
ENDIF ELSE BEGIN 
    overplot = 0
ENDELSE

  
wset, mainwindex

;mx = (mx > 0) < (imagesize(1)-mainxsize)
;my = (my > 0) < (imagesize(2)-mainysize)


; !d.table_size is the number of available colors
ncolors = !d.table_size

imageb = smartbytscl(image, min=imin, max=imax, top=ncolors)

dm_init_main_window, mx, my, mainxsize, mainysize, mainimage, image, imageb, imin, imax, ncolors, sz, colorscale, scale, overplot=overplot

outstr = ' '
IF keyword_set(file) THEN BEGIN
    openw, ifile, file, /get_lun
    print, 'Writing to file ', file
    mark = 1
ENDIF

FOR i=1, 25000 DO BEGIN

    IF keyword_set(exit_flag) EQ 0 THEN $
      tvrdc, x, y, 2, /device
    mscode = !err

    IF (mscode EQ 5) OR keyword_set(exit_flag) THEN BEGIN 
                                ; EXIT  (rt + left buttons)
        IF ver THEN print
        IF ver THEN print, 'EXIT     '
        IF keyword_set(file) THEN BEGIN
            close, ifile
            print, 'Closing output file ', file
        ENDIF
        return
    ENDIF
    
    ; store in blink ring
    IF (mscode EQ 6) THEN BEGIN 
        IF keyword_set(mark) THEN BEGIN
            print
            IF keyword_set(file) THEN printf, ifile, outstr
        ENDIF ELSE BEGIN 
            print, 'storing in blink ring                                            '
            window, /free, /pixmap, xsize=mainxsize, ysize=mainysize
            tv, mainimage
            blinkring = [blinkring, !d.window]
            print, blinkring
            wset, mainwindex
        ENDELSE 

        ; clear mouse buffer
        REPEAT BEGIN
            tvrdc, x, y, 0, /device
            mscode = !err
        ENDREP UNTIL (mscode EQ 0)
    ENDIF
    
    
    
    IF (mscode EQ 2) THEN BEGIN ; PAN / ZOOM
                                ; the idea is to catch user in this
                                ; loop until right mouse is released
        yold = y
        
        REPEAT BEGIN
            
            tvrdc, x, y, 0, /device
            mscode = !err
            
        ENDREP UNTIL (mscode EQ 0)
        newscale = scale
        IF y GT yold+5 THEN BEGIN 
            newscale = scale*2
            IF ver THEN print, 'ZOOM in to scale', newscale, format='(A,F14.10)'
        ENDIF
        
        IF y LT yold-5 THEN BEGIN
            newscale = scale/2.
            IF ver THEN print, 'ZOOM out to scale', newscale, format='(A,F14.10)'
        ENDIF
        
        mx = ((mx+(magx-mainxsize/2)/scale) > 0) < (imagesize(1))
        my = ((my+(magy-mainysize/2)/scale) > 0) < (imagesize(2))
        IF ver THEN print, 'panning to ', mx, my
        scale = newscale
        dm_init_main_window, mx, my, mainxsize, mainysize, mainimage, image, imageb, imin, imax, ncolors, sz, colorscale, scale, overplot=overplot
        
        
    ENDIF
    
    magx = x
    magy = y
    x = x*ncolors/mainxsize &  y=y*ncolors/mainxsize - 16
    
    
    IF ((mainxsize NE !d.x_size) OR (mainysize NE !d.y_size)) THEN BEGIN
                                ; change of main window size
        erase
        dm_init_main_window, mx, my, mainxsize, mainysize, mainimage, image, imageb, imin, imax, ncolors, sz, colorscale, scale, overplot=overplot
        blinkring = -1   
    ENDIF
    
    IF mscode AND 1 THEN stretch, x-y, x+y
    
    IF (mscode AND 3) EQ 3 THEN BEGIN
        
                                ;re-display image with different scale
        IF ver THEN print, 'Reassign color map to levels...'
        nextmin = imin+(imax-imin)/ncolors *(x-y)
        nextmax = imin+(imax-imin)/ncolors *(x+y)
        
        IF ver THEN print, imin, imax, nextmin, nextmax
        
        imin = nextmin
        imax = nextmax
        stretchmin = imin
        stretchmax = imax
        
        imageb = smartbytscl(image, min=imin, max=imax, top=ncolors)
        
        dm_init_main_window, mx, my, mainxsize, mainysize, mainimage, image, imageb, imin, imax, ncolors, sz, colorscale, scale, overplot=overplot
        
        stretch, 0, ncolors
        
    ENDIF
    
    IF (mscode EQ 0) THEN BEGIN 
        getcr, cr
;        cr = '                '+cr
        num_decimal_places = (5-round(alog10(abs(imax-imin)))) < 8
        strdec = string(num_decimal_places, format='(I1)')
        IF imagesize(3) GE 4 THEN valform = 'F10.'+strdec ELSE valform = 'I'
        IF imagesize(3) EQ 6 THEN valform = '2G12.7'
        form="($,'x=',i5,', y=',i5,', value=',"+valform+",a)"

        ptx = mx+(magx-mainxsize/2)/scale ;??? why is this -1 here?
        pty = my+(magy-mainysize/2)/scale ; ok, so it's not now...
        IF (ptx GE 0) AND (pty GE 0) AND $
          (ptx LE (imagesize(1)-1)) AND (pty LE (imagesize(2)-1)) THEN BEGIN
            
            ptval = image(ptx, pty)
                                ; if map is a recognized projection
            IF projsize NE 0 THEN BEGIN
                polar_inv_proj,ptx,pty,projsize,projedge,hemi,gall,galb,good
                IF keyword_set(hispec) THEN BEGIN 
                    hispeclb = hispec(*, (round(gall*2) MOD 720), $
                                      round(galb*2+180))
                    wset, specwindex
;                    yrangelow = 0.01
                    yrangelow = -1
                    hispeclb = (hispeclb/100.) > yrangelow
;                    plot, hirange, hispeclb, /ylog, yrange=[yrangelow, 200], $
;                      xtitle='V_LSR (km/s)', ytitle='K', ystyle=1, xstyle=1
                    plot, hirange, hispeclb, yrange=[-1, 1], $
                      xtitle='V_LSR (km/s)', ytitle='K', ystyle=1, xstyle=1
                    wset, mainwindex
                ENDIF

                CASE csys OF 
                    1: BEGIN    ; l,b
                        form="($,'x=',i5,', y=',i5,', l=',f7.3,' " + $
                          "b=',f7.3,', value=',"+valform+",a)"
                        outstr = string(form=form, ptx, pty, $
                                        gall, galb, ptval, cr)
                    ENDCASE
                    2: BEGIN    ; ra,dec or eclip
                        form="($,'x=',i5,', y=',i5,', RA=',a,' " + $
                          "dec=',a,', value=',"+valform+",a)"
                        IF keyword_set(eclip) THEN $
                          euler_2000, gall, galb, ra, dec, 6 $
                        ELSE $
                          euler_2000, gall, galb, ra, dec, 2
                        IF (epoch NE 2000) THEN $
                          precess, ra, dec, 2000., epoch
                        rastr = dec2hms(ra/15.)
                        decstr = strmid(dec2hms(dec),0,9)         
                        outstr = string(form=form, ptx, pty, $
                                        rastr, decstr, ptval, cr)
                    ENDCASE
                    3: BEGIN    ; ISSA
                        form="($,'l=',f7.3,' b=',f7.3,', plate=',i3,', "+ $
                          "x=',i4,' y=',i4,' value=',"+valform+",a)"
                        issa_coord2best,gall,galb,iplate,xpix,ypix
                        outstr = string(form=form, gall, galb, iplate, $
                                        xpix, ypix, ptval, cr)
                    ENDCASE
                    ELSE: BEGIN 
                        outstr = string(form=form, ptx, pty, ptval, cr)
                        
                    ENDCASE 
                ENDCASE
                
            ENDIF ELSE BEGIN 
;                form="($,'x=',i5,', y=',i5,', value=     <edge>               ',a)"
                outstr = string(form=form, ptx, pty, ptval, cr)
            ENDELSE
            
            print, outstr, form='($,a)'

            tmag = 2
            wset, zoomwindex
            tsize = (!d.x_size > !d.y_size)
            
            tsize = (tsize/tmag/2 + 1)*(tmag*2)
            
            toff = tsize/(tmag*2)
            x1 = ((magx-toff) > 0) < (sz(1)-toff*2)
            x2 = ((magx+toff-1) > (toff*2-1)) < (sz(1)-1)
            y1 = ((magy-toff) > 0) < (sz(2)-toff*2)
            y2 = ((magy+toff-1) > (toff*2-1)) < (sz(2)-1)
            
            
            IF ((x1 GT 0) AND (y1 GT 0)) THEN BEGIN 
                thumbnail = mainimage(x1:x2, y1:y2)
                tv, rebin(thumbnail, tsize, tsize, sample=1)
            ENDIF
            
            wset, mainwindex
        ENDIF
        
    ENDIF
    
    
ENDFOR
    
return

END

PRO dm_help_window

  charsizeold = !p.charsize
  !p.charsize = 1.2
  window, /free, xsize=460, ysize=400, title='DM Help Window'
  plot,[0,100],[1,30],/nodata, xmargin=[0, 0], ymargin=[0, 0], xstyle=5, ystyle=5
  
  xyouts,0,29, '--------------------------------------------------------------------'
  xyouts,0,28, ' Display Manager  20-Jan-1998 version                D. Finkbeiner  '
  xyouts,0,27, '--------------------------------------------------------------------'
  xyouts,0,26, ' MOUSE controls:  L = LEFT, M = MIDDLE, R = RIGHT                   '
  xyouts,0,25, '   L   - hold down and move for brightness and contrast             '
  xyouts,0,24, '   L+M - reassign color map (press L first)                         '
  xyouts,0,23, '   M   - click for recenter, flick up zoom in,  flick down zoom out '
  xyouts,0,22, '   R+L - exit (press R first to preserve stretch)                   '
  xyouts,0,21, '   R+M - add image to blink ring (press R first)                    '
  xyouts,0,20, ' (flick means press button and move mouse before releasing button)  '
  xyouts,0,19, '                                                                    '
  xyouts,0,18, ' KEYWORDS:                                                          '
  xyouts,0,17, '   min, max  - minimum and maximum values for color stretch         '
  xyouts,0,16, '   same      - same color stretch as last display                   '
  xyouts,0,15, '   align     - align image with previous display                    '
  xyouts,0,14, '   exit_flag - exit immediately after displaying image              '
  xyouts,0,13, '   wdelete   - delete dm windows                                    '
  xyouts,0,12, '   help      - this page                                            '
  xyouts,0,11, '--------------------------------------------------------------------'
  xyouts,0,10, ' COORDINATES                                                        '
  xyouts,0,9, ' Lambert projections used by Schlegel, Finkbeiner & Davis are       '
  xyouts,0,8, ' recognized.  for such projections the following keywords are used  '
  xyouts,0,7, '   hemi      - 1 for NGP, 0 for SGP (NGP default)                   '
  xyouts,0,6, '   epoch     - set to epoch for RA, dec display                     '
  xyouts,0,5, '   issa      - set flag for ISSA plate number and position on plate '
  xyouts,0,4, '   mark      - cause mouse R+M to xyouts,0,10 CR, to keep coord info'
  xyouts,0,3, '   file      - set to filename for list of coordinates from mark    '
  xyouts,0,2, '   overplot  - [n, 2] array is (l, b); [n,4] is (l,b,pol,theta)     '
  xyouts,0,1, '--------------------------------------------------------------------'

  !p.charsize = charsizeold

return
END
