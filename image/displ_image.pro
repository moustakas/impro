;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

PRO DISPL_IMAGE,image,HEADER=header,INDEX=ind,BOUNDS=limpix,LOOK_UP=lut,      $
  CUTS=lowhigh,TEXT=text,TIMEOUT=tout,LOG=log,SCALE=scale,PRT=prt

;_Purpose:      Display 2-dim images or subimages
;
;_Language:     IDL V5
;
;_Ext Calls:    EX_KEYWORD
;
;_Preparations: Input of image using READ_FITS is assumed complete
;
;_Par Input:    image       Name of image array
;
;_Keywords:     HEADER      Header for current image, no default
;               INDEX       IDL index of window to be opened. Default: 0
;               BOUNDS      Boundaries of subimage to be displayed,
;                           arranged as [xlow,ylow,xhigh,yhigh]
;                           Default: [0,NAXIS1-1,0,NAXIS2-1]
;               LOOK_UP     Colour look-up table index, from 0 to 37
;                           Default: previous setting
;               CUTS        Low and high cuts of data, blowing up the
;                           magnitude difference scale, arranged as
;                           [low_ADU,high_ADU]
;                           Default: [200,maximum of central median]
;               TEXT        Text to be written on window margin
;                           Default: read from header keyword 'OBJECT'
;                           if HEADER is specified
;               TIMEOUT     Timeout in seconds
;                           Default: infinite
;               LOG         Plot logarithmic image to increase contrast
;               SCALE       Scale image by factor (1: do not scale)
;                           Defaults to screen size
;               PRT         Writes image to PS file
;
;_History:      TG          11-Dec-91   Initial programming
;               TG          08-Apr-98   Included LOG keyword
;               TG          21-Apr-98   Included image scaling
;
;_Copyright:    (c) FOCES project, University Observatory Munich
;
;------------------------------------------------------------------------------

tmp = SIZE(image)                                        ; Get size information
nx_pix = tmp(1)                  ; ... and set subimage limits to their default
ny_pix = tmp(2)
tvpix = DATA_INPUT('TVPIX',[1280,1024])-[20,60]

b = image
IF KEYWORD_SET(scale) THEN scl = scale ELSE $               ; Scaling parameter
  scl = MIN([FLOAT(tvpix(0))/nx_pix,FLOAT(tvpix(1))/ny_pix])
IF scl NE 1. THEN $
  b = BILINEAR(image,FINDGEN(nx_pix*scl+0.5)/scl,FINDGEN(ny_pix*scl+0.5)/scl)
ix = FIX(nx_pix*scl)  &  iy = FIX(ny_pix*scl)

IF NOT KEYWORD_SET(limpix) THEN limpix = FIX([1,1,ix,iy])-1                   $
  ELSE limpix = FIX(scl*limpix)
xdif = limpix(2)-limpix(0)
ydif = limpix(3)-limpix(1)
IF KEYWORD_SET(lowhigh) THEN b = b>lowhigh(0)<lowhigh(1)   ; Set intensity cuts
IF KEYWORD_SET(log) THEN b = ALOG(b>1.)
b = b(limpix(0):limpix(2),limpix(1):limpix(3))

IF NOT KEYWORD_SET(text) THEN BEGIN
  IF NOT KEYWORD_SET(header) THEN text=' '                                    $
    ELSE text=EX_KEYWORD('OBJECT',header)
ENDIF
IF NOT KEYWORD_SET(ind) THEN ind = 0
IF KEYWORD_SET(lut) THEN LOADCT,lut                  ; Load colour lookup table
WINDOW,ind,RETAIN=2,TITLE=text,XPOS=tvpix(0)-xdif-10-20*ind>0,YPOS=5,         $
  XSIZE=xdif,YSIZE=ydif
TVSCL,b                                               ; Display scaled subimage
WSHOW,ind,1

IF (KEYWORD_SET(tout)) THEN BEGIN             ; Check timeout and delete window
  t = 1  &  tstart = SYSTIME(t)
  WHILE t-tstart LT ABS(tout) DO t = SYSTIME(t)
  IF (tout LT 0) THEN WSHOW,ind,0 ELSE WDELETE,ind
ENDIF
IF KEYWORD_SET(prt) THEN BEGIN
  PL,18.,18.,FN='c:/user/idl/image.ps',YOFF=3.,/COLOR
  TVSCL,b
  ENDPL
ENDIF

b = 0
RETURN
END