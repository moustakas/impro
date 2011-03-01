
pro plot2dh, xval, yval, xrange, yrange, figname, _extra = extra, $
    nx=nx, ny=ny, pixthresh=pixthresh, xysz = xysz, noplot = noplot, $
    maxpix = maxpix, img = img, notused = notused, over = over, weight = weight

if not keyword_set(xysz) then xysz=[7.5, 7.5]
if not keyword_set(pixthresh) then pixthresh = 5
if not keyword_set(nx) then nx = 150
if not keyword_set(ny) then ny = 150

if keyword_set(weight) then zval = weight else zval = yval * 0 + 1

;*******************************************************************************
; fix color table
;*******************************************************************************

if not keyword_set(noplot) then set_plot, 'ps'
loadct, 0, /silent
tvlct, r, g, b, /get
rr = shift(reverse(r), 16)
gg = shift(reverse(g), 16)
bb = shift(reverse(b), 16)
rr[0]=0 & gg[0]=0 & bb[0]=0 ;foreground
rr[2]=255 & gg[2]=255 & bb[2]=255 ;background
rr[16]=240 & gg[16]=240 & bb[16]=240 ;empty bin color
tvlct, rr, gg, bb
!p.background=2
!p.color = 0
!P.FONT = 0

if not keyword_set(noplot) and not keyword_set(over) then $
device, filename = figname + '.eps', /inches, /iso, /times, $
        color = 1, encapsul=1, xsize = xysz[0], ysize = xysz[1], scale = 1, $
        xoffset = 0.5, yoffset = 0.5

;*******************************************************************************
; Make points into an image
;*******************************************************************************

bin2d, nx, ny, xrange[0], xrange[1],  yrange[0], yrange[1],  xval, $
       yval, zval, map, num, pixthresh, mask=mask, tot=0, notused = notused

if keyword_set(maxpix) then num = num < maxpix
if keyword_set(weight) then cdd = alog10(num * map) else cdd = sqrt(num)
cdd[where(num eq 0)] = 0
cdd = cdd / max(cdd)

img = bytscl(cdd,top=195, max=1.2, min=0, /nan) + byte(16)
img[where(mask ne 0)] = byte(16)

;*******************************************************************************
; Plot
;*******************************************************************************

if keyword_set(noplot) then return

plotimage,img, xrange = xrange, yrange=yrange,  $
  imgxrange=xrange, imgyrange=yrange, interp = 0, _extra = extra

if n_elements(notused) ge 2 then $
oplot, xval[notused], yval[notused], psym=3, symsize = 0.1

setplotcolors

end

