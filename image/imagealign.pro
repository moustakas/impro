pro imagealign, im1, im2, xr

;+
;NAME:
;	IMAGEALIGN

;PURPOSE:
;	This procedure take two images of the same field and aligns them
;exactly. 

;CALLING SEQUENCE:
;	IMAGEALIGN, Im1, Im2, Boxsize

;INPUTS:
;	Im1: The image used as the 'standard' to which Im2 will be
;aligned. 

;	Im2: The image that should be moved to be aligned with Im1. 

;OPTIONAL INPUTS:

;	Boxsize: This is the size of the area used in the cross
;correlation calculation.  A larger Boxsize means a larger area is
;included in the calculation, which includes more stars in the cross
;correlation, which usually leads to more accurate alignment. 
;
;	However, the computation time goes as Boxsize SQUARED.  The
;default boxsize is 128. 
;
;OUTPUTS:
;	There are no outputs.  The amount by which Im2 needs to be moved
;is printed on the screen and you are given the choice as to whether to
;move Im2 by this amount or not. 

;RESTRICTIONS:
;	The input images must be FLOATING POINT (i.e., REAL) numbers. 

;	You are asked to move the cursor onto a star and click; this
;defines the center of the box used in the cross correlation.  This
;center doesn't HAVE to be a star; it can be any point.  But this point
;must be at least as far from the edge of the image as Boxsize. 

;PROCEDURE:
;	This routine first removes the sky background by convolving each
;image with a small second-derivative type of kernal (6 X 6 pixels). 
;Then it cross-correlates the two images and finds the position where the
;cross-correlation is maximum.  Then it prints out the amount by which
;the Im2 needs to be shifted and asks you if you want to actually move
;Im2 by this amount.  When doing the moving, we use the IDL SHIFT
;procedure, in which the edges of Im2 are 'wrapped around' so that no
;information is lost; you can move it back to the original position with
;SHIFT and no information is lost. 

;EXAMPLE:
;	Get Im1 into a window.  Then type IMAGEALIGN, Im1, Im2, Boxsize. 

;HISTORY:
;	Written by Carl Heiles.  Documented 13 Dec 1997.  Modified 13
;Dec 1997 to check if images are real numbers. 
;Modified 9 Nov 1998 to do images of arbitrary size (it was 512)
;Updated documentation 19 Oct 99.
;-

; begin with two 512 by 512 images called im1, im2;
; end with im2 being shifted so that it is aligned with im1

if (n_params() eq 0) then begin
print, 'BEGIN with two images of identical size, say called im1, im2'
print, 'end with im2 being shifted so that it is aligned with im1'
print, 'you should have the image im1 displayed before calling this'
print, ' '
print, 'call this by: imagealign, im1, im2'
print, 'there is an optional third parameter, the box size: default is 128'
print, 'you will be prompted for the window nr of im1, and'
print, 'you will be asked to click on a star to define the shifting'
print, ' '
print, '***** IMPORTANT ***** the input images MUST be FLOAT!!!!!'
return
endif

testsize1 = size( im1)
testsize2 = size(testsize1)
testsize1im1 = testsize1
tsts = testsize1( testsize2(3)-2)
if ( (tsts ne 4) and (tsts ne 5) ) then begin
	print, 'Im1 image is not real numbers! Returning.'
	print, 'REDEFINE IT AS REAL NUMBERS AND TRY AGAIN.'
;	stop
	return
endif

testsize1 = size( im2)
testsize2 = size(testsize1)
testsize1im2 = testsize1
tsts = testsize1( testsize2(3)-2)
if ( (tsts ne 4) and (tsts ne 5) ) then begin
	print, 'Im2 image is not real numbers! Returning.'
	print, 'REDEFINE IT AS REAL NUMBERS AND TRY AGAIN.'
;	stop
	return
endif

if ( (testsize1im2[1] ne testsize1im1[1]) or $
	(testsize1im2[2] ne testsize1im1[2]) ) then begin
	print, 'The images are of DIFFERENT SIZE! Returning.'
;	stop
	return
endif

x512 = testsize1im1[1]
y512 = testsize1im1[2]

if (n_params() ne 3) then begin
 xr = 128l
print, 'box size is the default, namely 128.'
print, 'if you want a different box size, enter it as the third parameter.'
print, 'for example, imagealign, im1, im2, 64 will run much faster...'
print, 'but give possibly unacceptable results.'
endif
;xr = 256l
;xr=64l

print, ' '
print, '***** IMPORTANT ***** the input images MUST be FLOAT!!!!!'
print, ' '

read, wndw, prompt = 'type the window nr of the first image: ' 
wset, wndw

print, 'box size is ', xr
print, 'click with the mouse on a nice bright star, not too near the edge'
print, 'THAT IS: FURTHER FROM THE EDGE THAN THE BOX SIZE!'
cursor, xcntr, ycntr, /device
print, 'got the cursor! coordinates are ', xcntr, ycntr
print, ' '
print, 'sit back and relax...this will take a while............'
print, 'Getting impatient???? Next time, use a smaller Boxsize!'

;FIRST, define a small kernal to get rid of sky background...
kernal = fltarr(6,6) - 4./32.
kernal(2:3, 2:3) = 1.

;APPLY THIS KERNAL TO BOTH IMAGES...GENERATE NEW ONES CALLED 'TEST'
test1 = convol(im1, kernal)
test1 = test1 > 0.
test2 = convol(im2, kernal)
test2 = test2 > 0.

;DEFINE A NEW, SMALLER IMAGE FROM IM1, CENTERED ON THE CHOSEN STAR
xb = xcntr - xr
xe = xcntr + xr - 1
yb = ycntr - xr
ye = ycntr + xr - 1
;print, 'xb, xe = ', xb, xe
if (xe gt x512) or (xb lt 0) or (ye gt y512) or (yb lt 0) then begin 
print, 'THE STAR YOU SELECTED IS TOO NEAR THE EDGE!' 
print, 'USE EITHER A DIFFERENT STAR OR A SMALLER BOX SIZE!!'
;STOP
return
endif

im1mod = test1(xb:xe, yb:ye)
;im1mod = im1mod - (total(im1mod1)/(float(xe-xb+1)^2))
;help, im1mod

;DEFINE THE SECOND IMAGE AS THE KERNAL; MAKE IT HALF THE SIZE OF THE 
;NEW, SMALLER IM1:

kb = xr/2
ke = kb + xr - 1
;print, 'kb, ke = ', kb, ke

im2mod = test2(xb:xe, yb:ye)
krnl = im2mod( kb:ke, kb:ke)
;help, krnl

;GET THE RESULT OF THE CROSS CORRELATION; IT TAKES A WHILE, SO BEEP WHEN DONE:
cv11 = convol(im1mod, krnl)
print, string(7b)

print, 'max and indx nr of max are: ', max(cv11, indx), indx

;THE RULE IS...
; if a = fltarr( nrcols, nrrows)
; then
; one-d index = rownr*nrcols + colnr*nrrows

indxy = indx/(2*xr)
indxx = indx - indxy*2*xr
print, 'indxx, indxy = ', indxx, indxy
print, indxx, indxy

shiftx = indxx - xr
shifty = indxy - xr

print, 'to align these two images, you must:
print, 'shift image number 2 to RIGHT by: ', shiftx
print, 'shift image number 2 to UP : ', shifty

print, 'the most convenient way to do this is...'
print, 'image2 = shift(image2,', shiftx, ',', shifty,')'

print, ' '
read, choice, prompt = 'if you wish me to do this, type 1 now, otherwise 0: '

if (choice eq 1) then im2 = shift(im2, shiftx, shifty)

return
end
