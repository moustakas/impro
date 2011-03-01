function imcut,image,xc,yc,imsize,help=help

;+
; Function to cut out a square subsection out of an image,
; with the user specifying the center of the subsection.
; Will return a non-square image if the specified center is
; too close to the edges.
;
; INPUTS
;	image	original image
;	xc,yc	center location of the subsection 
;	imsize	size of the subsection (256 default)
;
; RETURNS
;	the image subsection
;
; HISTORY
; Written by MCL(UCB): 10/01/95
;
;-

if keyword_set(help) then begin
	print,'function imcut(image,xc,yc,imsize,help=help)'
	retall
endif

if n_elements(imsize) eq 0 then imsize=256

dd = fix(imsize/2)
sz = size(image)
x0 = round(xc) - dd > 0
x1 = round(xc) + dd < (sz(1)-1)
y0 = round(yc) - dd > 0
y1 = round(yc) + dd < (sz(2)-1)
if not(odd(imsize)) then begin
	y1 = y1-1
	x1 = x1-1
endif
if (x1-x0+1) lt imsize or (y1-y0+1) lt imsize then begin
;	print,imsize,dd
;	print,x1,x0,x1-x0+1
;	print,y1,y0,y1-y0+1
	print,'* too close to edge *'
endif

return,image(x0:x1,y0:y1)

end
