;+
; NAME:
;	Disk_Region
; PURPOSE:
;	Function returns a matrix of 0 and 1, 1 indicating the inside of circle.
; CALLING:
;	disk_mask = Disk_Region( center, SIZE=sz, RADIUS=radius, COUNT=npix )
; INPUTS:
;	center =
; KEYWORDS:
;	SIZE =
;	RADIUS =
;	COUNT =
; OUTPUTS:
;
; PROCEDURE:
;
; HISTORY:
;	Written: Frank Varosi, HSTX @ NASA/GSFC, 1997.
;-

function Disk_Region, center, SIZE=sz, RADIUS=radius, COUNT=npix

	if N_elements( sz ) LE 0 then begin
		npix = 0
		return,(-1)
	   endif

	if N_elements( center ) ne sz[0] then begin
		npix = 0
		return,(-1)
	   endif

	d2 = [(findgen( sz[1] ) - center[0])^2] # replicate( 1, sz[2] ) $
		+ replicate( 1, sz[1] ) # [(findgen( sz[2] ) - center[1])^2]

return, where( d2 LE float(radius)^2, npix )
end
