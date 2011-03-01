;------------------------------------------------------------------------
; polar.pro         [IDL]                        20-JUL-95     03-MAY-96
;------------------------------------------------------------------------
; WRITTEN BY: Doug Finkbeiner                                  UCB - DPF 
;                                                                     
;        FOR: The University of California, Berkeley
;             Department of Astronomy
;             Campbell Hall
;             Berkeley, CA  94720
;------------------------------------------------------------------------  
; PURPOSE:
;	Create Polar projections of DIRBE maps
; 
;	Call Create_DIRBE_pixmap
;------------------------------------------------------------------------  

;------------------------------------------------------------------------  
; Procedure POLAR_INV_PROJ
;------------------------------------------------------------------------  
; INPUTS:
;	coordx	: x coordinate 
;	coordy	: y coordinate 
;	size	: diameter of "good" region (default = 1024)
; 				for big projection, 4096
;       hemi    : +1 = ngp    -1 = sgp
;------------------------------------------------------------------------  
; KEYWORDS:
;------------------------------------------------------------------------  
; OUTPUTS:
;	l       : galactic longitude
;       b	: galactic latitude
;    	good	: 1=pixel in projection, 0=pixel off edge of projection
;------------------------------------------------------------------------  

pro polar_proj_ns,l,b,xn,yn,xs,ys

sinb=sin(b/!radeg)
; NGP
rho=sqrt(1.-sinb)
xn=rho*cos(l/!radeg)
yn=-rho*sin(l/!radeg)
; SGP
rho=sqrt(1.+sinb)
xs=rho*cos(l/!radeg)
ys=+rho*sin(l/!radeg)

return
end



;------------------------------------------------------------------------  
; POLAR.pro                *** END  OF  FILE ***               UCB - DPF
;------------------------------------------------------------------------  
