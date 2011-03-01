;------------------------------------------------------------------------
; polar_inv_proj.pro      [IDL]                  20-JUL-95     03-MAY-96
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

pro polar_inv_proj,coordx,coordy,size,edge,hemi,l,b,good

  x = (coordx-edge+0.5)*(2./size)-1.
  y = (coordy-edge+0.5)*(2./size)-1.

  rho = sqrt(x*x+y*y)
  good = (rho LT sqrt(2.0))
  rho = good*rho
  theta = acos(1.0-rho*rho)
  rho = 0
  b = hemi*( !pi/2.-theta)  
  l = -hemi*atan(y/x)
  y = 0
  l = (x GT 0.)*l+(x LT 0.)*(l+!pi)

; Convert from radians to degrees: !radeg =  57.2958

  b = b*!radeg
  l = l*!radeg
  
  l = (l > 0.) + (l LT 0.)*(l+360.)

return
end

