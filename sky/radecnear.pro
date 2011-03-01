; J. Moustakas
; March 1999

; determines if catalog stars are near object centers
; inputs:  angle (degrees), ra/dec of objects and starlist
; outputs: byte array of ones and zeros


PRO radecnear, theta, raobj, decobj, racat, deccat, near

   theta = theta*!dtor

   uvobj = ll2uv([[raobj],[decobj]]) ; (n,3) array
   uvcat = ll2uv([[racat],[deccat]]) ; (n,3) array
   dot = double(uvobj#transpose(uvcat))

   near = dot GT cos(theta)
   near = total(near, 1) NE 0

END 

