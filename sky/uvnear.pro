; J. Moustakas
; March 1999

; determines if catalog stars are near object centers
; inputs:  angle (degrees), unit vectors of objects and starlist
; outputs: byte array of ones and zeros


PRO uvnear, theta, uvobj, uvcat, near

   theta = theta*!dtor

   dot = double(uvobj#transpose(uvcat))

   near = dot GT cos(theta)
   near = total(near, 1) NE 0

END 

