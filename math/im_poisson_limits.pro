;+
; NAME:
;   IM_POISSON_LIMITS()
;
; PURPOSE:
;   For a given number of events, calculate the Poisson upper and
;   lower limits for a given 1-sided confidence interval. 
;
; INPUTS: 
;   n_i  - number of (counting) events or samples (must be a scalar) 
;   conf - desired confidence interval; for example, 0.8413, 0.9772,
;     and 0.9987 correspond to the usual 1-sided intervals for 1-, 2-, 
;     and 3-sigma Gaussian probabilities (must be a scalar)

; OPTIONAL INPUTS: 
;
; KEYWORD PARAMETERS: 
;
; OUTPUTS: 
;   limit - two-element array corresponding to the lower- and
;     upper- confidence limit
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Formulae taken from Gehrels (1986).
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2008 Oct 17, NYU - written based on a Perl script
;     kindly provided by G. Rudnick 
;   jm09mar18nyu - use a default value (0.8413) for CONF
;   jm11mar31ucsd - updated to use the more accurate eq (9) for
;     computing upper limits
;
; Copyright (C) 2008-2010, John Moustakas
; 
; This program is free software; you can redistribute it and/or modify 
; it under the terms of the GNU General Public License as published by 
; the Free Software Foundation; either version 2 of the License, or
; (at your option) any later version. 
; 
; This program is distributed in the hope that it will be useful, but 
; WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
; General Public License for more details. 
;-

function im_poisson_limits, n_i, conf

    if (n_elements(n_i) ne 1L) then begin
       splog, 'N_I must be a scalar'
       return, -1.0
    endif

    if (n_elements(conf) eq 0L) then conf = 0.8413D else $
      if (n_elements(conf) ne 1L) then begin
       splog, 'CONF must be a scalar'
       return, -1.0
    endif
    
    if (n_i lt 0.0) then begin
       splog, 'N_I must be positive!'
       return, -1.0
    endif
    
; define a grid of confidence intervals
    conf_grid = [0.8413D, 0.90D, 0.95D, 0.975D, 0.9772D, $
      0.990D, 0.995D, 0.9987D, 0.999D, 0.9995D]
    if (conf lt min(conf_grid)) or (conf gt max(conf_grid)) then begin
       splog, 'CONF must be in the interval ['+string(min(conf_grid),$
         format='(G0.0)')+'-'+string(max(conf_grid),format='(G0.0)')+']'
       return, -1.0
    endif

; 'S' parameter corresponding to CONF_GRID
    SS = [1.0D, 1.282D, 1.645D, 1.960D, 2.00D, 2.326D, $
      2.576D, 3.000D, 3.090D, 3.291D]

; additional parameters needed to compute the lower limit
    beta = [0.0D, 0.010D, 0.031D, 0.058D, 0.062D, 0.103D, $
      0.141D, 0.222D, 0.241D, 0.287D]
    gamma = [0.0D, -4.0D, -2.50D, -2.22D, -2.19D, -2.07D, $
      -2.0D, -1.88D, -1.85D, -1.80D]

; equations (9) and (14); note that eq (9) becomes slightly less
; inaccurate than eq (7) for large CONF (>0.99) and small n
;   limit_up_grid = n_i + SS*sqrt(n_i+1.0D) + (SS^2.0+2.0D)/3.0D ; eq (10) is inaccurate!!
    limit_up_grid = (n_i+1)*(1-1D/(9*(n_i+1))+SS/(3*sqrt(n_i+1)))^3
    if (n_i lt 1.0) then limit_dn_grid = 0.0D else begin
       limit_dn_grid = n_i*(1.0D - 1.0D/(9.0D*n_i) - SS/(3.0D*sqrt(n_i)) + $
         beta*n_i^gamma)^3.0D
    endelse
    
    limit_up = interpol(limit_up_grid,conf_grid,conf)
    if (limit_dn_grid[0] eq 0.0) then limit_dn = limit_dn_grid else $
      limit_dn = interpol(limit_dn_grid,conf_grid,conf)
    
return, [limit_dn,limit_up]
end
