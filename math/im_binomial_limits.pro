;+
; NAME:
;   IM_BINOMIAL_LIMITS()
;
; PURPOSE:
;   For a given number of events, calculate the Binomial upper and
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
;   Formulae taken from Gehrels 1986.
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2008 Oct 17, NYU - written based on a Perl script
;     kindly provided by G. Rudnick 
;
; Copyright (C) 2008, John Moustakas
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

function im_binomial_limits, n_i, conf

; for two given numbers of events calculate the binomial upper and
; lower limits for the 1-sided confidence interval of their ratio
; (n1/(n1 + n2).  Formulas are (18), (21), (26) from Gehrels, 1986, ApJ, 303,
; 336. 0.8413, 0.9772, and 0.9987 correspond to the usual 1-sided
; intervals for 1, 2, and 3 sigma gaussian probabilities

 my $n1 = shift;
 my $n2 = shift;
 my $conflim_d = shift;

 my ($imatch, $iconf);

 my ($diff_d, $conf_up, $conf_dn, $w_d, $h_d, $lambda_d, $epsilon_d);

 my ($p1u_d, $p1l_d, $p2u_d);

 ; allow only these limits
 my @conf_d = (0.8413, 0.90, 0.95, 0.975, 0.9772, 0.990, 0.995,
               0.9987, 0.999, 0.9995);

 ; 'S' parameter for these limits
 my @S_d = (1.0, 1.282, 1.645, 1.960, 2.00, 2.326, 2.576, 3.000,
            3.090, 3.291);

 my $ntot = $n1 + $n2;

 ; in this version just use the value of @conf_d that's closest to the
 ; input confidence limit
 ; find closest value
 for($iconf = 0, $diff_d = 1.e2; $iconf < @conf_d; $iconf++) {
   if(abs($conf_d[$iconf] - $conflim_d) < $diff_d) {
     $diff_d = abs($conf_d[$iconf] - $conflim_d);
     $imatch = $iconf;
   }
 }

 ; upper limits
 ; special cases
 if($n2 == 1) {
   $p1u_d = $conflim_d**(1.0/$ntot);
 } elsif ($n1 == 0) {
   $p1u_d = 1 - (1 - $conflim_d)**(1.0/$n2);
 } else {

   ; general case
   $lambda_d = ($S_d[$imatch]**2 - 3.0) / 6.0;

   $h_d = 2.0 * (1.0 / (2.0 * $n2 - 1.0) + 1.0 / (2.0 * $n1 + 1))**(-1);

   $w_d = $S_d[$imatch] * ($h_d + $lambda_d)**0.5 / $h_d +
     (1.0 / (2.0 * $n2 - 1.0) - 1.0 / (2.0 * $n1 + 1)) *
       ($lambda_d + 5.0/6.0 - 2.0 / (3.0 * $h_d));

   $epsilon_d = 0.64 * (1 - $S_d[$imatch]) * exp(-$n2);

   $p1u_d = ( ($n1 + 1.0) * exp(2.0 * $w_d) + $epsilon_d * $n2) /
     ( ($n1 + 1.0) * exp(2.0 * $w_d) + $n2);
 }

 ; lower limits
 ; recalculate but switch n2 and n1
 ; special cases
 if($n1 == 1) {
   $p2u_d = $conflim_d**(1.0/$ntot);
 } elsif ($n2 == 0) {
   $p2u_d = 1 - (1 - $conflim_d)**(1.0/$n1);
 } else {

   ; general case
   $lambda_d = ($S_d[$imatch]**2 - 3.0) / 6.0;

   $h_d = 2.0 * (1.0 / (2.0 * $n1 - 1.0) + 1.0 / (2.0 * $n2 + 1))**(-1);

   $w_d = $S_d[$imatch] * ($h_d + $lambda_d)**0.5 / $h_d +
     (1.0 / (2.0 * $n1 - 1.0) - 1.0 / (2.0 * $n2 + 1)) *
       ($lambda_d + 5.0/6.0 - 2.0 / (3.0 * $h_d));

   $epsilon_d = 0.64 * (1 - $S_d[$imatch]) * exp(-$n1);

   $p2u_d = ( ($n2 + 1.0) * exp(2.0 * $w_d) + $epsilon_d * $n1) /
     ( ($n2 + 1.0) * exp(2.0 * $w_d) + $n1);
 }
 $p1l_d = 1 - $p2u_d;

 $conf_dn = $p1l_d;
 $conf_up = $p1u_d;

 return $conf_dn, $conf_up;
    
    if (n_elements(n_i) ne 1L) or (n_elements(conf) ne 1L) then begin
       splog, 'N_I and CONF must be scalars'
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

; eqs. (10) and (14) from Geherls (1986)    
    limit_up_grid = n_i + SS*sqrt(n_i+1.0D) + (SS^2.0+2.0D)/3.0D
    if (n_i lt 1.0) then limit_dn_grid = 0.0D else begin
       limit_dn_grid = n_i * (1.0D - 1.0D/(9.0D*n_i) - SS/(3.0D*sqrt(n_i)) + $
         beta*n_i^(gamma^3.0D))
    endelse

    limit_up = interpol(limit_up_grid,conf_grid,conf)
    limit_dn = interpol(limit_dn_grid,conf_grid,conf)
    
return, [limit_dn,limit_up]
end
