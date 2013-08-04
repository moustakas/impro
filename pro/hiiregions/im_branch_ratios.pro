;+
; NAME:
;       IM_BRANCH_RATIOS()
;
; PURPOSE:
;       Compute and return the branching ratios among various
;       collisionally excited emission lines, in the low-density
;       limit. 
;
; INPUTS: 
;
; OPTIONAL INPUTS: 
;       dens - electron density (default 100 cm^-3); note that if the
;              density is made too large this program does not account
;              for collisional de-excitation 
;       temp - electron temperature (default 10,000 K)
;
; KEYWORD PARAMETERS: 
;
; OUTPUTS: 
;       branch - structure with the various branching ratios (see code
;                for details)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;       From B. Moore: "Many of the temperature diagnostic
;       ratios follow the form 
;
;         [F(Lower line #1) + F(Lower line #2)] / F(Upper line)
;
;       It is possible to calculate a diagnostic ratio with only one
;       of the lower lines present, but only if the program can be
;       made aware of the missing line.  Since the diagnostic ratios
;       are initialized, the alternate ratios are not easily accounted
;       for.  Instead, this procedure makes it possible to account for
;       the missing line by assuming its flux from the flux of the
;       other lower line scaled by the ratio of the lines given by the
;       atomic data.  This is done by calling the nlevel program for
;       default values of the density and temperature, and passing
;       these values back out in a structure.
;
;       Note that the fluxes are pulled out of the results from
;       IM_NLEVEL() by index, NOT wavelength.  Therefore, the user
;       must be diligent in checking the results if they change any
;       aspect of IM_NLEVEL()."
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2007-Nov-26, NYU: based almost entirely on
;          B. Moore's BRANCHCALC.
;
; Copyright (C) 2007, John Moustakas
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

function im_branch_ratios, dens=dens, temp=temp

    if (n_elements(temp) eq 0L) then temp = 1D4 ; 10^4 [K]
    if (n_elements(dens) eq 0L) then dens = 1D2 ; [cm^-3]

    branch = {$
      o_i:    0.0,$
      n_ii:   0.0,$
      si_iii: 0.0,$
      s_iii:  0.0,$
      ar_iii: 0.0,$
      o_iii:  0.0,$
      cl_iv:  0.0,$
      ne_iii: 0.0,$
      ar_v:   0.0,$
      ne_v:   0.0}

    level = im_nlevel('o_i',dens=dens,temp=temp)
    branch.o_i = level.emissivity[3,0]/level.emissivity[3,1] ;    F(6300)/F6363)

    level = im_nlevel('n_ii',dens=dens,temp=temp)
    branch.n_ii = level.emissivity[3,2]/level.emissivity[3,1] ;   F(6584)/F(6548)

    level = im_nlevel('si_iii',dens=dens,temp=temp)
    branch.si_iii = level.emissivity[3,0]/level.emissivity[2,0] ;   F(1883)/F(1892)

    level = im_nlevel('s_iii',dens=dens,temp=temp)
    branch.s_iii = level.emissivity[3,2]/level.emissivity[3,1] ;   F(9532)/F(9069)

    level = im_nlevel('ar_iii',dens=dens,temp=temp)
    branch.ar_iii = level.emissivity[3,0]/level.emissivity[3,1] ;   F(7135)/F(7751)

    level = im_nlevel('o_iii',dens=dens,temp=temp)
    branch.o_iii = level.emissivity[3,2]/level.emissivity[3,1] ;   F(5007)/F(4959)

    level = im_nlevel('cl_iv',dens=dens,temp=temp)
    branch.cl_iv = level.emissivity[3,2]/level.emissivity[3,1] ;   F(8045)/F(7531)

    level = im_nlevel('ne_iii',dens=dens,temp=temp)
    branch.ne_iii = level.emissivity[3,0]/level.emissivity[3,1] ;   F(3869)/F(3968)

    level = im_nlevel('ar_V',dens=dens,temp=temp)
    branch.ar_v = level.emissivity[3,2]/level.emissivity[3,1] ;   F(7005)/F(6435)

    level = im_nlevel('ne_V',dens=dens,temp=temp)
    branch.ne_v = level.emissivity[3,2]/level.emissivity[3,1] ;   F(3426)/F(3346)

return, branch
end
