;+
; NAME:
;       WRITE_TEMDEN_LOOKUP_TABLE
;
; PURPOSE:
;       Write a look-up table for use with IM_TEMDEN().
;
; INPUTS:
;
; OPTIONAL INPUTS:
;       max_temp   - maximum temperature to consider (default 25,000 K)
;       min_temp   - minimum temperature to consider (default 5,000 K)
;       step_temp  - temperature step size
;       max_dens   - maximum density to consider (default 10,000 cm^-3)
;       min_dens   - minimum density to consider (default 1 cm^-3)
;       step_dens  - density step size
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       Binary FITS lookup table.
;
; COMMENTS:
;       Currently only the following ratios are written out.
;
;       Temperature-sensitive line-ratios:
;          O_III = (4959+5007)/4363
;          S_III = (9069+9532)/6312
;          O_II  = (3726+3729)/7325
;          N_II  = (6548+6584)/5755
; 
;       Density-sensitive line-ratios:
;          O_II = 3726/3729
;          S_II = 6716/6731
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2007 Nov 27, NYU - written
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

pro write_temden_lookup_table, table, max_temp=max_temp, min_temp=min_temp, $
  step_temp=step_temp, max_dens=max_dens, min_dens=min_dens, $
  step_dens=step_dens, force=force

; current grid of ions; see B. Moore's routine DIAGNOSTICS to add more
; ions; also note that the temp grid, below, is optimized for the
; [OIII] 4959+5007/4363, which varies rapidly with decreasing temp;
; other ions might behave differently

    ionlist = ['s_ii','o_ii','n_ii','s_iii','o_iii']
    nion = n_elements(ionlist)
    
; set up the temp and dens grid, where the current limits correspond
; to reasonable limits for the atomic data

    if (n_elements(min_temp) eq 0L) then min_temp = 3D3 else min_temp = min_temp>3D3 ; [K]
    if (n_elements(max_temp) eq 0L) then max_temp = 2.5D4 else max_temp = max_temp<2.5D4 ; [K]
    if (n_elements(step_temp) eq 0L) then step_temp = 0.005D ; 0.002D

    if (n_elements(min_dens) eq 0L) then min_dens = 1.0D else min_dens = min_dens>1.0D ; [cm^-3]
    if (n_elements(max_dens) eq 0L) then max_dens = 1D5 else max_dens = max_dens<1D5
    if (n_elements(step_dens) eq 0L) then step_dens = 0.05D ; 0.02D ; 0.01D

    grid_temp = 10.0D^(dindgen((alog10(max_temp)-alog10(min_temp))/step_temp+1L)*step_temp+alog10(min_temp))
    grid_dens = 10.0D^(dindgen((alog10(max_dens)-alog10(min_dens))/step_dens+1L)*step_dens+alog10(min_dens))

    ndens = n_elements(grid_dens)
    ntemp = n_elements(grid_temp)
    splog, 'NTEMP = '+string(ntemp,format='(I0.0)')
    splog, 'NDENS = '+string(ndens,format='(I0.0)')
    
; initialize the lookup table

    table = {$
      date:                        '', $
      grid_temp:            grid_temp, $
      grid_dens:            grid_dens, $
      sii_ratio:      fltarr(ntemp,ndens), $ ; 6716/6731
      oii_dens_ratio: fltarr(ntemp,ndens), $ ; 3726/3729
      oii_temp_ratio: fltarr(ntemp,ndens), $ ; (3726+3729)/7325
      nii_ratio:      fltarr(ntemp,ndens), $ ; (6548+6584)/5755
      siii_ratio:     fltarr(ntemp,ndens), $ ; (9069+9532)/6312
      oiii_ratio:     fltarr(ntemp,ndens)}   ; (4959+5007)/4363
    
; loop on each ION/RATIO; note that I'm using INTERPOLATE in
; combination with FINDEX because when the ratio hits against the
; limits of GRID_TEMP and/or GRID_DENS the last physically meaningful
; temperature and/or density is returned, rather than INTERPOL's
; crappy extrapolation ; when adding more ions, see B. Moore's
; DIAGINIT for the appropriate emissivity indices
    
    t0 = systime(1)
    for ii = 0L, nion-1L do begin
       splog, 'Computing line-ratios for ion '+strupcase(ionlist[ii])
       level = im_nlevel(ionlist[ii],dens=grid_dens,temp=grid_temp,matrix=matrix)
       case strlowcase(ionlist[ii]) of
          's_ii': table.sii_ratio = level.emissivity[2,0]/level.emissivity[1,0]     ; 6716/6731
          'o_ii': begin
             table.oii_dens_ratio = level.emissivity[2,0]/level.emissivity[1,0]     ; 3726/3729
             table.oii_temp_ratio = (level.emissivity[2,0]+level.emissivity[1,0])/$ ; (3726+3729)/7325
               ((level.emissivity[4,1]+level.emissivity[3,1])+(level.emissivity[4,2]+level.emissivity[3,2]))
          end
          'n_ii' : table.nii_ratio  = (level.emissivity[3,1]+level.emissivity[3,2])/level.emissivity[4,3]  ; (6548+6584)/5755
          's_iii': table.siii_ratio = (level.emissivity[3,1]+level.emissivity[3,2])/level.emissivity[4,3]  ; (9069+9532)/6312
          'o_iii': table.oiii_ratio = (level.emissivity[3,1]+level.emissivity[3,2])/level.emissivity[4,3]  ; (4959+5007)/4363
          's_iii': table.siii_ratio = (level.emissivity[3,2]+level.emissivity[3,1])/level.emissivity[4,3]  ; (9069+9532)/6312
          else: splog, 'No entry for ion '+ionlist[ii]
       endcase
;      plot, table.grid_temp, table.oiii_ratio[*,2], charsize=2                                 
    endfor 
    splog, format='("Total time to compute lookup table = ",G0," minutes.")', (systime(1)-t0)/60.0

; write out    
    
    tablefile = getenv('IMPRO_DIR')+'/etc/temden_table.fits'
    if file_test(tablefile,/regular) and (not keyword_set(force)) then begin
       splog, 'Overwrite '+tablefile+' [Y/N]? '
       cc = get_kbrd(1)
       if (strupcase(cc) eq 'Y') then begin
          splog, 'Writing '+tablefile
          table.date = hogg_iso_date()
          mwrfits, table, tablefile, /create
       endif else splog, 'Table not over-written.'
    endif else begin
       splog, 'Writing '+tablefile
       table.date = hogg_iso_date()
       mwrfits, table, tablefile, /create
    endelse

return
end
