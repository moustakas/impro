;+
; NAME:
;       BALMER_TEMPERATURE()
;
; PURPOSE:
;       Compute the temperature dependence of the Balmer decrement, as
;       implemented in RETURN_TBALMER().
;
; CALLING SEQUENCE:
;       data - balmer_temperature(npoly=,coeff_data=,/doplot)
;
; INPUTS:
;
; OPTIONAL INPUTS:
;       npoly - order of the polynomial fit between temperature and
;               the various Balmer decrements (default 2, quadratic) 
;
; KEYWORD PARAMETERS:
;       doplot - generate debugging plots and wait for a keystroke 
;
; OUTPUTS:
;       data - output data structure giving the interpolated Balmer
;              decrements at the grid of temperatures
;
; OPTIONAL OUTPUTS:
;       coeff_data - coefficients of the polynomial fits
;
; PROCEDURES USED:
;       PLOTSYM, LEGEND, 
;
; COMMENTS:
;       References: Storey & Hummer (1995); Hummer & Storey (1987);
;       Brocklehurst (1971); Caplan & Deharveng (1986)
;    
;       Coefficients from Brocklehurst (1971):    
;
;          temp = [0.5,1.0,2.0]*1E4
;          log_temp = alog(temp/1E4)
;          
;          hahb = [3.032,2.859,2.744]    ; [Ha/Hb]
;          hghb = [0.4585,0.4685,0.4755] ; [Hg/Hb]
;          hbhg = 1.0/hghb               ; [Hb/Hg]
;          hahg = hahb/hghb              ; [Ha/Hg]
;          hdhb = [0.2516,0.2591,0.2643] ; [Hd/Hb]
;          hahd = hahb/hdhb              ; [Ha/Hd]
;
;       Coefficients are from Storey & Hummer (1995) for n_e =
;       100/cm3. 
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2003 Aug 06 - written
;       jm04jan22uofa - updated and generalized, converted to a
;                       function
;       jm04nov03uofa - documented and cleaned up
;
; Copyright (C) 2003-2004, John Moustakas
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

function balmer_temperature, npoly=npoly, coeff_data=coeff_data, doplot=doplot

    if n_elements(npoly) eq 0L then npoly = 2 ; quadratic
    
; initialize the output data structure
    
    temp = [0.30,0.50,0.75,1.00,1.25,1.50,2.00]*1E4
    log_temp = alog(temp/1E4)

    data = {$
      temperature: 0.0, $
      hahb:        0.0, $
      hghb:        0.0, $
      hbhg:        0.0, $
      hahg:        0.0, $
      hdhb:        0.0, $
      hahd:        0.0, $
      hghd:        0.0}
    data = replicate(data,n_elements(temp))
    data.temperature = temp

    ha = [1.048E-24,6.687E-25,4.625E-25,3.536E-25,2.860E-25,2.398E-25,1.807E-25]
    hb = [3.265E-25,2.199E-25,1.579E-25,1.235E-25,1.014E-25,8.600E-26,6.579E-26]
    hg = [1.468E-25,1.008E-25,7.334E-26,5.784E-26,4.777E-26,4.067E-26,3.128E-26]
    hd = [7.997E-26,5.530E-26,4.043E-26,3.198E-26,2.647E-26,2.257E-26,1.739E-26]
    
    data.hahb = ha/hb ; [Ha/Hb]
    data.hghb = hg/hb ; [Hg/Hb]
    data.hbhg = hb/hg ; [Hb/Hg]
    data.hahg = ha/hg ; [Ha/Hg]
    data.hdhb = hd/hb ; [Hd/Hb]
    data.hahd = ha/hd ; [Ha/Hd]
    data.hghd = hg/hd ; [Hg/Hd]

    if keyword_set(doplot) then begin
       
;   struct_print, data

       coeff_data = {$
         coeff_hahb:  fltarr(npoly+1), $
         coeff_hghb:  fltarr(npoly+1), $
         coeff_hbhg:  fltarr(npoly+1), $
         coeff_hahg:  fltarr(npoly+1), $
         coeff_hdhb:  fltarr(npoly+1), $
         coeff_hahd:  fltarr(npoly+1)}

       bigtemp = findgen((25000.0-1000.0)/1000.0+1)*1000.0+1000.0

       window, 0, xsize=500, ysize=500
       plotsym, 0, 1, /fill

; ---------------------------------------------------------------------------    
; fit the Ha/Hb temperature dependence
; ---------------------------------------------------------------------------    

       coeff_data.coeff_hahb = poly_fit(alog10(data.temperature),alog10(data.hahb),npoly)
       yfit = poly(alog10(bigtemp),coeff_data.coeff_hahb)

       plot, alog10(data.temperature), alog10(data.hahb), ps=8, xsty=3, ysty=3, $
         charsize=2.0, charthick=2.0, xtitle='log Temperature [K]', $
         ytitle='log Balmer Decrement'
       oplot, alog10(bigtemp), yfit, line=0, thick=2.0
       legend, textoidl('H\alpha/H\beta'), /right, /top, box=0, charsize=2.0, charthick=2.0
       cc = get_kbrd(1)
       
; ---------------------------------------------------------------------------    
; fit the Ha/Hg temperature dependence
; ---------------------------------------------------------------------------    

       coeff_data.coeff_hahg = poly_fit(alog10(data.temperature),alog10(data.hahg),npoly)
       yfit = poly(alog10(bigtemp),coeff_data.coeff_hahg)

       plot, alog10(data.temperature), alog10(data.hahg), ps=8, xsty=3, ysty=3, $
         charsize=2.0, charthick=2.0, xtitle='log Temperature [K]', $
         ytitle='log Balmer Decrement'
       oplot, alog10(bigtemp), yfit, line=0, thick=2.0
       legend, textoidl('H\alpha/H\gamma'), /right, /top, box=0, charsize=2.0, charthick=2.0
       cc = get_kbrd(1)

; ---------------------------------------------------------------------------    
; fit the Hb/Hg temperature dependence
; ---------------------------------------------------------------------------    

       coeff_data.coeff_hbhg = poly_fit(alog10(data.temperature),alog10(data.hbhg),npoly)
       yfit = poly(alog10(bigtemp),coeff_data.coeff_hbhg)

       plot, alog10(data.temperature), alog10(data.hbhg), ps=8, xsty=3, ysty=3, $
         charsize=2.0, charthick=2.0, xtitle='log Temperature [K]', $
         ytitle='log Balmer Decrement'
       oplot, alog10(bigtemp), yfit, line=0, thick=2.0
       legend, textoidl('H\beta/H\gamma'), /right, /top, box=0, charsize=2.0, charthick=2.0
       cc = get_kbrd(1)
       
; ---------------------------------------------------------------------------    
; fit the Ha/Hd temperature dependence
; ---------------------------------------------------------------------------    

       coeff_data.coeff_hahd = poly_fit(alog10(data.temperature),alog10(data.hahd),npoly)
       yfit = poly(alog10(bigtemp),coeff_data.coeff_hahd)

       plot, alog10(data.temperature), alog10(data.hahd), ps=8, xsty=3, ysty=3, $
         charsize=2.0, charthick=2.0, xtitle='log Temperature [K]', $
         ytitle='log Balmer Decrement'
       oplot, alog10(bigtemp), yfit, line=0, thick=2.0
       legend, textoidl('H\alpha/H\delta'), /right, /top, box=0, charsize=2.0, charthick=2.0
       cc = get_kbrd(1)

    endif
       
return, data
end    
    
