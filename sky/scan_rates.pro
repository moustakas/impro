;+
; NAME:
;	SCAN_RATES()
;
; PURPOSE:
;	Compute drift scan rates and offsets for integrated long-slit
;	spectroscopy.  
;
; CALLING SEQUENCE:
;       scan = scan_rates(pa,dec,scanlen,exptime,object=,$
;          offset_along_slit=,offset_perp_slit=,nscan=,$
;          /verbose,/doplot)
;
; INPUTS:
;	pa      - slit position angle (measured from North to East,
;                 where North is up and East is to the left)
;                 [degrees].  nominally PA should correspond to the
;                 major axis of the galaxy
;	dec     - object declination [string degrees-minutes-seconds);
;                 used to compute the correction to the RA scan rate 
;	scanlen - total scan length perpendicular to the slit [arcsec] 
;	exptime - total exposure time [s]
;
; OPTIONAL INPUTS:
;       object      - object name
;	offset_along_slit - offset the scan *along* the slit by OFFSET_ALONG_SLIT
;                     [arcsec] (default 0.0) 
;       nscan       - number of scans ("passes") of the slit over the
;                     galaxy (default 10) equal to EXPTIME / SCANTIME
;
; KEYWORD PARAMETERS:
;	verbose     - print the SCAN data structure to STDOUT 
;
; OUTPUTS:
;	scan        - output data structure with the following fields
;          object     - input
;          ra_offset  - slit offset in RA (East is positive) [arcsec]
;          dec_offset - slit offset in DEC (North is positive) [arcsec]
;          ra_rate    - scan rate in RA (East is positive) [arcsec/s]
;          dec_rate   - scan rate in DEC (North is positive) [arcsec/s]
;          scantime   - duration of one scan ("pass") [s]
;          pa         - input
;          exptime    - input
;          scanlen    - input
;          scanlen    - input
;          nscan      - input
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;       The cardinal directions are given by (N,E) such that North is
;       up and positive (y axis) and East is measured positive **along
;       the x axis**.  The coordinate system of the galaxy is (x',y')
;       such that y' always coincides with the galaxy's major axis.  
;
;       The scan length is the length of the scan *perpendicular* to
;       the slit (and, presumably, the galaxy's major axis).  The slit
;       offset is *along* the slit (along the galaxy major axis).  
;
; EXAMPLE:
;	Compute the drift rates and offsets for an object at 12 hours
;	declination, requiring a slit position angle of 135 degrees,
;	and an offset along the slit of 80 arcsec EAST.  Compute for
;	an exposure time of 1200 seconds and 20 drift scan passes.
;	
;	IDL> scan = scan_rates(135,'12:00:00',55.0,1200.0,$
;       IDL>    offset_along_slit=80.0,nscan=20,/verbose)
;
; PROCEDURES USED:
;	IM_OFFSET_AND_ROTATE(), IM_HMS2DEC(), PRINT_STRUCT
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2002 February 9, U of A
;       jm03may22uofa - rewritten completely and more transparently 
;       jm03may27uofa - added OBJECT optional input
;       jm05mar09uofa - added DOPLOT keyword; changed OFFSET_ALONG_SLIT to
;                       OFFSET_ALONG_SLIT; added OFFSET_PERP_SLIT
;                       optional input 
;-

function scan_rates, pa, dec, scanlen, exptime, object=object, $
  offset_along_slit=offset_along_slit, offset_perp_slit=offset_perp_slit, $
  nscan=nscan, verbose=verbose, doplot=doplot

    npa = n_elements(pa)
    ndec = n_elements(dec)
    nscanlen = n_elements(scanlen)
    nexptime = n_elements(exptime)
    
    if (npa*ndec*nscanlen*nexptime)[0] eq 0L then begin
       print, 'Syntax - scan = scan_rates(pa,dec,scanlen,exptime,object=,$'
       print, '   offset_along_slit=,offset_perp_slit=,nscan=,/verbose,/doplot)'
       return, -1L
    endif

    if (npa ne ndec) and (npa ne nscanlen) and (npa ne nexptime) then begin
       print, 'Dimensions of PA, DEC, SCANLEN, and EXPTIME are not the same.'
       return, -1
    endif    

    noffset_along_slit = n_elements(offset_along_slit)
    if noffset_along_slit eq 0L then offset_along_slit = replicate(0.0,npa) else begin
       if noffset_along_slit ne npa then begin
          print, 'Dimensions of PA and OFFSET_ALONG_SLIT are not the same.'
          return, -1L
       endif
    endelse

    noffset_perp_slit = n_elements(offset_perp_slit)
    if noffset_perp_slit eq 0L then offset_perp_slit = replicate(0.0,npa) else begin
       if noffset_perp_slit ne npa then begin
          print, 'Dimensions of PA and OFFSET_PERP_SLIT are not the same.'
          return, -1L
       endif
    endelse

    nnscan = n_elements(nscan) 
    if nnscan eq 0L then nscan = replicate(10.0,npa) else begin
       nscan = float(nscan)
       if nnscan ne npa then begin
          print, 'Dimensions of PA and NSCAN are not the same.'
          return, -1L
       endif
    endelse
    
    nobject = n_elements(object) 
    if nobject eq 0L then begin
       object = strarr(npa)
       for k = 0L, npa-1L do object[k] = 'Object_'+string(k,format='(I3.3)')
    endif else begin
       if nobject ne npa then begin
          print, 'Dimensions of PA and OBJECT are not the same.'
          return, -1L
       endif
    endelse 
    
    negtime = where(exptime le 0.0,nneg)
    if nneg ne 0L then begin
       print, 'EXPTIME must be greater than zero.'
       return, -1L
    endif
    
; call this routine recursively
    
    if npa gt 1L then begin 

       for i = 0L, npa-1L do begin
          
          scan1 = scan_rates(pa[i],dec[i],scanlen[i],exptime[i],object=object[i],$
            offset_along_slit=offset_along_slit[i],offset_perp_slit=offset_perp_slit[i],$
            nscan=nscan[i],verbose=0,doplot=doplot)
          if (i eq 0L) then scan = scan1 else scan = [ [scan], [scan1] ]
          
       endfor

       scan = reform(scan)
       if keyword_set(verbose) then print_struct, scan
       return, scan

    endif
    
; compute the offsets and the scan rates
    
    scantime = exptime / float(nscan)

    offsets = im_offset_and_rotate([-scanlen/2.0+offset_perp_slit,+offset_along_slit],pa)
;   offsets = im_offset_and_rotate([-scanlen/2.0,+offset_along_slit],pa)
    rates = im_offset_and_rotate([scanlen/scantime,0.0],pa)

    decra = im_hms2dec(dec)*!dtor
    rates[0] = rates[0] / cos(decra)
    
    EPS = (machar()).eps ; avoid round-off errors

    offsets = offsets * (abs(offsets) gt 1E-4)
    rates   = rates * (abs(rates) gt 1E-4)

    scan = {$
      object:            object,           $
      ra_offset:        -offsets[0],       $ ; RA is positive along the -x axis
      dec_offset:        offsets[1],       $
      ra_rate:          -rates[0],         $
      dec_rate:          rates[1],         $
      scantime:          scantime,         $
      offset_along_slit: offset_along_slit,$
      offset_perp_slit:  offset_perp_slit, $
      pa:                pa,               $
      exptime:           exptime,          $
      scanlen:           scanlen,          $
      nscan:             nscan}

    if keyword_set(verbose) then print_struct, scan

    if keyword_set(doplot) then begin

       if !d.window ne 0L then im_window, 0, xratio=0.5, /square
       
       width = fix(((abs(scan.ra_offset)>abs(scan.dec_offset))+scanlen)>8*60.0) ; [arcsec]
       axis = findgen(10*width+1)/10.0-width/2.0

       plot, [0], [0], /nodata, xrange=minmax(axis), yrange=minmax(axis), $
         xthick=2.0, ythick=2.0, charsize=2.0, charthick=2.0, xsty=3, ysty=3, $
         xtitle=textoidl('\Delta\alpha [arcsec]'), ytitle=textoidl('\Delta\delta [arcsec]'), $
         title=scan.object, xgridstyle=1, ygridstyle=1, yticklen=1, $
         xticklen=1
;      oplot, !x.crange, [0,0], line=2, thick=1
;      oplot, [0,0], !y.crange, line=2, thick=1
       plots, 0.0, 0.0, ps=7, syms=2, thick=3

; visualize the slits

       im_oplot_box, scan.scanlen, 3.5*60.0, scan.pa, line=0, thick=2, $
         xoffset=scan.offset_perp_slit[0], yoffset=scan.offset_along_slit[0], $
         corners=corners
       djs_oplot, [corners[0,0],corners[0,3]], [corners[1,0],corners[1,3]], line=0, color='red', thick=2
       djs_oplot, [corners[0,1],corners[0,2]], [corners[1,1],corners[1,2]], line=0, color='red', thick=2

       plotsym, 8, 2, /fill
       plots, offsets[0], offsets[1], ps=8

       scanlabel = [$
         'Scan Length = '+strtrim(string(scan.scanlen,format='(I3)'),2)+'"', $
         'Slit PA = '+strtrim(string(scan.pa,format='(I3)'),2), $
         'Perp Slit Offset = '+strtrim(string(scan.offset_perp_slit,format='(I5)'),2)+'"', $
         'Along Slit Offset = '+strtrim(string(scan.offset_along_slit,format='(I5)'),2)+'"']
       legend, textoidl(scanlabel), /left, /top, box=0, charsize=1.3, $
         charthick=2.0, textcolor=djs_icolor('black'), /clear

       cc = get_kbrd(1)

    endif

return, scan
end
