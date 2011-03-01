;------------------------------------------------------------------
;+
; NAME:
;    rjc_overscan
;
; PURPOSE:
;    apply overscan correction to given image
;    
;
;
;  CALLING SEQUENCE:
;
;      rjc_overscan, image, biassec, order=order,
;           functype=functype
;
;  INPUTS:
; 
;     image : [nx, ny] image to be overscan corrected
;     biassec : image section to be used for overscan statistics
;               [4] = [x0,x1,y0,y1]
;
;  OPTIONAL INPUTS:
;     
;     order : order of fitting function to use (default = 3)
;     functype : type of fitting function to use
;     -right now i only fit a polynomial
;     niter : number of rejection iterations (default == 3)
;     doplot : plot the overscan and the fit
;  -
;-------------------------------------------------------------

pro im_overscan, image, biassec=biassec, order=order, functype=functype, $
  niter=niter, doplot=doplot, hdr=hdr, median=median, silent=silent
   
   IF NOT keyword_set(biassec) THEN BEGIN
      IF NOT keyword_set(hdr) THEN BEGIN
         splog, 'You must specify either a hdr or a biassec'
         stop
      ENDIF ELSE BEGIN
         bs = sxpar(hdr, 'BIASSEC')
         biassec = strsplit(bs, '[]:,', /extract)
         ;;Now correct for the fact that IRAF uses 1 indexing and
         ;;IDL uses zero
         biassec = biassec - 1.0
      ENDELSE
   ENDIF
   
   IF n_elements(image) EQ 0 THEN BEGIN
      splog, 'Please input an image with more than 0 pixels'
      stop
   endIF
   
   if (not keyword_set(silent)) then splog, 'Subtracting the overscan region.'

   IF NOT keyword_set(order) THEN order = 8
   IF NOT keyword_set(maxiter) THEN maxiter = 3
   
   nx = n_elements(image(*,0))
   ny = n_elements(image(0,*))
   
   biasimage = image[biassec[0]:biassec[1],biassec[2]:biassec[3]]

   if keyword_set(median) then begin

      djs_iterstat, biasimage, sigrej=3.0, mean=median
      overscan_image = median

   endif else begin
      
; Now collapse this is the x direction

      biasvector = djs_avsigclip(biasimage, 2)
;     biasvector = djs_avsigclip(biasimage, 1)
      pixels = findgen(n_elements(biasvector))
      
; Fit the bias section

      fit = poly_fit(pixels, biasvector, order, yfit)
      
      FOR iter = 0, maxiter -1 DO BEGIN
         residual = biasvector - yfit
         djs_iterstat, residual, median=median, sigma=sigma
         k = where(abs(residual-median) LT 5*sigma)
         mask = where(abs(residual-median) GT 5*sigma)
         fit = poly_fit(pixels(k), biasvector(k), order, yfit)
      endFOR
      
      djs_iterstat, biasvector, median=median, sigma=sigma1

      IF keyword_set(doplot) THEN BEGIN
         window, 0, xsize=600, ysize=600
         !P.multi = [0, 1, 2]
         !P.position=[0.2, 0.5, 0.9, 0.9]
         
         djs_plot, pixels, biasvector, /xstyle, /ystyle, $
           ytitle='Overscan Counts', $
           charsize=1.5, yrange=[median-8*sigma1, median+8*sigma1]
         djs_oplot, pixels, yfit, color='red', thick=2
         IF mask(0) GT -1 THEN $
           djs_oplot, pixels(mask), yfit(mask), ps=6, color='red'
         !P.position=[0.2, 0.1, 0.9, 0.4]
         djs_plot, pixels, biasvector-yfit, /xstyle, /ystyle, $
           xtitle='Pixel Number', ytitle='Residual Counts', charsize=1.5, $
           yrange=[-8*sigma, 8*sigma]
         IF mask(0) GT -1 THEN $
           djs_oplot, pixels(mask), residual(mask), ps=6, color='red'
         !P.multi=0
         !P.position=0 
      endIF
      
; Now that we have a good polynomial fit, we reconstruct the full image

      fitvector = biasvector * 0.0
      fit = reform(fit)
      FOR icoeff = 0, order DO $
        fitvector = fitvector + fit(icoeff) * pixels^float(icoeff)
      
; Now create an image the size of the original

      overscan_image = rebin(rebin(fitvector, 1, ny), nx, ny)

   endelse
      
; Now do the subtraction

   image = image - overscan_image

; if the hdr was given to the program, update it
   
   if keyword_set(hdr) then begin
      sxaddpar, hdr, 'OVERSCAN', '['+string(biassec[0]+1,format='(I0)')+':'+string(biassec[1]+1,format='(I0)')+','+$
        string(biassec[2]+1,format='(I0)')+':'+string(biassec[3]+1,format='(I0)')+']', ' overscan region (mean='+$
        strtrim(string(median),2)+')', before='HISTORY'
      sxaddhist, "'Overscan region subtracted "+hogg_iso_date()+"'", hdr
   endif
   
return
end
