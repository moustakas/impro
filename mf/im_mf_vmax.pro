;+
; NAME:
;   PRIMUS_MF_VMAX()
;
; PURPOSE:
;   Compute a Vmax-weighted mass function.
;
; INPUTS: 
;   mass - stellar mass for each galaxy [NGAL]
;   oneovervmax - 1/Vmax for each object (h^3 Mpc^-3) [NGAL] 
;
; OPTIONAL INPUTS: 
;   minmass - minimum stellar mass; see FULLBIN (default min(mass)) 
;   binsize - logarithmic width of each stellar mass bin 
;     (default 0.1 dex)
;
; KEYWORD PARAMETERS: 
;
; OUTPUTS: 
;   mf_data - the results
;     ngal - total number of galaxies
;     nbins - number of mass bins
;     binmass - stellar mass bin [NBINS]
;     phi - number density (h^3 Mpc^-3) [NBINS] 
;     errphi - error on PHI (h^3 Mpc^-3) [NBINS]
;     fullbin - Boolean array indicating whether all the galaxies in
;       each bin have MASS>MINMASS [NBINS]
;     number - number of objects in each bin [NBINS]
;
; COMMENTS:
;   Basically a glorified wrapper on IM_HIST1D().
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Aug 25, UCSD - based on MF_VMAX()
;-

function im_mf_vmax, mass, oneovervmax, minmass=minmass, $
  binsize=binsize, histmin=histmin, histmax=histmax, $
  nolog=nolog, debug=debug

    ngal = n_elements(mass)
    if (ngal eq 0L) then begin
       doc_library, 'im_mf_vmax'
       return, -1
    endif
    if (ngal ne n_elements(oneovervmax)) then begin
       splog, 'Dimensions of MASS and ONEOVERVMAX do not agree'
       return, -1
    endif
    zero = where(oneovervmax le 0.0,nzero)
    if (nzero ne 0L) then message, 'Some ONEOVERVMAX are <= 0!'

; round the bin centers to two decimal points
    if (n_elements(binsize) eq 0) then binsize = 0.1 ; dex
    if (n_elements(minmass) eq 0) then minmass = min(mass)
    if (n_elements(histmin) eq 0) then $
      histmin = ceil((min(mass)+binsize/2.0)*100.0)/100.0-binsize/2.0

; generate weighted histogram    
    phi = im_hist1d(mass,oneovervmax,binsize=binsize,obin=binmass,$
      binedge=0,omin=omin,omax=omax,h_err=phierr,histmin=histmin,$
      histmax=histmax,reverse_indices=rev)

; compute FULLBIN and NUMBER
    nbins = n_elements(phi)
    fullbin = intarr(nbins)+1 ; start with all full
    number = intarr(nbins)
    for ii = 0L, nbins-1 do begin
       number[ii] = rev[ii+1]-rev[ii] ; number of objects per bin
       if (rev[ii] eq rev[ii+1]) then fullbin[ii] = 0 else $
         fullbin[ii] = total(mass[rev[rev[ii]:rev[ii+1]-1]] lt minmass) eq 0.0
    endfor

; pack into a structure
    null = fltarr(50>n_elements(phi))-999.0
    mf_data = {$
      ngal:         ngal,$
      nbins:       nbins,$
      fullbin: fix(null),$
      number:  fix(null),$
      mass:         null,$
      phi:          null,$
      phierr:       null}
    mf_data.fullbin[0:nbins-1] = fullbin
    mf_data.number[0:nbins-1] = number
    mf_data.mass[0:nbins-1] = binmass

    if keyword_set(nolog) then begin
       mf_data.phi[0:nbins-1] = phi
       mf_data.phierr[0:nbins-1] = phierr
    endif else begin
       mf_data.phi[0:nbins-1] = alog10(phi+(fullbin eq 0))*(fullbin ne 0)
       mf_data.phierr[0:nbins-1] = phierr/(phi+(fullbin eq 0))/alog(10)*(fullbin ne 0)
    endelse
;   if total(phi le 0) gt 0 then stop
    
; QAplot    
    if keyword_set(debug) then begin
       ploterror, binmass, alog10(phi), phierr/phi/alog(10.0), $
         ps=10, xsty=3, ysty=3
       niceprint, binmass, phi
       cc = get_kbrd(1)
    endif

return, mf_data
end
