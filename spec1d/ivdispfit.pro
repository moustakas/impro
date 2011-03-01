;+
; NAME:
;   ivdispfit
;
; PURPOSE:
;   Compute velocity dispersions for galaxy spectra.
;
; CALLING SEQUENCE:
;   vdans = vdispfit(objflux, objivar, [ objloglam, hdr=, zobj=, npoly=, $
;    eigenfile=, eigendir=, columns=, sigma=, sigerr=, yfit=, $
;    plottitle=, /doplot, /debug ])
;
; INPUTS:
;   objflux    - Galaxy spectrum (spectra); array of [NPIX,NGALAXY].
;   objivar    - Galaxy inverse variance; array of [NPIX,NGALAXY].
;
; OPTIONAL INPUTS:
;   objloglam  - Log-10 wavelengths; this can be either an NPIX vector
;                if all the galaxy spectra have the same wavelength mapping,
;                or an array with the same dimensions as OBJFLUX.
;                Either OBJLOGLAM or HDR must be specified.
;   hdr        - FITS header from which to read COEFF0, COEFF1 for the
;                wavelength mapping.
;                Either OBJLOGLAM or HDR must be specified.
;   zobj       - Redshift for each galaxy; default to 0.
;   npoly      - Number of polynomial terms to append to eigenspectra;
;                default to 5.
;   eigenfile  - File name for eigenvectors; default to 'spEigenVdisp*.fits'
;   eigendir   - Directory name for EIGENFILE; default to
;                '$IDLSPEC2D_DIR/templates'
;   columns    - Column numbers of the eigenspectra image to use in the
;                PCA fit; default to all columns.
;   plottitle  - Title of plot (if /DOPLOT is set).
;   doplot     - If set, then make plots.
;   debug      - If set, then wait for keystroke after plot.
;
; OUTPUTS:
;   vdans      - Output structure [NGALAXY] with the following elements:
;                vdisp : Velocity dispersion in km/sec.
;                vdisp_err : Error for VDISP in km/sec.
;                vdispchi2 : Minimum chi^2
;                vdispnpix : Number of pixels overlapping the templates
;                            and used in the fits
;                vdispdof : Degrees of freedom = the number of pixels
;                           overlapping the templates minus the number of
;                           templates minus the number of polynomial terms
;                           minus 1 (the last 1 is for the velocity dispersion)
;                vdisptheta : ???
;
; OPTIONAL OUTPUTS:
;   yfit       - Best-fit template (actually, the one with the closest
;                velocity dispersion to the best-fit sigma); wavelengths
;                outside of those available with the templates have their
;                values set to zero [NPIX,NGALAXY]
;
; COMMENTS:
;   Note that the wavelength spacing in the galaxy and stellar template spectra
;   must be the same.
;
;   We currently mask within +/- 280 km/sec of the following wavelengths
;   that could have emission lines:
;     linelist = [3725.94, 3727.24, 3970.072, 4101.73, 4340.46, $
;      4861.3632, 4958.911, 5006.843, 6300.32, 6548.05, 6562.801, $
;      6583.45, 6716.44, 6730.82]
;
;   The constructed over-sampled and smoothed eigenspectra are stored
;   in a common block between calls if the eigen-vector file is the same.
;
; EXAMPLES:
;
; BUGS:
;
; DATA FILES:
;   $IDLSPEC2D_DIR/templates/spEigenVdisp*.fits
;
; PROCEDURES CALLED:
;   airtovac
;   combine1fiber
;   computechi2()
;   djs_filepath()
;   find_nminima
;   mrdfits()
;   poly_array()
;   splog
;   sxpar()
;
; INTERNAL SUPPORT ROUTINES:
;   create_vdans()
;   vdisp_gconv()
;
; REVISION HISTORY:
;   13-Mar-2001  Written by D. Schlegel, Princeton
;------------------------------------------------------------------------------
; Create output structure
function create_vdans, nstar

   vdans = create_struct( $
    name = 'VDANS'+strtrim(string(nstar),1), $
    'vdisp'      , 0.0, $
    'vdisp_err'  , 0.0, $
    'vdispchi2'  , 0.0, $
    'vdispnpix'  ,  0L, $
    'vdispdof'   ,  0L, $
    'vdisptheta' , fltarr(nstar) )

   return, vdans
end
;------------------------------------------------------------------------------
function vdisp_gconv, x, sigma

   ; Special case for no smoothing
   if (sigma EQ 0) then return, x

   khalfsz = round(4*sigma+1)
   xx = findgen(khalfsz*2+1) - khalfsz

   kernel = exp(-xx^2 / (2*sigma^2))
   kernel = kernel / total(kernel)

   return, convol(x, kernel, /center, /edge_truncate)
end

;------------------------------------------------------------------------------
function ivdispfit, objflux, objivar, objloglam, $
 hdr=hdr, zobj=zobj, npoly=npoly, $
 eigenfile=eigenfile, eigendir=eigendir, columns=columns, $
 sigma=sigma, sigerr=sigerr, yfit=yfit, $
 plottitle=plottitle, doplot=doplot1, debug=debug

   common com_vdispfit, bigflux, bigloglam, bigmask, nsamp, bigsig, $
    nbigpix, nsig, dsig, nstar, lastfile

   if (NOT keyword_set(objloglam) AND NOT keyword_set(hdr)) then $
    message, 'Must specify either OBJLOGLAM or HDR!'

   if (n_elements(npoly) EQ 0) then npoly = 5
   dims = size(objflux, /dimens)
   npixobj = dims[0]
   if (size(objflux, /n_dimen) EQ 1) then nobj = 1 $
    else nobj = dims[1]
   if (NOT keyword_set(zobj)) then zobj = fltarr(nobj)
   if (NOT keyword_set(lastfile)) then lastfile = ''
   if (NOT keyword_set(plottitle)) then plottitle = ''

   ; Plot if either /DOPLOT or /DEBUG is set.
   if (keyword_set(doplot1)) then doplot = doplot1
   if (keyword_set(debug)) then doplot = 1

   ;---------------------------------------------------------------------------
   ; If multiple object flux vectors exist, then call this routine recursively.

   if (nobj GT 1) then begin
      sigma = fltarr(nobj)
      sigerr = fltarr(nobj)
      if (arg_present(yfit)) then yfit = fltarr(npixobj,nobj)
      lamdims = size(objloglam, /n_dimens)
      for iobj=0, nobj-1 do begin
         if (lamdims EQ 1) then thisloglam = objloglam $
          else if (lamdims EQ 2) then thisloglam = objloglam[*,iobj]
         vdans1 = vdispfit(objflux[*,iobj], objivar[*,iobj], $
          thisloglam, hdr=hdr, zobj=zobj[iobj], npoly=npoly, $
          eigenfile=eigenfile, eigendir=eigendir, columns=columns, yfit=yfit1)
         if (iobj EQ 0) then vdans = vdans1 $
          else vdans = [[vdans], [vdans1]]
         if (keyword_set(yfit)) then yfit[*,iobj] = yfit1
      endfor
      return, vdans
   endif

   ;----------
   ; Determine the wavelength mapping for the object spectra,
   ; which are the same for all of them.

   if (NOT keyword_set(objloglam)) then begin
      objloglam0 = sxpar(hdr, 'COEFF0')
      objdloglam = sxpar(hdr, 'COEFF1')
      objloglam = objloglam0 + dindgen(npixobj) * objdloglam
   endif else begin
      objdloglam = objloglam[1] - objloglam[0]
   endelse
   restloglam = objloglam - alog10(1 + zobj[0]) ; De-redshift this!

   ;---------------------------------------------------------------------------
   ; Generate the over-sampled eigen-templates for the stellar spectra.
   ; This is saved in a common block between calls.

   ;----------
   ; Find the template file matching EIGENFILE

   if (NOT keyword_set(eigenfile)) then $
    eigenfile = 'spEigenVdisp*.fits'
   if (NOT keyword_set(eigendir)) then $
    eigendir = concat_dir(getenv('IDLSPEC2D_DIR'), 'templates')
   allfiles = findfile(djs_filepath(eigenfile, root_dir=eigendir), count=ct)
   if (ct EQ 0) then $
    message, 'Unable to find EIGENFILE matching '+eigenfile
   thisfile = allfiles[ (reverse(sort(allfiles)))[0] ]

   if (NOT keyword_set(bigflux) OR (thisfile NE lastfile)) then begin
      splog, 'Computing grid of dispersion templates from ' + thisfile

      nsamp = 10
      nsig = 35
      dsig = 25.0
      bigsig = findgen(nsig) * dsig ; in km/sec

      ;----------
      ; Read the stellar templates

      eflux = mrdfits(thisfile, 0, ehdr)
      naxis1 = sxpar(ehdr,'NAXIS1')
      nstar = sxpar(ehdr,'NAXIS2') > 1
      eloglam0 = sxpar(ehdr, 'COEFF0')
      edloglam = sxpar(ehdr, 'COEFF1')
      eloglam = eloglam0 + dindgen(naxis1) * edloglam

      ; Pixel size in km/sec for these oversampled (smaller) pixels
      cspeed = 2.99792458e5
      pixsz = (10.^(edloglam)-1) * cspeed / nsamp

      ;----------
      ; Re-samples to higher resolution by a factor of NSAMP

      nbigpix = (naxis1 - 1) * nsamp + 1
      bigloglam = eloglam0 + dindgen(nbigpix) * edloglam / nsamp
      bigflux = fltarr(nbigpix, nstar, nsig)

      splog, 'Oversampling eigentemplates'
      for istar=0, nstar-1 do begin
         ; Burles counter...
         print, format='("Template ",i5," of ",i5,a1,$)', $
          istar+1, nstar, string(13b)

         combine1fiber, eloglam, eflux[*,istar], $
          newloglam=bigloglam, newflux=tmpflux, maxiter=0
         bigflux[*,istar,0] = tmpflux
         if (istar EQ 0) then bigmask = tmpflux NE 0
      endfor

      ;----------
      ; Generate array of broadened templates

      splog, 'Broadening eigentemplates'
      for isig=1, nsig-1 do begin
         for istar=0, nstar-1 do begin
            ; Burles counter...
            print, format='("Template ",i5," of ",i5,a1,$)', $
             isig*nstar+istar+1, nsig*nstar, string(13b)

            bigflux[*,istar,isig] = $
             vdisp_gconv(bigflux[*,istar,0], bigsig[isig]/pixsz)
         endfor
      endfor

      ;----------
      ; Mask out the first few and last few pixels, since those would not
      ; have been smoothed properly.  The masked region is 350 km/sec
      ; for a native binning of 70 km/sec.

      bigmask[0:4*nsamp-1] = 0
      bigmask[nbigpix-4*nsamp:nbigpix-1] = 0

      ;----------
      ; Mask out emission lines, setting bigmask=0 near these wavelengths.

      linelist = [3725.94, 3727.24, 3970.072, 4101.73, 4340.46, $
       4861.3632, 4958.911, 5006.843, 6300.32, 6548.05, 6562.801, $
       6583.45, 6716.44, 6730.82]
      vaclist = linelist
      airtovac, vaclist
      vaclist = alog10(vaclist)
      mwidth = 6.e-4 ; Mask out any pixels within +/- 420 km/s
      for iline=0, n_elements(vaclist)-1 do $
       bigmask = bigmask AND (bigloglam LT vaclist[iline] - mwidth $
        OR bigloglam GT vaclist[iline] + mwidth)

      lastfile = thisfile
   endif else begin
      print, 'Using previously cached velocity dispersion templates'
   endelse

   ;----------
   ; Create the output structure

   vdans = create_vdans(nstar+npoly)

   ;----------
   ; Find the pixel numbers to use from the object and the templates

   ; Find the sub-pixel shifts in the object
;   subshift = round(((bigloglam[0]-restloglam[0]) / objdloglam MOD 1) * nsamp)
   ; Bug fix since the IDL MOD function does the wrong thing for negatives.
   subdum = (bigloglam[0]-restloglam[0]) / objdloglam
   subshift = round((subdum-floor(subdum)) * nsamp)
   indx = subshift + nsamp * lindgen(nbigpix/nsamp)

   if (max(restloglam) LT min(bigloglam[indx]) $
    OR min(restloglam) GT max(bigloglam[indx])) then begin
;      splog, 'No wavelength overlap with template'
      vdans.vdisp = 0.0
      vdans.vdisp_err = -4L
      yfit = fltarr(npixobj)
      return, vdans
   endif

   if (restloglam[0] LT bigloglam[indx[0]]) then begin
      ipixt0 = 0L
      junk = min(abs(restloglam - bigloglam[indx[0]]), ipixo0)
   endif else begin
      ipixo0 = 0L
      junk = min(abs(bigloglam[indx] - restloglam[0]), ipixt0)
   endelse

   npixcomp = (npixobj - ipixo0 + 1) < (n_elements(indx) - ipixt0)

   indxo = ipixo0 + lindgen(npixcomp) ; Indices for object spectrum
   indxt = indx[ipixt0 + lindgen(npixcomp)] ; Indices for template spectra

   ;----------
   ; Add more eigen-templates that represent polynomial terms.

   if (keyword_set(npoly)) then $
    polyflux = poly_array(npixcomp,npoly)

   ;----------
   ; Select which eigenvectors to use (default to all)

   if (n_elements(columns) EQ 0) then iuse = lindgen(nstar) $
    else iuse = columns
   nuse = n_elements(iuse)

   ;----------
   ; Fit for chi^2 at each possible velocity dispersion

   chi2arr = fltarr(nsig)

   objsmall = objflux[indxo]
   sqivar = sqrt( objivar[indxo] ) * bigmask[indxt]

   acoeffarr = fltarr(nuse+npoly,nsig)
   for isig=0, nsig-1 do begin
      eigenflux = bigflux[indxt,iuse,isig]
      if (keyword_set(npoly)) then eigenflux = [[eigenflux], [polyflux]]
      chi2arr[isig] = computechi2(objsmall, sqivar, eigenflux, acoeff=acoeff)
      acoeffarr[*,isig] = acoeff
   endfor

   ;----------
   ; Fit for the dispersion value at the minimum in chi^2

;   findchi2min, bigsig, chi2arr, minchi2, sigma, sigerr, $
;    plottitle=plottitle, doplot=doplot, debug=debug
   ; Use only the 3 points nearest the minimum for the fit.
   ; If the minimum is at a dispersion of zero, then duplicate the
   ; next point as a negative dispersion value simply for the benefit
   ; of computing an error, and to prevent an error code from being
   ; generated in the call to FIND_NMINIMA.
   junk = min(chi2arr, imin)
   if (imin GT 0) then $
    sigma = find_nminima(chi2arr, bigsig, $
     width=1.5*dsig, ypeak=minchi2, xerr=sigerr, $
     errcode=errcode, plottitle=plottitle, xtitle=textoidl('\sigma [km/s]'), $
     doplot=doplot, debug=debug) $
   else $
    sigma = find_nminima([chi2arr[1],chi2arr], [-bigsig[1],bigsig], $
     width=1.5*dsig, ypeak=minchi2, xerr=sigerr, $
     errcode=errcode, plottitle=plottitle, doplot=doplot, debug=debug)
   vdans.vdisp = sigma > 0 ; Numerical round-off can push this negative
   ; Set VDISP_ERR to the error-code if it is non-zero
   vdans.vdisp_err = sigerr * (errcode EQ 0) + errcode
   vdans.vdispchi2 = minchi2
   vdans.vdispnpix = npixcomp
   vdans.vdispdof = npixcomp - nstar - npoly - 1 ; One dof is for the vel. disp.

   ;----------
   ; Return the best-fit template (actually, the one with the closest
   ; velocity dispersion to the best-fit sigma).

   junk = min(abs(bigsig - sigma), isig)
   if (arg_present(yfit)) then begin
      eigenflux = bigflux[indxt,iuse,isig]
      if (keyword_set(npoly)) then eigenflux = [[eigenflux], [polyflux]]
      yfit = fltarr(npixobj)
      yfit[indxo] = acoeffarr[*,isig] ## eigenflux
   endif
   vdans.vdisptheta = acoeffarr[*,isig]

   ;----------
   ; If the best-fit value is at the maximum dispersion value tested,
   ; then we don't really know the answer and should set the error
   ; to a large value.

   if (sigma GE max(bigsig)) then begin
      vdans.vdisp_err = -3L
   endif

   return, vdans
end
;------------------------------------------------------------------------------
