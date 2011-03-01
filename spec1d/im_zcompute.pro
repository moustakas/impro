;+
; NAME:
;   zcompute
;
; PURPOSE:
;   Compute relative redshift of object(s) vs. eigen-templates.
;
; CALLING SEQUENCE:
;   zans = zcompute(objflux, objivar, starflux, [starmask, nfind=, $
;    poffset=, pspace=, pmin=, pmax=, mindof=, width=, minsep=, $
;    plottitle=, /doplot, /debug, /verbose ]
;
; INPUTS:
;   objflux    - Object fluxes [NPIXOBJ,NOBJ]
;   objivar    - Object inverse variances [NPIXOBJ,NOBJ]
;   starflux   - Eigen-template fluxes [NPIXSTAR,NTEMPLATE]
;
; OPTIONAL INPUTS:
;   starmask   - Eigen-template mask; 0=bad, 1=good [NPIXSTAR]
;   nfind      - Number of solutions to find per object; default to 1.
;   poffset    - Offset between all objects and templates, in pixels.
;                A value of 10 indicates that STARFLUX begins ten pixels
;                after OBJFLUX, i.e. OBJFLUX[i+10] = STARFLUX[i] for the
;                case when the relative redshift should be zero.  If the
;                wavelength coverage of the templates is larger, then the
;                value of ZOFFSET will always be negative.
;                [Scalar or vector of size NOBJ]
;   pspace     - The spacing in redshifts to consider; default to 1 [pixels];
;                [Scalar or vector of size NOBJ]
;   pmin       - The smallest redshift to consider [pixels].
;                [Scalar or vector of size NOBJ]
;   pmax       - The largest redshift to consider [pixels].
;                [Scalar or vector of size NOBJ]
;   mindof     - Minimum number of degrees of freedom for a valid fit;
;                default to 10.
;   width      - Parameter for FIND_NMINIMA(); default to 3 * PSPACE.
;   minsep     - Parameter for FIND_NMINIMA(); default to the same as WIDTH.
;   plottitle  - Title of plot (if /DOPLOT is set).
;   doplot     - If set, then make plots.
;   debug      - If set, then wait for keystroke after plot.
;   verbose    - If set, then log using SPLOG instead of PRINT.
;
; OUTPUTS:
;   zans       - Output structure [NOBJECT,NFIND] with the following elements:
;                z : The relative redshift.
;                z_err : Error in the redshift, based upon the local quadratic
;                        fit to the chi^2 minimum. 
;                chi2 : Fit value for the best (minimum) chi^2
;                dof : Number of degrees of freedom, equal to the number of
;                      pixels in common between the object and templates
;                      minus the number of templates.
;                theta : Mixing angles [NTEMPLATE].  These are computed at the
;                        nearest integral redshift, e.g. at ROUND(ZOFFSET).
;
; COMMENTS:
;   Fits are done to chi^2/DOF, not to chi^2.
;
; EXAMPLES:
;
; BUGS:
;
; DATA FILES:
;   $IDLSPEC2D_DIR/etc/TEMPLATEFILES
;
; PROCEDURES CALLED:
;   computechi2()
;   find_nminima()
;   splog
;
; INTERNAL SUPPORT ROUTINES:
;   create_zans()
;
; REVISION HISTORY:
;   10-Jul-2000  Written by D. Schlegel, Princeton
;------------------------------------------------------------------------------
; Create output structure
function create_zans, nstar, nfind

   zans1 = create_struct( $
    name = 'ZANS'+strtrim(string(nstar),1), $
    'z'     , 0.0, $
    'z_err' , 0.0, $
    'chi2'  , 0.0, $
    'dof'   ,  0L, $
    'theta' , fltarr(nstar) )

   return, replicate(zans1, nfind)
end

;------------------------------------------------------------------------------
function im_zcompute, objflux, objivar, starflux, starmask, nfind=nfind, $
 poffset=poffset, pspace=pspace, pmin=pmin, pmax=pmax, $
 mindof=mindof, width=width, minsep=minsep, $
 plottitle=plottitle, doplot=doplot1, debug=debug, verbose=verbose, $
 bestlag=bestlag

   if (NOT keyword_set(nfind)) then nfind = 1
   if (NOT keyword_set(pspace)) then pspace = 1
   if (NOT keyword_set(width)) then width = 3 * pspace
   if (NOT keyword_set(plottitle)) then plottitle = ''

   ; Plot if either /DOPLOT or /DEBUG is set.
   if (keyword_set(doplot1)) then doplot = doplot1
   if (keyword_set(debug)) then doplot = 1

   ;---------------------------------------------------------------------------
   ; Check dimensions of object vectors

   ndim = size(objflux, /n_dimen)
   dims = size(objflux, /dimens)
   npixobj = dims[0]
   if (ndim EQ 1) then nobj = 1 $
    else if (ndim EQ 2) then nobj = dims[1] $
    else message, 'OBJFLUX is neither 1-D or 2-D'

   if total(abs(size(objflux, /dimens)-size(objivar, /dimens))) NE 0 $
    OR size(objflux, /n_dimen) NE size(objivar, /n_dimen) THEN  $
    message, 'Dimensions of OBJFLUX and OBJIVAR do not match'

   ;---------------------------------------------------------------------------
   ; If multiple object vectors exist, then call this routine recursively.

   if (nobj GT 1) then begin
      t0 = systime(1)
      for iobj=0L, nobj-1 do begin
         if (n_elements(poffset) EQ 1) then poffset1 = poffset $
          else poffset1 = poffset[iobj]
         if (n_elements(pspace) EQ 1) then pspace1 = pspace $
          else pspace1 = pspace[iobj]
         if (n_elements(pmin) EQ 1) then pmin1 = pmin $
          else pmin1 = pmin[iobj]
         if (n_elements(pmax) EQ 1) then pmax1 = pmax $
          else pmax1 = pmax[iobj]

         zans1 = zcompute(objflux[*,iobj], objivar[*,iobj], $
          starflux, starmask, nfind=nfind, $
          poffset=poffset1, pspace=pspace1, pmin=pmin1, pmax=pmax1, $
          mindof=mindof, width=width, minsep=minsep, $
          plottitle=plottitle+ ' -- Object #'+strtrim(string(iobj+1),2), $
          doplot=doplot, debug=debug)
         if (iobj EQ 0) then zans = zans1 $
          else zans = [[zans], [zans1]]
         if (keyword_set(verbose)) then $
          splog, 'Object #', iobj, '  Elap time=', systime(1)-t0, $
           ' (sec)  z=', zans1[0].z, ' (pix)' $
         else $
          print, 'Object #', iobj, '  Elap time=', systime(1)-t0, $
           ' (sec)  z=', zans1[0].z, ' (pix)'
      endfor
      return, zans
   endif

   ;---------------------------------------------------------------------------

   if (NOT keyword_set(mindof)) then mindof = 10
   if (NOT keyword_set(width)) then width = 3 * pspace
   if (NOT keyword_set(minsep)) then minsep = width

   ndim = size(starflux, /n_dimen)
   dims = size(starflux, /dimens)
   npixstar = dims[0]
   if (ndim EQ 1) then nstar = 1 $
    else if (ndim EQ 2) then nstar = dims[1] $
    else message, 'STARFLUX is neither 1-D or 2-D'

   if (NOT keyword_set(starmask)) then begin
      starmask = bytarr(npixstar) + 1
   endif else begin
      if (n_elements(starmask) NE npixstar) then $
       message, 'Dimensions of STARFLUX and STARMASK do not match'
   endelse

   if (n_elements(poffset) EQ 0) then poffset = 0
   if (n_elements(poffset) GT 1) then $
    message, 'ZOFFSET must be a scalar'
   pixoffset = round(poffset)

   if (n_elements(pmin) EQ 0) then $
    pmin = -2 * ((npixobj < npixstar) + 1) + pixoffset
   if (n_elements(pmax) EQ 0) then $
    pmax = pmin + 2 * ((npixobj < npixstar) - 1)

   if (n_elements(pmin) GT 1) then $
    message, 'PMIN must be a scalar'
   if (n_elements(pmax) GT 1) then $
    message, 'PMAX must be a scalar'
   if (pmin GT pmax) then $
    message, 'PMAX must be greater than or equal to PMIN'

   nlag = ((pmax - pmin + 1) / pspace) > 1
   lags = - lindgen(nlag) * pspace + pixoffset - long(pmin) ; must be integers

   chi2arr = fltarr(nlag)
   dofarr = fltarr(nlag)
   thetaarr = fltarr(nstar,nlag)
   zans = create_zans(nstar, nfind)

   ;---------------------------------------------------------------------------

   sqivar = sqrt(objivar)
   objmask = objivar NE 0

   for ilag=0L, nlag-1 do begin
      j1 = lags[ilag] < (npixstar-1L)
      if (j1 LT 0) then i1 = -j1 $
       else i1 = 0L
      j1 = j1 > 0L
      j2 = npixstar-1 < (npixobj+j1-i1-1L)
      i2 = i1 + j2 - j1
      chi2arr[ilag] = computechi2( objflux[i1:i2], $
       sqivar[i1:i2] * starmask[j1:j2], starflux[j1:j2,*], $
       acoeff=acoeff, dof=dof)
      dofarr[ilag] = dof
      thetaarr[*,ilag] = acoeff

      ; Burles counter of lag number...
;     print, format='("Lag ",i5," of ",i5,a1,$)', ilag, nlag, string(13b)
   endfor

   ;-----
   ; Fit this chi2/DOF minimum

   indx = where(dofarr GT mindof, ct)
   if (ct GT width) then begin
;      xpeak1 = find_npeaks(-chi2arr[indx]/dofarr[indx], lags[indx], $
;       nfind=nfind, minsep=minsep, width=width, $
;       ypeak=ypeak1, xerr=xerr1, npeak=npeak)
;      zans[0:npeak-1].chi2 = -ypeak1
      xpeak1 = find_nminima(chi2arr[indx], lags[indx], $
       dofarr=dofarr[indx], nfind=nfind, minsep=minsep, width=width, $
       ypeak=ypeak1, xerr=xerr1, npeak=npeak, errcode=errcode, $
       plottitle=plottitle, xtitle='Lag [pixels]', doplot=doplot)
      bestlag = round(xpeak1) ; jm05jun28uofa
      zans[0:npeak-1].z = poffset - xpeak1
      ; Set Z_ERR equal to the error-code if it is non-zero
      zans[0:npeak-1].z_err = xerr1 * (errcode EQ 0) + errcode
      zans[0:npeak-1].chi2 = ypeak1
      for ipeak=0L, npeak-1 do begin
         junk = min(abs(lags-xpeak1[ipeak]), ilag)
         zans[ipeak].dof = dofarr[ilag]
         zans[ipeak].theta = thetaarr[*,ilag]
      endfor
      zans[0:npeak-1].chi2 = zans[0:npeak-1].chi2 * zans[0:npeak-1].dof

      ; Wait for a keystroke...
      if (keyword_set(debug)) then begin
         print, 'Press any key...'
         cc = strupcase(get_kbrd(1))
      endif

   endif else if (ct GE 1) then begin
      zans[0].chi2 = -min(-chi2arr[indx]/dofarr[indx], ilag)
      zans[0].z = poffset - lags[indx[ilag]]
      zans[0].z_err = 0
      zans[0].dof = dofarr[indx[ilag]]
      zans[0].chi2 = zans[0].chi2 * zans[0].dof
      zans[0].theta = thetaarr[*,indx[ilag]]
   endif

   return, zans
end
;------------------------------------------------------------------------------
