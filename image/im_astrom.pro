;+
; NAME:
;   dfastrom
;
; PURPOSE:
;   Compute astrometric solution for arbitrary images.
;
; CALLING SEQUENCE:
;   dfastrom, input, [hdr=, starcat=, maxnstar=, /ps, $
;    angle=, /sdss, success=, infostr= ]
;
; INPUTS:
;   input      - Either an image or an input FITS file name with an image.
;                If a second HDU is present in this file, it is assumed to
;                contain the inverse variance, which we use for masking.
;
; OPTIONAL INPUTS:
;   hdr        - Input header if INPUT is an image; required for the 2.5-m
;   starcat    - Optional star catalogue.  If one is not given, USNO SA2.0
;                is used. 
;   maxnstar   - If DAO find finds more than maxnstar stars, cull the
;                list.  Never use fewer than 1/4 of the found stars,
;                however.
;
; KEYWORDS:
;   ps         - Make PS plots of QA info.
;   angle      - Rotation angle (initial guess)
; 
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   hdr        - Modified FITS header with astrometric solution
;   infostr    - String for log
;   success    - set to unity if the attempt to compute and stuff the
;                astrometry was successful
;
; COMMENTS:
;   We assume that we know the scale and rotation well enough, then solve
;   for the X,Y offsets by correlating with Tycho stars.
;
;   If INPUT is the name of a FITS file, then update its header to
;   contain the astrometric solution.
;
; BUGS:
;
; PROCEDURES CALLED:
;   astr_25m
;   djs_iterstat
;   djs_modfits
;   find
;   gsssadxy
;   gsssputast
;   gsssxyad
;   mkhdr
;   offset_from_pairs
;   pt_errlog
;   pt_tweak_astr
;   rdss_fits()
;
; REVISION HISTORY:
;   17-May-2000  Written by D. Schlegel, Princeton.
;   13-Sep-2000  Now adds contents of log_info structure to FITS
;                       header -DPF
;   13-Sep-2000  Only call DAO find once - DPF
;   13-Sep-2000  added maxnstar keyword - DPF (in Schlegel's account)
;   14-Sep-2000  modified to keep ncull brightest stars if list is
;                  trimmed
;   02-Nov-2000  added seeing estimate and FWHM_* FITS keywords
;   04-Nov-2000  GSSS implemented; 3rd order fit; new
;           offset_from_pairs criterion; new gsainit - DPF
;   05-Nov-2000  added angle, sdss keywords - DPF
;   26-Nov-2000  Made findsig and binsz configurable.  Also fixed
;        astr_25m so 2.5m astrometry works. 
;-
;------------------------------------------------------------------------------

pro im_astrom, image, starcat=starcat, datasec=datasec, maxnstar=maxnstar, ps=ps, $
  gsa=gsa, success=success, infostr=infostr, $
  searchrad=searchrad, maxsep=maxsep, radial=radial, binsz=binsz, $
  findsig=findsig, dmax=dmax


; Start with success bit set to false
  success= 0B


  if (size(gsa,/tname) NE 'STRUCT') then message, 'Must pass gsa structure with initial guess'

; -------- Set defaults
  if NOT keyword_set(searchrad) then  searchrad = .5  ; radius [deg] of offset space to search
  if NOT keyword_set(maxsep)    then  maxsep    = 10. ; criterion for match on first pass [pix]
  if NOT arg_present(radial)    then  radial    = 0B  ; fit radial distortion terms
  if NOT keyword_set(binsz)     then  binsz     = 5   ; pixels - binsize for offset_from_pairs
  if NOT keyword_set(findsig)   then  findsig   = 5.  ; DAO find sigma cut
  if NOT keyword_set(dmax)      then  dmax      = 200
  if NOT keyword_set(infostr)   then  infostr=''

  nstar_min = 10                ; minimum number of stars for acceptable fit
  sdss = 1
  xvec0 = gsa.amdx

   ; start with central image

   imsize = size(image, /dimens)
   cut0 = imsize/4     ; 512
   cut1 = imsize/4*3-1 ; 1535
   if (keyword_set(sdss)) then begin
      cut0 = [0,0]
      cut1 = imsize - 1
   endif
   im_small = image[cut0[0]:cut1[0], cut0[1]:cut1[1]]

; bin down image for statistics
   bf = 8
   isb = imsize/bf
   imsam = rebin(image[0:isb[0]*bf-1, 0:isb[1]*bf-1], isb[0], isb[1], /sample)
   imsig = djsig(imsam)


; -------- Set params for DAO find  
  sharplim = [.2, 1.0]
  roundlim = [-1.0, 1.0]*1.5
  hmin = imsig*findsig

  find, image, imx, imy, flux, sharp, round, hmin, 3., roundlim, $
    sharplim, /silent

; when there are no good stars returned, an array of zero flux stars
; may be returned (I think this is a DAO find bug).  -DPF
   count= 0
   if n_elements(flux) GE 1 then goodstar = where(flux GT 0, count)
   IF count LT 20 THEN BEGIN 
      errstr = string('---> DAO Find:', count, ' stars found -- giving up!', $
                      format='(A,I4,A)')
      pt_errlog, infostr, errstr
      return
   ENDIF 

   medround = median(round)
   print, 'DAO FIND roundness parameter : ', medround
   IF abs(medround) GT 0.5 THEN BEGIN 
      print, 'WARNING: Extremely distorted PSF'
      pt_errlog, 'WARNING: Extremely distorted PSF'
   ENDIF 

; check if too many stars
   IF keyword_set(maxnstar) THEN BEGIN 
      nstar = n_elements(imx)
      IF (nstar GT maxnstar) THEN BEGIN 
         ncull = (nstar/4) > maxnstar    ; never use fewer than 1/4 of stars
         print, 'PTASTROM: found', nstar,' stars - culling to', ncull
;         ind = long(findgen(ncull)*float(nstar)/float(ncull))
         sind = sort(flux)
         ind = sind[nstar-ncull-1:nstar-1] ; keep brightest ncull
         imx = imx[ind]
         imy = imy[ind]
         flux = flux[ind]
         sharp = 0  ; throw away
         round = 0  ; throw away
      END
   END 


   cat = starcat ; so we don't overwrite starcat
   gsssadxy, gsa, cat.ra, cat.dec, catx, caty
   cat.x = catx
   cat.y = caty

; display full image

;   display, image, min=immean-100, max=immean+201
;   oplot, imx, imy, ps=6, syms=2
;   oplot, catx, caty, ps=7, syms=2

; cut out center section
   wsmall = where((imx GT cut0[0]) AND (imx LT cut1[0]) AND $
                  (imy GT cut0[1]) AND (imy LT cut1[1]), count_small)

   IF count_small LT 5 THEN BEGIN 
      errstr = string('---> DAO Find:', count_small, $
                      ' stars found -- giving up!', format='(A,I4,A)')
      pt_errlog, infostr, errstr
      return
   ENDIF 

   ;----------
   ; Search for the angle between the catalogue stars and image stars

   ang = angle_from_pairs(cat.x, cat.y, imx[wsmall], imy[wsmall], $
    dmax=dmax, binsz=binsz, bestsig=bestsig) ; 4 arcmin
   if bestsig gt 12 then begin 
      print, 'Best Angle: ', ang, '  sigma: ', bestsig
      cd = fltarr(2, 2)
      cd[0, 0] = gsa.amdx[0]
      cd[0, 1] = gsa.amdx[1]
      cd[1, 1] = gsa.amdy[0]
      cd[1, 0] = gsa.amdy[1]
      angrad = ang * !pi / 180.
      mm = [[cos(angrad), sin(angrad)], [-sin(angrad), cos(angrad)]]
      cd = cd # mm
      print, cd
      print
      print, mm
      gsa.amdx[0] = cd[0, 0]
      gsa.amdx[1] = cd[0, 1]  
      gsa.amdy[0] = cd[1, 1]
      gsa.amdy[1] = cd[1, 0]
      gsssadxy, gsa, cat.ra, cat.dec, catx, caty
      cat.x = catx
      cat.y = caty
   endif else begin
      print, 'WARNING:  I think I am lost, but I will try anyway...'
   endelse

   ;----------
   ; Search for the offset between the catalogue stars and image stars

   xyshift = offset_from_pairs(cat.x, cat.y, imx[wsmall], imy[wsmall], $
       dmax=dmax, binsz=binsz, errflag=errflag, bestsig=bestsig)


   IF errflag THEN BEGIN 
      print, 'XY shift FAILED on first attempt...'
      print, '  Trying again with larger search radius'
      print
      xyshift = offset_from_pairs(cat.x, cat.y, imx[wsmall], imy[wsmall], $
                   dmax=4*dmax, binsz=binsz*2, errflag=errflag, bestsig=bestsig)
      IF errflag THEN BEGIN 
         errstr = string('XY shift FAILED:', bestsig, ' sigma -- giving up!',$
                        format='(A,F6.2,A)')
         pt_errlog, infostr, errstr
         return
      ENDIF 
   ENDIF 

   print, '---> XYSHIFT: ', xyshift

;=============================== Begin Fink

   xcen = gsa.ppo3/gsa.xsz-.5d   ; 1023.5
   ycen = gsa.ppo6/gsa.ysz-.5d   ; 1023.5
;  NOTE FITS crpix is 1-indexed but argument of xy2ad is 0-indexed

   refpix = [xcen, ycen] - xyshift
   gsssxyad, gsa, refpix[0], refpix[1], racen, deccen

;  update astrometry structure with new CRVALs
   gsa.crval = [racen, deccen]
;  update catalogue .x and .y fields
   gsssadxy, gsa, cat.ra, cat.dec, catx, caty
   cat.x = catx & cat.y=caty

   im_template = {im_specs, $
                  x:   0.0, $
                  y:   0.0, $
                  ra:  0.0, $
                  dec: 0.0, $
                  adu: 0.0}


   im = replicate(im_template, n_elements(imx))
   im.x = imx
   im.y = imy

   IF keyword_set(display) THEN BEGIN 
       display, im_small, min=immean-5*imsig, max=immean+10*imsig
       oplot, cat.x-cut0, cat.y-cut0, ps=6
       oplot, imx-cut0, imy-cut0, ps=4
   ENDIF 

; find matches between USNO catalogue and image stars

   gsa1 = gsa

;----------
; Tweak astrometry structure with cat (ra,dec) and im (x,y) comparison

; first pass

   pt_tweak_astr, cat, im, maxsep, gsa1, errflag=errflag, infostr=infostr
   gsa2 = gsa1


   maxrad = [1000, 1200, 2000]

   pt_tweak_astr, cat, im, maxsep/2, gsa2, radial=radial, $
     maxrad=maxrad[0], infostr=infostr, errflag=errflag
   pt_tweak_astr, cat, im, maxsep/4., gsa2, radial=radial, $
     maxrad=maxrad[1], infostr=infostr, errflag=errflag
   pt_tweak_astr, cat, im, maxsep/8., gsa2, radial=radial, $
     infostr=infostr, errflag=errflag
   pt_tweak_astr, cat, im, maxsep/8., gsa2, radial=radial, $
     infostr=infostr, errflag=errflag, /qa, ps=ps, pt_logstr=pt_logstr, $
     log_info=log_info, catmatch=catmatch

; if anything set the errflag, DO NOT update 
   IF keyword_set(errflag) THEN return   

   success = 1B
   gsa = gsa2

   xvec1 = gsa.amdx

   angerr = acos((transpose(xvec0[0:1]) # xvec1[0:1])/sqrt(total(xvec0[0:1]^2)*total(xvec1[0:1]^2)))/!dtor
   print, 'Initial guess rotated by: ', angerr
   if angerr[0] gt 1. then print, 'WARNING - initial guess very poor!'

   return
end
;------------------------------------------------------------------------------
