;#############################################################################
;
; Copyright (C) 2001-2009, Michele Cappellari
; E-mail: cappellari_at_astro.ox.ac.uk
; 
; Updated versions of the software are available from my web page
; http://purl.org/cappellari/idl
;
; If you have found this software useful for your research,
; I would appreciate an acknowledgment to the use of the
; "Penalized Pixel-Fitting method by Cappellari & Emsellem (2004)".
;
; This software is provided as is without any warranty whatsoever.
; Permission to use, for non-commercial purposes is granted.
; Permission to modify for personal or internal use is granted,
; provided this copyright and disclaimer are included unchanged
; at the beginning of the file. All other rights are reserved.
;
;#############################################################################
;+
; NAME:
;   PPXF
;
; PURPOSE:
;   Extract galaxy stellar kinematics (V, sigma, h3, h4, h5, h6)
;   by fitting a template to an observed spectrum in pixel space,
;   using the Penalized Pixel-Fitting (pPXF) method described in
;   Cappellari M., & Emsellem E., 2004, PASP, 116, 138.
;
;   The following key optional features are also available:
;   1) An optimal template, positive linear combination of different
;      input templates, can be fitted together with the kinematics.
;   2) One can enforce smoothness on the template weights during the fit. 
;      This is useful to attach a physical meaning to the weights e.g. 
;      in term of the star formation history of a galaxy.
;   3) Additive and/or multiplicative polynomials can be included to adjust
;      the continuum shape of the template to the observed spectrum.
;   4) Iterative sigma clipping can be used to clean the spectrum.
;   5) It is possible to fit a mirror-symmetric LOSVD to two spectra at
;      the same time. This is useful for spectra taken at point-symmetric
;      spatial positions with respect to the center of an equilibrium
;      stellar system.
;   6) One can include sky spectra in the fit, to deal with cases
;      where the sky dominates the observed spectrum and an accurate
;      sky subtraction is critical.
;   7) One can derive an estimate of the reddening in the spectrum.
;
; CALLING SEQUENCE:
;   PPXF, star, galaxy, noise, velScale, start, sol, $
;       BESTFIT=bestFit, BIAS=bias, /CLEAN, DEGREE=degree, ERROR=error, $
;       GOODPIXELS=goodPixels, LAMBDA=lambda, MDEGREE=mdegree, MOMENTS=moments, $
;       /OVERSAMPLE, /PLOT, POLYWEIGHTS=polyWeights, /QUIET, REDDENING=reddening, $
;       REGUL=regul, SKY=sky, VSYST=vsyst, WEIGHTS=weights
;
; INPUT PARAMETERS:
;   STAR: vector containing the spectrum of a single template star or array of
;       dimensions STAR[nPixels,nTemplates] containing different templates to
;       be optimized during the fit of the kinematics. nPixels has to be >= the
;       number of galaxy pixels.
;     - To apply linear regularization to the WEIGHTS via the keyword REGUL,
;       STAR should be an array of two STAR[nPixels,nAge] or three dimensions 
;       STAR[nPixels,nAge,nMetal]. This can be useful e.g. to try to attach a 
;       physical meaning to the output WEIGHTS, in term of the galaxy star 
;       formation history. In that case the templates can represent e.g. 
;       single stellar population SSP models and will be arranged in sequence 
;       of increasing age and metallicity along the second or third dimension 
;       of the array. 
;     - STAR and GALAXY do not need to span the same wavelength range. However
;       an error will be returned by PPXF, if the velocity shift in pixels,
;       required to match the galaxy with the templates, becomes larger than
;       nPixels. In that case one has to truncate either the galaxy or the
;       templates to make the two rest-frame spectral ranges more similar.
;   GALAXY: vector containing the spectrum of the galaxy to be measured. The
;       star and the galaxy spectra have to be logarithmically rebinned but the
;       continuum does *not* have to be subtracted. The rebinning may be
;       performed with the LOG_REBIN routine that is distributed with PPXF.
;     - For high redshift galaxies, one should bring the spectra close to the
;       restframe wavelength, before doing the PPXF fit, to prevent too large 
;       velocity shifts of the templates. This can be done by dividing the 
;       observed wavelenghts by (1+z), where z is a rough estimate of the 
;       galaxy redshift, before the logarithmic rebinning.
;     - GALAXY can also be an array of dimensions GALAXY[nGalPixels,2] containing
;       two spectra to be fitted, at the same time, with a reflection-symmetric
;       LOSVD. This is useful for spectra taken at point-symmetric spatial
;       positions with respect to the center of an equilibrium stellar system.
;       For a discussion of the usefulness of this two-sided fitting
;       see e.g. Section 3.6 of Rix & White (1992, MNRAS, 254, 389).
;     - IMPORTANT: 1) For the two-sided fitting the VSYST keyword has to be used.
;       2) Make sure the spectra are rescaled to be not too many order of
;       magnitude different from unity, to avoid over or underflow problems
;       in the calculation. E.g. units of erg/(s cm^2 A) may cause problems!
;   NOISE: vector containing the 1*sigma error (per pixel) in the galaxy spectrum.
;       If GALAXY is a Nx2 array, NOISE has to be an array with the same dimensions.
;     - IMPORTANT: the penalty term of the pPXF method is based on the *relative*
;       change of the fit residuals. For this reason the penalty will work as
;       expected even if no reliable estimate of the NOISE is available
;       (see Cappellari & Emsellem [2004] for details).
;       If no reliable noise is available this keyword can just be set to:
;           NOISE = galaxy*0+1 ; Same weight for all pixels
;   VELSCALE: velocity scale of the spectra in km/s per pixel. It has to be the
;       same for both the galaxy and the template spectra.
;   START: two elements vector [velStart,sigmaStart] with the initial estimate
;       for the velocity and the velocity dispersion in km/s.
;     - Unless a good initial guess is available, it is recommended to set the starting
;       sigma >= 3*velScale in km/s (i.e. 3 pixels). In fact when the LOSVD is severely
;       undersampled, and far from the true solution, the chi^2 of the fit becomes weakly
;       sensitive to small variations in sigma (see pPXF paper). In some instances the
;       near-constancy of chi^2 may cause premature convergence of the optimization.
;     - In the case of two-sided fitting a good starting value for the
;       velocity is velStart=0.0 (in this case VSYST will generally be nonzero).
;       Alternatively on should keep in mind that velStart refers to the first
;       input galaxy spectrum, while the second will have velocity -velStart.
;
; KEYWORDS:
;   BESTFIT: a named variable to receive a vector with the best fitting
;       template: this is a linear combination of the templates, convolved with
;       the best fitting LOSVD, with added polynomial continuum terms.
;   BIAS: This parameter biases the (h3,h4,...) measurements towards zero
;       (Gaussian LOSVD) unless their inclusion significantly decreses the
;       error in the fit. Set this to BIAS=0.0 not to bias the fit: the
;       solution (including [V,sigma]) will be noisier in that case. The
;       default BIAS should provide acceptable results in most cases, but it
;       would be safe to test it with Monte Carlo simulations. This keyword
;       precisely corresponds to the parameter \lambda in the Cappellari &
;       Emsellem (2004) paper. Note that the penalty depends on the *relative*
;       change of the fit residuals, so it is insensitive to proper scaling
;       of the NOISE vector. A nonzero BIAS can be safely used even without a
;       reliable error spectrum, or with equal weighting for all pixels.
;   /CLEAN: set this keyword to use the iterative sigma clipping method
;       described in Section 2.1 of Cappellari et al. (2002, ApJ, 578, 787).
;       This is useful to remove from the fit unmasked bad pixels, residual
;       gas emissions or cosmic rays.
;     - IMPORTANT: This is recommended *only* if a reliable estimate of the
;       error spectrum is available. See also note below for SOL.
;   DEGREE: degree of the *additive* Legendre polynomial used to correct
;       the template continuum shape during the fit (default: 4).
;       Set DEGREE = -1 not to include any additive polynomial.
;   ERROR: a named variable that will contain a vector of *formal* errors
;       (1*sigma) for the parameters in the output vector SOL. This option
;       can be used when speed is essential, to obtain an order of magnitude
;       estimate of the uncertainties, but we *strongly* recommend to run Monte
;       Carlo simulations to obtain more reliable errors. In fact these errors
;       can be severely underestimated in the region where the penalty effect
;       is most important (sigma < 2*velScale).
;     - These errors are meaningless unless Chi^2/DOF~1 (see parameter SOL below).
;       However if one *assume* that the fit is good, a corrected estimate of the 
;       errors is: errorCorr = error*sqrt(chi^2/DOF) = error*sqrt(sol[6]).
;     - IMPORTANT: when running Monte Carlo simulations to determine the error,
;       the penalty (BIAS) should be set to zero, or better to a very small value.
;       See Section 3.4 of Cappellari & Emsellem (2004) for an explanation.
;   GOODPIXELS: integer vector containing the indices of the good pixels in the
;       GALAXY spectrum (in increasing order). Only these pixels are included
;       in the fit. If the /CLEAN keyword is set, in output this vector will
;       contain the pixels that were actually used in the fit.
;     - IMPORTANT: in all likely situations this keyword *has* to be specified.
;   LAMBDA: When the keyword REDDENING is used, the user has to pass in this
;       keyword a vector with the same dimensions of GALAXY, giving the restframe
;       wavelength in Angstrom of every pixel in the input galaxy spectrum.
;       If one uses my LOG_REBIN routine to rebin the spectrum before the PPXF fit:
;           LOG_REBIN, lamRange, galaxy, galaxyNew, logLam
;       the wavelength can be obtained as EXP(logLam)
;   MDEGREE: degree of the *multiplicative* Legendre polynomial (with mean of 1)
;       used to correct the continuum shape during the fit (default: 0). The
;       zero degree multiplicative polynomial is always included in the fit as
;       it corresponds to the weights assigned to the templates.
;       Note that the computation time is longer with multiplicative
;       polynomials than with the same number of additive polynomials.
;     - IMPORTANT: Multiplicative polynomials cannot be used when 
;       the REDDENING keyword is set.  
;   MOMENTS: Order of the Gauss-Hermite moments to fit. Set this keyword to 4
;       to fit [h3,h4] and to 6 to fit [h3,h4,h5,h6]. Note that in all cases
;       the G-H moments are fitted (nonlinearly) *together* with [V,sigma]. If
;       this keyword is not set (or set to 2) then only [V,sigma] are fitted
;       and the other parameters are returned as zero. If MOMENTS=0 then only
;       the templates and the continuum polynomials are fitted.
;   /OVERSAMPLE: Set this keyword to oversample the template by a factor 3x
;       before convolving it with a well sampled LOSVD. This is mainly useful
;       for testing: one has to be careful if the results depend significantly
;       from this keyword.
;   /PLOT: set this keyword to plot the best fitting solution and the residuals
;       at the end of the fit.
;   POLYWEIGHTS: vector with the weights of the additive Legendre polynomials.
;       The best fitting additive polynomial can be explicitly evaluated as
;           x = range(-1d,1d,n_elements(galaxy))
;           apoly = 0d ; Additive polynomial
;           for j=0,DEGREE do apoly += legendre(x,j)*polyWeights[j]
;     - When doing a two-sided fitting (see help for GALAXY parameter), the additive
;       polynomials are allowed to be different for the left and right spectrum.
;       In that case the output weights of the additive polynomials alternate between
;       the first (left) spectrum and the second (right) spectrum.
;   /QUIET: set this keyword to suppress verbose output of the best fitting
;       parameters at the end of the fit.
;   REDDENING: Set this keyword to an initail estimate of the reddening E(B-V)>=0 
;       to fit a positive reddening together with the kinematics and the templates.
;       After the fit the input estimate is replaced with the best fitting E(B-V) value.
;       The fit assumes the exctinction curve of Calzetti et al. (2000, ApJ, 533, 682)
;       but any other prescriptions could be trivially implemented by modifying the
;       function REDDENING_CURVE below.
;     - IMPORTANT: The MDEGREE keyword cannot be used when REDDENING is set.  
;   REGUL: If this keyword is nonzero, the program applies second-degree 
;       linear regularization to the WEIGHTS during the PPXF fit.
;       Regularization is done in one or two dimensions depending on whether
;       the array of templates STAR has two or three dimensions respectively. 
;       Large REGUL values correspond to smoother WEIGHTS output. The WEIGHTS 
;       tend to a linear trend for REGUL >> 1. When this keyword is nonzero the 
;       solution will be a trade-off between smoothness of WEIGHTS and 
;       goodness of fit. Typical values for this parameter are 0 < REGUL < 1.
;     - For a detailed explanation see Section 18.5 of Press et al. (1992, 
;       Numerical Recipes 2nd ed.) available here http://www.nrbook.com/a/bookfpdf.php
;       The adopted implementation corresponds to their equation (18.5.10).
;   SKY: vector containing the spectrum of the sky to be included in the fit, or array
;       of dimensions SKY[nPixels,nSky] containing different sky spectra to add to
;       the model of the observed GALAXY spectrum. The SKY has to be log-rebinned as
;       the GALAXY spectrum and needs to have the same number of pixels.
;     - The sky is generally subtracted from the data before the PPXF fit. However,
;       for oservations very heavily dominated by the sky spectrum, where a very
;       accurate sky subtraction is critical, it may be useful *not* to subtract
;       the sky from the spectrum, but to include it in the fit using this keyword.
;   VSYST: galaxy systemic velocity (zero by default). The input initial guess
;       and the output velocities are measured with respect to this velocity.
;       The value assigned to this keyword is *crucial* for the two-sided fitting.
;       In this case VSYST can be determined from a previous normal one-sided
;       fit to the galaxy velocity profile. After that initial fit, VSYST
;       can be defined as the measured velocity at the galaxy center.
;       More accurately VSYST is the value which has to be subtracted to obtain
;       a nearly anti-symmetric velocity profile at the two opposite sides of
;       the galaxy nucleus.
;     - IMPORTANT: this value is generally *different* from the systemic
;       velocity one can get from the literature. Do not try to use that!
;   WEIGHTS: a named variable to receive the value of the weights by which each
;       template was multiplied to best fit the galaxy spectrum. The optimal
;       template can be computed with an array-vector multiplication:
;           TEMPLATE = STAR # WEIGHTS (in IDL syntax)
;     - When the SKY keyword is used WEIGHTS[0:nTemplates-1] contains the weights
;       for the templates, while WEIGHTS[nTemplates:*] gives the ones for the sky.
;       In that case the best fitting galaxy template and sky are given by:
;           TEMPLATE = STAR # WEIGHTS[0:nTemplates-1]
;           BESTSKY = SKY # WEIGHTS[nTemplates:*]
;     - When doing a two-sided fitting (see help for GALAXY parameter) *together*
;       with the SKY keyword, the sky weights are allowed to be different for the
;       left and right spectrum. In that case the output sky weights alternate
;       between the first (left) spectrum and the second (right) spectrum.
;
; OUTPUT PARAMETER:
;   SOL: seven elements vector containing in output the values of
;       [Vel,Sigma,h3,h4,h5,h6,Chi^2/DOF] of the best fitting solution.
;       Vel is the velocity, Sigma is the velocity dispersion, h3-h6 are the
;       Gauss-Hermite coefficients. The model parameter are fitted simultaneously.
;     - I hardcoded the following safety limits on the fitting parameters:
;         a) Vel is constrained to be +/-2000 km/s from the first input guess
;         b) velScale/10 < Sigma < 1000 km/s
;         c) -0.3 < [h3,h4,...] < 0.3 (limits are extreme value for real galaxies)
;     - In the case of two-sided LOSVD fitting the output values refer
;       to the first input galaxy spectrum, while the second spectrum will
;       have by construction kinematics parameters [-Vel,Sigma,-h3,h4,-h5,h6].
;       If VSYST is nonzero (as required for two-sided fitting), then the
;       output velocity is measured with respect to VSIST.
;     - IMPORTANT: if Chi^2/DOF is not ~1 it means that the errors are not
;       properly estimated, or that the template is bad and it is *not* safe
;       to set the /CLEAN keyword.
;     - When MDEGREE > 1 then SOL contains in output the 7+MDEGREE elements
;       [Vel,Sigma,h3,h4,h5,h6,Chi^2/DOF,cx1,cx2,...,cxn], where cx1,cx2,...,cxn
;       are the coefficients of the multiplicative Legendre polynomials
;       of order 1,2,...,n. The polynomial can be explicitly evaluated as:
;           x = range(-1d,1d,n_elements(galaxy))
;           mpoly = 1d ; Multiplicative polynomial
;           for j=1,MDEGREE do mpoly += legendre(x,j)*sol[6+j]
;
;--------------------------------
; IMPORTANT: Proper usage of pPXF
;--------------------------------
;
; The PPXF routine can give sensible quick results with the default BIAS
; parameter, however, like in any penalized/filtered/regularized method, the
; optimal amount of penalization generally depends on the problem under study.
;
; The general rule here is that the penalty should leave the line-of-sight
; velocity-distribution (LOSVD) virtually unaffected, when it is well
; sampled and the signal-to-noise ratio (S/N) is sufficiently high.
;
; EXAMPLE: If you expect an LOSVD with up to a high h4 ~ 0.2 and your
; adopted penalty (BIAS) biases the solution towards a much lower h4 ~ 0.1,
; even when the measured sigma > 3*velScale and the S/N is high, then you
; are *misusing* the pPXF method!
;
; THE RECIPE: The following is a simple practical recipe for a sensible
; determination of the penalty in pPXF:
;
; 1. Choose a minimum (S/N)_min level for your kinematics extraction and
;    spatially bin your data so that there are no spectra below (S/N)_min;
;
; 2. Perform a fit of your kinematics *without* penalty (PPXF keyword BIAS=0).
;    The solution will be noisy and may be affected by spurious solutions,
;    however this step will allow you to check the expected mean ranges in
;    the Gauss-Hermite parameters [h3,h4] for the galaxy under study;
;
; 3. Perform a Monte Carlo simulation of your spectra, following e.g. the
;    included ppxf_simulation.pro routine. Adopt as S/N in the simulation the
;    chosen value (S/N)_min and as input [h3,h4] the maximum representative
;    values measured in the non-penalized pPXF fit of the previous step;
;
; 4. Choose as penalty (BIAS) the *largest* value such that, for sigma > 3*velScale,
;    the mean difference between the output [h3,h4] and the input [h3,h4]
;    is well within the rms scatter of the simulated values
;    (see e.g. Fig.2 of Emsellem et al. 2004, MNRAS, 352, 721).
;
;--------------------------------
;
; REQUIRED ROUTINES:
;       RANGE: by M. Cappellari (included in the PPXF distribution)
;       BVLS: by M. Cappellari from http://www-astro.physics.ox.ac.uk/~mxc/idl/
;       MPFIT: by C.B. Markwardt http://purl.com/net/mpfit
;       ROBUST_SIGMA: by H. Freudenreich from http://idlastro.gfsc.nasa.gov/
;
; MODIFICATION HISTORY:
;   V1.0 -- Created by Michele Cappellari, Leiden, 10 October 2001.
;   V3.47 -- First released version. MC, Leiden, 8 December 2003
;   V3.5 -- Included /OVERSAMPLE option. MC, Leiden, 11 December 2003
;   V3.6 -- Added MDEGREE option for multiplicative polynomials.
;           Linear implementation: fast, works well in most cases, but
;           can fail in certain cases. MC, Leiden, 19 March 2004
;   V3.7 -- Revised implementation of MDEGREE option. Nonlinear implementation:
;           straightforward, robust, but slow. MC, Leiden, 23 March 2004
;   V3.71 -- Updated documentation. MC, Leiden, 31 March 2004
;   V3.72 -- Corrected program stop after fit when MOMENTS=2.
;           Bug was introduced in V3.7. MC, Leiden, 28 April 2004
;   V3.73 -- Corrected bug: keyword ERROR was returned in pixels
;           instead of km/s. Decreased lower limit on fitted dispersion.
;           Thanks to Igor V. Chilingarian. MC, Leiden, 7 August 2004
;   V4.0 -- Introduced optional two-sided fitting assuming a reflection-symmetric
;           LOSVD for two input spectra. MC, Vicenza, 16 August 2004
;   V4.1 -- Corrected implementation of two-sided fitting of the LOSVD.
;           Thanks to Stefan van Dongen for reporting problems.
;           MC, Leiden, 3 September 2004
;   V4.11 -- Increased maximum number of iterations ITMAX in BVLS.
;           Thanks to Jesus Falcon-Barroso for reporting problems.
;           Introduced error message when velocity shift is too big.
;           Corrected output when MOMENTS=0. MC, Leiden, 21 September 2004
;   V4.12 -- Handle special case where a single template without additive
;           polynomials is fitted to the galaxy. MC, Leiden, 11 November 2004
;   V4.13 -- Updated documentation. MC, Vicenza, 30 December 2004
;   V4.14 -- Make sure input NOISE is a positive vector. MC, Leiden, 12 January 2005
;   V4.15 -- Verify that GOODPIXELS is monotonic and does not contain duplicated values.
;           After feedback from Richard McDermid. MC, Leiden, 10 February 2005
;   V4.16 -- Print number of nonzero templates. Do not print outliers in /QUIET mode.
;           MC, Leiden, 20 January 2006
;   V4.17 -- Updated documentation with important note on penalty determination.
;           MC, Oxford, 6 October 2007
;   V4.2 -- Introduced optional fitting of SKY spectrum. Many thanks to
;           Anne-Marie Weijmans for testing. MC, Oxford, 15 March 2008
;   V4.21 -- Use LA_LEAST_SQUARES (IDL 5.6) instead of SVDC when fitting
;           a single template. Please let me know if you need to use PPXF
;           with an older IDL version. MC, Oxford, 17 May 2008
;   V4.22 -- Added keyword POLYWEIGHTS. MC, Windhoek, 3 July 2008
;   V4.23 -- Corrected error message for too big velocity shift. 
;           MC, Oxford, 27 November 2008
;   V4.3 -- Introduced REGUL keyword to perform linear regularization of WEIGHTS 
;           in one or two dimensions. MC, Oxford, 4 Mach 2009
;   V4.4 -- Introduced Calzetti et al. (2000) REDDENING_CURVE function to
;           estimate the reddening from the fit. MC, Oxford, 18 September 2009
;
;
; J. Moustakas, 2009-Nov-13, UCSD - generalized reddening curve K_LAMBDA()
; jm09nov18ucsd - added support for simultaneously fitting photometry 
; jm09nov18ucsd - added support to fit for the reddening when MOMENTS=0
;-
;----------------------------------------------------------------------------
;;FUNCTION reddening_curve, lambda, ebv
;;;
;;; Reddening curve of Calzetti et al. (2000, ApJ, 533, 682; here C+00).
;;; This is reliable between 0.12 and 2.2 micrometres.
;;; - LAMBDA is the restframe wavelength in Angstrom of each pixel in the 
;;; input galaxy spectrum (1 Angstrom = 1d-4 micrometres)
;;; - EBV is the assumed E(B-V) colour excess to redden the spectrum.
;;; In output the vector FRAC gives the fraction by which the flux at each
;;; wavelength has to be multiplied, to model the dust reddening effect.
;;
;;k1 = lambda*0 
;;lam = 1e4/lambda ; Convert Angstrom to micrometres and take 1/lambda          
;;rv = 4.05d ; C+00 equation (5) 
;;
;;w1 = where(lambda ge 6300d, m1, COMPLEMENT=w2, NCOMPLEMENT=m2) 
;;; C+00 equation (3) but extrapolate for lam>2.2
;;if m1 gt 0 then k1[w1] = rv + 2.659d*(1.040d*lam[w1] - 1.857d)
;;; C+00 equation (4) but extrapolate for lam<0.12
;;if m2 gt 0 then k1[w2] = rv + $
;;    2.659d*(1.509d*lam[w2] - 0.198d*lam[w2]^2 + 0.011d*lam[w2]^3 - 2.156d)
;;fact = 10d^(-0.4*ebv*(k1>0))  ; Calzetti+00 equation (2) with opposite sign
;;
;;return, fact ; The model spectrum has to be multiplied by this vector
;;END
;----------------------------------------------------------------------------
PRO im_ppxf, star, galaxy, noise, velScale, start, sol, $
    BESTFIT=bestFit, BIAS=bias, CLEAN=clean, DEGREE=degree, ERROR=error, $
    GOODPIXELS=goodPixels, MDEGREE=mdegree, MOMENTS=moments, OVERSAMPLE=oversample, $
    POLYWEIGHTS=polyweights, PLOT=plot, QUIET=quiet, SKY=sky, VSYST=vsyst, WEIGHTS=weights, $
    REGUL=regul, LAMBDA=lambda, REDDENING=reddening, losvd=losvd, photo=photo, $
  outphoto=outphoto, vmaxshift=vmaxshift, sigmamax=sigmamax
compile_opt idl2
;on_error, 2

; Do extensive checking of possible input errors
;
s1 = size(star)
if s1[0] eq 3 then begin
    reg_dim = s1[2:3]
    star = reform(star,s1[1],s1[2]*s1[3])
    s1 = size(star)
endif else reg_dim = s1[2]
if n_elements(regul) eq 0 then regul = 0    
s2 = size(galaxy)
s3 = size(noise)
s4 = size(sky)
if (s1[0] gt 2 || s2[0] gt 2 || s3[0] gt 2) then message, 'Wrong input dimensions'
if ~array_equal(s2,s3) then message, 'GALAXY and NOISE must have the same size/type'
if (s1[1] lt s2[1]) then message, 'STAR length cannot be smaller than GALAXY'
if n_elements(reddening) gt 0 then begin
    if ~array_equal(size(lambda),s2) then $
        message, 'LAMBDA and GALAXY must have the same size/type'
    if n_elements(mdegree) gt 0 then $
        message, 'MDEGREE cannot be used with REDDENING keyword'
endif else lambda = 0
if ~array_equal(noise gt 0, 1) then message, 'NOISE must be a positive vector'
if s4[0] gt 0 && s4[1] ne s2[1] then message, 'SKY must have the same size as GALAXY'
if n_elements(degree) le 0 then degree = 4 else degree = degree > (-1)
if n_elements(mdegree) le 0 then mdegree = 0 else mdegree = mdegree > 0
if keyword_set(oversample) then factor = 3 else factor = 1
nGood = n_elements(goodPixels)
if nGood le 0 then begin
    nGood = s2[1]
    goodPixels = indgen(nGood)
endif else begin
    if ~array_equal((goodPixels-shift(goodPixels,1))[1:*] gt 0, 1) then $
        message, 'GOODPIXELS is not monotonic or contains duplicated values'
    if goodPixels[0] lt 0 || goodPixels[nGood-1] gt s2[1]-1 then $
        message, 'GOODPIXELS are outside the data range'
endelse
if n_elements(bias) le 0 then bias = 0.7d*sqrt(500d/n_elements(goodPixels)) ; pPXF paper pg.144 left
if n_elements(moments) eq 0 then moments = 2 else $
    if total(moments eq [0,2,4,6]) eq 0 then message, 'MOMENTS should be 0, 2, 4 or 6'
if moments ge 2 && n_elements(start) ne 2 then message, 'START must have two elements [V,sigma]'
if s2[0] eq 2 then goodPixels = [goodPixels,s2[1]-1+goodPixels]  ; two-sided fitting of LOSVD
if n_elements(vsyst) eq 0 then begin
    if s2[0] eq 2 then message, 'VSYST must be defined for two-sided fitting'
    vsyst = 0d
endif

; jm09nov16ucsd - deal with the photometry
if (n_elements(photo) gt 0) then splog, 'Using photometry!'

; Parameters to be passed to the fitting function, converted to DOUBLE
;
if n_elements(sky) eq 0 then begin
   functArgs = {BIAS:bias, DEGREE:degree, FACTOR:factor, GALAXY:double(galaxy), $
     GOODPIXELS:goodPixels, MDEGREE:mdegree, NOISE:double(noise), $
     STAR:double(star), VSYST:vsyst/velScale, $
     REGUL:regul, REG_DIM:reg_dim, LAMBDA:lambda}
   if (n_elements(photo) gt 0) then functargs = create_struct(functargs,photo) ; jm09nov16ucsd - 
endif else begin
   functArgs = {BIAS:bias, DEGREE:degree, FACTOR:factor, GALAXY:double(galaxy), $
     GOODPIXELS:goodPixels, MDEGREE:mdegree, NOISE:double(noise), $
     SKY:double(sky), STAR:double(star), VSYST:vsyst/velScale, $
     REGUL:regul, REG_DIM:reg_dim, LAMBDA:lambda}
endelse

; jm09dec08ucsd
if (n_elements(vmaxshift) eq 0) then vmaxshift = 2000D ; maximum velocity shift (km/s)
if (n_elements(sigmamax) eq 0) then sigmamax = 1000D ; maximum velocity dispersion (km/s)

if moments gt 0 then begin

    ; Explicitly specify the step for the numerical derivatives (pixel/100)
    ; in MPFIT routine and force safety limits on the fitting parameters.
    ;
    if mdegree gt 0 then begin
        start1 = dblarr(moments+mdegree*s2[0])  ; Set [h3,h4,...] to zero as initial guess
        parinfo = replicate({step:1d-2,limits:[0d,0d],limited:[1,1]}, moments+mdegree*s2[0])
        parinfo[moments:*].limits = [-1d,1d] ; force <100% corrections
        parinfo[moments:*].step = 1d-3
    endif else if n_elements(reddening) gt 0 then begin
        start1 = dblarr(moments+1)  ; Set [h3,h4,...] to zero as initial guess
        start1[moments] = reddening
        parinfo = replicate({step:1d-2,limits:[0d,0d],limited:[1,1]}, moments+1)
        parinfo[moments].limits = [0d,10d] ; force positive E(B-V) < 10 mag
        parinfo[moments].step = 1d-3
    endif else begin
        start1 = dblarr(moments)  ; Set [h3,h4,...] to zero as initial guess
        parinfo = replicate({step:1d-2,limits:[0d,0d],limited:[1,1]}, moments)
    endelse
    start1[0] = start/velScale  ; Convert velocity scale to pixels
    parinfo[0].limits = start1[0] + [-vmaxshift,vmaxshift]/velScale ; +/-VMAXSHIFT km/s from first guess
    parinfo[1].limits = [0.1d,sigmamax/velScale] ; hard-coded velScale/10<sigma<SIGMAMAX km/s
;   parinfo[0].limits = start1[0] + [-2d3,2d3]/velScale ; +/-2000 km/s from first guess
;   parinfo[1].limits = [0.1d,1d3/velScale] ; hard-coded velScale/10<sigma<1000 km/s
    if moments gt 2 then begin
        parinfo[2:moments-1].limits = [-0.3d,0.3d] ; -0.3<[h3,h4,...]<0.3
        parinfo[2:moments-1].step = 1d-3
    endif 
    if s1[1] le 2*(abs(vsyst/velScale)+abs(start1[0])+5d*start1[1]) then $
        message, 'Velocity shift too big: Adjust wavelength ranges of spectrum and templates'

    ; Here the actual calculation starts.
    ; If required, once the minimum is found, clean the pixels deviating
    ; more than 3*sigma from the best fit and repeat the minimization
    ; until the set of cleaned pixels does not change any more.
    ;
    good = goodPixels
    for j=0,4 do begin ; Do at most five cleaning iterations
        res = mpfit('fitfunc_pxf_optimal_template', start1, ERRMSG=errmsg, $
            FTOL=1d-4, FUNCTARGS=functArgs, NFEV=ncalls, PARINFO=parinfo, $
            PERROR=error, /QUIET)
        if errmsg ne '' then message, errmsg
        if ~keyword_set(clean) then break
        goodOld = goodPixels
        goodPixels = good
; jm09nov13ucsd - photometry
        if (n_elements(photo) gt 0) then begin
           tmp = fitfunc_pxf_optimal_template(res, BIAS=bias, /CLEAN, $
             DEGREE=degree, FACTOR=factor, GALAXY=galaxy, $
             GOODPIXELS=goodPixels, MDEGREE=mdegree, NOISE=noise, $
             QUIET=quiet, SKY=sky, STAR=star, VSYST=vsyst/velScale, $
             REGUL=regul, REG_DIM=reg_dim, LAMBDA=lambda, losvd=losvd, $
             weff=photo.weff, maggies=photo.maggies, $
             errmaggies=photo.errmaggies, modelmaggies=photo.modelmaggies)
        endif else begin
           tmp = fitfunc_pxf_optimal_template(res, BIAS=bias, /CLEAN, $
             DEGREE=degree, FACTOR=factor, GALAXY=galaxy, $
             GOODPIXELS=goodPixels, MDEGREE=mdegree, NOISE=noise, $
             QUIET=quiet, SKY=sky, STAR=star, VSYST=vsyst/velScale, $
             REGUL=regul, REG_DIM=reg_dim, LAMBDA=lambda, losvd=losvd)
        endelse
        if array_equal(goodOld,goodPixels) then break
    endfor

endif else begin

; --------------------------------------------------
; jm09nov18ucsd - added support for MOMENTS = 0, but still solve for
; the reddening!   
   if n_elements(reddening) gt 0 then begin
      moments = 2
      start1 = [start/velscale,double(reddening)]
      parinfo = replicate({fixed: 1, step:1d-2,limits:[0d,0d],limited:[0,0]},3)

      parinfo[2].fixed = 0
      parinfo[2].limited = [1,1]
      parinfo[2].limits = [0d,10d] ; force positive E(B-V) < 10 mag
      parinfo[2].step = 1d-3

      good = goodPixels
      for j=0,4 do begin        ; Do at most five cleaning iterations
         res = mpfit('fitfunc_pxf_optimal_template', start1, ERRMSG=errmsg, $
           FTOL=1d-4, FUNCTARGS=functArgs, NFEV=ncalls, PARINFO=parinfo, $
           PERROR=error, /QUIET)
         if errmsg ne '' then message, errmsg
         if ~keyword_set(clean) then break
         goodOld = goodPixels
         goodPixels = good
         tmp = fitfunc_pxf_optimal_template(res, BIAS=bias, /CLEAN, $
           DEGREE=degree, FACTOR=factor, GALAXY=galaxy, $
           GOODPIXELS=goodPixels, MDEGREE=mdegree, NOISE=noise, $
           QUIET=quiet, SKY=sky, STAR=star, VSYST=vsyst/velScale, $
           REGUL=regul, REG_DIM=reg_dim, LAMBDA=lambda, losvd=losvd)
         if array_equal(goodOld,goodPixels) then break
      endfor
   endif else begin
; --------------------------------------------------
; jm09nov18ucsd - original code, which does not deal with reddening
; when moments=0
       if mdegree gt 0 then message, 'MDEGREE>0 not implemented with MOMENTS=0'
       res = start
       res[0:1] = res[0:1]/velScale ; Convert velocity scale to pixels
       error = start*0d             ; No output error on parameters
       ncalls = 1

    endelse

endelse

; Evaluate scatter at the bestfit (with BIAS=0)
; and also get the output BESTFIT and WEIGHTS.
;

if (n_elements(photo) gt 0) then begin ; jm09nov13ucsd
   err1 = fitfunc_pxf_optimal_template(res, BESTFIT=bestFit1, BIAS=0, DEGREE=degree, $
     FACTOR=factor, GALAXY=galaxy, GOODPIXELS=goodPixels, MDEGREE=mdegree, $
     NOISE=noise, SKY=sky, STAR=star, VSYST=vsyst/velScale, WEIGHTS=weights, $
     REGUL=regul, REG_DIM=reg_dim, LAMBDA=lambda, losvd=losvd, $
     weff=photo.weff, maggies=photo.maggies, $
     errmaggies=photo.errmaggies, modelmaggies=photo.modelmaggies)

   ngoodpix = n_elements(goodpixels)
   npix = n_elements(galaxy)
   nphoto = n_elements(photo.weff)
   
   bestfit = bestfit1[0:npix-1]
   err = err1[0:ngoodpix-1]
   chi2 = im_robust_sigma(err, /ZERO)^2 ; Robust computation of Chi^2/DOF.

   photo_bestfit = bestfit1[npix:npix+nphoto-1]
   photo_err = err1[ngoodpix:ngoodpix+nphoto-1]
   photo_chi2 = im_robust_sigma(photo_err, /ZERO)^2 ; Robust computation of Chi^2/DOF.
   outphoto = {bestmaggies: photo_bestfit, photo_chi2: photo_chi2}
endif else begin
   err = fitfunc_pxf_optimal_template(res, BESTFIT=bestFit, BIAS=0, DEGREE=degree, $
     FACTOR=factor, GALAXY=galaxy, GOODPIXELS=goodPixels, MDEGREE=mdegree, $
     NOISE=noise, SKY=sky, STAR=star, VSYST=vsyst/velScale, WEIGHTS=weights, $
     REGUL=regul, REG_DIM=reg_dim, LAMBDA=lambda, losvd=losvd)
   chi2 = im_robust_sigma(err, /ZERO)^2 ; Robust computation of Chi^2/DOF.
endelse

sol = dblarr(7+mdegree*s2[0])
sol[0] = res[0:(n_elements(start)>moments)-1]
sol[0] = sol[0:1]*velScale ; Bring velocity scale back to km/s
error[0:1] = error[0:1]*velScale ; Convert errors to km/s
sol[6] = chi2
if mdegree ge 1 then sol[7] = res[moments:*]
if n_elements(reddening) gt 0 then reddening = res[moments] ; Replace input with best fit

if degree ge 0 then polyweights = weights[0:(degree+1)*s2[0]-1] ; output weights for the additive polynomials
weights = weights[(degree+1)*s2[0]:*] ; output weights for the templates (or sky) only

; Print final results on the screen.
;
if ~keyword_set(quiet) then begin
    print, 'feval', 'V', 'sigma', 'h3', 'h4', 'h5', 'h6', 'Chi2/DOF', FORMAT='(8A10)'
    print, ncalls, sol[0:6], FORMAT='(i10,2f10.1,5f10.3)'
    nw = n_elements(weights)
    if n_elements(reddening) gt 0 then print, 'Reddening E(B-V): ', reddening, FORMAT='(a,g0.3)'
    print, 'Nonzero Templates: ', total(weights gt 0), ' / ', nw, FORMAT='(a,g0,a,g0)'
    if n_elements(weights) le 20 then begin
        print, 'Template weights:'
        print, weights, FORMAT='(10g11.3)'
    endif
endif

; Plot final data-model comparison if required.
; The colors below were chosen for the colormap n. 12 (LOADCT, 12)
; loadct, 12 ; 16 levels. green=40, blue=100, magenta=130, red=200
;
if keyword_set(plot) then begin
    mn = min(bestfit[goodPixels], MAX=mx)
    resid = mn + galaxy - bestfit
    plot, galaxy, XTITLE='pixels', YTITLE='counts', /XSTYLE, /YNOZERO, $
        YRANGE=[min(resid[goodPixels]),mx], YSTYLE=2, XRANGE=[-0.02,1.02]*s2[1]*s2[0]
    oplot, bestfit, COLOR=200, THICK=2
    n = n_elements(goodPixels)
    oplot, goodPixels, replicate(mn,n*s2[0]), PSYM=3, COLOR=40
    oplot, goodPixels, resid[goodPixels], PSYM=4, COLOR=40, SYMSIZE=0.3
    w = where(goodPixels - shift(goodPixels,1) gt 1, m)
    for j=0,m-1 do begin
        x = range(goodPixels[w[j]-1],goodPixels[w[j]])
        oplot, x, resid[x], COLOR=100
    endfor
    if m gt 0 then w = [0,w-1,w,n-1] else w = [0,n-1]  ; Add first and last point
    for j=0,n_elements(w)-1 do $
        oplot, goodPixels[w[[j,j]]], [mn,bestfit[goodPixels[w[j]]]], COLOR=40
endif

END
;----------------------------------------------------------------------------
