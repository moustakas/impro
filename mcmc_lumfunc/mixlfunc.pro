;+
; NAME:
;
;     MIXLFUNC
;
; PURPOSE:
;
;     ESTIMATE A LUMINOSITY FUNCTION ASSUMING A MIXTURE OF GAUSSIAN
;     FUNCTIONS MODEL.
;
; EXPLANATION:
;
;     THIS ROUTINE PERFORMS BAYESIAN INFERENCE ON A LUMINOSITY
;     FUNCTION ASSUMING A MIXTURE OF GAUSSIAN FUNCTIONS MODEL. THE
;     METROPLIS-HASTINGS ALGORITHM IS USED TO SIMULATE RANDOM DRAWS
;     FROM THE POSTERIOR PROBABILITY DISTRIBUTION OF THE PARAMETERS
;     FOR THE MIXTURE OF GAUSSIAN FUNCTIONS MODEL. THESE RANDOM DRAWS
;     CAN THEN BE USED TO COMPUTE ESTIMATES AND UNCERTAINTIES ON THE
;     LUMINOSITY FUNCTION, AS WELL AS ANY QUANTITIES DERIVED FROM
;     IT. SEE THE REFERENCE BELOW FOR MORE DETAILS.
;
;     THE SELECTION FUNCTION MUST BE CONSTRUCTED BY THE USER IN THE
;     SUBROUTINE SELECTION_FUNCTION. THE SELECTION FUNCTION SHOULD
;     CONSTRUCTED AS A STRUCTURE CALLED 'SELFUNC'. THIS STRUCTURE MUST
;     CONTAIN THE FOLLOWING TAGS:
;
;           ZGRID - A VECTOR CONTAINING THE GRID OF VALUES FOR
;                   LOG-REDSHIFTS THAT THE SELECTION FUNCTION IS
;                   EVALUATED AT. ZGRID MUST BE REGULAR.
;           LGRID - A VECTOR CONTAINING THE GRID OF VALUES FOR
;                   LOG-LUMINOSITY THAT THE SELECTION FUNCTION IS
;                   EVAULATED AT. LGRID MUST BE REGULAR AT A GIVEN
;                   REDSHIFT.
;         DETPROB - A VECTOR CONTAINING THE VALUES OF THE DETECTION
;                   PROBABILITY AT ZGRID AND LGRID. DETPROB[i] SHOULD
;                   BE THE DETECTION PROBABILITY FOR A SOURCE AT
;                   REDSHIFT ZGRID[i] AND WITH LUMINOSITY LGRID[i].
;
;     A SIMPLE EXAMPLE IS GIVEN IN THE SUBROUTINE SELECTION_FUNCTION
;     FOR GUIDANCE REGARDING THE FORMAT.
;
; INPUT PARAMETERS
;
;      LUM - A VECTOR CONTAINING THE LOGARITHM OF THE LUMINOSITIES FOR
;            THE SOURCES IN ONE'S SURVEY.
;     ZETA - A VECTOR CONTAINING THE LOGARITHM OF THE REDSHIFTS FOR
;            THE SOURCES IN ONE'S SURVEY.
;   NGAUSS - THE NUMBER OF GAUSSIAN FUNCTIONS TO USE FOR FITTING THE
;            LUMINOSITY FUNCTION. OFTEN, 2-4 IS SUFFICIENT.
;  SKYFRAC - THE FRACTION OF THE SKY COVERED BY THE SURVEY.
;
; OUTPUT PARAMETERS
;
;     POST - A STRUCTURE CONTAINING THE RANDOM DRAWS FROM THE
;            PROBABILITY DISTRIBUTION OF THE MIXTURE OF GAUSSIAN
;            FUNCTIONS PARAMETERS, GIVEN THE OBSERVED DATA AND
;            SELECTION FUNCTION. THE VALUES OF POST ARE SIMULATED
;            USING THE METROPOLIS-HASTINGS ALGORITHM. THE TAGS OF POST
;            ARE:
;
;         PI - AN NGAUSS ARRAY CONTAINING THE RANDOM DRAWS FOR THE
;              GAUSSIAN WEIGHTS.
;         MU - AN 2 x NGAUSS ARRAY CONTAINING THE RANDOM DRAWS FOR THE
;              GAUSSIAN CENTROIDS. POST.MU[0,K] IS THE LUM-COMPONENT
;              OF THE MEAN FOR THE K-TH GAUSSIAN, AND POST.MU[1,K] IS
;              THE ZETA COMPONENT.
;      COVAR - AN 2 x 2 x NGAUSS ARRAY CONTAINING THE RANDOM DRAWS FOR
;              THE GAUSSIAN COVARIANCE MATRICES. POST.COVAR[0,0,K] IS
;              THE VARIANCE OF THE K-TH GAUSSIAN IN THE LUMINOSITY
;              DIRECTION, POST.COVAR[1,1,K] IS THE VARIANCE OF THE
;              H-TH GAUSSIAN IN THE ZETA DIRECTION, AND
;              POST.COVAR[0,1,K] IS THE COVARIANCE BETWEEN LUM AND
;              ZETA FOR THE K-TH GAUSSIAN.
;          N - A SCALAR CONTAINING THE RANDOM DRAWS FOR THE LUMINOSITY
;              FUNCTION NORMALIZATION, I.E., THE TOTAL NUMBER OF
;              OBJECTS.
;      GMEAN - THE RANDOM DRAWS FOR THE PRIOR MEAN ON THE GAUSSIAN
;              CENTROIDS. IN GENERAL THIS IS UNINTERESTING AND SHOULD
;              BE IGNORED. ONLY PRESENT IF UNIFORM = 0
;       AMAT - THE RANDOM DRAWS FOR THE PRIOR SCALE MATRIX ON THE
;              GAUSSIAN COVARIANCE MATRICES. AS WITH GMEAN, THIS
;              SHOULD USUALLY BE IGNORED. ONLY PRESENT IF UNIFORM =
;              0.
;    LOGPOST - THE VALUE OF THE LOGARITHM OF THE POSTERIOR PROBABILITY
;              DISTRIBUTION FOR EACH RANDOM SET OF LUMINOSITY FUNCTION
;              PARAMETERS.
;    DETPROB - THE DETECTION PROBABILITY FOR EACH RANDOM SET OF
;              LUMINOSITY FUNCTION PARAMETERS.
;
; OPTIONAL INPUT PARAMETERS
;
;     PLIM - A STRUCTURE CONTAINING VALUES FOR PRIOR CONSTRAINTS ON
;            THE MIXTURE OF GAUSSIAN FUNCTIONS PARAMETERS. THE TAGS
;            SHOULD BE:
;
;      MLLOW - THE LOWER LIMIT ON THE LUMINOSITY COMPONENT OF THE
;              GAUSSIAN CENTROIDS. THE DEFAULT IS -1D30.
;     MLHIGH - THE UPPER LIMIT ON THE LUMINOSITY COMPONENT OF THE
;              GAUSSIAN CENTROIDS. THE DEFAULT IS +1D30.
;      MZLOW - THE LOWER LIMIT ON THE ZETA COMPONENT OF THE
;              GAUSSIAN CENTROIDS. THE DEFAULT IS -1D30.
;     MZHIGH - THE UPPER LIMIT ON THE ZETA COMPONENT OF THE
;              GAUSSIAN CENTROIDS. THE DEFAULT IS +1D30.
;      VLLOW - THE LOWER LIMIT ON THE GAUSSIAN FUNCTION VARIANCES IN
;              THE LUMINOSITY DIRECTION. THIS SHOULD ALWAYS BE GREATER
;              THAN THE LUMINOSITY GRID SPACING USED FOR COMPUTING THE
;              SELECTION FUNCTION. THE DEFAULT IS 0.01.
;     VLHIGH - THE UPPER LIMIT ON THE GAUSSIAN FUNCTION VARIANCES IN
;              THE LUMINOSITY DIRECTION. THE DEFAULT IS +1D30.
;      VZLOW - THE LOWER LIMIT ON THE GAUSSIAN FUNCTION VARIANCES IN
;              THE ZETA DIRECTION. THIS SHOULD ALWAYS BE GREATER THAN
;              THE ZETA GRID SPACING USED FOR COMPUTING THE SELECTION
;              FUNCTION. THE DEFAULT IS 0.05^2.
;     VZHIGH - THE UPPER LIMIT ON THE GAUSSIAN FUNCTION VARIANCES IN
;              THE ZETA DIRECTION. THE DEFAULT IS +1D30.
;       CLOW - THE LOWER LIMIT ON THE GAUSSIAN FUNCTION
;              CORRELATIONS. THE DEFAULT IS -0.98.
;      CHIGH - THE UPPER LIMIT ON THE GAUSSIAN FUNCTION
;              CORRELATIONS. THE DEFAULT IS +0.98.
;
;  UNIFORM - IF UNIFORM = 1, THEN ASSUME A UNIFORM PRIOR ON ALL
;            PARAMETERS. OTHERWISE, ASSUME THE PRIOR DESCRIBED IN
;            Kelly, Fan, & Vestergaard (2008). THE DEFAULT IS UNIFORM
;            = 0, I.E., USE THE KFV08 PRIOR.
;  NCHAINS - THE NUMBER OF CHAINS TO USE IN THE MARKOV CHAIN MONTE
;            CARLO. THE DEFAULT IS 1.
;  MAXITER - THE NUMBER OF MCMC ITERATIONS TO PERFORM AFTER BURN-IN
;            BEFORE. THE DEFAULT IS MAXITER = 5E4.
; BURNITER - THE NUMBER OF BURN-IN ITERATIONS TO PERFORM. THE BURN-IN
;            STAGE IS NECESSARY IN ORDER FOR THE MCMC TO 'FORGET'
;            ABOUT ITS INITIAL CONDITIONS AND CONVERGE TO THE
;            POSTERIOR DISTRIBUTION. THE DEFAULT IS BURNITER = 2E4.
;   SILENT - IF SILENT = 1, THEN BE QUIET. THE DEFAULT IS SILENT = 0.
;
; NOTES:
;
;   PLEASE REPORT ANY BUGS TO B.KELLY (bkelly@as.arizona.edu). 
;   ANY QUESTIONS SHOULD ALSO BE DIRECTED TO B.KELLY.
;
; ROUTINES CALLED:
;
;   MRANDOMN, RANDOMDIR, RANDOMWISH, MRANDOMT, NEGBIN, IDL ASTRO
;   LIBRARY ROUTINES
;
; REVISION HISTORY
;     Written by Brandon C. Kelly, Steward Observatory, May 2008
;
; REFERENCE:
;
;     "A FLEXIBLE METHOD OF ESTIMATING LUMINOSITY FUNCTIONS", KELLY,
;        B.C., FAN, X., & VESTERGAARD, M. 2008, ACCEPTED BY ApJ,
;        (arXiv:0805.2946)
;-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                                                                      ;
;                   ROUTINE TO CONSTRUCT THE SELECTION FUNCTION.                       ;
;                                                                                      ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function selection_function

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Simple example for the selection function for a quasar survey over z
; < 6 and an nu * F_nu > 1d-13. K-corrections are ignored, as this
; example is only included to illustrate how to construct a selection
; function.
;
;fmin = -13.0
;fmax = -8.0
;zmin = -1.0
;zmax = alog10(6.0)
;
;ngrid = 128L
;zgrid = zmin + (zmax - zmin) * dindgen(ngrid) / (ngrid - 1)
;fgrid = fmin + (fmax - fmin) * dindgen(ngrid) / (ngrid - 1)
;
;                                ;get luminosity distance, assume
;                                ;standard cosmology
;lumdist = lumdist( 10d^zgrid, h0=71.0, omega_m=0.30, lambda0=0.70 )
;
;pc = 3.086d18                   ;parsec, in cm
;clight = 3d18                   ;speed of light in angstroms / sec
;
;lumdist = lumdist * 1d6 * pc ;luminosity distance is in Mpc, convert to cm
;
;lgrid = dblarr(ngrid * ngrid)
;                                ;convert flux densities to
;                                ;log-luminosities: nu * L_nu
;for i = 0, ngrid - 1 do begin
;    
;    log_flux = alog10(fgrid) + (1d - alpha) * alog10(2500 * (1d + 10d^zgrid[i]) / 7481) ;log(nu * f_nu)
;
;    lgrid[ngrid*i] = log_flux + alog10(4d * !pi) + 2d *
;    alog10(lumdist[i])
;
;endfor
;
;                                ;make zgrid into an array, and then a
;                                ;vector. now the values of zgrid
;                                ;correspond to the values of lgrid
;zgrid = (zgrid ## replicate(1d, ngrid))[*]
;
;                                ;detection probability is uniform for
;                                ;fmin < alog10(nu * F_nu) < fmax
;detprob = replicate(1d, ngrid * ngrid)
;
;selfunc = {detprob:detprob, zgrid:zgrid, lgrid:lgrid}
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

return, selfunc
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                                                      ;
; Routine to compute the inverse of a 2 x 2 matrix using Cramer's Rule ;
;                                                                      ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function cramer_invert, A

determ = A[0,0] * A[1,1] - A[0,1] * A[1,0]

A_inv = [[A[1,1], -A[0,1]], [-A[1,0], A[0,0]]]

return, A_inv / determ
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                                               ;
; Routine to do Metropolis-Hastings update                      ;
;                                                               ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function metro_update, logpost_old, logpost_new, jump_oldnew, jump_newold, $
                  seed = seed

if n_elements(jump_oldnew) eq 0 then jump_oldnew = 0d
if n_elements(jump_newold) eq 0 then jump_newold = 0d

lograt = logpost_new - jump_oldnew - (logpost_old - jump_newold)

accept = 0

if lograt gt 0 then accept = 1 else begin
    
    u = randomu(seed)

    if alog(u) le lograt then accept = 1
    
endelse

return, accept
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                                            ;
; Routine to get the detection probability, given THETA.     ;
;                                                            ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function detprob_theta, theta, selfunc, gind

nobs = n_elements(lum)
ngauss = n_elements(theta) / 6

pind = indgen(ngauss) ;induces for pi in theta
mlind = ngauss + indgen(ngauss) ;induces for mu_l in theta
mzind = 2 * ngauss + indgen(ngauss) ;induces for mu_z in theta
vlind = 3 * ngauss + indgen(ngauss) ; induces for var_l in theta
vzind = 4 * ngauss + indgen(ngauss) ; induces for var_z in theta
cind = 5 * ngauss + indgen(ngauss) ;induces for covar_lz in theta

pi = theta[pind]
lmu = theta[mlind]
zmu = theta[mzind]
lvar = theta[vlind]
zvar = theta[vzind]
lzcov = theta[cind]
corrlz = lzcov / sqrt(lvar * zvar)

lgrid = selfunc.lgrid
zgrid = selfunc.zgrid

sorted = sort(lgrid)
unique = uniq(lgrid[sorted])
dl = lgrid[sorted[unique[1]]] - lgrid[sorted[unique[0]]]

sorted = sort(zgrid)
unique = uniq(zgrid[sorted])
dz = zgrid[sorted[unique[1]]] - zgrid[sorted[unique[0]]]

detprob = dblarr(n_elements(gind))

for k = 0, n_elements(gind) - 1 do begin

    log_zmarg = -0.5 * alog(2d * !pi * zvar[gind[k]]) - $
      0.5 * (zgrid - zmu[gind[k]])^2 / zvar[gind[k]]
    
    clmean = lmu[gind[k]] + lzcov[gind[k]] / zvar[gind[k]] * (zgrid - zmu[gind[k]])

    clvar = lvar[gind[k]] - lzcov[gind[k]]^2 / zvar[gind[k]]

                                ;get detection probability for
                                ;each mixture component

    determ = zvar[gind[k]] * lvar[gind[k]] - lzcov[gind[k]]^2

    arg = (zgrid - zmu[gind[k]])^2 / zvar[gind[k]] + (lgrid - lmu[gind[k]])^2 / lvar[gind[k]] - $
      2d * lzcov[gind[k]] * (lgrid - lmu[gind[k]]) * (zgrid - zmu[gind[k]]) / $
      (lvar[gind[k]] * zvar[gind[k]])

    prob = 1d / (2d * !pi * sqrt(determ)) * exp( -0.5 * arg / (1d - corrlz[gind[k]]^2) )

    zero = where(finite(prob) eq 0, nzero)
    if nzero gt 0 then prob[zero] = 0d

    detprob[k] = total( prob * selfunc.detprob, /nan ) * dl * dz

endfor

return, detprob
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                                                  ;
;          ROUTINE TO GET THE LOG-POSTERIOR                        ;
                                                                   ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro logpost_theta, theta, lum, zeta, fhatj, detprob, gind, selfunc, gmean, $
                   Amat, logpost, uniform

nobs = n_elements(lum)

ngauss = n_elements(theta) / 6

pind = indgen(ngauss) ;induces for pi in theta
mlind = ngauss + indgen(ngauss) ;induces for mu_m in theta
mzind = 2 * ngauss + indgen(ngauss) ;induces for mu_z in theta
vlind = 3 * ngauss + indgen(ngauss) ; induces for var_m in theta
vzind = 4 * ngauss + indgen(ngauss) ; induces for var_z in theta
cind = 5 * ngauss + indgen(ngauss) ;induces for corr_mz in theta

pi = theta[pind]
lmu = theta[mlind]
zmu = theta[mzind]
lvar = theta[vlind]
zvar = theta[vzind]
lzcov = theta[cind]
corrlz = lzcov / sqrt(lvar * zvar)

                                ;get estimated density of the observed
                                ;data

for j = 0, n_elements(gind) - 1 do begin

    determ = zvar[gind[j]] * lvar[gind[j]] - lzcov[gind[j]]^2

    arg = (zeta - zmu[gind[j]])^2 / zvar[gind[j]] + (lum - lmu[gind[j]])^2 / lvar[gind[j]] - $
      2d * lzcov[gind[j]] * (lum - lmu[gind[j]]) * (zeta - zmu[gind[j]]) / $
      (lvar[gind[j]] * zvar[gind[j]])

    logprob = -0.5 * alog(4d * !pi^2 * determ) - 0.5 * arg / (1d - corrlz[gind[j]]^2)

    fhatj[0,gind[j]] = pi[gind[j]] * exp(logprob)

    detprob[gind[j]] = detprob_theta( theta, selfunc, gind[j] )
    
endfor

fhat = ngauss gt 1 ? total(fhatj, 2, /nan) : fhatj

detprob0 = total(pi * detprob)

if detprob0 le 1d-33 then begin

    logpost = -1d300
    return

endif

loglik = total( alog(fhat) ) - nobs * alog( detprob0 )

if finite(loglik) eq 0 then loglik = -1d300

if not uniform then begin
                                ;log of the prior for the mixture
                                ;model parameters

;; now get log-prior density

    lzvar = dblarr(2, 2, ngauss)
    lzvar[0,0,*] = lvar
    lzvar[1,1,*] = zvar
    lzvar[0,1,*] = lzcov
    lzvar[1,0,*] = lzcov
    
    lzvar_inv = lzvar
    for k = 0, ngauss - 1 do lzvar_inv[*,*,k] = cramer_invert(lzvar[*,*,k])
    
;Gaussian means have Cauchy prior
    
    logprior_mu = 0d
    scale = 1d^2
    
    gcovar_inv = dblarr(2, 2)
    for k = 0, ngauss - 1 do gcovar_inv = gcovar_inv + lzvar_inv[*,*,k] / ngauss
    gcovar = scale * cramer_invert(gcovar_inv)
    gcovar_inv = gcovar_inv / scale
    
    gcovdet = gcovar[0,0] * gcovar[1,1] - gcovar[0,1] * gcovar[1,0]
    
    for k = 0, ngauss - 1 do logprior_mu = logprior_mu - 0.5 * alog(gcovdet) - $
      1.5 * alog( 1d + ([lmu[k], zmu[k]] - gmean) ## gcovar_inv ## $
                  transpose([lmu[k], zmu[k]] - gmean) )
    
    logprior_mu = logprior_mu[0] ;make sure logprior_mu is scalar
    
;Gaussian covariances have independent inverse-wishart prior with
;common scale matrix Amat and DOF = 1
    
    logprior_sigma = 0d
    
    Amat_det = Amat[0,0] * Amat[1,1] - Amat[0,1] * Amat[1,0]
    
    for k = 0, ngauss - 1 do begin
        
        lzvardet = lzvar[0,0,k] * lzvar[1,1,k] - lzvar[0,1,k] * lzvar[1,0,k]
        
        logprior_sigma = logprior_sigma + $
          0.5 * alog(Amat_det) - 2d * alog(lzvardet) - 0.5 * trace(Amat ## lzvar_inv[*,*,k])
        
    endfor
    
    logprior = logprior_mu + logprior_sigma

endif else logprior = 0d

logpost = logprior + loglik

return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                                                      ;
; Routine to compute conditional posterior for prior hyper-parameters  ;
;                                                                      ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function hyper_prior, zgrid, lgrid, dz, dl, lmu, zmu, gmean, gcovar_inv, zeroind

ngauss = n_elements(lmu)

gcovar = cramer_invert(gcovar_inv)

lmu0 = gmean[0]
zmu0 = gmean[1]
lvar0 = gcovar[0,0]
lzcov0 = gcovar[0,1]
zvar0 = gcovar[1,1]
lzcorr0 = lzcov0 / sqrt(lvar0 * zvar0)

zsqr = (lmu - lmu0)^2 / lvar0 + (zmu - zmu0)^2 / zvar0 - $
  2d * lzcorr0 * (lmu - lmu0) * (zmu - zmu0) / sqrt(lvar0 * zvar0)

gcovdet = gcovar[0,0] * gcovar[1,1] - gcovar[0,1] * gcovar[1,0]

logprior = -0.5 * alog( gcovdet ) - 1.5d * alog( 1d + zsqr / (1d - lzcorr0^2) )

logprior = total(logprior)

zsqr = (lgrid - lmu0)^2 / lvar0 + (zgrid - zmu0)^2 / zvar0 - $
  2d * lzcorr0 * (lgrid - lmu0) * (zgrid - zmu0) / sqrt(lvar0 * zvar0)

gcovinv_det = 1d / gcovdet

norm = sqrt( gcovinv_det ) * ( 1d + zsqr / (1d - lzcorr0^2) )^(-1.5d)
norm[zeroind] = 0d
norm = total(norm * dl * dz, /nan) > 1d-300

return, logprior - ngauss * alog(norm)
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                                                      ;
; Routine to report the acceptance rates for the Metropolis-Hastings   ;
; Algorithm.                                                           ;
;                                                                      ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro mha_report, naccept, ntries, ngauss, uniform

print, 'Acceptance Rates Since Last Report Are: '
print, ''
print, 'Gaussian Weights, PI:'
print, float(naccept[0]) / ntries

print, 'Gaussian Means (Gaussian 1, Gaussian 2, etc.)'
print, float(naccept[1:ngauss]) / ntries

print, 'Gaussian Covariance Matrices (Gaussian 1, Gaussian 2, etc.)'
print, float(reform(naccept[1+ngauss:2*ngauss])) / ntries

if not uniform then begin
    
    print, 'Cauchy Prior Mean:'
    print, float(naccept[1+2*ngauss]) / ntries
    
endif

print, ''

return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                                                                      ;
;                                    MAIN ROUTINE                                      ;
;                                                                                      ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro mixlfunc, lum, zeta, ngauss, skyfrac, post, plim=plim, uniform=uniform, $
              nchains=nchains, maxiter=maxiter, burniter=burniter, silent=silent

if n_params() lt 4 then begin
    print, 'Syntax- mixlfunc, lum, zeta, ngauss, skyfrac, post, plim=plim, /uniform,'
    print, '                  nchains=nchains, maxiter=maxiter, burniter=burniter, /silent'
    return
endif

if not keyword_set(uniform) then uniform = 0
if not keyword_set(silent) then silent = 0
if n_elements(plim) eq 0 then plim = {mllow:-1d30, mlhigh:1d30, mzlow:-1d30, mzhigh:1d30, vllow:0.1^2, $
                                      vlhigh:1d30, vzlow:0.05^2, vzhigh:1d30, clow:-0.98, chigh:0.98}

;check PLIM structure

ntags = 0
ntags = ntags + tag_exist(plim, 'mllow')
ntags = ntags + tag_exist(plim, 'mlhigh')
ntags = ntags + tag_exist(plim, 'vllow')
ntags = ntags + tag_exist(plim, 'vlhigh')
ntags = ntags + tag_exist(plim, 'mzlow')
ntags = ntags + tag_exist(plim, 'mzhigh')
ntags = ntags + tag_exist(plim, 'vzlow')
ntags = ntags + tag_exist(plim, 'vzhigh')
ntags = ntags + tag_exist(plim, 'clow')
ntags = ntags + tag_exist(plim, 'chigh')

if ntags lt 10 then begin
    print, 'PLIM is missing some tags.'
    return
endif

nobs = n_elements(lum)

if n_elements(zeta) ne nobs then begin
    print, 'ZETA and LUM must be of the same size.'
    return
endif

selfunc = selection_function() ;get survey selection function

;get initial parameters to start the Markov chain

if n_elements(nchains) eq 0 then nchains = 1L
if n_elements(maxiter) eq 0 then maxiter = 50000L
if n_elements(burniter) eq 0 then burniter = 20000L

npar = 6 * ngauss

theta = dblarr(npar, nchains)
pind = indgen(ngauss) ;induces for pi in theta
mlind = ngauss + indgen(ngauss) ;induces for mu_l in theta
mzind = 2 * ngauss + indgen(ngauss) ;induces for mu_z in theta
vlind = 3 * ngauss + indgen(ngauss) ; induces for var_l in theta
vzind = 4 * ngauss + indgen(ngauss) ; induces for var_z in theta
cind = 5 * ngauss + indgen(ngauss) ;induces for corr_lz in theta

Amat = fltarr(2, 2, nchains) ;scale matrix of prior on Gaussian covariances
gmean = fltarr(2, nchains) ;mean for gaussian centroids

scale = 1d^2 ;scale for the prior Cauchy prior covariance matrix, ignored if UNIFORM = 1

;get initial values for mixture parameters

pimin = 1d-5                    ;minimum value of gaussian weights, PI
coef = linfit(zeta, lum)
ulimit = 2d * stddev(lum - coef[0] - coef[1] * zeta)
anchor = [plim.mzlow, coef[0] + coef[1] * plim.mzlow]

zmu = dblarr(ngauss)
lmu = dblarr(ngauss)
covar = dblarr(2,2,ngauss)

for i = 0, nchains - 1 do begin

    repeat begin

        nu = randomu(seed, ngauss, gamma = 1)
        pi = nu / total(nu)

        bad = where(pi le pimin, nbad)
        
    endrep until nbad eq 0

    ind = (permute(nobs))[0:ngauss-1]

    for k = 0, ngauss - 1 do begin
        
        zmu[k] = zeta[ind[k]]
        lhigh = plim.mlhigh < (ulimit + coef[0] + coef[1] * zmu[k])
        llow = plim.mllow > (coef[0] + coef[1] * zmu[k] - 2.5 * ulimit)
        lmu[k] = (lhigh - llow) * randomu(seed) + llow

    endfor

    zsig = stddev(zeta)
    lsig = stddev(lum)
    covar0 = dblarr(2,2)

    covar0[0,0] = lsig^2
    covar0[1,1] = zsig^2
    covar0[0,1] = 0.5 * lsig * zsig
    covar0[1,0] = 0.5 * lsig * zsig

    for k = 0, ngauss - 1 do begin

        repeat begin
            
            covar[*,*,k] = randomwish(seed, 10, covar0 / 10)
            lvar = covar[0,0,k]
            zvar = covar[1,1,k]
            lzcorr = covar[0,1,k] / sqrt(zvar * lvar)

            move_on = 1
            if lvar gt plim.vlhigh or lvar lt plim.vllow or $
              zvar gt plim.vzhigh or zvar lt plim.vzlow or $
              lzcorr gt plim.chigh or lzcorr lt plim.clow $
              then move_on = 0

        endrep until move_on

    endfor

    dist = abs(zmu - anchor[0]) + abs(lmu - anchor[1])
    sorted = sort(dist)
    zmu = zmu[sorted]
    lmu = lmu[sorted]

    theta[pind,i] = pi
    theta[mlind,i] = lmu
    theta[mzind,i] = zmu
    theta[vlind,i] = covar[0,0,*]
    theta[vzind,i] = covar[1,1,*]
    theta[cind,i] = covar[0,1,*]

endfor

if not uniform then begin
    ;get initial guess for prior parameters
    for i = 0, nchains - 1 do begin
        
        gcovar = dblarr(2,2)
        for k = 0, ngauss - 1 do begin
            
            lzvar = [[theta[vlind[k],i], theta[cind[k],i]], $
                     [theta[cind[k],i], theta[vzind[k],i]]]
            gcovar = gcovar + cramer_invert(lzvar) / ngauss
            
        endfor
        
        gcovar = cramer_invert(gcovar)
        
        gmean[*,i] = [mean(theta[mlind,i]), mean(theta[mzind,i])]
        Amat[*,*,i] = randomwish(seed, 4, gcovar / 4)
        
    endfor

    hmllow = min(lum) - (max(lum) - min(lum)) > plim.mllow
    hmlhigh = max(lum) + (max(lum) - min(lum)) < plim.mlhigh
    hmzlow = min(zeta) - (max(zeta) - min(zeta)) > plim.mzlow
    hmzhigh = max(zeta) + (max(zeta) - min(zeta)) < plim.mzhigh

    zgrid_hyper = hmzlow + (hmzhigh - hmzlow) * dindgen(32) / 31.0
    dzhyper = zgrid_hyper[1] - zgrid_hyper[0]
    lgrid_hyper = hmllow + (hmlhigh - hmllow) * dindgen(32) / 31.0
    dlhyper = lgrid_hyper[1] - lgrid_hyper[0]
    
    zgrid_hyper = (transpose(zgrid_hyper ## replicate(1, 32)))[*]
    lgrid_hyper = (lgrid_hyper ## replicate(1, 32))[*]
    zeroind = where(lgrid_hyper ge (coef[0] + coef[1] * zgrid_hyper + ulimit))

endif

jvar_mu = dblarr(2, 2, ngauss)
A_mu = dblarr(2, 2, ngauss)
diag = [0,3]

for k = 0, ngauss - 1 do begin
    
    jvar_mu[0,0,k] = 1d-3
    jvar_mu[1,1,k] = 1d-3

    A_mu0 = jvar_mu[*,*,k]
        
    choldc, A_mu0, P, /double   ;do Cholesky decomposition
    
    for j = 0, 2 - 1 do for l = j, 2 - 1 do A_mu0[l,j] = 0d
    
    A_mu0[diag] = P
    A_mu[*,*,k] = A_mu0
    
endfor

jumpdof = replicate(long(nobs) / ngauss, ngauss)
rfact = 10d

niter = 0L
nburn = 2500L
nskip = 10L                  ;only save every 10 iterations of the MHA
checkiter = 2500L
convergence = 0
ndraws = 0L
ntries = 0L
burnin = 1
naccept = lonarr(2 + 2 * ngauss)

logpost_draw = dblarr(maxiter / nskip, nchains)
theta_draw = dblarr(maxiter / nskip, nchains, npar)

mudraw = dblarr(2, nburn, ngauss, nchains)
nmdraw = 0L

logpost = dblarr(nchains)
fhatj = dblarr(nobs, ngauss, nchains)
detprob = dblarr(ngauss, nchains)

if not uniform then begin

    gmean_draw = dblarr(maxiter / nskip, nchains, 2)
    Amat_draw = dblarr(maxiter / nskip, nchains, 4)

endif

for i = 0, nchains - 1 do begin

    fhatj0 = fhatj[*,*,i]
    detprob0 = detprob[*,i]
    gind = indgen(ngauss)
    
    logpost_theta, theta[*,i], lum, zeta, fhatj0, detprob0, gind, selfunc, $
      gmean[*,i], Amat[*,*,i], logpost0, uniform
    
    logpost[i] = logpost0
    fhatj[*,*,i] = fhatj0
    detprob[0,i] = detprob0

endfor

if not silent then print, 'Performing Markov Chain Simulation...'

repeat begin
                                ;first draw a proposal for pi from
                                ;its jumping density
    if not silent and niter gt 0 and (niter mod 1000) eq 0 then print, niter

    for i = 0, nchains - 1 do begin

                                ;update posterior before starting
                                ;iteration
        gind = indgen(ngauss)
        fhatj0 = fhatj[*,*,i]
        detprob0 = detprob[*,i]

        logpost_theta, theta[*,i], lum, zeta, fhatj0, detprob0, gind, selfunc, $
          gmean[*,i], Amat[*,*,i], logpost0, uniform

        fhatj[*,*,i] = fhatj0
        detprob[0,i] = detprob0            
        logpost[i] = logpost0
        
        if ngauss gt 1 then begin
                                ;First update Gaussian weights,
                                ;PI. Draw proposed values of PI from
                                ;Dirichlet distribution

            pi = reform(theta[pind,i])
            alpha = nobs * pi / rfact + 1
            pi_prop = reform(randomdir( seed, alpha ))
            
            bad = where(pi_prop le pimin, nbad)
            
            if nbad eq 0 then begin
                
                alpha_prop = nobs * pi_prop / rfact + 1
                
                log_jfor = lngamma(total(alpha)) - total(lngamma(alpha)) + $
                  total( (alpha - 1d) * alog(pi_prop) )
                
                log_jback = lngamma(total(alpha_prop)) - total(lngamma(alpha_prop)) + $
                  total( (alpha_prop - 1d) * alog(pi) )
                
                tprop = theta[*,i]
                tprop[pind] = pi_prop
                
                fhatj0 = fhatj[*,*,i]
                detprob0 = detprob[*,i]
                gind = indgen(ngauss)
                
                logpost_theta, tprop, lum, zeta, fhatj0, detprob0, gind, selfunc, $
                  gmean[*,i], Amat[*,*,i], logpost_prop, uniform
                
                accept = metro_update( logpost[i], logpost_prop, log_jfor, log_jback, seed=seed )
                
                if accept eq 1 then begin
                                ;keep this proposal
                    naccept[0] = naccept[0] + 1
                    logpost[i] = logpost_prop
                    theta[pind,i] = pi_prop
                    fhatj[*,*,i] = fhatj0
                    detprob[*,i] = detprob0
                    
                endif

            endif

        endif
                                ;now update the gaussian means and
                                ;variances, one gaussian at a time
        for k = 0, ngauss - 1 do begin    

                                ;first do gaussian means
            
            lzmu = [theta[mlind[k],i], theta[mzind[k],i]]
            
                                ;draw proposal for (mu_l,mu_z) from
                                ;truncated normal

            lzprop = lzmu + mrandomt(seed, 4, jvar_mu[*,*,k], A=A_mu[*,*,k])

            if lzprop[1] ge plim.mzlow and lzprop[1] le plim.mzhigh and $
              lzprop[0] ge plim.mllow and lzprop[0] le plim.mlhigh then begin
                
                tprop = theta[*,i]
                tprop[mlind[k]] = lzprop[0]
                tprop[mzind[k]] = lzprop[1]
                
                fhatj0 = fhatj[*,*,i]
                detprob0 = detprob[*,i]

                gind = k
                
                logpost_theta, tprop, lum, zeta, fhatj0, detprob0, gind, selfunc, $
                  gmean[*,i], Amat[*,*,i], logpost_prop, uniform
                
                accept = metro_update( logpost[i], logpost_prop, seed=seed )
            
                if accept eq 1 then begin
                                ;keep this proposal
                    naccept[1 + k] = naccept[1 + k] + 1
                    logpost[i] = logpost_prop
                    theta[mlind[k],i] = lzprop[0]
                    theta[mzind[k],i] = lzprop[1]
                    fhatj[*,*,i] = fhatj0
                    detprob[0,i] = detprob0

                endif

            endif
                                ;Now get proposal for gaussian
                                ;covariance matrix from Wishart
                                ;distribution
            lzvar = [[theta[vlind[k],i], theta[cind[k],i]], $
                     [theta[cind[k],i], theta[vzind[k],i]]]
            lzvar_inv = cramer_invert(lzvar)

            vprop = randomwish(seed, jumpdof[k], lzvar / jumpdof[k])

            cprop = vprop[0,1] / sqrt(vprop[0,0] * vprop[1,1])

            if vprop[0,0] ge plim.vllow and vprop[0,0] le plim.vlhigh and $
              vprop[1,1] ge plim.vzlow and vprop[1,1] le plim.vzhigh and $
              cprop ge plim.clow and cprop le plim.chigh then begin
 
                vprop_inv = cramer_invert(vprop)
                
                vjumpdet = (lzvar[0,0] * lzvar[1,1] - lzvar[0,1] * lzvar[1,0]) / jumpdof[k]^2
                vpropdet = vprop[0,0] * vprop[1,1] - vprop[0,1] * vprop[1,0]

                log_jfor = -jumpdof[k] / 2d * alog( vjumpdet ) + $
                  (jumpdof[k] - 3d) / 2d * alog( vpropdet ) - $
                  0.5 * trace( (lzvar_inv * jumpdof[k]) ## vprop )
                
                vjumpdet2 = vpropdet / jumpdof[k]^2

                lzvardet = lzvar[0,0] * lzvar[1,1] - lzvar[0,1] * lzvar[1,0]

                log_jback = -jumpdof[k] / 2d * alog( vjumpdet2 ) + $
                  (jumpdof[k] - 3d) / 2d * alog( lzvardet ) - $
                  0.5 * trace( (vprop_inv * jumpdof[k]) ## lzvar )
                
                tprop = theta[*,i]
                tprop[vlind[k]] = vprop[0,0]
                tprop[vzind[k]] = vprop[1,1]
                tprop[cind[k]] = vprop[0,1]
                
                gind = k

                fhatj0 = fhatj[*,*,i]
                detprob0 = detprob[*,i]
                
                logpost_theta, tprop, lum, zeta, fhatj0, detprob0, gind, selfunc, $
                  gmean[*,i], Amat[*,*,i], logpost_prop, uniform
                
                accept = metro_update( logpost[i], logpost_prop, log_jback, log_jfor, seed=seed )

                if accept eq 1 then begin
                                ;keep this proposal
                    naccept[1 + k + ngauss] = naccept[1 + k + ngauss] + 1
                    logpost[i] = logpost_prop
                    theta[vlind[k],i] = vprop[0,0]
                    theta[vzind[k],i] = vprop[1,1]
                    theta[cind[k],i] = vprop[0,1]
                    fhatj[*,*,i] = fhatj0
                    detprob[0,i] = detprob0
                    
                endif

            endif
        
        endfor
    
        if not uniform then begin

;; update prior parameters, do gibbs sampler
        
            lzmu = [[theta[mlind,i]], [theta[mzind,i]]]
            
            gcovar_inv = dblarr(2,2)
            for k = 0, ngauss - 1 do begin
                
                lzvar = [[theta[vlind[k],i], theta[cind[k],i]], $
                         [theta[cind[k],i], theta[vzind[k],i]]]
                gcovar_inv = gcovar_inv + cramer_invert(lzvar) / ngauss
                
            endfor
            
            gcovar_inv = gcovar_inv / scale
            
            logprior_old = hyper_prior(zgrid_hyper, lgrid_hyper, dzhyper, dlhyper, lzmu[*,0], $
                                       lzmu[*,1], gmean[*,i], gcovar_inv, zeroind)
            
                                ;draw new value of GMEAN from
                                ;multivariate normal density
            jvar_gmean = correlate(transpose(lzmu), /covar)
            ridge = 0.1 ;add ridge coefficient to matrix is not ill-conditioned
            jvar_gmean = jvar_gmean + ridge * [[jvar_gmean[0,0], 0], [0, jvar_gmean[1,1]]]
            
            gmprop = gmean[*,i] + mrandomn(seed, jvar_gmean / ngauss)
    
            bad = where(gmprop[0] le plim.mllow or gmprop[0] ge plim.mlhigh, nbad1)
            bad = where(gmprop[1] le plim.mzlow or gmprop[1] ge plim.mzhigh, nbad2)
 
            if nbad1 eq 0 and nbad2 eq 0 then begin
                
                logprior_new = hyper_prior(zgrid_hyper, lgrid_hyper, dzhyper, dlhyper, $
                                           lzmu[*,0], lzmu[*,1], gmprop, gcovar_inv, zeroind)
                
                lograt = logprior_new - logprior_old
                
                accept = 0
                
                if lograt gt 0 then accept = 1 else begin
                    
                    u = randomu(seed)
                    
                    if alog(u) le lograt then accept = 1
                    
                endelse
                
                if accept eq 1 then begin
                                ;keep this proposal
                    naccept[1 + 2 * ngauss] = naccept[1 + 2 * ngauss] + 1
                    gmean[*,i] = gmprop
                    logprior_old = logprior_new
                    
                endif

            endif
            
                                ;finally, draw new values of Amat
            Ahat = 0d
            dofA = ngauss + 3
            
            for k = 0, ngauss - 1 do begin
                
                lzvar = [[theta[vlind[k],i], theta[cind[k],i]], $
                         [theta[cind[k],i], theta[vzind[k],i]]]
                
                Ahat = Ahat + cramer_invert(lzvar)
                
            endfor
            
            Ahat = cramer_invert(Ahat)
                                ;draw new value of Amat from Wishart
                                ;distribution
            Amat[*,*,i] = randomwish(seed, dofA, Ahat)
            
        endif

        if burnin gt 0 then begin
            
            for k = 0, ngauss - 1 do begin
                
                mind = [mlind[k], mzind[k]]
                mudraw[*, nmdraw, k, i] = theta[mind,i]
                
            endfor
            
            if n_elements(logpost_arr) eq 0 then logpost_arr = logpost $
            else logpost_arr = [logpost_arr, logpost]
            
            detprobi = total(theta[pind,i] * detprob[*,i])

            if n_elements(detprob_arr) eq 0 then detprob_arr = detprobi $
            else detprob_arr = [detprob_arr, detprobi]
            

        endif

    endfor

    if burnin then nmdraw = nmdraw + 1L

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; markov chain iteration complete ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    niter = niter + 1
    ntries = ntries + nchains

    if not burnin and (niter mod nskip) eq 0 then begin

        for i = 0, nchains - 1 do begin

            theta_draw[ndraws,i,*] = theta[*,i]

            if not uniform then begin

                gmean_draw[ndraws,i,*] = gmean[*,i]
                Amat_draw[ndraws,i,*] = Amat[*,*,i]

            endif

            logpost_draw[ndraws,i] = logpost[i]

        endfor
    
        ndraws = ndraws + 1
        
    endif

    if not burnin and niter eq checkiter then begin
                                ;report acceptance rates for MHA
        if not silent then mha_report, naccept, ntries, ngauss, uniform

        ntries = 0L
        naccept = lonarr(2 + 2 * ngauss, nchains)
    
        checkiter = checkiter + 2500L

    endif

    if niter eq burniter and burnin ge 1 then begin
                                ;burnin complete
        ndraws = 0L
        niter = 0L
        burnin = 0
        
        if not silent then begin

            print, ''
            print, '-------------------------------------------------'
            print, '-                BURN-IN COMPLETE               -'
            print, '-------------------------------------------------'
            print, ''

        endif

        delvarx, mudraw, logpost_arr
        burniter = 0L

    endif

    if niter eq nburn and burnin ge 1 then begin

        if not silent then mha_report, naccept, ntries, ngauss, uniform
        
;; update MHA jumping distribution parameters automatically

                                ;first update covariance matrices for
                                ;gaussian mean jumps

        for k = 0, ngauss - 1 do begin
            
            if naccept[1 + k] / float(ntries) le 0.05 then $
              jvar_mu[*,*,k] = jvar_mu[*,*,k] / 1d2 $
            else begin

                jvar_mu1 = 0d
                for i = 0, nchains - 1 do jvar_mu1 = jvar_mu1 + $
                  cramer_invert( correlate(mudraw[*,*,k,i], /cov) )
                jvar_mu1 = nchains * cramer_invert(jvar_mu1)
                
                if la_determ(jvar_mu1, /double, /check) gt 0d then $
                  jvar_mu[*,*,k] = jvar_mu1 else $
                  jvar_mu[*,*,k] = jvar_mu[*,*,k] / 1d4
                
            endelse
            
                                ;now update degrees of freedom for
                                ;Wishart jumping distribution

            arate = naccept[1 + k + ngauss] / float(ntries)
            
            case 1 of
                
                arate le 0.01 : jumpdof[k] = jumpdof[k] * 6L
                arate gt 0.01 and arate le 0.05 : jumpdof[k] = jumpdof[k] * 4L
                arate gt 0.05 and arate le 0.10 : jumpdof[k] = jumpdof[k] * 3L
                arate gt 0.10 and arate le 0.20 : jumpdof[k] = jumpdof[k] * 2L
                arate ge 0.70 : jumpdof[k] = jumpdof[k] / 6L
                arate lt 0.70 and arate ge 0.50 : jumpdof[k] = jumpdof[k] / 4L
                arate lt 0.50 and arate ge 0.40 : jumpdof[k] = jumpdof[k] / 2L
                else : jumpdof[k] = jumpdof[k]
                
            endcase
            
        endfor
                                ;finally, update reduction factor for
                                ;controlling the Dirichlet jumping
                                ;distribution

        arate = naccept[0] / float(ntries)
        
        case 1 of
            
            arate le 0.01 : rfact = rfact / 10L
            arate gt 0.01 and arate le 0.05 : rfact = rfact / 6L
            arate gt 0.05 and arate le 0.10 : rfact = rfact / 3L
            arate gt 0.10 and arate le 0.20 : rfact = rfact / 2L
            arate ge 0.70 : rfact = rfact * 10L
            arate lt 0.70 and arate ge 0.50 : rfact = rfact * 4L
            arate lt 0.50 and arate ge 0.40 : rfact = rfact * 2L
            else : rfact = rfact
            
        endcase

        burnwait = 2500L

        nburn = niter + burnwait

        ntries = 0L
        naccept = lonarr(2 + 2 * ngauss)

        A_mu = dblarr(2, 2, ngauss)
        diag = [0,3]

        for k = 0, ngauss - 1 do begin
                                ;store covariance into dummy variable
                                ;for input into CHOLDC
            A_mu0 = jvar_mu[*,*,k]
            
            choldc, A_mu0, P, /double ;do Cholesky decomposition
            
            for j = 0, 2 - 1 do for l = j, 2 - 1 do A_mu0[l,j] = 0d
            
            A_mu0[diag] = P
            A_mu[*,*,k] = A_mu0
            
        endfor
                
        mudraw = dblarr(2, burnwait, ngauss, nchains)
        nmdraw = 0L
        
    endif 
                                ;maximum iterations reached       
    if niter ge maxiter then convergence = 1

endrep until convergence

theta_draw = reform(theta_draw, nchains * ndraws, npar)
logpost_draw = logpost_draw[*]

if not uniform then begin

    gmean_draw = reform(gmean_draw, nchains * ndraws, 2)
    Amat_draw = reform(Amat_draw, nchains * ndraws, 4)

endif

nsamp = n_elements(theta_draw[*,0])

if not uniform then $
  post = {pi:dblarr(ngauss), mu:dblarr(2, ngauss), covar:dblarr(2, 2, ngauss), n:0L, $
          gmean:[0d,0d], Amat:dblarr(2,2), logpost:0d, detprob:0d} $
else $
  post = {pi:dblarr(ngauss), mu:dblarr(2, ngauss), covar:dblarr(2, 2, ngauss), n:0L, $
          logpost:0d, detprob:0d}

post = replicate(post, nsamp)

if not silent then print, 'Getting posterior draws of N...'

gind = indgen(ngauss)

for i = 0L, nsamp - 1 do begin

    theta0 = reform(theta_draw[i,*])

    pi = theta0[pind]
    lmu = theta0[mlind]
    zmu = theta0[mzind]
    lvar = theta0[vlind]
    zvar = theta0[vzind]
    corr = theta0[cind] / sqrt(theta0[vlind] * theta0[vzind])

    post[i].pi = pi
    post[i].mu[0,*] = lmu
    post[i].mu[1,*] = zmu
    post[i].covar[0,0,*] = lvar
    post[i].covar[1,1,*] = zvar
    post[i].covar[0,1,*] = theta0[cind]
    post[i].covar[1,0,*] = post[i].covar[0,1,*]
    
    post[i].logpost = logpost_draw[i]

    if not uniform then begin

        post[i].gmean = reform(gmean_draw[i,*])
        post[i].amat = reform(Amat_draw[i,*], 2, 2)

    endif

    detprob = detprob_theta(theta0, selfunc, gind)
    detprob = total(pi * detprob) > 1d-10

    detprob = detprob * skyfrac

    nmis = negbin(nobs, detprob, seed=seed)
    post[i].n = nmis + nobs
    post[i].detprob = detprob

endfor

return
end
