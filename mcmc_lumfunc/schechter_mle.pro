;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; FIT A SCHECHTER (GAMMA) LUMINOSITY FUNCTION VIA MAXIMUM LIKELIHOOD.
;
; AUTHOR : BRANDON C. KELLY, STEWARD OBSERVATORY, SEPT. 2007
;
; INPUTS :
;
;    LUM - NDATA-ELEMENT VECTOR CONTAINING THE LOGARITHM OF THE
;          LUMINOSITIES FOR THE NDATA GALAXIES IN ONE'S SAMPLE.
;
;    LLIM - DETECTION LIMIT, IN LOG-LUMINOSITY.
;
; OPTIONAL INPUTS :
;
;    POISSON - USE THE (INCORRECT) POISSON LIKELHOOD FUNCTION INSTEAD
;              OF THE BINOMIAL LIKELIHOOD FUNCTION.
;
; OUTPUTS :
;
;    THETA - 3-ELEMENT VECTOR CONTAINING THE MAXIMUM LIKELIHOOD
;            ESTIMATES (MLE). THETA[0] CONTAINS THE LUMINOSITY
;            FUNCTION NORMALIZATION, I.E., THETA[0] IS THE TOTAL
;            NUMBER OF GALAXIES IN ONE'S SURVEY. THETA[1] IS THE
;            SCHECHTER FUNCTION SHAPE PARAMETER, AND THETA[2] IS
;            LSTAR.
;
; OPTIONAL OUTPUTS :
;
;    COVAR - THE ASYMPTOTIC COVARIANCE MATRIX OF THETA.
;
; ROUTINES CALLED :
;
;    HESSIAN
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function schechter_loglik, theta

common lumdata, lum, llim, skyfrac

if total(finite(theta)) lt 3 then begin
    loglik = -1d-33
    return, loglik
endif

nlum = n_elements(lum)

logntot = theta[0]
ntot = exp(logntot)
k = exp(theta[1])
loglstar = theta[2]
lstar = exp(loglstar)

lognorm = lngamma(ntot + 1) - lngamma(ntot - nlum + 1)

detprob = (1d - igamma(k, llim / lstar, /double, itmax=1d6)) * skyfrac > 1d-33

det_loglik = (k - 1) * (alog(lum) - loglstar) - loglstar - lngamma(k) - lum / lstar
det_loglik = total( det_loglik )

loglik = lognorm + (ntot - nlum) * alog( 1d - detprob ) + det_loglik

return, loglik
end

function poisson_loglik, theta

common lumdata, lum, llim, skyfrac

nlum = n_elements(lum)

if total(finite(theta)) lt 3 then begin
    loglik = -1d-33
    return, loglik
endif

logntot = theta[0]
ntot = exp(logntot)
k = exp(theta[1])
loglstar = theta[2]
lstar = exp(loglstar)

det_loglik = (k - 1) * (alog(lum) - loglstar) - loglstar - lngamma(k) - lum / lstar
det_loglik = total( det_loglik )

detprob = ( 1d - igamma(k, llim / lstar, /double, itmax=1d6) ) * skyfrac > 1d-33

loglik = nlum * logntot + det_loglik - ntot * detprob

return, loglik
end

pro schechter_mle, lum, llim, skyfrac, theta, covar=covar, poisson=poisson, silent=silent

if n_params() lt 4 then begin
    print, 'Syntax- schechter_mle, lum, llim, skyfrac, theta, covar=covar, /poisson, /silent'
    return
endif

common lumdata, lumc, llimc, skyfracc

lumc = lum
llimc = llim
skyfracc = skyfrac

if not keyword_set(poisson) then poisson = 0
if not keyword_set(silent) then silent = 0
ndata = n_elements(lum)

naccept = 0
nreject = 0

niter = 0

k0 = 1.0
lstar0 = 1d44
ntot0 = ndata / (1d - igamma(k0, llim / lstar0, /double))

theta = [alog(ntot0), alog(k0), alog(lstar0)]
                                ;place constraints on parameters
tbnd = transpose([[alog(ndata),1d30], [-1d6,1d6], [alog(1d35),alog(1d50)]])
gbnd = transpose([0,0])

if poisson then likname = 'poisson_loglik' else likname = 'schechter_loglik'

constrained_min, theta, tbnd, gbnd, 0, likname, inform, /maximize

theta[0] = theta[0] / alog(10)
theta[1] = theta[1] / alog(10)
theta[2] = theta[2] / alog(10)

infomat = -1d * hessian( theta * alog(10), likname )

covar = la_invert( infomat, /double )

return
end
