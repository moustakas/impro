;+
; NAME:
;  mcc_test
;
; PURPOSE:
;  Test the MCMC code.  Currently tests a constant and a power law.
;  The default mode compares the recovered likelihood distribution to the known
;  analytic likelihood.  If ntest= is sent, then the distribution of means for
;  a series of chains is compared to the expected.
;
; CALLING SEQUENCE:
;  test, type=, ntest=
;
; OPTIONAL INPUTS:
;  type: A type of function to test.  Possible valures are 'constant'
;     and 'power'.  Default is both.
;  ntest: Run this many chains and compare the resulting distribution of mean
;         parameters against the expected distribution.  Currently only
;         supported for "constant" type.
;
; MODIFICATION HISTORY:
;  Created: 2006-10-14: Erin Sheldon, NYU
;-


function mcmc_test_const_like, p
  common mcmc_test_block, x, y, ivar, psigma, npars
  chi2 = total( (y-p[0])^2*ivar )
  return, -0.5*chi2
;  return, -0.5*(chi2-n_elements(y))
end 

function mcmc_test_const_step, seed, pars
  common mcmc_test_block, x, y, ivar, psigma, npars
  return, pars + psigma*randomn(seed)
end 
function mcmc_test_const_truepars, sigma=sigma
  sigma = 1.0
  return, 10.0
end 
pro mcmc_test_const_setup
  common mcmc_test_block, x, y, ivar, psigma, npars
  ny = 10
  pars = mcmc_test_const_truepars(sigma=sigma)
  y = replicate(pars, ny)
  yerr = replicate(sigma, ny)

  ivar = 1.0/yerr^2
  y[*] = y[*] + yerr*randomn(seed, ny)
  psigma = yerr[0]
  npars = 1
end 


function mcmc_test_pow_like, p

  common mcmc_test_block, x, y, ivar, psigma, npars

  ;; is ivar a matrix?
  n_dim = size(ivar, /n_dim)

  model = p[0]*x^p[1]

  diff = y-model
  if n_dim eq 1 then begin 
      chi2 = total(diff^2*ivar )
  endif else begin 
      chi2 = transpose(diff)#(reform(ivar##diff))
  endelse 

  return, -0.5*chi2

end 

function mcmc_test_pow_step, seed, pars
  common mcmc_test_block, x, y, ivar, psigma, npars

  newpars = pars + psigma*randomn(seed, npars)
  return, newpars
end 

function mcmc_test_pow_truepars
  return, [0.85,-0.4]
end 

pro mcmc_test_pow_setup
  common mcmc_test_block, x, y, ivar, psigma, npars

  nx = 50
  minx = 0.025
  maxx = 10.0
  p = mcmc_test_pow_truepars()
  npars = n_elements(p)

  logbin_spacing, minx, maxx, nx, bmin, bmax, dim=2, bin_mean=x
  
  y = p[0]*x^p[1]
  yerr = y*0.2
  y = y + yerr*randomn(seed, nx)

  ivar = 1d/yerr^2

  ;; sigma of random numbers to generate
  psigma = [0.1,0.1]

;  help, x, y, yerr, ivar
;  colprint, x, y, yerr, ivar

end 

pro mcmc_test, type=typearr, ntest=ntest, outst=outst, sigfrac=sigfrac

  common mcmc_test_block, x, y, ivar, psigma, npars

  ;; Default, test both
  if n_elements(typearr) eq 0 then typearr=['constant','power']

  mcmc=obj_new('mcmc')

  for pi=0L, n_elements(typearr)-1 do begin 
      type = typearr[pi]

      print
      print,'Testing type: '+ntostr(type)

      tm = systime(1)

      case type of
          'constant': begin 

              nstep = 10000
              burnin = 100

              ;; We are testing distribution of means against expected error
              if n_elements(ntest) ne 0 then begin 

                  meanarray = fltarr(ntest)
                  sdevarray = meanarray

                  for i=0L, ntest-1 do begin 
                      print,backspace(14)+'i = ',i+1,format='(A,I10,$)'
                      
                      mcmc_test_const_setup
                      parguess = 5.0
                      pars = mcmc->run('mcmc_test_const_step', 'mcmc_test_const_like', $
                                       nstep, parguess, seed=seed, /log)
                      
                      mom = moment(pars[burnin-1:nstep-1], sdev=sd)
                      meanarray[i] = mom[0]
                      sdevarray[i] = sd
                  endfor 
                  
                  outst= {means:meanarray, sdevs:sdevarray}
                  
                  mom = moment(meanarray, sdev=sd)
                  if n_elements(sigfrac) eq 0 then sigfrac=4.0
                  plothist, meanarray, bin=sd/sigfrac, xhist, yhist, /noplot
                  ;; normalize histogram
                  yhist = yhist

                  pplot, xhist, yhist/qgauss(yhist, xhist, 100), psym=10, $
                    aspect=!gratio
                  pars = mcmc_test_const_truepars(sigma=sigma)
                  meanerr = sigma/sqrt(n_elements(y))
                  oplot, xhist, gaussprob(xhist, pars, meanerr), color=!red

                  legend, 'Nchain = '+ntostr(ntest),/left,box=0,charsize=1
                  legend, ['MCMC Means','Expected Distribution'], $
                    line=0, color=[!p.color, !red], $
                    /right, box=0, charsize=1


              endif else begin ;; ntest != 0
                  ;; Testing if the likelihood function has the right shape
                  nstep = 1000000

                  mcmc_test_const_setup
                  parguess = 5.0

                  pars = mcmc->run('mcmc_test_const_step', 'mcmc_test_const_like', $
                                   nstep, parguess, seed=seed, printstep=10000,/log)

                  mom = moment(pars[burnin-1:nstep-1], sdev=sd)
                  mn = mom[0]

                  if n_elements(sigfrac) eq 0 then sigfrac=10.0
                  plothist, pars[burnin-1:nstep-1], bin=sd/sigfrac, xhist, yhist, $
                    /noplot
                  yhist = yhist/qgauss(yhist, xhist, 100)

                  tpars=mcmc_test_const_truepars(sigma=sigma)
                  mom = moment(pars[burnin-1:nstep-1],sdev=sd)
                  xrange = mom[0] + sd*[-5,5]
                  pplot, xhist, yhist, psym=10, xrange=xrange, xstyle=3, $
                    aspect=!gratio, xtitle='Value', ytitle='Likelihood'

                  oplot,[tpars,tpars], [0,1.e10], line=0, color=!DarkGreen

                  nh=n_elements(xhist)
                  plike = fltarr(nh)
                  for pari=0L, nh-1 do begin 
                      plike[pari] = exp(mcmc_test_const_like(xhist[pari]))
                  endfor 
                  plike = plike/qgauss(plike, xhist, 100)
                  oplot, xhist, plike, color=!red

                  legend,'Fit Type: Constant',/left,box=0,charsize=1
                  legend, ['MCMC','Analytic likelihood','Truth'], line=0, $
                    color=[!p.color, !red,!DarkGreen], $
                    /right, box=0, charsize=1

              endelse 
          end 
          'power': begin 

              pold = !p.multi
              !p.multi = [0,0,2]

              nstep = 1000000
              burnin = 1000
              
              mcmc_test_pow_setup
              parguess = [0.5, -1.0]
              
              pars = mcmc->run('mcmc_test_pow_step', 'mcmc_test_pow_like', $
                               nstep, parguess, printstep=10000,/log)
              
              mom = moment(pars[0,burnin-1:nstep-1], sdev=sd)
              mean0 = mom[0]
              sd0 = sd
              mom = moment(pars[1,burnin-1:nstep-1], sdev=sd)
              mean1 = mom[0]
              sd1 = sd

              normrange = mean0 + 5*sd0*[-1,1]
              powrange = mean1 + 5*sd1*[-1,1]

              yerr=sqrt(1d/ivar)
              pow_chisq_conf_gen, x, y, yerr, powrange, normrange, $
                100, 100, powvals=powvals, normvals=normvals, $
                pow_like=powlike, norm_like=normlike, /nodisplay


              if n_elements(sigfrac) eq 0 then sigfrac=10.0
              tpars=mcmc_test_pow_truepars()              
              
              plothist, pars[0,burnin-1:nstep-1], bin=sd/sigfrac, xhist, yhist,$
                min=normrange[0], max=normrange[1], $
                /noplot
              plot, xhist, yhist/qgauss(yhist,xhist,100), psym=10, $
                xtitle='Norm', ytitle='Likelihood'
              oplot, normvals, normlike/qgauss(normlike, normvals, 100), color=!red
              oplot,[tpars[0],tpars[0]], [0,1.e10], line=0, color=!DarkGreen


              legend,'Fit Type: Power Law',/left,box=0,charsize=1
              legend, ['MCMC','Analytic Likelihood','Truth'], line=0, $
                color=[!p.color, !red,!DarkGreen], $
                /right, box=0, charsize=1


              plothist, pars[1,burnin-1:nstep-1], bin=sd/sigfrac, xhist, yhist,$
                min=powrange[0], max=powrange[1], $
                /noplot
              plot, xhist, yhist/qgauss(yhist,xhist,100), psym=10, $
                xtitle='Index', ytitle='Likelihood'
              oplot, powvals, powlike/qgauss(powlike, powvals, 100), color=!red
              oplot,[tpars[1],tpars[1]], [0,1.e10], line=0, color=!DarkGreen

              !p.multi=pold
              
          end 
          else: message,'Unknown type: '+ntostr(type)
      endcase
      
      ptime,systime(1)-tm

  endfor 
  
end 
