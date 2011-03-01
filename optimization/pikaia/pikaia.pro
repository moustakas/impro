;********************************************************************
	function urand
; 	Common block to make iseed visible to urand_init (and to save
; 	it between calls)
	common share1,iseed
	urand = randomu( iseed )
	return,urand
	end

;******************************************************************


;**********************************************************************
;                       REPRODUCTION MODULE
;**********************************************************************
;
;     SELECT:   Parent selection by roulette wheel algorithm
;               called by: PIKAIA
;
;     RNKPOP:   Ranks initial population
;               called by: PIKAIA, NEWPOP
;
;     GENREP:   Inserts offspring into population, for full
;               generational replacement
;               called by: PIKAIA
;
;     STDREP:   Inserts offspring into population, for steady-state
;               reproduction
;               called by: PIKAIA
;               calls:     FF
;
;     NEWPOP:   Replaces old generation with new generation
;               called by: PIKAIA
;               calls:     FF, RNKPOP
;
;**********************************************************************
	pro select,np,jfit,fdif,idad,nparents
	compile_opt IDL2
	common share1,iseed
;======================================================================
;     Selects nparents parents from the population, using roulette wheel
;     algorithm with the relative fitnesses of the phenotypes as
;     the "hit" probabilities [see Davis 1991, chap. 1].
;======================================================================
;
;     implicit none
;
;     Input:
;     integer        np, jfit(np), nparents
;     real           fdif
;
;     Output:
;     integer        idad
;
	np1 = np+1
	dice = rebin(randomu(iseed,nparents)*np*np1,nparents,np+1)
	rtfit = make_array(1,np,value=np1) + fdif*(np1-2*(jfit+1))
	rtfit = [[fltarr(nparents)],[rebin(total(rtfit,/cumulative),nparents,np)]]
	diff = (rtfit - dice) > 0
	w = where( diff - shift(diff,0,1) eq diff and diff gt 0)
	cols = w mod nparents
	rows = w/nparents

	idad = intarr(nparents)
	idad[cols] = rows-1

	return
	end
;******************************************************************


;******************************************************************
	pro rnkpop,n,arrin,indx,rank
	compile_opt IDL2
;======================================================================
;     Calls external sort routine to produce key index and rank order
;     of input array arrin (which is not altered).
;======================================================================
;
; 	implicit none

; 	Input
; 	integer    n
; 	real       arrin(n)

; 	Output
; 	integer    indx(n),rank(n)
	rank = lonarr(n)

; 	Compute the key index
	indx = sort(arrin[0:n-1])
; 	and the rank order
	rank[indx] = n - indgen(n) - 1
	return
	end
;******************************************************************


;******************************************************************
	pro stdrep,ndim,n,np,irep,ielite,phs1,phs2,oldph,fitns,ifit,jfit,nnew,_extra=extra
	compile_opt IDL2
	common func,fname
;======================================================================
;     steady-state reproduction: insert offspring pair into population
;     only if they are fit enough (replace-random if irep=1 or
;     replace-worst if irep=2).
;======================================================================
;
;     implicit none
;
;     Input:
;     integer        ndim, n, np, irep, ielite
;     real           ff, ph(ndim,2)
;     external       ff
;
;     Input/Output:
;     real           oldph(ndim,np), fitns(np)
;     integer        ifit(np), jfit(np)
;
;     Output:
;     integer        nnew
;
	nnew = 0

	npeeps = n_elements(ifit)
	n1 = n_elements(phs1[0,*])
	n2 = n_elements(phs2[0,*])
	phnew = fltarr(n,npeeps + n1 + n2)

	; calculate fitness for first and second offspring sets
	if n_elements(extra) gt 0 then begin
		fit1 = call_function(fname,n,phs1,_extra=extra)
		fit2 = call_function(fname,n,phs2,_extra=extra)
	endif else begin
		fit1 = call_function(fname,n,phs1)
		fit2 = call_function(fname,n,phs2)
	endelse

	; put fitness measurements together
	totfit = [fitns,fit1,fit2]

	; sort, and sort the sort (for reverse-lookup)
	sortind = sort(totfit)
	sortsort = sort(sortind)

	; get the indexes for everything
	i1 = sortsort[0:npeeps-1]
	i2 = sortsort[npeeps:npeeps+n1-1]
	i3 = sortsort[npeeps+n1:npeeps+n1+n2-1]

	; fill the new array in with the old data, from the proper locations
	phnew[*,i1] = oldph
	phnew[*,i2] = phs1
	phnew[*,i3] = phs2

	; count the number of new people
	w = where( i1 gt n1+n2, c )
	nnew = npeeps-c

	; store the result
	oldph = phnew[*,n1+n2:*]
	fitns = totfit[sortind[n1+n2:*]]
	;oldph = phnew[*,0:npeeps-1]
	;fitns = totfit[sortind[0:npeeps-1]]

; 	Compute the key index
	ifit = sort(fitns)
	jfit = intarr(npeeps)
; 	and the rank order
	jfit[ifit] = npeeps - indgen(npeeps) - 1

return
end
;**********************************************************************


;**********************************************************************
;                         GENETICS MODULE
;**********************************************************************
;
;     ENCODE:    encodes phenotype into genotype
;                called by: PIKAIA
;
;     DECODE:    decodes genotype into phenotype
;                called by: PIKAIA
;
;     CROSS:     Breeds two offspring from two parents
;                called by: PIKAIA
;
;     MUTATE:    Introduces random mutation in a genotype
;                called by: PIKAIA
;
;     ADJMUT:    Implements variable mutation rate
;                called by: PIKAIA
;
;**********************************************************************


;**********************************************************************
	pro encode,n,nd,ph,gn
	compile_opt IDL2
;======================================================================
;     encode phenotype parameters into integer genotype
;     ph(k) are x,y coordinates [ 0 < x,y < 1 ]
;======================================================================
;
;
;     Inputs:
;     integer   n, nd
;     real      ph(n)
;
;     Output:
;     integer   gn(n*nd)
;
;     Local:
;     integer   ip, i, j, ii
;     real      z
;

	z=10.^nd
	npeeps = n_elements(ph[0,*])
	ndims = n*nd

; 	let the dimensional juggling begin!
	tens = rebin(reform(rebin(10L^(lindgen(nd)),nd,n),ndims),ndims,npeeps)
	temp = reform(rebin(reform(ph,1,n,npeeps),nd,n,npeeps),nd*n,npeeps)*tens
	gn = floor( 10.0*( temp - floor(temp) ) )

	return
	end
;******************************************************************


;**********************************************************************
	pro decode,n,nd,gn,ph
	compile_opt IDL2
;======================================================================
;     decode genotype into phenotype parameters
;     ph(k) are x,y coordinates [ 0 < x,y < 1 ]
;======================================================================
;
;
;     Inputs:
;     integer   n, nd, gn(n*nd)
;
;     Output:
;     real      ph(n)

	z=10.^(-nd)
	ndims = n*nd
	npeeps = n_elements(gn[0,*])
	tens = rebin(reform(rebin(10L^(make_array(nd,/integer,value=nd-1)-indgen(nd)),nd,n),ndims),ndims,npeeps)
	ph = rebin(gn*tens,n,npeeps)*nd*z

	return
	end
;**********************************************************************


;**********************************************************************
	pro cross,n,nd,pcross,gns1,gns2
	compile_opt IDL2
	common share1,iseed
;======================================================================
;     breeds two parent chromosomes into two offspring chromosomes
;     breeding occurs through crossover starting at position ispl
;======================================================================
;
;     Inputs:
;     integer        n, nd
;     real           pcross
;
;     Input/Output:
;     integer        gn1(n*nd), gn2(n*nd)


;     Use crossover probability to decide whether a crossover occurs
	rands = randomu(iseed,n_elements(gns1[0,*]))
	w = where( rands lt pcross, c )
	isps = long(randomu(iseed,c)*n*nd)

	for i=0,c-1 do begin
		t = gns2[isps[i]:*,w[i]]
		gns2[isps[i]:*,w[i]] = gns1[isps[i]:*,w[i]]
		gns1[isps[i]:*,w[i]] = t
	endfor

	return
	end
;******************************************************************


;******************************************************************
	pro mutate,n,nd,pmut,gn
	compile_opt IDL2
	common share1,iseed
;======================================================================
;     Mutations occur at rate pmut at all gene loci
;======================================================================
;
;
;     Input:
;     integer        n, nd
;     real           pmut
;
;     Input/Output:
;     integer        gn(n*nd)
;
;     Subject each locus to mutation at the rate pmut

;       for i=0,n*nd-1 do begin
;          if urand() lt pmut then begin
;             gn(i)=long(urand()*10.)
;          endif
;       endfor

	w = where(randomu(iseed,size(gn,/dimensions)) lt pmut,c)
	if c gt 0 then gn[w] = long(randomu(iseed,c)*10)

      return
      end
;******************************************************************


;**********************************************************************
	pro adjmut,np,fitns,ifit,pmutmn,pmutmx,pmut
	compile_opt IDL2
;======================================================================
;     dynamical adjustment of mutation rate; criterion is relative
;     difference in absolute fitnesses of best and median individuals
;======================================================================
;
;     implicit none
;
;     Input:
;     integer        np, ifit(np)
;     real           fitns(np), pmutmn, pmutmx
;
;     Input/Output:
;     real           pmut
;
;     Local:
;     real           rdif, rdiflo, rdifhi, delta
;     parameter      (rdiflo=0.05, rdifhi=0.25, delta=1.5)
           
      rdiflo=0.05
      rdifhi=0.25  
      delta=1.5
 
      rdif=abs(fitns[ifit[np-1]]-fitns[ifit[np/2-1]])/(fitns[ifit[np-1]]+fitns[ifit[np/2-1]])
      if rdif le rdiflo then pmut=pmutmx<pmut*delta
      if rdif ge rdifhi then pmut=pmutmn>pmut/delta
 
      return
      end
;******************************************************************


;******************************************************************
	pro setctl,ctrl,n,np,ngen,nd,pcross,pmutmn,pmutmx,pmut,imut,fdif,irep,ielite,ivrb,status
	compile_opt IDL2

;     Set control variables and flags from input and defaults
;
;
;     Input
;     integer  n
;
;     Input/Output
;     real     ctrl(12)
;
;     Output
;     integer  np, ngen, nd, imut, irep, ielite, ivrb, status
;     real     pcross, pmutmn, pmutmx, pmut, fdif
;
;     Local
;     integer  i
;     real     DFAULT(12)
;     save     DFAULT
;     data     DFAULT /100,500,5,.85,2,.005,.0005,.25,1,1,1,0/

;    DFAULT = [100,20,4,.85,2,.005,.0005,.25,1,1,1,1]
      DFAULT = [100,200,4,.85,2,.005,.0005,.25,1,3,1,0]
;
      if n_elements(ctrl) eq 0 then ctrl = dfault
      test = where(ctrl lt 0, c)
      if c gt 0 then ctrl[test] = DFAULT[test]
; 
      np = long(ctrl[0])
      ngen = long(ctrl[1])
      nd = long(ctrl[2])
      pcross = ctrl[3]
      imut = long(ctrl[4])
      pmut = ctrl[5]
      pmutmn = ctrl[6]
      pmutmx = ctrl[7]
      fdif = ctrl[8]
      irep = long(ctrl[9])
      ielite = long(ctrl[10])
      ivrb = long(ctrl[11])
      status = 0
;
;     Print a header
;
      if (ivrb gt 0) then begin
 
         f2 = '(/1x,60("*"),/'
         f2 = f2+'" *",13x,"PIKAIA Genetic Algorithm Report ",13x,"*",/'
	 f2 = f2 + '1x,60("*"),//,"   Number of Generations evolving: ",i4,/'
         f2 = f2 + '"       Individuals per generation: ",i4,/'
         f2 = f2 +'"    Number of Chromosome segments: ",i4,/'
	 f2 = f2 + '"    Length of Chromosome segments: ",i4,/'
	 f2 = f2 + '"            Crossover probability: ",f9.4,/'
	 f2 = f2 + '"            Initial mutation rate: ",f9.4,/'
	 f2 = f2 + '"            Minimum mutation rate: ",f9.4,/'
	 f2 = f2 + '"            Maximum mutation rate: ",f9.4,/'
	 f2 = f2 + '"    Relative fitness differential: ",f9.4)'
	print,format=f2,ngen,np,n,nd,pcross,pmut,pmutmn,pmutmx,fdif
;
        f3='("                    Mutation Mode: ",A)'
         if (imut eq 1) then print,format=f3, 'Constant'
         if (imut eq 2) then print,format=f3, 'Variable'
;
        f4 = '("                Reproduction Plan: ",A)'
         if (irep eq 1) then print,format=f4, 'Full generational replacement'
         if (irep eq 2) then print,format=f4, 'Steady-state-replace-random'
         if (irep eq 3) then print,format=f4, 'Steady-state-replace-worst'
;
      endif
 
;     Check some control values
      f10  = '" ERROR: illegal value for imut (ctrl(5))"'
      if (imut ne 1 and imut ne 2) then begin
         print,f10 
         status = 5
      endif
 
      f11 = '" ERROR: illegal value for fdif (ctrl(9))"'
      if (fdif gt 1) then begin
         print,f11
         status = 9
      endif
 
      f12 = '" ERROR: illegal value for irep (ctrl(10))"'
      if (irep ne 1  and  irep ne 2  and  irep ne 3) then begin
	 print,f12
         status = 10
      endif
 
      f13 = ' ERROR: illegal value for pcross (ctrl(4))'
      if (pcross gt 1.0  or  pcross lt 0 ) then begin
	 print,f13
         status = 4
      endif
 
      f14 = ' ERROR: illegal value for pcross (ctrl(4))'
      if (ielite ne 0  and  ielite ne 1) then begin
 	 print,f14
         status = 11
      endif
 
      f15 = '" WARNING: dangerously high value for pmut (ctrl(6));",/" (Should enforce elitism with ctrl(11)=1.)"'
      if (irep eq 1  and  imut eq 1  and  pmut gt 0.5  and ielite eq 0) then begin
	 print,f15
      endif
 
      f16 = '" WARNING: dangerously high value for pmutmx (ctrl(8));",/" (Should enforce elitism with ctrl(11)=1.)"'
      if (irep eq 1 and imut eq 2 and pmutmx gt 0.5 and ielite eq 0) then begin
	 print,f16
      endif
 
      f17 = ' WARNING: dangerously low value of fdif (ctrl(9))'
      if (fdif lt 0.33) then begin
	 print,f17
      endif
 
      return
      end
;*************************************************************************


;********************************************************************
	pro report,ivrb,ndim,n,np,nd,oldph,fitns,ifit,pmut,ig,nnew
	compile_opt IDL2
;********************************************************************

	common datstore,bestft,pmutpv
 
;     Write generation report to standard output
;
;     Input:
;     integer ifit(np),ivrb,ndim,n,np,nd,ig,nnew
;     real oldph(ndim,np),fitns(np),pmut
;
;     Output: none
;
;     Local
;     real bestft,pmutpv
;     save bestft,pmutpv
;     integer ndpwr,k
;     logical rpt
;     data bestft,pmutpv /0,0/
;
      rpt=0.
 
      if (pmut ne pmutpv) then begin
         pmutpv=pmut
         rpt=1.
      endif
 
      if (fitns(ifit(np-1)) ne bestft) then begin
         bestft=fitns(ifit(np-1))
         rpt=1.
      endif
 
      if (rpt ne 0.  or  ivrb ge 2) then begin
 
;        Power of 10 to make integer genotypes for display

         ndpwr = round(10.^nd)
 
         print,format='(/i6,i6,f10.6,4f10.6)', ig+1,nnew,pmut,fitns(ifit(np-1)), fitns(ifit(np-2)), fitns(ifit((np+1)/2 -1))
	 for k=0,n-1 do begin
            print,format='(22x,3i10)',round(ndpwr*oldph(k,ifit(np-1))),round(ndpwr*oldph(k,ifit(np-2))),round(ndpwr*oldph(k,ifit((np+1)/2 -1)))
	endfor
 
      endif
      end
;*************************************************************************


;*************************************************************************
pro pikaia,n,ctrl,x,f,status,fname=functionname,_extra=extra
compile_opt IDL2
;
common share1,iseed
common datstore,bestft,pmutpv
common func,fname

if n_elements(functionname) eq 0 then fname = 'func' else fname = functionname
;
;     Optimization (maximization) of user-supplied "fitness" function
;     fname over n-dimensional parameter space  x  using a basic genetic
;     algorithm method.
;
;
;     Paul Charbonneau & Barry Knapp
;     High Altitude Observatory
;     National Center for Atmospheric Research
;     Boulder CO 80307-3000
;     <paulchar@hao.ucar.edu>
;     <knapp@hao.ucar.edu>
;
;     Version of 1995 April 13
;
;	(IDL version 2007 July - comments largely kept intact)
;      Vectorized and optomized for IDL by Conor Mancone, Astronomy Department, University of Florida
;
;     Genetic algorithms are heuristic search techniques that
;     incorporate in a computational setting, the biological notion
;     of evolution by means of natural selection.  This subroutine
;     implements the three basic operations of selection, crossover,
;     and mutation, operating on "genotypes" encoded as strings.
;
;     References:
;
;        Charbonneau, Paul.  "Genetic Algorithms in Astronomy and
;           Astrophysics."  Astrophysical J. (Supplement), vol 101,
;           in press (December 1995).
;
;        Goldberg, David E.  Genetic Algorithms in Search, Optimization,
;           & Machine Learning.  Addison-Wesley, 1989.
;
;        Davis, Lawrence, ed.  Handbook of Genetic Algorithms.
;           Van Nostrand Reinhold, 1991.
;
;
;     implicit none
;
;     Input:
;     integer   n
;     real      ff
;     external  ff
;
;      o Integer  n  is the parameter space dimension, i.e., the number
;        of adjustable parameters.  (Also the number of chromosomes.)
;
;      o Function  ff  is a user-supplied scalar function of n vari-
;        ables, which must have the calling sequence f = ff(n,x), where
;        x is a real parameter array of length n.  This function must
;        be written so as to bound all parameters to the interval [0,1];
;        that is, the user must determine a priori bounds for the para-
;        meter space, and ff must use these bounds to perform the appro-
;        priate scalings to recover true parameter values in the
;        a priori ranges.
;
;        By convention, ff should return higher values for more optimal
;        parameter values (i.e., individuals which are more "fit").
;        For example, in fitting a function through data points, ff
;        could return the inverse of chi**2.
;
;        In most cases initialization code will have to be written
;        (either in a driver or in a separate subroutine) which loads
;        in data values and communicates with ff via one or more labeled
;        common blocks.  An example exercise driver and fitness function
;        are provided in the accompanying file, xpikaia.f.
;
;
;      Input/Output:
;      real ctrl(12)
;
;      o Array  ctrl  is an array of control flags and parameters, to
;        control the genetic behavior of the algorithm, and also printed
;        output.  A default value will be used for any control variable
;        which is supplied with a value less than zero.  On exit, ctrl
;        contains the actual values used as control variables.  The
;        elements of ctrl and their defaults are:
;
;           ctrl( 1) - number of individuals in a population (default
;                      is 100)
;           ctrl( 2) - number of generations over which solution is
;                      to evolve (default is 500)
;           ctrl( 3) - number of significant digits (i.e., number of
;                      genes) retained in chromosomal encoding (default
;                      is 6)  (Note: This number is limited by the
;                      machine floating point precision.  Most 32-bit
;                      floating point representations have only 6 full
;                      digits of precision.  To achieve greater preci-
;                      sion this routine could be converted to double
;                      precision, but note that this would also require
;                      a double precision random number generator, which
;                      likely would not have more than 9 digits of
;                      precision if it used 4-byte integers internally.)
;           ctrl( 4) - crossover probability; must be  <= 1.0 (default
;                      is 0.85)
;           ctrl( 5) - mutation mode; 1/2=steady/variable (default is 2)
;           ctrl( 6) - initial mutation rate; should be small (default
;                      is 0.005) (Note: the mutation rate is the proba-
;                      bility that any one gene locus will mutate in
;                      any one generation.)
;           ctrl( 7) - minimum mutation rate; must be >= 0.0 (default
;                      is 0.0005)
;           ctrl( 8) - maximum mutation rate; must be <= 1.0 (default
;                      is 0.25)
;           ctrl( 9) - relative fitness differential; range from 0
;                      (none) to 1 (maximum).  (default is 1.)
;           ctrl(10) - reproduction plan; 1/2/3=Full generational
;                      replacement/Steady-state-replace-random/Steady-
;                      state-replace-worst (default is 3)
;           ctrl(11) - elitism flag; 0/1=off/on (default is 0)
;                      (Applies only to reproduction plans 1 and 2)
;           ctrl(12) - printed output 0/1/2=None/Minimal/Verbose
;                      (default is 0)
;
;
;     Output:
;     real      x(n), f

	x = fltarr(n)

;     integer   status
;
;      o Array  x(1:n)  is the "fittest" (optimal) solution found,
;         i.e., the solution which maximizes fitness function ff
;
;      o Scalar  f  is the value of the fitness function at x
;
;      o Integer  status  is an indicator of the success or failure
;         of the call to pikaia (0=success; non-zero=failure)
;
;
;     Constants
;     integer   NMAX, PMAX, DMAX
;

	NMAX = 32
	PMAX = 512 
	DMAX = 8
;
;      o NMAX is the maximum number of adjustable parameters
;        (n <= NMAX)
;
;      o PMAX is the maximum population (ctrl(1) <= PMAX)
;
;      o DMAX is the maximum number of Genes (digits) per Chromosome
;        segement (parameter) (ctrl(3) <= DMAX)
;
;
;     Local variables
;     integer        np, nd, ngen, imut, irep, ielite, ivrb, k, ip, ig,
;    +               ip1, ip2, new, newtot
;     real           pcross, pmut, pmutmn, pmutmx, fdif
;
;     real           ph(NMAX,2), oldph(NMAX,PMAX), newph(NMAX,PMAX)

;
;
;     User-supplied uniform random number generator
;     real           urand
;     external       urand
;
;     Function urand should not take any arguments.  If the user wishes
;     to be able to initialize urand, so that the same sequence of
;     random numbers can be repeated, this capability could be imple-
;     mented with a separate subroutine, and called from the user's
;     driver program.  An example urand function (and initialization
;     subroutine) which uses the function ran2 from Press, et al,
;     Numerical Recipes, 2nd ed., Cambridge Univ Press, 1992, is
;     provided in the accompanying file, xpikaia.f.
;
;
;  rather than do a data block, I will zero these out now,
;  and keep them in a common block shared between pro pikaia
;  and pro report
;
	pmutpv = 0.
	bestft = 0.
;
;     Set control variables from input and defaults


	setctl,ctrl,n,np,ngen,nd,pcross,pmutmn,pmutmx,pmut,imut, fdif,irep,ielite,ivrb,status

	if (status  ne  0) then begin
		print," Control vector (ctrl) argument(s) invalid"
		return
	endif
 
;     Make sure locally-dimensioned arrays are big enough

	if (n gt NMAX  or  np gt PMAX  or  nd gt DMAX) then begin
		print," Number of parameters, population, or genes too large"
		status = -1
		return
	endif


	ph = fltarr(n,2)
	oldph = randomu(iseed,n,np)
	newph = fltarr(n,np)


; 	integer        gn1(NMAX*DMAX), gn2(NMAX*DMAX)
; 	integer        ifit(PMAX), jfit(PMAX)

	gn1 = lonarr(n*nd)
	gn2 = lonarr(n*nd)
	ifit = lonarr(np)
	jfit = lonarr(np)

; 	real           fitns(PMAX)

; 	compute fitness of initial population
	if n_elements(extra) gt 0 then $
		fitns = call_function(fname,n,oldph,_extra=extra) $
	else $
		fitns = call_function(fname,n,oldph)
; 	wascii,transpose(fitns),'fitns_new'
 
; 	Rank initial population by fitness order
	rnkpop,np,fitns,ifit,jfit

; 	Main Generation Loop

	for ig=0,ngen-1 do begin
		newtot=0

; 		select first and second parents
		select,np,jfit,fdif,ips1,np/2
		select,np,jfit,fdif,ips2,np/2

; 		2. encode parent phenotypes
		encode,n,nd,oldph[*,ips1],gns1
		encode,n,nd,oldph[*,ips2],gns2

; 		3. breed
		cross,n,nd,pcross,gns1,gns2
		mutate,n,nd,pmut,gns1
		mutate,n,nd,pmut,gns2

; 		4. decode offspring genotypes
		decode,n,nd,gns1,phs1
		decode,n,nd,gns2,phs2

; 		5. insert into population
; 		I've only implemented one replacement technique, it's the only one I need
		stdrep,NMAX,n,np,irep,ielite,phs1,phs2,oldph,fitns,ifit,jfit,new,_extra=extra
		newtot = newtot+new

; 		if running full generational replacement: swap populations
; 		if (irep eq 1) then newpop,ielite,NMAX,n,np,oldph,newph,ifit,jfit,fitns,newtot

; 		adjust mutation rate?
		if (imut eq 2) then adjmut,np,fitns,ifit,pmutmn,pmutmx,pmut

; 		print generation report to standard output?

		if (ivrb gt 0) then report,ivrb,NMAX,n,np,nd,oldph,fitns,ifit,pmut,ig,newtot

; 	End of Main Generation Loop
	endfor

; 	Return best phenotype and its fitness

	x = oldph[*,ifit[np-1]]
	f = fitns[ifit[np-1]]

end