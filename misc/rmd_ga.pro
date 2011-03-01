; +
; NAME:
;  RMD_GA
;
; PURPOSE:
;  Function wrapper for RMD_GA object class that attempts to find the
;  global minimum using a genetic algorithm search strategy.  This is
;  a *simple* GA (SGA) in the sense that it uses roulette-wheel reproduction,
;  simple head-to-head and tail-to-tail swapping crossover and bit-flip
;  mutation.  Parameter encoding is binary.  Termination criteria are
;  either the maximum number of iterations ITMAX has been reached or
;  the average fitness of the whole population is within FTOL of the best
;  fit in the population.
;
;  In general the SGA attempts to find a global maximum in the fitness
;  function.  However we wish to find a global minimum of some objective
;  function.  In order to perform the mapping from the objective function
;  to be minimized, FUNC, to fitness, FITNESS, we implement the following
;  relationship:
;
;     FITNESS = MAX(FUNC/STRETCH_FACTOR)-FUNC/STRETCH_FACTOR
;
;  where STRETCH_FACTOR is a number that can be tuned to accept a narrow
;  or wide range of fitness values. Note that you can override this mapping
;  and implement your own, specifying it using the OBJECTIVE_FUNCTION
;  keyword.
;
;  For more detailed information on the simple genetic algorithm please
;  consult the text: "Genetic Algorithms in Search, Optimization, and Machine
;  Learning" by D.E.Goldberg (Addison-Wesley, 1989).
;
; CATEGORY:
;  OPTIMIZATION
;
; CALLING SEQUENCE:
;  Result = RMD_GA(FUNCTION_VALUE = function_value,...)
;
; RETURN VALUE:
;  Returns the best parameters found that minimize the function.
;
; PARAMETERS:
;  FTOL: The stopping criteria in terms of the change in the average
;        fitness relative to the maximum fitness value.  (default: 0.1)
;
; KEYWORDS:
;  FUNCTION_NAME:       name (string) of the function to be minimized
;
;     The requirements on the function to be minimized, FUNCTION_NAME,
;     are that it be written as a FUNCTION and it must accept a parameter
;     vector P and keywords passed in via the EXTRA keyword inheritance
;     mechanism.  An example of a function to be minimized is shown
;     below.  Keywords are passed into the function using FUNCTARG
;     described below.
;
;         function test_rmd_fun,p,_Extra = extra
;         y = 4.0+sin(4.*p)*exp(-0.5*((p-2.)/4.0)^2)
;         return,y
;         end
;
;  OBJECTIVE_FUNCTION:  User-supplied function (optional) that replaces
;                       RMD_GA_OBJ_FUN as described above.  The requirements
;                       on the function are that it convert a maximization
;                       problem into a minimization problem.  Note that keywords
;                       should be allowed via the _Extra keyword inheritance
;                       mechanism.
;  BOLTZMANN:           If set then the objective function is given by the
;                       Boltzmann distribution.  Here the fitness F is given
;                       by F = exp(-E/(max(E)-min(E))) where E is the function
;                       evaluation.  Note that setting this keyword overrides
;                       a user-specified OBJECTIVE_FUNCTION.
;  FUNCTION_VALUE:      Function evaluated at the best parameter set found
;  PRANGE:              A 2 x N vector defining the lower and upper bounds
;                       for the N parameters to be refined in the minimization.
;  QUIET:               If set, no intermediate output will be provided
;                       (i.e. ITERPROC will not be called).  Default: 1B
;  FUNCTARGS:           Structure containing exogenous information necessary
;                       to evaluate the function to be minimized.  This is
;                       a convenience to the user and is not required.
;  ITERPROC:            Name (string) of the procedure to be called at
;                       each evolutionary step.  An example of
;                       such a procedure is TEST_GA_ITERPROC, shown below.
;                       (Passed into this function via _EXTRA).
;
;     The requirements on the intermediate reporting procedure are
;     that it be written as a PROCEDURE and it must accept the
;     following arguments: FUNC, P (the current best parameter set),
;     ITER (the current iteration, and INTERRUPT (a byte variable equal 1B or 0B
;     indicating if the algorithm should cease at the next update--note that
;     interrupt must be a variable since it is expected to be passed by reference).
;     Optional keywords are FUNCTARGS, ITERARGS, and OREF.  OREF is provided
;     so that you can extract information on the GA object (its state)
;     at the current iteration.  This can be useful if you want to get
;     all of the current members of the population and display them, for instance.
;     An example is provided at the end of this code listing.
;
;  ITERARGS:            Structure containing exogenous information necessary
;                       to perform the intermediate reporting procedure
;                       ITERPROC. (Passed into this function via _EXTRA).
;  NCALLS:              Number of function calls over the course of the
;                       minimization procedure.
;  PCROSS:              Probability that two genes will be crossed (default: 0.9)
;  PMUTATE:             Probability that a single gene will have one of its
;                       chromosomes bit-flipped (default:0.01)
;  ITMAX:               Maximum number of generations to proceed in the evolution (default:10)
;  NPOP:                Number of chromosomes in the population (default: 20)
;  GENE_LENGTH:         Length in bits of each member of the population.  The more
;                       bits the better the numerical precision of the result (default: 5)
;  STRETCH_FACTOR:      Amount by which to *stretch* the range of accepted fitness
;                       values.
;
; EXAMPLE USAGE:
;  See the example below TEST_RMD_GA, TEST_GA_ITERPROC, TEST_FUNC_1
;
; REQUIREMENTS:
;  IDL 6.0 and higher
;
; REQUIRED PROGRAMS:
;  RMD_GA__DEFINE
;
; COMMON BLOCKS:
;  NONE
;
; SIDE EFFECTS:
;  NONE
;
; AUTHOR:
;  Robert Dimeo
;  NIST Center for Neutron Research
;  National Institute of Standards and Technology
;  100 Bureau Drive-Stop 8562
;  Gaithersburg, MD 20899
;
; DISCLAIMER
;   This software is provided as is without any warranty whatsoever.
;   Permission to use, copy, modify, and distribute modified or
;   unmodified copies is granted, provided this disclaimer
;   is included unchanged.
;
; ACKNOWLEDGMENT:
;  This work is part of the DAVE development effort at the NCNR and
;  was supported in part by the National Science Foundation
;  under Agreement No. DMR-0086210.
;
; MODIFICATION HISTORY:
;  Written by RMD 11/05/04
;
; -
; *********************************** ;
function rmd_ga,  ftol,                                     $
                  function_value = function_value,          $
                  function_name = function_name,            $
                  prange = prange,                          $
                  ncalls = ncalls,                          $
                  quiet = quiet,                            $
                  pcross = pcross,                          $
                  pmutate = pmutate,                        $
                  boltzmann = boltzmann,                    $
                  itmax = itmax,                            $
                  objective_function = ofun,                $
                  gene_length = gene_length,                $
                  stretch_factor = stretch_factor,          $
                  npop = npop,                              $
                  functargs = ft,                           $
                  _Extra = extra

compile_opt idl2,hidden
if n_params() eq 0 then ftol = 1e-1
if n_elements(boltzmann) ne 0 then boltzmann = 1B else boltzmann = 0B
if n_elements(ofun) eq 0 then ofun = 'rmd_ga_obj_fun'
keep_best = 1B ; always keep the best solution in the population
               ; from one generation to the next.
if n_elements(stretch_factor) eq 0 then stretch_factor = 1.0
if n_elements(quiet) eq 0 then quiet = 1B
if n_elements(function_name) eq 0 then return,0
if n_elements(prange) eq 0 then return,0
if n_elements(ft) eq 0 then functargs = {iterstop:0} else functargs = ft
if n_elements(pcross) eq 0 then pcross = 0.9
pcross = 0.0 > pcross < 1.0
if n_elements(pmutate) eq 0 then pmutate = 0.1
pmutate = 0.0 > pmutate < 1.0
if n_elements(itmax) eq 0 then itmax = 10
itmax = 1 > itmax
if n_elements(npop) eq 0 then npop = 20
npop = 2 > npop
; Make sure that the number in the population is even
if npop mod 2 ne 0 then npop++
if n_elements(gene_length) eq 0 then gene_length = 5
gene_length = 2 > gene_length

prange_size = size(prange)
;if prange_size[0] ne 2 then return,0
if prange_size[1] ne 2 then return,0
o = obj_new('rmd_ga',   ftol,                                     $
                        function_name = function_name,            $
                        boltzmann = boltzmann,                    $
                        pcross = pcross,                          $
                        pmutate = pmutate,                        $
                        itmax = itmax,                            $
                        npop = npop,                              $
                        gene_length = gene_length,                $
                        prange = prange,                          $
                        keep_best = keep_best,                    $
                        stretch_factor = stretch_factor,          $
                        objective_function = ofun,                $
                        functargs = functargs,                    $
                        quiet = quiet,                            $
                        _Extra = extra                            )

ret = o->create_population()
fitness = o->evaluate_fitness()
ret = o->evolve()

o->get_property,best_parms = best_parms,ncalls = ncalls
;print,'Best parms: '
;print,best_parms
function_value = call_function(function_name,best_parms,_Extra = functargs)
;print,'FEVAL: ',function_value
obj_destroy,o
return,best_parms
end
; *********************************** ;
; Begin example usage.  The following
; three routines are necessary for the
; example: TEST_RMD_GA, TEST_GA_ITERPROC,
; and TEST_FUNC_1.
; *********************************** ;
function test_func_1,p,_EXTRA = extra
x = p[0] & y = p[1]
z = 3.*(1.-x)^2 * exp(-x^2 - (y+1.)^2) - $
    10.*(x/5 - x^3 - y^5)*exp(-x^2-y^2) - $
    (1./3)*exp(-(x+1)^2 - y^2)
return,z
end
; *********************************** ;
pro test_ga_iterproc,func,                   $
                     p,                      $
                     iter,                   $
                     interrupt,              $
                     functargs = functargs,  $
                     oref = oga,             $
                     _Extra = extra
compile_opt hidden,idl2
oga->get_property,ave_fitness = ave_fitness
x = 1+indgen(iter+1)
y = ave_fitness[0:iter]
tvlct,r,g,b,/get
rnew = reverse(r) & gnew = reverse(g) & bnew = reverse(b)
tvlct,rnew,gnew,bnew
wset,extra.winpix
plot,[x],[y],psym = -4,title = 'Function evaluation',xtitle = 'Generation', $
   ytitle = '<F(p)>'
wset,extra.winvis
device,copy = [0,0,!d.x_size,!d.y_size,0,0,extra.winpix]
tvlct,r,g,b
end
; *********************************** ;
function my_obj_fun,z,_Extra = extra
return,(extra.iter+1.) * (max(z)-z)
end
; *********************************** ;
pro test_rmd_ga
; Example implementation of the SGA to minimize a function of
; two variables.  Here the range of the parameters is defined
; as follows:
;  -9.<P[0]<9. and -9.<P[1]<9.
; The actual global maximum (within the specified range) is
; (x,y) = (0.228,-1.626)

prange = [[-9.0,9.0],[-9.,9.]]
ofun = 'my_obj_fun'
func = 'test_func_1'
quiet = 0B

if ~quiet then begin
   xsize = 400 & ysize = 400
   window,0,xsize = xsize,ysize = ysize
   winvis = 0
   window,/free,/pixmap,xsize = xsize,ysize = ysize
   winpix = !d.window
   iterargs = {winvis:winvis,winpix:winpix}
   iterproc = 'test_ga_iterproc'
endif
ftol = 1.e-2
p = rmd_ga(       ftol,                               $
                  function_value = function_value,    $
                  function_name = func,               $
                  prange = prange,                    $
;                  /boltzmann,                         $
                  ncalls = ncalls,                    $
                  quiet = quiet,                      $
                  objective_function = ofun,          $
                  pcross = 0.95,                      $
                  gene_length = 25,                   $
                  pmutate = 0.01,                     $
                  stretch_factor = 1.,                $
                  itmax = 100,                        $
                  iterproc = iterproc,                $
                  iterargs = iterargs,                $
                  npop = 250)

if ~quiet then wdelete,winpix
this_format = '(f15.3)'
print,'Best parameters: '+strtrim(string(p[0],format = this_format),2)+ $
      ','+strtrim(string(p[1],format = this_format),2)
print,'Function value: '+strtrim(string(function_value,format = this_format),2)
end