; +
; NAME:
;  RMD_GA__DEFINE
;
; PURPOSE:
;  Class definition for minimizing a function using the
;  genetic algorithm search strategy that attempts to find the
;  global minimum using a genetic algorithm search strategy.  This is
;  a *simple* GA (SGA) in the sense that it uses roulette-wheel reproduction,
;  simple head-to-head and tail-to-tail swapping crossover and bit-flip
;  mutation.  Parameter encoding is binary.
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
;  or wide range of fitness values.  Note that you can override this mapping
;  and implement your own and specify it using the OBJECTIVE_FUNCTION
;  keyword.
;
;  For more detailed information on the simple genetic algorithm please
;  consult the text: "Genetic Algorithms in Search, Optimization, and Machine
;  Learning" by D.E.Goldberg (Addison-Wesley, 1989).
;
; CATEGORY:
;  OPTIMIZATION, OBJECTS
;
; CALLING SEQUENCE:
;  O = OBJ_NEW('RMD_GA',FTOL,FUNCTION_NAME = function_name,....)
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
;     below.  Keywords are passed into the function using FUNCTARGS
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
;                       (i.e. ITERPROC will not be called).  Default: 0B
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
;     An example of this is provided at the end of this code listing.
;
;  ITERARGS:            Structure containing exogenous information necessary
;                       to perform the intermediate reporting procedure
;                       ITERPROC. (Passed into this function via _EXTRA).
;  NCALLS:              Number of function calls over the course of the
;                       minimization procedure.
;  PCROSS:              Probability that two genes will be crossed (default: 0.9)
;  PMUTATE:             Probability that a single gene will have one of its
;                       chromosomes bit-flipped (default:0.01)
;  ITMAX:               Maximum number of iterations/generations to proceed
;                       in the evolution (default: 20)
;  NPOP:                Number of chromosomes in the population (default: 20)
;  STRETCH_FACTOR:      Amount by which to *stretch* the range of accepted fitness
;                       values.
;
; METHODS:
;  INIT:                (F) Object creation (ret 1)
;  CLEANUP:             (P) Object destruction
;  CREATE_POPULATION:   (F) Creates initial population (ret 1)
;  DECODE_POPULATION:   (F) Returns parameters in population
;  REPRODUCE:           (F) Roulette wheel reproduction (ret 1)
;  CROSSOVER:           (F) Head-to-head and tail-to-tail swapping (ret 1)
;  MUTATE:              (F) Bit flip (ret 1)
;  EVOLVE:              (F) Cycles through reproduction, crossover, and mutation (ret 1)
;  EVALUATE_FITNESS:    (F) Returns the fitness values of each member of the population
;  GET_PROPERTY:        (P) Accessor method
;  SET_PROPERTY:        (P) Accessor method
;
; INSTANCE DATA:
;  PCROSS:        Crossover probability (0.<PCROSS<1.)
;  PMUTATE:       Mutation probability (0.<PMUTATE<1.)
;  ITMAX:         Maximum number of generations for population to evolve
;  NPOP:          Number of members in the population
;  FUNC:          Objective function to be minimized
;  NPARMS:        Number of parametes in FUNC
;  KEEP_BEST:     Keyword when set results in a *greedy* GA
;                 (i.e. always keeps the best chromosome from
;                 one generation to the next)
;  GENE_LENGTH:   Length of the binary encoding for a single parameter
;  STRETCH_FACTOR:Amount by which to *stretch* the range of accepted fitness
;                 values.
;  INTERRUPT:     Byte value indicating whether the evolution should cease
;                 at the current generation.
;  AVE_FIT_PTR:   Heap variable that keeps a running tally at each generation
;                 of the changing average fitness in the population.
;  FIT_PTR:       Heap variable containing the fitness values in the current
;                 generation for all members of the population.
;  BEST_IN_POP_PTR:
;                 Best fit chromosome up until this generation.
;  QUIET:         Byte variable that instructs the algorithm whether or not to
;                 display intermediate results using ITERPROC.
;  ITERPROC:      Name (string) of the procedure to execute at each evolutionary
;                 step if QUIET is not set.
;  ITERARGS_PTR:  Heap variable containing the structure of values to be
;                 passed onto ITERPROC.
;  FUNCTARGS_PTR: Heap variable that points to a structure of information
;                 required by FUNC.  It is passed into FUNC via EXTRA keyword
;                 inheritance.
;  NCALLS:        Number of function evaluations during the course of the
;                 minimization.
;  POP_PTR:       Heap variable that contains a string array representing the
;                 binary-encoded population.
;  EXTRA_PTR:     Heap variable that points to a structure of exogenous information.
;
; REQUIREMENTS:
;  IDL 6.0 and higher
;
; REQUIRED PROGRAMS:
;  NONE
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
function rmd_ga_obj_fun,z,_Extra = extra
return,max(z)-z
end
; *********************************** ;
function rmd_ga_boltzmann,z,_Extra = extra
zmax = max(z,min = zmin)
return,exp(-z/(zmax-zmin))
end
; *********************************** ;
pro rmd_ga::cleanup
compile_opt idl2,hidden
ptr_free,self.prange_ptr,self.extra_ptr,self.pop_ptr
ptr_free,self.seed_ptr,self.fit_ptr,self.functargs_ptr
ptr_free,self.iterargs_ptr,self.ave_fit_ptr,self.best_fit_ptr
ptr_free,self.best_in_pop_ptr,self.dave_fit_ptr
ptr_free,self.feval_ptr,self.best_so_far_ptr
ptr_free,self.raw_fit_ptr
end
; *********************************** ;
pro rmd_ga_iterproc, func,                   $
                     p,                      $
                     iter,                   $
                     interrupt,              $
                     functargs = functargs,  $
                     oref = oga,             $
                     _Extra = iterargs

compile_opt hidden,idl2
; Default intermediate reporting procedure
print,'Iteration: '+strtrim(string(iter),2)
print,' *********************** '
for i = 0,n_elements(p)-1 do begin
   strout = 'p['+strtrim(string(i),2)+']='+strtrim(string(p[i]),2)
   print,strout
endfor
print
end
; *********************************** ;
function rmd_ga::decode_one,pop
compile_opt idl2,hidden
self->get_property,  gene_length=gene_length,      $
                     prange=prange
two_vec = 2^indgen(gene_length)
du = (prange[1,*]-prange[0,*])/(2^gene_length - 1.)
nparms = n_elements(du)
parms = dblarr(1,nparms)
i = 0 & pop = transpose(pop)
for j = 0,nparms-1 do begin
   gene_string = reform(pop[i,j*gene_length:(j+1)*gene_length-1])
   parms[i,j] = prange[0,j]+du[j]*(total((fix(gene_string))*reverse(two_vec)))
endfor
return,parms
end
; *********************************** ;
function rmd_ga::decode_population
compile_opt idl2,hidden
self->get_property,  pop=pop,                      $
                     gene_length=gene_length,      $
                     prange=prange,npop=npop
two_vec = 2^indgen(gene_length)
du = (prange[1,*]-prange[0,*])/(2^gene_length - 1.)
nparms = n_elements(du)
parms = dblarr(npop,nparms)
for i = 0,npop-1 do begin
;   parms[i,*] = self->decode_one(transpose(pop[i,*]))
   for j = 0,nparms-1 do begin
      gene_string = reform(pop[i,j*gene_length:(j+1)*gene_length-1])
      parms[i,j] = prange[0,j]+du[j]*(total((fix(gene_string))*reverse(two_vec)))
   endfor
endfor
return,parms
end
; *********************************** ;
pro rmd_ga::set_property,  prange = prange,                 $
                           pmutate = pmutate,               $
                           pcross = pcross,                 $
                           function_name = function_name,   $
                           gene_length = gene_length,       $
                           pop = pop,                       $
                           itmax = itmax,                   $
                           ftol = ftol,                     $
                           npop = npop

compile_opt idl2,hidden
if n_elements(ftol) ne 0 then self.ftol = ftol
if n_elements(itmax) ne 0 then self.itmax = itmax
if n_elements(pop) ne 0 then *self.pop_ptr = pop
if n_elements(prange) ne 0 then *self.prange_ptr = prange
if n_elements(pmutate) ne 0 then self.pmutate = pmutate
if n_elements(function_name) ne 0 then self.func = function_name
if n_elements(pcross) ne 0 then self.pcross = pcross
if n_elements(npop) ne 0 then self.npop = npop
if n_elements(gene_length) ne 0 then self.gene_length = gene_length
end
; *********************************** ;
pro rmd_ga::get_property,  prange = prange,                 $
                           pmutate = pmutate,               $
                           pcross = pcross,                 $
                           function_name = function_name,   $
                           gene_length = gene_length,       $
                           pop = pop,                       $
                           nparms = nparms,                 $
                           parms = parms,                   $
                           best_parms = best_parms,         $
                           itmax = itmax,                   $
                           ncalls = ncalls,                 $
                           fitness = fitness,               $
                           feval = feval,                   $
                           ave_fitness = ave_fitness,       $
                           dave_fitness = dave_fitness,     $
                           best_fitness = best_fitness,     $
                           raw_fitness = raw_fitness,       $
                           ftol = ftol,                     $
                           npop = npop
compile_opt idl2,hidden
ftol = self.ftol
itmax = self.itmax
ncalls = self.ncalls
prange = *self.prange_ptr
pmutate = self.pmutate
pcross = self.pcross
npop = self.npop
gene_length = self.gene_length
function_name = self.func
if arg_present(raw_fitness) then raw_fitness = *self.raw_fit_ptr
if arg_present(nparms) then nparms = (size(prange))[2]
if arg_present(parms) then begin
   if n_elements(*self.pop_ptr) eq 0 then parms = 0B
   parms = self->decode_population()
endif
if arg_present(feval) then feval = *self.feval_ptr
if arg_present(pop) and (n_elements(*self.pop_ptr)) ne 0    $
   then pop = *self.pop_ptr else pop = 0B
if arg_present(fitness) then begin
   if n_elements(*self.fit_ptr) eq 0 then fitness = 0B else $
      fitness = *self.fit_ptr
endif
if arg_present(best_parms) then $;begin
;   self->get_property,fitness = fitness,parms = parms
;   max_fit = max(fitness,imax)
;   best_parms = reform(parms[imax,*])
   best_parms = self->decode_one(reform(*self.best_so_far_ptr))
;endif
if arg_present(ave_fitness) then begin
   if n_elements(self.ave_fit_ptr) ne 0 then $
      ave_fitness = *self.ave_fit_ptr else $
      ave_fitness = 0B
endif
if arg_present(dave_fitness) then begin
   if n_elements(self.dave_fit_ptr) ne 0 then $
      dave_fitness = *self.dave_fit_ptr else $
      dave_fitness = 0B
endif
if arg_present(best_fitness) then best_fitness = *self.best_fit_ptr
end
; *********************************** ;
function rmd_ga::evaluate_fitness
compile_opt idl2,hidden
parms = self->decode_population()
raw_fitness = fltarr(self.npop)
for i = 0,self.npop-1 do $
   raw_fitness[i] = call_function(self.func,parms[i,*],_Extra = *self.functargs_ptr)
*self.raw_fit_ptr = raw_fitness
*self.fit_ptr = call_function(self.objective_function,   $
   raw_fitness/self.stretch_factor,iter = self.iter)
self.ncalls = self.ncalls+self.npop
return,*self.fit_ptr
end
; *********************************** ;
function rmd_ga::reproduce
compile_opt idl2,hidden
; Spin a *weighted* roulette wheel to determine who
; gets propagated into the next generation.
regions = [0.,total(*self.fit_ptr/total(*self.fit_ptr),/cumulative)]
s = *self.seed_ptr
spin = randomu(s,self.npop)
*self.seed_ptr = s
old_pop = *self.pop_ptr
*self.pop_ptr = old_pop[value_locate(regions,spin),*]
if self.keep_best then (*self.pop_ptr)[0,*] = *self.best_so_far_ptr;*self.best_in_pop_ptr
return,1
end
; *********************************** ;
function rmd_ga::crossover
compile_opt idl2,hidden
;list_1 = indgen(self.npop)
s = *self.seed_ptr
list_2 = sort(randomu(s,self.npop))
chromosome_length = self.gene_length * self.nparms
num_div = chromosome_length - 1
crossPoints = num_div*randomu(s,self.npop)
rcross = randomu(s,self.npop)

*self.seed_ptr = s
; Shuffle the population
old_pop = (*self.pop_ptr)[list_2,*]

for i = 0,self.npop-1,2 do begin
   start = crossPoints[i]
   mate_1 = reform(old_pop[i,*])
   mate_2 = reform(old_pop[i+1,*])
   child_1 = [mate_1[0:start],mate_2[start+1:chromosome_length-1]]
   child_2 = [mate_2[0:start],mate_1[start+1:chromosome_length-1]]
   if rcross[i] le self.pcross then begin
      (*self.pop_ptr)[i,*] = child_1
      (*self.pop_ptr)[i+1,*] = child_2
   endif else begin
      (*self.pop_ptr)[i,*] = mate_1
      (*self.pop_ptr)[i+1,*] = mate_2
   endelse
endfor

return,1
end
; *********************************** ;
function rmd_ga::mutate
compile_opt idl2,hidden
s = *self.seed_ptr
toss = randomu(s,self.npop)
mutate_index = where(toss le self.pmutate,count_mutate)
if count_mutate gt 0 then begin
   pop = *self.pop_ptr
   total_length = self.gene_length * self.nparms
   rand_pos = fix(total_length*randomu(s,count_mutate))
   for i = 0,count_mutate-1 do   $
      (*self.pop_ptr)[mutate_index[i],rand_pos[i]] = string(fix(~byte(fix('1'))))
endif
*self.seed_ptr = s
if self.keep_best then (*self.pop_ptr)[0,*] = *self.best_so_far_ptr;*self.best_in_pop_ptr
return,1
end
; *********************************** ;
function rmd_ga::evolve
compile_opt idl2,hidden
*self.ave_fit_ptr = fltarr(self.itmax)
*self.dave_fit_ptr = fltarr(self.itmax)
*self.best_fit_ptr = fltarr(self.itmax)
*self.feval_ptr = fltarr(self.itmax)
   *self.best_so_far_ptr = (*self.pop_ptr)[0,*]
i = (self.iter = 0L)
self.interrupt = 0B
termination = 0B
while (i lt self.itmax) and ~termination do begin
   self.iter = i
   ; Evaluate the fitness of the current population
   ret = self->evaluate_fitness()
   ; Get the raw fitness values
   raw_fitmin = min(*self.raw_fit_ptr,imin)
   *self.best_in_pop_ptr = (*self.pop_ptr)[imin,*]

   (*self.ave_fit_ptr)[i] = (moment((*self.raw_fit_ptr)))[0]
   (*self.dave_fit_ptr)[i] = sqrt((moment((*self.raw_fit_ptr)))[1])

   if self.iter eq 0 then begin
      self.best_raw_fitness = raw_fitmin
      *self.best_so_far_ptr = (*self.pop_ptr)[imin,*]
   endif else begin
      if min(*self.raw_fit_ptr) lt self.best_raw_fitness then begin
         self.best_raw_fitness = raw_fitmin
         *self.best_so_far_ptr = (*self.pop_ptr)[imin,*]
      endif
   endelse
   (*self.best_fit_ptr)[self.iter] = self.best_raw_fitness
   (*self.feval_ptr)[self.iter] = self.best_raw_fitness

   ; Now execute the three main GA operations
   ret = self->reproduce()
   ret = self->crossover()
   ret = self->mutate()

   termination = ((abs(raw_fitmin-(*self.ave_fit_ptr)[i])/abs(raw_fitmin)) le self.ftol)

   if ~self.quiet then begin
      interrupt = self.interrupt
      self->get_property,best_parms = p
      call_procedure, self.iterproc,self.func,p,self.iter,  $
                      interrupt,                            $
                      oref = self,                          $
                      functargs = *self.functargs_ptr,      $
                      _Extra = *self.iterargs_ptr
      self.interrupt = interrupt
   endif
   i++
   if self.interrupt then i = self.itmax
endwhile

return,1
end
; *********************************** ;
function rmd_ga::create_population
compile_opt idl2,hidden
; Create the initial population
self->get_property,  npop=npop,prange=prange,      $
                     gene_length=gene_length

total_length = (size(prange))[2] * gene_length
if n_elements(*self.seed_ptr) ne 0 then s = *self.seed_ptr
*self.pop_ptr = strtrim(string(round(randomu(s,npop,total_length))),2)
*self.seed_ptr = s
return,1
end
; *********************************** ;
function rmd_ga::init,  ftol,                            $
                        function_name = function_name,   $
                        pcross = pcross,                 $
                        itmax = itmax,                   $
                        pmutate = pmutate,               $
                        npop = npop,                     $
                        gene_length = gene_length,       $
                        prange = prange,                 $
                        keep_best = keep_best,           $
                        functargs = functargs,           $
                        boltzmann = boltzmann,           $
                        objective_function = ofun,       $
                        iterproc = iterproc,             $
                        iterargs = iterargs,             $
                        stretch_factor = stretch_factor, $
                        quiet = quiet,                   $
                        _Extra = extra
compile_opt idl2,hidden
if n_params() eq 0 then ftol = 0.1
self.ftol = ftol
if n_elements(ofun) eq 0 then ofun = 'rmd_ga_obj_fun'
if n_elements(boltzmann) ne 0 then ofun = 'rmd_ga_boltzmann'
self.objective_function = ofun
if n_elements(itmax) eq 0 then itmax = 20
self.itmax = itmax
if n_elements(stretch_factor) eq 0 then stretch_factor = 1.0
self.stretch_factor = stretch_factor
if n_elements(iterproc) eq 0 then iterproc = 'rmd_ga_iterproc'
if (strupcase(iterproc) eq 'RMD_GA_ITERPROC') and $
   (n_elements(iterargs) eq 0) then begin
   iterargs = {iterstop:1}
   interrupt = 0B
endif
self.iterproc = iterproc
self.iterargs_ptr = ptr_new(iterargs)
if n_elements(quiet) eq 0 then quiet = 1B
self.quiet = quiet
if n_elements(keep_best) eq 0 then keep_best = 0B
self.keep_best = keep_best
if n_elements(function_name) eq 0 then return,0
self.func = function_name
if n_elements(pcross) eq 0 then pcross = 0.9
self.pcross = 0.0 > pcross < 1.0
if n_elements(pmutate) eq 0 then pmutate = 0.1
self.pmutate = 0.0 > pmutate < 1.0
if n_elements(npop) eq 0 then npop = 20
self.npop = 2 > npop
; Make sure that the number in the population is even
if self.npop mod 2 ne 0 then self.npop++
if n_elements(gene_length) eq 0 then gene_length = 5
self.gene_length = 2 > gene_length
if n_elements(prange) eq 0 then return,0
prange_size = size(prange)
if prange_size[1] ne 2 then return,0
self.nparms = (size(prange))[2]
self.prange_ptr = ptr_new(prange,/no_copy)
self.pop_ptr = ptr_new(/allocate_heap)
self.fit_ptr = ptr_new(/allocate_heap)
self.dave_fit_ptr = ptr_new(/allocate_heap)
if n_elements(extra) ne 0 then $
   self.extra_ptr = ptr_new(extra) else $
   self.extra_ptr = ptr_new(/allocate_heap)
self.seed_ptr = ptr_new(/allocate_heap)
self.functargs_ptr = ptr_new(/allocate_heap)
if n_elements(functargs) eq 0 then functargs = {iterstop:0}
*self.functargs_ptr = functargs
self.ave_fit_ptr = ptr_new(/allocate_heap)
self.ncalls = 0L
self.best_in_pop_ptr = ptr_new(/allocate_heap)
self.best_fit_ptr = ptr_new(/allocate_heap)
self.feval_ptr = ptr_new(/allocate_heap)
*self.feval_ptr = fltarr(self.itmax)
self.best_so_far_ptr = ptr_new(/allocate_heap)
self.raw_fit_ptr = ptr_new(/allocate_heap)
self.best_raw_fitness = 0.0
self.iter = 0L
return,1
end
; *********************************** ;
pro rmd_ga__define
compile_opt idl2,hidden
define = {  rmd_ga,                    $
            pcross:0.0,                $
            pmutate:0.0,               $
            npop:0L,                   $
            nparms:0L,                 $
            keep_best:0B,              $
            gene_length:0L,            $
            stretch_factor:0.,         $
            interrupt:0B,              $
            ave_fit_ptr:ptr_new(),     $
            ftol:0.,                   $
            iter:0L,                   $
            fit_ptr:ptr_new(),         $
            functargs_ptr:ptr_new(),   $
            best_in_pop_ptr:ptr_new(), $
            best_raw_fitness:0.0,      $
            feval_ptr:ptr_new(),       $
            raw_fit_ptr:ptr_new(),     $
            best_fit_ptr:ptr_new(),    $
            dave_fit_ptr:ptr_new(),    $
            iterargs_ptr:ptr_new(),    $
            best_so_far_ptr:ptr_new(), $
            objective_function:'',     $
            iterproc:'',               $
            itmax:0L,                  $
            ncalls:0L,                 $
            quiet:0B,                  $
            func:'',                   $  ; function to be maximized (the fitness function)
            pop_ptr:ptr_new(),         $
            extra_ptr:ptr_new(),       $
            seed_ptr:ptr_new(),        $
            prange_ptr:ptr_new()       $
         }
end
; *********************************** ;