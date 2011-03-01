; +
; NAME:
;  RMD_GA_ISING_DEFINE
;
; PURPOSE:
;  Demonstration of the use of a Genetic Algorithm for
;  finding the minimum energy configuration of a random
;  Ising spin system in one-dimension (an Ising glass).
;  This is a simple modification of the RMD_GA__DEFINE
;  class definition for illustration purposes.
;
; CATEGORY:
;  OPTIMIZATION, OBJECTS
;
; AUTHOR:
;  Robert Dimeo
;  NIST Center for Neutron Research
;  National Institute of Standards and Technology
;  100 Bureau Drive-Stop 8562
;  Gaithersburg, MD 20899
;
; CALLING SEQUENCE:
;  RMD_GA_ISING,NSPINS
;
; REQUIRED PROGRAMS:
;  TVIMAGE
;
; COMMON BLOCKS:
;  NONE
;
; SIDE EFFECTS:
;  NONE
;
; MODIFICATION HISTORY:
;  Written by RMD 11/12/04
;
; -
; *********************************** ;
function rmd_ga_ising_obj_fun,z,_Extra = extra
return,max(z)-z
end
; *********************************** ;
function rmd_ga_boltzmann,z,_Extra = extra
zmax = max(z,min = zmin)
dz = zmax-zmin
return,exp(-sqrt(extra.iter+1)*((z-zmin)/(zmax-zmin)))
end
; *********************************** ;
pro rmd_ga_ising::cleanup
compile_opt idl2,hidden
ptr_free,self.extra_ptr,self.pop_ptr
ptr_free,self.seed_ptr,self.fit_ptr,self.functargs_ptr
ptr_free,self.iterargs_ptr,self.ave_fit_ptr,self.best_fit_ptr
ptr_free,self.best_in_pop_ptr,self.dave_fit_ptr
ptr_free,self.feval_ptr,self.best_so_far_ptr
ptr_free,self.raw_fit_ptr
end
; *********************************** ;
pro rmd_ga_iterproc, func,                   $
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
function rmd_ga_ising::decode_one,pop
compile_opt idl2,hidden
return,-1+2*round(fix(strtrim(string(round(fix(pop))),2)))
end
; *********************************** ;
function rmd_ga_ising::decode_population
compile_opt idl2,hidden
return,-1+2*round(fix(string(round(fix(*self.pop_ptr)))))
end
; *********************************** ;
pro rmd_ga_ising::set_property,                             $
                           pmutate = pmutate,               $
                           pcross = pcross,                 $
                           function_name = function_name,   $
                           pop = pop,                       $
                           itmax = itmax,                   $
                           ftol = ftol,                     $
                           npop = npop

compile_opt idl2,hidden
if n_elements(ftol) ne 0 then self.ftol = ftol
if n_elements(itmax) ne 0 then self.itmax = itmax
if n_elements(pop) ne 0 then *self.pop_ptr = pop
if n_elements(pmutate) ne 0 then self.pmutate = pmutate
if n_elements(function_name) ne 0 then self.func = function_name
if n_elements(pcross) ne 0 then self.pcross = pcross
if n_elements(npop) ne 0 then self.npop = npop
end
; *********************************** ;
pro rmd_ga_ising::get_property,                             $
                           pmutate = pmutate,               $
                           pcross = pcross,                 $
                           function_name = function_name,   $
                           pop = pop,                       $
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
pmutate = self.pmutate
pcross = self.pcross
npop = self.npop

function_name = self.func
if arg_present(raw_fitness) then raw_fitness = *self.raw_fit_ptr
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
if arg_present(best_parms) then $
   best_parms = self->decode_one(reform(*self.best_so_far_ptr))
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
function rmd_ga_ising::evaluate_fitness
compile_opt idl2,hidden
parms = self->decode_population()
raw_fitness = fltarr(self.npop)
;for i = 0,self.npop-1 do $
;   raw_fitness[i] = call_function(self.func,parms[i,*],_Extra = *self.functargs_ptr)
raw_fitness = call_function(self.func,parms,_Extra = *self.functargs_ptr)
*self.raw_fit_ptr = raw_fitness
*self.fit_ptr = call_function(self.objective_function,   $
   raw_fitness/self.stretch_factor,iter = self.iter)
self.ncalls = self.ncalls+self.npop
return,*self.fit_ptr
end
; *********************************** ;
function rmd_ga_ising::reproduce
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
function rmd_ga_ising::crossover
compile_opt idl2,hidden
list_1 = indgen(self.npop)
s = *self.seed_ptr
list_2 = list_1[sort(randomu(s,self.npop))]
chromosome_length = self.nspins
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
function rmd_ga_ising::mutate
compile_opt idl2,hidden
s = *self.seed_ptr
toss = randomu(s,self.npop)
mutate_index = where(toss le self.pmutate,count_mutate)
if count_mutate gt 0 then begin
   pop = *self.pop_ptr
   total_length = self.nspins
   rand_pos = fix(total_length*randomu(s,count_mutate))
   for i = 0,count_mutate-1 do   $
      (*self.pop_ptr)[mutate_index[i],rand_pos[i]] = string(fix(~byte(fix('1'))))
endif
*self.seed_ptr = s
if self.keep_best then (*self.pop_ptr)[0,*] = *self.best_so_far_ptr;*self.best_in_pop_ptr
return,1
end
; *********************************** ;
function rmd_ga_ising::evolve
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
function rmd_ga_ising::create_population
compile_opt idl2,hidden
; Create the initial population
if n_elements(*self.seed_ptr) ne 0 then s = *self.seed_ptr
*self.pop_ptr = strtrim(string(round(randomu(s,self.npop,self.nspins))),2)
*self.seed_ptr = s
return,1
end
; *********************************** ;
function rmd_ga_ising::go
if self.micro_ga then begin
   ; Get the current parameters out
   self->get_property,  npop = npop,         $
                        pcross = pcross,     $
                        pmutate = pmutate,   $
                        itmax = itmax
   ; Set the micro-GA parameters
   self.pcross = 1.0
   self.npop = 8
   self.itmax = 10
   self.pmutate = 0.0
   ntimes = 5
   ; Evolve the micro-GA
   for i = 0,ntimes-1 do begin
      ret = self->create_population()
      ; Keep the best from the micro-GA
      if i gt 0 then $
         (*self.pop_ptr)[0,*] = *self.best_so_far_ptr
      ret = self->evaluate_fitness()
      ret = self->evolve()
   endfor
   self.pcross = pcross
   self.npop = npop
   self.itmax = itmax
   self.pmutate = pmutate
   ret = self->create_population()
   ; Keep the best from the micro-GA
   (*self.pop_ptr)[0,*] = *self.best_so_far_ptr
   ; Evolve the new population
   fitness = self->evaluate_fitness()
   ret = self->evolve()
endif else begin
   ret = self->create_population()
   fitness = self->evaluate_fitness()
   ret = self->evolve()
endelse
return,1
end
; *********************************** ;
function rmd_ga_ising::init,                             $
                        ftol,                            $
                        nspins = nspins,                 $
                        function_name = function_name,   $
                        pcross = pcross,                 $
                        itmax = itmax,                   $
                        pmutate = pmutate,               $
                        npop = npop,                     $
                        keep_best = keep_best,           $
                        functargs = functargs,           $
                        objective_function = ofun,       $
                        iterproc = iterproc,             $
                        iterargs = iterargs,             $
                        stretch_factor = stretch_factor, $
                        boltzmann = boltzmann,           $
                        micro_ga = micro_ga,             $
                        quiet = quiet,                   $
                        _Extra = extra
compile_opt idl2,hidden
if n_elements(micro_ga) eq 0 then micro_ga = 0B
self.micro_ga = micro_ga
if n_params() eq 0 then ftol = 0.1
self.ftol = ftol
if n_elements(nspins) eq 0 then nspins = 10
self.nspins = nspins
if n_elements(ofun) eq 0 then ofun = 'rmd_ga_ising_obj_fun'
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
pro rmd_ga_ising__define
compile_opt idl2,hidden
define = {  rmd_ga_ising,              $
            pcross:0.0,                $
            pmutate:0.0,               $
            npop:0L,                   $
            keep_best:0B,              $
            micro_ga:0B,               $
            nspins:0L,                 $
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
            seed_ptr:ptr_new()         $
         }
end
; *********************************** ;
