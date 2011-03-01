; +
; NAME:
;  RMD_PSO
;
; PURPOSE:
;  Function that attempts to find the global minimum of a function
;  using a particle swarm optimization (PSO) strategy.
;
;  In this algorithm the individual agents (swarm members) "fly" through
;  parameter space seeking the global function minimum.  At each time
;  step in the search every agent determines the value of the function
;  at its location in parameter space and compares this value with its
;  own previous best value as well as the overall best parameter found
;  among all of the agents.  The velocity of each agent is determined by
;  an attraction between its current position, the best solution for that
;  particular agent, and the best solution for all agents found.
;  The velocity of each agent is updated at each time step according to
;  the equation:
;
;     v[t+dt] = b*v[t] + dv[t] where
;     dv[t] = w1*r*(pbest[t] - p[t]) + w2*r*(gbest[t] - p[t])
;
;  w1 and w2 are weights (also known as learning factors) with
;  typical values ranging between 0 and 4, b is a velocity damping factor
;  (also called an inertia weight), r is a uniform random deviate,
;  p is the current parameter set, pbest is the best parameter set found
;  for this particular agent, and gbest is the best parameter set found
;  for all particular agents.  The inertia weight decreases from 1 at each
;  time step so that after k updates, the value is b^k.  One points to
;  understand about the learning factors: (1) w1 governs the "cognitive"
;  aspect of the motion (that which depends on its own memory) and w2
;  governs the social aspect of the motion (that which depends on the
;  discovery of the entire swarm).  Note that the scales of the parameters,
;  as defined by PRANGE, is absorbed into the learning factors.
;
;  The position in parameter space is updated according to the equation:
;     p[t+dt] = p[t] + dt*v[t]
;
;  Note that "elastic" boundary conditions are implemented so that if an
;  agent exceeds the limits of the parameter range (PRANGE) then it will
;  be "reflected" (or bounce) from the boundary.
;
;  For additional details on this optimization technique please consult
;  the following web site: http://www.swarmintelligence.org/
;
; RETURN_VALUE:
;  Best parameters found in the minimization.
;
; CATEGORY:
;  OPTIMIZATION, OBJECTS, SWARM INTELLIGENCE
;
; CALLING SEQUENCE:
;  RESULT = RMD_PSO(PRANGE = prange,FUNCTION_NAME = function_name,....)
;
; PARAMETERS:
;  NONE
;
; KEYWORDS:
;  FTOL:          Fractional tolerance defined as the average fitness
;                 divided by the best fitness.  FTOL is the termination
;                 criteria for the quantity defined above (DEF: 0.1)
;  PRANGE:        Range of parameters (2 by N array where N is the number
;                 of parameters in the optimization).
;  WEIGHTS:       Learning factors.  (w1 and w2 in equation above)  (DEF: [1.,1.])
;  FUNCTION_NAME: Name of the function to be used in the minimization
;  FUNCTION_VALUE:Value of FUNCTION_NAME at the result found in the minimization.
;  FUNCTARGS:     Exogenous inputs required for function evaluation above
;                 and beyond the parameters in the optimization.
;  NUM_PARTICLES: Number of particles in the swarm (DEF: 10)
;  ITMAX:         Maximum number of iterations to perform if FTOL termination
;                 criteria not met.
;  VEL_FRACTION:  The fraction of the parameter range by which the initial
;                 velocity is created. (DEF: 0.1)
;  VEL_DECREMENT: The fractional amount by which the velocity will be decremented
;                 at each iteration.  Also called velocity damping and inertia weight.
;                 This is the parameter designated b in the equation above.  (DEF: 0.95)
;  ITERPROC:      Name (string) of the procedure to be called at
;                 each iteration step.  An example of such a procedure is
;                 PSO_ITERPROC, shown below.
;
;     The requirements on the intermediate reporting procedure are
;     that it be written as a PROCEDURE and it must accept the
;     following arguments: FUNC, P (the current best parameter set),
;     ITER (the current iteration, and INTERRUPT (a byte variable equal 1B or 0B
;     indicating if the algorithm should cease at the next update--note that
;     interrupt must be a variable since it is expected to be passed by reference).
;     Optional keywords are FUNCTARGS, ITERARGS, and OREF.  OREF is provided
;     so that you can extract information on the PSO object (its state)
;     at the current iteration.  This can be useful if you want to get
;     all of the current members of the swarm and display them, for instance.
;
;  ITERARGS:      Structure containing exogenous information necessary
;                 to perform the intermediate reporting procedure
;                 ITERPROC.
;  NCALLS:        Number of function calls over the course of the
;                 optimization procedure.
;  QUIET:         If set, no intermediate output will be provided
;                 (i.e. ITERPROC will not be called).  Default: 1B
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
; EXAMPLE:
;  Example usage found at end of this source listing namded
;  PSO_MIN_EXAMPLE.  To run this code, compile this file and
;  type PSO_MIN_EXAMPLE at the IDL prompt.
;
;     IDL> PSO_MIN_EXAMPLE
;
; AUTHOR:
;  Robert Dimeo
;  NIST Center for Neutron Research
;  National Institute of Standards and Technology
;  100 Bureau Drive-Stop 8562
;  Gaithersburg, MD 20899
;
; DISCLAIMER
;  This software is provided as is without any warranty whatsoever.
;  Permission to use, copy, modify, and distribute modified or
;  unmodified copies is granted, provided this disclaimer
;  is included unchanged.
;
; ACKNOWLEDGMENT:
;  This work is part of the DAVE development effort at the NCNR and
;  was supported in part by the National Science Foundation
;  under Agreement No. DMR-0086210.
;
; MODIFICATION HISTORY:
;  Written by RMD (12/20/04)
;  Updated header documentation, tightened up some
;  of the code, and corrected a few errors in the implementation
;  of the intermediate reporting (ITERPROC) -- RMD (12/22/04)
;  Replaced the periodic boundary conditions with "elastic"
;  boundaries -- RMD (12/22/04)
;  Removed the velocity tolerance and redefined FTOL so that
;  it is not problem nor scale-dependent -- RMD (12/22/04)
;  Removed the scale dependence of the learning factors
;     -- RMD (12/22/04)
;
; -
; *************************************************** ;
pro pso::cleanup
compile_opt idl2,hidden
ptr_free,self.v_ptr,self.p_ptr,self.pbest_ptr
ptr_free,self.gbest_ptr,self.prange_ptr,self.f_ptr
ptr_free,self.fbest_ptr,self.gfbest_ptr,self.iterargs_ptr
ptr_free,self.functargs_ptr,self.seed_ptr
end
; *************************************************** ;
pro pso_iterproc,    func,                   $
                     p,                      $
                     iter,                   $
                     interrupt,              $
                     functargs = functargs,  $
                     oref = opso,            $
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
; *************************************************** ;
function pso::evaluate_function
compile_opt idl2,hidden
p = *self.p_ptr
psize = size(p)
nparms = psize[1]
f = fltarr(self.num_particles)
for i = 0,self.num_particles-1 do   $
   f[i] = call_function(self.func,reform(p[*,i]),_Extra = *self.functargs_ptr)
*self.f_ptr = f
self.ncalls = self.ncalls + self.num_particles
return,1
end
; *************************************************** ;
pro pso::get_property,  ncalls = ncalls,           $
                        presult = presult,         $
                        fresult = fresult,         $
                        params = params,           $
                        feval_ave = feval_ave
compile_opt idl2,hidden
params = *self.p_ptr
ncalls = self.ncalls
presult = *self.gbest_ptr
fresult = *self.gfbest_ptr
if arg_present(feval_ave) then begin
   pop_stats = moment(*self.f_ptr)
   f = abs((pop_stats[0])-(*self.gfbest_ptr))
   feval_ave = abs(0.5-0.5*(1.+(2./!pi)*atan(f)))
endif
end
; *************************************************** ;
function pso::go_swarm
compile_opt idl2,hidden
; This is the main driver that determines the trajectories
; of the individual swarm members.
factor = self.decrement^(dindgen(self.niter))
prange = *self.prange_ptr
;for i = 0,self.niter-1 do begin
i = 0
ftol = 1.e5
; Get some numbers which don't change and that we'll
; need within the WHILE loop
prange = *self.prange_ptr
plo = reform(prange[0,*])
phi = reform(prange[1,*])
nparams = n_elements(plo)
nparticles = self.num_particles
uparams = 1+bytarr(nparams)
uparticles = 1+bytarr(nparticles)

dp = (phi-plo)#uparticles
w = self.weights
; Scale the weights to reflect the scale of the
; parameter ranges
w1 = ((uparams#uparticles)*(0.5/dp))*w[0]
w2 = ((uparams#uparticles)*(0.5/dp))*w[1]
;w1 = ((uparams#uparticles))*w[0]
;w2 = ((uparams#uparticles))*w[1]
s = *self.seed_ptr

while (ftol gt self.ftol) and (i lt self.niter) do begin
;while (i lt self.niter) do begin
   p = *self.p_ptr
   ; Evaluate the function and store the best parameters
   ; found in the current population and those found so far
   ; overall.
   ret = self->evaluate_function()
   if i eq 0 then begin
   ; Store the best parameters so far of each individual
      *self.pbest_ptr = *self.p_ptr ; best fitness-based parameters
      *self.fbest_ptr = *self.f_ptr ; best fitness values
      fbest = min(*self.f_ptr,index)
      *self.gbest_ptr = (*self.p_ptr)[*,index]  ; best parameters so far
      *self.gfbest_ptr = fbest   ; best fitness value so far
   endif else begin
      for j = 0,self.num_particles-1 do begin
         if (*self.f_ptr)[j] lt (*self.fbest_ptr)[j] then begin
            (*self.fbest_ptr)[j] = (*self.fbest_ptr)[j]
            (*self.pbest_ptr)[*,j] = (*self.p_ptr)[*,j]
         endif
      endfor
      fbest_current = min(*self.f_ptr,index)
      if fbest_current lt (*self.gfbest_ptr) then begin
         *self.gfbest_ptr = fbest_current
         *self.gbest_ptr = (*self.p_ptr)[*,index]
      endif
   endelse

   ; Calculate the termination criteria (tolerance)
   f = (moment(*self.f_ptr))[0]
   ftol = abs(0.5*(1.+(2./!pi)*atan(f)))

   if i lt 5 then ftol = 1.e4
   ; Update the particle velocity
   v = (vold = *self.v_ptr)
   p = *self.p_ptr

   gbest = *self.gbest_ptr
   pbest = *self.pbest_ptr
   dv =  ((w1*randomu(s,nparams,nparticles))*(uparams#uparticles))*(pbest-p)+ $
         ((w2*randomu(s,nparams,nparticles))*(uparams#uparticles))*(gbest#uparticles-p)
   *self.v_ptr = factor[i]*vold+dv

   ; Use "elastic" boundary conditions so that, if an agent
   ; exceeds one of the boundaries, its velocity will be
   ; reversed.
   ptest = p + (self.dt) * (*self.v_ptr)
   for j = 0,nparams-1 do begin
      too_big = where(ptest[j,*] gt prange[1,j],count_big)
      too_small = where(ptest[j,*] lt prange[0,j],count_small)
      if count_small gt 0 then begin
         (*self.v_ptr)[j,too_small] = -(*self.v_ptr)[j,too_small]
      endif
      if count_big gt 0 then begin
         (*self.v_ptr)[j,too_big] = -(*self.v_ptr)[j,too_big]
      endif
   endfor

   *self.p_ptr = p + self.dt * (*self.v_ptr)
   if ~self.quiet then begin

      interrupt = self.interrupt

      call_procedure,   self.iterproc,self.func,            $
                        gbest,i,                            $
                        interrupt,                          $
                        oref = self,                        $
                        functargs = *self.functargs_ptr,    $
                        _Extra = *self.iterargs_ptr
      self.interrupt = interrupt
   endif
   i++
   if self.interrupt then i = self.niter
endwhile
return,1
end
; *************************************************** ;
function pso::randomize_positions
compile_opt idl2,hidden
; Pick random locations for the number of particles in
; the swarm.
plo = reform((*self.prange_ptr)[0,*]) & phi = reform((*self.prange_ptr)[1,*])
dp = phi-plo
nparticles = self.num_particles
psize = size(*self.prange_ptr)
nparams = (psize[0] gt 1) ? psize[2]:psize[0]
uparticles = rebin([1B],nparticles,/sample)

*self.p_ptr = (plo#uparticles)+(dp#uparticles)*randomu(s,nparams,nparticles)
*self.v_ptr = (dp#uparticles)*self.velocity_fraction*(-1.+2.*(randomu(s,nparams,nparticles) gt 0.0))
*self.seed_ptr = s

return,1
end
; *************************************************** ;
function pso::init,  prange = prange,                 $
                     ftol = ftol,                     $
                     function_name = function_name,   $
                     num_particles = num_particles,   $
                     itmax = itmax,                   $
                     weights = weights,               $
                     vel_decrement = vel_decrement,   $
                     vel_fraction = vel_fraction,     $
                     iterproc = iterproc,             $
                     iterargs = iterargs,             $
                     quiet = quiet,                   $
                     functargs = functargs,           $
                     _Extra = extra

compile_opt idl2,hidden
if n_elements(iterproc) eq 0 then iterproc = 'pso_iterproc'
if (strupcase(iterproc) eq 'PSO_ITERPROC') and $
   (n_elements(iterargs) eq 0) then begin
   iterargs = {iterstop:1}
   interrupt = 0B
endif
self.iterproc = iterproc
self.iterargs_ptr = ptr_new(iterargs)
if n_elements(quiet) eq 0 then quiet = 1B
self.quiet = quiet

self.dt = 1.
self.ncalls = 0L
if n_elements(ftol) eq 0 then ftol = 0.1
self.ftol = ftol
if n_elements(vel_fraction) eq 0 then vel_fraction = 0.1
if n_elements(weights) eq 0 then self.weights = [1.,1.] else $
   self.weights = weights
if n_elements(itmax) eq 0 then self.niter = 20 else $
   self.niter = itmax
if n_elements(prange) eq 0 then return,0
if n_elements(function_name) eq 0 then return,0
if n_elements(num_particles) eq 0 then num_particles = 10
self.num_particles = num_particles
self.prange_ptr = ptr_new(prange,/no_copy)
self.pbest_ptr = ptr_new(/allocate_heap)
self.gbest_ptr = ptr_new(/allocate_heap)
self.v_ptr = ptr_new(/allocate_heap)
self.p_ptr = ptr_new(/allocate_heap)
self.f_ptr = ptr_new(/allocate_heap)
self.gfbest_ptr = ptr_new(/allocate_heap)
self.fbest_ptr = ptr_new(/allocate_heap)
self.velocity_fraction = vel_fraction
self.func = function_name
self.seed_ptr = ptr_new(/allocate_heap)
if n_elements(vel_decrement) eq 0 then vel_decrement = 0.95
self.decrement = vel_decrement < 1.0
self.functargs_ptr = ptr_new(/allocate_heap)
if n_elements(functargs) eq 0 then functargs = {iterstop:0}
*self.functargs_ptr = functargs
return,1
end
; *************************************************** ;
pro pso__define
compile_opt idl2,hidden
void =   {  pso,                                   $
            func:'',                               $
            dt:0.0,                                $
            num_particles:0L,                      $
            velocity_fraction:0.,                  $
            ftol:0.0,                              $
            weights:fltarr(2),                     $
            decrement:0.0,                         $
            interrupt:0B,                          $
            ncalls:0L,                             $
            niter:0L,                              $
            iterargs_ptr:ptr_new(),                $
            iterproc:'',                           $
            quiet:0B,                              $
            seed_ptr:ptr_new(),                    $
            functargs_ptr:ptr_new(),               $
            v_ptr:ptr_new(),                       $
            p_ptr:ptr_new(),                       $
            f_ptr:ptr_new(),                       $
            pbest_ptr:ptr_new(),                   $
            gbest_ptr:ptr_new(),                   $
            fbest_ptr:ptr_new(),                   $
            gfbest_ptr:ptr_new(),                  $
            prange_ptr:ptr_new()                   $
         }
end
; *************************************************** ;
function rmd_pso, ftol = ftol,                              $
                  function_name = function_name,            $
                  function_value = function_value,          $
                  prange = prange,                          $
                  quiet = quiet,                            $
                  itmax = itmax,                            $
                  ncalls = ncalls,                          $
                  functargs = functargs,                    $
                  weights = weights,                        $
                  vel_decrement = vel_decrement,            $
                  vel_fraction = vel_fraction,              $
                  num_particles = num_particles,            $
                  _Extra = extra
compile_opt idl2,hidden
if n_elements(iterproc) eq 0 then iterproc = 'pso_iterproc'
if (strupcase(iterproc) eq 'PSO_ITERPROC') and $
   (n_elements(iterargs) eq 0) then begin
   iterargs = {iterstop:1}
   interrupt = 0B
endif
if n_elements(functargs) eq 0 then functargs = {iterstop:0}
if n_elements(quiet) eq 0 then quiet = 1B
if n_elements(ftol) eq 0 then ftol = 0.1
if n_elements(vel_fraction) eq 0 then vel_fraction = 0.1
if n_elements(weights) eq 0 then weights = [1.,1.]
if n_elements(itmax) eq 0 then itmax = 20
if n_elements(prange) eq 0 then return,0
if n_elements(function_name) eq 0 then return,0
if n_elements(num_particles) eq 0 then num_particles = 10
if n_elements(vel_decrement) eq 0 then vel_decrement = 0.95

o = obj_new('pso',   prange = prange,                 $
                     ftol = ftol,                     $
                     function_name = function_name,   $
                     num_particles = num_particles,   $
                     itmax = itmax,                   $
                     weights = weights,               $
                     vel_decrement = vel_decrement,   $
                     vel_fraction = vel_fraction,     $
                     iterproc = iterproc,             $
                     iterargs = iterargs,             $
                     functargs = functargs,           $
                     quiet = quiet,                   $
                     _Extra = extra                   )

ret = o->randomize_positions()
ret = o->go_swarm()
o->get_property,ncalls = ncalls,presult = presult,fresult = fresult
obj_destroy,o
function_value = fresult
return,presult
end
; *************************************************** ;
; Begin example usage
; *************************************************** ;
function linspace,xlo,xhi,nx
; Function that creates nx regularly-spaced set of values
; that vary between xlo and xhi
return,xlo+((xhi-xlo)/(nx-1.0))*findgen(nx)
end
; *************************************************** ;
function f1_eval,p,_EXTRA = extra
; This function has a couple of minima and maxima.  This
; particular function is the one that the algorithm will
; call and evaluate agents positions, one at a time.
x = 1.*p[0] & y = 1.*p[1]
z = 3.*(1.-x)^2 * exp(-x^2 - (y+1.)^2) - $
    10.*(x/5 - x^3 - y^5)*exp(-x^2-y^2) - $
    (1./3)*exp(-(x+1)^2 - y^2)
return,z
end
; *********************************** ;
function f1_plot,x1,y1,_Extra = extra
; This is the same function as above except that this
; one is used for calculating the image for display
; purposes.
ux = 1+bytarr(n_elements(x1))
uy = 1+bytarr(n_elements(y1))
x = 1.*x1#uy & y = 1.*ux#y1
z = 3.*(1.-x)^2 * exp(-x^2 - (y+1.)^2) - $
    10.*(x/5 - x^3 - y^5)*exp(-x^2-y^2) - $
    (1./3)*exp(-(x+1)^2 - y^2)
return,z
end
; *************************************************** ;
function pso_test_convert,xr,yr,p,xsize,ysize
; Convert the parameter values (data) to pixels

; xsize,ysize: size in pixels of the display window
; p: the parameter (2 by nparticles)
; x,y: data coordinates for the 2-d function
xd = reform(p[0,*]) & yd = reform(p[1,*])
xp = fix((-xr[0]+xd)*xsize/(xr[1]-xr[0]))
yp = fix((-yr[0]+yd)*ysize/(yr[1]-yr[0]))
return,transpose([[xp],[yp]])
end
; *************************************************** ;
pro pso_test_iterproc,  func,                   $
                        p,                      $
                        iter,                   $
                        interrupt,              $
                        functargs = functargs,  $
                        oref = opso,            $
                        _Extra = extra
compile_opt hidden,idl2
xlo = min(extra.x,max = xhi) & ylo = min(extra.y,max = yhi)
wset,extra.winpix
title = '!3Iteration: '+strtrim(string(iter+1),2)
device,copy = [0,0,!d.x_size,!d.y_size,0,0,extra.cleanpix]
; Get the full parameter set out of the PSO object
opso->get_property,params = params
pplot = pso_test_convert([xlo,xhi],[ylo,yhi],params,!d.x_size,!d.y_size)
for j = 0,extra.n-1 do plots,[pplot[0,j]],[pplot[1,j]], $
   psym = 2,symsize = 1.5,                      $
   color = 255,/device,/continue
xyouts,0.3,0.05,title,/normal,color = 255,charsize = 2.0,charthick = 2.
wset,extra.winvis
device,copy = [0,0,!d.x_size,!d.y_size,0,0,extra.winpix]
wait,extra.delay
end
; *************************************************** ;
pro pso_min_example

; Uses files listed above such as LINSPACE, F1_EVAL,
; F1_PLOT, PSO_TEST_CONVERT, PSO_TEST_ITERPROC

prange = [[-3,9.],[-9.,4.]]   ; Range of the two parameters
func = 'f1_eval'              ; Function to be minimized
plot_func = 'f1_plot'         ; Same as above but modified
                              ; so that it can be plotted.
device,decomposed = 0
loadct,1,/silent
n = 20                        ; number of agents in the swarm
delay = 5.e-5                 ; delay in seconds in updating the display
                              ; i.e. "animation frame rate".  Set this to
                              ; a larger number like 0.1 to see the swarm
                              ; motion in very deliberate steps.

; Set up the image plot parameters
nx = (ny = 100)
x = linspace(prange[0,0],prange[1,0],nx)
y = linspace(prange[0,1],prange[1,1],ny)
z = call_function(plot_func,x,y,_Extra = extra)
; Set up the pixmaps and plot windows to display the
; swarm motion.
winvis = 0
xsize = 400 & ysize = 400
window,winvis,xsize = xsize,ysize = ysize
window,/free,/pixmap,xsize = xsize,ysize = ysize
winpix = !d.window
window,/free,/pixmap,xsize = xsize,ysize = ysize
cleanpix = !d.window
wset,winpix
tv,bytscl(smooth(congrid(z,xsize,ysize),5))
wset,winvis
device,copy = [0,0,!d.x_size,!d.y_size,0,0,winpix]
; Copy the pristine image to the CLEANPIX pixmap
; to allow rapid updates as the swarm motion proceeds.
wset,cleanpix
device,copy = [0,0,!d.x_size,!d.y_size,0,0,winpix]
; Create a structure to be passed into the intermediate
; reporting procedure.  The information here will allow
; us to see the swarm in motion during the minimization.
iterargs =  {  winpix:winpix,    $
               winvis:winvis,    $
               cleanpix:cleanpix,$
               x:x,y:y,n:n,      $
               nx:nx,ny:ny,      $
               delay:delay       }
; Call the PSO minimization routine
p = rmd_pso(         ftol = 1.e-1,                          $
                     function_name = func,                  $
                     function_value = fval,                 $
                     ncalls = ncalls,                       $
                     weights = [1.,2.],                     $
                     itmax = 400,                           $
                     quiet = 0B,                            $
                     iterargs = iterargs,                   $
                     iterproc = 'pso_test_iterproc',        $
                     vel_fraction = 5.e-2,                  $
                     num_particles = n,                     $
                     vel_decrement = 0.98,                 $
                     prange = prange                        )
; Clean up the pixmaps and print the results.
wdelete,winpix
wdelete,cleanpix
print,'Result: ',p
print,'Value: ',fval
print,'# function calls: ',ncalls
end
