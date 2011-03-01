; *********************************** ;
pro ga_ising_iterproc,                       $
                     func,                   $
                     p,                      $
                     iter,                   $
                     interrupt,              $
                     functargs = functargs,  $
                     oref = oga,             $
                     _Extra = extra
compile_opt hidden,idl2
oga->get_property,fitness = fitness,itmax = itmax,ave_fitness = ave_fitness, $
   dave_fitness = dave_fitness
x = 1+indgen(iter+1)
y = ave_fitness[0:iter]
dy = dave_fitness[0:iter]
loadct,0,/silent & white = 255 & black = 0
wset,extra.winpix
plot,[x],[y],psym = -8,title = '!3Average Energy in Population',xtitle = 'Generation', $
   ytitle = '<E(S)>',xrange = [0,iter+2],/xsty,/nodata,color = black,background = white
oplot,[x],[y],psym = -8,color = black

plots,!x.crange,[extra.emin,extra.emin],linestyle = 2,thick = 2.0,/data, $
   color = black
wset,extra.winvis
device,copy = [0,0,!d.x_size,!d.y_size,0,0,extra.winpix]

tvlct,r,g,b,/get
loadct,1,/silent
j = *functargs.jptr
ns = n_elements(p)
nb = ns-1
im = fltarr(ns+nb)
im[0:ns+nb-1:2] = p[0:ns-1]
im[1:ns+nb-2:2] = j[0:nb-1]
wset,extra.impix
tv,bytscl(congrid(rebin(im,ns+nb,2,/sample),!d.x_size,!d.y_size))
;tv,bytscl(congrid(rebin(p,ns,2,/sample),!d.x_size,!d.y_size))
wset,extra.imvis
device,copy = [0,0,!d.x_size,!d.y_size,0,0,extra.impix]
tvlct,r,g,b

strout = strtrim(string(iter),2)
widget_control,extra.label_id,set_value = strout
event = widget_event(/nowait,extra.stop_but)
evname = tag_names(event,/structure_name)
if evname eq '' or event.id eq 0 then return

if evname eq 'WIDGET_BUTTON' and event.select then begin
   interrupt = 1B
endif

end
; *********************************** ;
function ising_energy,s,_Extra = extra
j = *extra.jptr
z = s
zsize = size(z)
if zsize[0] eq 2 then begin
   npop = zsize[1] & nspins = zsize[2]
   ; Sum over the second second dimension
   jmat = rebin(transpose(j),npop,nspins-1,/sample)
   energy = 1.*total(-jmat*z[*,0:nspins-2]*z[*,1:nspins-1],2)
endif else begin
   nspins = zsize[1]
   energy = 1.*(total(-j[0:nspins-2] * z[0:nspins-2] * z[1:nspins-1]))
endelse
return,energy
end
; *********************************** ;
pro rmd_ga_ising, nspins,                 $
                  npop = npop,            $
                  random = random,        $
                  pcross = pcross,        $
                  pmutate = pmutate,      $
                  micro_ga = micro_ga,    $
                  boltzmann = boltzmann,  $
                  ntimes = ntimes,        $
                  ftol = ftol

if n_elements(ntimes) eq 0 then ntimes = 1
if n_elements(boltzmann) eq 0 then boltzmann = 0B
if n_elements(micro_ga) eq 0 then micro_ga = 0B
if n_elements(pmutate) eq 0 then pmutate = 0.033
if n_elements(pcross) eq 0 then pcross = 0.9
device,decomposed = 0
loadct,0,/silent
if n_params() eq 0 then $
   nspins = 20L
if n_elements(npop) eq 0 then npop = fix(nspins^1.75)
nbonds = nspins - 1L
; Create a random bonding between neighboring spins
if ~keyword_set(random) then $
   j = replicate(-1.,nbonds) else $
   j = 0.25*randomn(s,nbonds)
;j[0:nbonds-1:2] = randomn(s,nbonds/2+1)
jptr = ptr_new(j,/no_copy)
; Determine the theoretical minimum energy
emin = -total(abs(*jptr))
; Now use a GA to minimize the energy of the spins
func = 'ising_energy'
iterproc = 'ga_ising_iterproc'
keep_best = 1B

functargs = {jptr:jptr}
quiet = 0B


title = 'Press to interrupt minimization'
stop_base = widget_base(/col,title = title,tlb_frame_attr = 2)
strout = 'Iteration '+strtrim(string(0),2)
label_id = widget_label(stop_base,value = strout,/dynamic_resize)
stop_but = widget_button(stop_base,value = 'Interrupt', xsize = 250)
widget_control,stop_base,/realize

if ~quiet then begin
   xsize = 400 & ysize = 400
   ypos = 200
   window,0,xsize = xsize,ysize = ysize,xpos = 500,ypos = ypos
   winvis = 0
   window,1,xsize = xsize,ysize = 50,xpos = 500,ypos = ypos - 75,title = 'Spin chain configuration'
   imvis = !d.window
   window,/free,/pixmap,xsize = xsize,ysize = ysize
   winpix = !d.window

   window,/free,/pixmap,xsize = xsize,ysize = 50
   impix = !d.window
   iterargs = {winvis:winvis,winpix:winpix,emin:emin,    $
               label_id:label_id,stop_but:stop_but,      $
               imvis:imvis,impix:impix}
   iterproc = 'ga_ising_iterproc'
endif

nth = 30
th = (2.0*!pi/nth)*findgen(nth)
xc = cos(th) & yc = sin(th)
usersym,xc,yc
if n_elements(ftol) eq 0 then ftol = 1.e-4
itmax = 200
success = bytarr(ntimes)
for i = 0,ntimes-1 do begin
o = obj_new('rmd_ga_ising',                              $
                        ftol,                            $
                        nspins = nspins,                 $
                        function_name = func,            $
                        pcross = pcross,                 $
                        itmax = itmax,                   $
                        pmutate = pmutate,               $
                        npop = npop,                     $
                        micro_ga = micro_ga,             $
                        keep_best = keep_best,           $
                        functargs = functargs,           $
                        boltzmann = boltzmann,           $
                        iterproc = iterproc,             $
                        iterargs = iterargs,             $
                        quiet = quiet)

ret = o->go()
o->get_property,best_parms = best_parms,ncalls = ncalls
fval = call_function(func,best_parms,jptr = jptr)
print,'*********************************************'
;print,best_parms
print,'Minimum energy (GA): '+strtrim(string(fval),2)
print,'Minimum energy (theory): '+strtrim(string(emin),2)
if fval eq emin then success[i] = 1B
print
nconfig = nspins* alog10(2.0) / alog10(10.0)
this_format = '(g15.3)'
print,'Number of possible configurations: '+strtrim(string((1.e1)^nconfig,format = this_format),2)
print,'Number of function calls: '+strtrim(string(ncalls),2)
print
print,'Percentage probed: '+strtrim(string(100.*ncalls/((1.e1)^nconfig)),2)

obj_destroy,o
endfor
ok = where(success,count_success)
print,'Success rate: '+strtrim(string(count_success),2)+'/'+strtrim(string(ntimes),2)
if ~quiet then wdelete,winpix
if ~quiet then wdelete,impix
ptr_free,jptr
widget_control,stop_base,/destroy
end