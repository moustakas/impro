;;;;;; HIIphot package 			  ;;;;;;
;;;;;; D. Thilker, R. Braun, R. Walterbos ;;;;;;
;;;;;; July 1998 - August 2001	 	  ;;;;;;

;;;;;; In order to run, do the following: ;;;;;;

;;;;;; 1) copy this file to the directory ;;;;;;
;;;;;;    of your choice, then type       ;;;;;;
;;;;;;    'setenv HIIphot $PWD'           ;;;;;;
;;;;;; 2) cd to your data directory	  ;;;;;;
;;;;;; 3) start IDL			  ;;;;;;
;;;;;; 4) type '.comp $HIIphot/HIIphot'   ;;;;;;
;;;;;; 5) type 'HIIphot'		  ;;;;;;

function index2coord, index, nx
;converts 'where' index to IDL pixel coord (for 2D images only!)
;inputs:
;       index   integer, 'where' index
;       nx      integer, X size of image
;output:
;       index2coord     lonarr(2), IDL pixel position @ index

index2coord=[(index mod nx),index/nx]

return,index2coord
end
;-----------------------------------------------------

function coord2index, coord, nx
;converts IDL pixel coord to 'where' index (for 2D images only!)
;inputs:
;       coord   lonarr(2), IDL pixel position
;       nx      integer, X size of image

;output:
;       coord2index     integer, 'where' index for pixel @ coord
 
coord2index=coord(0)+coord(1)*nx
 
return,coord2index
end
;-----------------------------------------------------

function coord2coord, coord_in, offset

coord2coord=[coord_in(0)+offset(0),coord_in(1)+offset(1)]

return,coord2coord
end

;-----------------------------------------------------
;-----------------------------------------------------

PRO HIIphot, NO_DISPLAY=no_display, BATCH=batch, HELP=help

; version 5.1
; David A. Thilker
; July 2001

; fixed bugs left in v5.0 / noticed by Marc Sauvage (thanks!)
; added PPP methodology, took out most IDL "saves" 
; added check on median surface brightness (before growth)

;;;;;;;;;;;;;;;;;;;;;;;; INITIALIZATION / SESSION STARTUP

npar=N_params() ; number of parameters used calling HIIphot

code_version='HIIphot v5.1'

if keyword_set(HELP) then begin
	echo_help
	return
endif

;this bit guarantees that if user has given something in quotes for
;one of the keywords, that the code detects the right keyword
keyword_set_batch=0L
size_batch=size(batch) & size_batch=size_batch(1)
if size_batch eq 7 then begin
	if strlen(batch) gt 0 then keyword_set_batch=1L
endif

;get some ugly messages out of the way
if not keyword_set_batch then proceed=read_key(0)
working_dir=getenv('PWD')
astrolib
junk=avg([0.,0.,0.])

if not keyword_set_batch then scr_erase

if keyword_set (NO_DISPLAY) then display=0L else display=1L

;setup constants to interpret (read_key) results
ynum=121 & cynum=89 & nnum=110 & cnnum=78
spcnum=32 & retnum=10 & qnum=113 & cqnum=81

if (not keyword_set_batch) then begin
	print,'Interactive syntax:'
	print,''
	print,'HIIphot [, /NO_DISPLAY, BATCH=''file'', /HELP ]'
	print,''
	print,'Everything in brackets is optional.'
	print,'To prepare for batch mode jobs, use HIIphot_batch.'
	print,''
	print,'Hit (q) to exit, otherwise any key to run interactively.'
	proceed=read_key(1)
	if (proceed eq qnum) or (proceed eq cqnum) then return
 	batch='NULL'
	root='HIIphot'
	footprint_lev=0.2
	seed_lev=0.5
	f_remain_crit=0.9
	rho_crit=0.25
	termgrad_arr=[1.5,0.,0.,0.,0.]
	multiple_cutoffs=0L
	docongrid=0L
	interact=1L 
endif else begin ; non-interactive session, read input from *.HIIphot_batch
	interact=0L
	display=0L
	;judge how many lines the shell prepends to spawn output
	spawn,'pwd',spawn_out
	spawn_off=n_elements(spawn_out)-1
	;read in the *.batch file
	spawncmd=strcompress('wc '+batch)
	spawn,spawncmd,spawn_out
	nlines=0L
	junk=0L
	reads,spawn_out(spawn_off),nlines,junk,junk
	openr,1,batch
	batchinfo=strarr(nlines-1)
	;stop if the *.batch file looks incorrect 
	if nlines lt 28 then begin
		close,1
		print,'Warning! Looks like an old/incomplete batch file...'
		return
	endif
	junk=''
	readf,1,junk 
	readf,1,batchinfo
	hs_filename=strmid(batchinfo(0),strpos(batchinfo(0),'=')+2,132)
	c_filename=strmid(batchinfo(1),strpos(batchinfo(1),'=')+2,132)
	h_filename=strmid(batchinfo(2),strpos(batchinfo(2),'=')+2,132)
	b_filename=strmid(batchinfo(3),strpos(batchinfo(3),'=')+2,132)
	gc=float(strmid(batchinfo(4),strpos(batchinfo(4),'=')+2,132))
	nc=float(strmid(batchinfo(5),strpos(batchinfo(5),'=')+2,132))
	scalec=float(strmid(batchinfo(6),strpos(batchinfo(6),'=')+2,132))
	ghc=float(strmid(batchinfo(7),strpos(batchinfo(7),'=')+2,132))
	nhc=float(strmid(batchinfo(8),strpos(batchinfo(8),'=')+2,132))
	scalehc=float(strmid(batchinfo(9),strpos(batchinfo(9),'=')+2,132))
	scaleh=float(strmid(batchinfo(10),strpos(batchinfo(10),'=')+2,132))
	scaleh_em=float(strmid(batchinfo(11),strpos(batchinfo(11),'=')+2,132))
	NIIpercent=float(strmid(batchinfo(12),strpos(batchinfo(12),'=')+2,132))
	dist_mpc=float(strmid(batchinfo(13),strpos(batchinfo(13),'=')+2,132))
	scale=float(strmid(batchinfo(14),strpos(batchinfo(14),'=')+2,132))
	sigma_find=float(strmid(batchinfo(15),strpos(batchinfo(15),'=')+2,132))
	PSFfwhm_pix=float(strmid(batchinfo(16),strpos(batchinfo(16),$
		'=')+2,132))
	size_max_pc=float(strmid(batchinfo(17),strpos(batchinfo(17),$
		'=')+2,132))
	bkgd_annulus_pc=float(strmid(batchinfo(18),strpos(batchinfo(18),$
		'=')+2,132))
	docongrid=round(float(strmid(batchinfo(19),strpos(batchinfo(19),$
		'=')+2,132)))
	root=strmid(batchinfo(20),strpos(batchinfo(20),'=')+2,132)
	dump_stamps=round(float(strmid(batchinfo(21),strpos(batchinfo(21),$
		'=')+2,132)))
	footprint_lev=float(strmid(batchinfo(22),strpos(batchinfo(22),$
		'=')+2,132))
	seed_lev=float(strmid(batchinfo(23),strpos(batchinfo(23),'=')+2,132))
	f_remain_crit=float(strmid(batchinfo(24),strpos(batchinfo(24),$
		'=')+2,132))
	rho_crit=float(strmid(batchinfo(25),strpos(batchinfo(25),'=')+2,132))
	termgrad_str=strmid(batchinfo(26),strpos(batchinfo(26),'=')+2,132)
	termgrad_arr=fltarr(5)
	termgrad_arr(0)=float(strmid(termgrad_str,0,9))
	termgrad_arr(1)=float(strmid(termgrad_str,10,9))
	termgrad_arr(2)=float(strmid(termgrad_str,20,9))
	termgrad_arr(3)=float(strmid(termgrad_str,30,9))
	termgrad_arr(4)=float(strmid(termgrad_str,40,9))
	if max(termgrad_arr(1:4)) gt 1.e-6 then multiple_cutoffs=1L else $
		multiple_cutoffs=0L
	close,1

endelse

;echo which mode we have entered
wait,0.5
if not keyword_set_batch then scr_erase
print,'-------------------------------------------------------------------------------'
if ((interact eq 1) and (display eq 0)) then print,$
	'Starting interactive, text-only HIIphot session.'
if ((interact eq 1) and (display eq 1)) then print,$
	'Starting interactive HIIphot session.'
if (interact ne 1) then print,$
	'Starting non-interactive HIIphot session.'
print,'-------------------------------------------------------------------------------'
wait,0.5
if not keyword_set_batch then scr_erase

;;;;;;;;;;;;;;;;;;;;;; SELECT & READ IMAGES

;get all 4 image names, if they are to be used
select_images,display,interact,hs_filename,c_filename,h_filename,b_filename,$
	use_cimage,use_himage,use_bimage

;read images, extract header info when needed
;use astrometry from cont-sub image if available, else check continuum img
if not keyword_set_batch then scr_erase
print,'Reading LINE image...'
fxread, hs_filename, hsimage, hshead	;in supplied units
index_bad=where(finite(hsimage) eq 0,count_bad)
if count_bad gt 0 then hsimage(index_bad)=0.
if interact then begin
	ghc=fxpar(hshead,'GAIN')
	if ghc eq 0. then ghc=3.
	nhc=fxpar(hshead,'NCOMBINE')
	if nhc eq 0 then nhc=1
endif
naxis=fxpar(hshead,'NAXIS*')
naxis_orig=naxis
extast,hshead,astrometry,noparams
equinox_img=fxpar(hshead,'EQUINOX')
if equinox_img eq 0 then equinox_img=fxpar(hshead,'EPOCH')
equinox_img=float(equinox_img)
use_coords=1L
;demand that the header contains CDn_m astrometry in order to use WCS info
;it can also contain ROTA + CDELT (AIPS-type) astrometry, as long as the
;CD matrix is filled in a consistent manner.
if noparams lt 2 then use_coords=0L
crota=fltarr(2)
if use_coords then begin
	getrot,hshead,rot,cdelt
	crota(0)=0.
	crota(1)=rot
	print,'Found WCS info.'
	if equinox_img lt 1.e-6 then begin
		print,'EQUINOX keyword is missing! Assuming J2000 coords.'
		equinox_img=2000.
	endif
	astrometry_big=astrometry
	astrometry_big.cd=astrometry_big.cd/2.
	astrometry_big.crpix=astrometry_big.crpix*2.
endif
print,''

if use_cimage then begin
print,'Reading CONTINUUM image...'
fxread,c_filename,cimage,chead ; in supplied units 
index_bad=where(finite(cimage) eq 0,count_bad)
if count_bad gt 0 then cimage(index_bad)=0.
if interact then begin
	gc=fxpar(chead,'GAIN')
	if gc eq 0. then gc=3.
	nc=fxpar(chead,'NCOMBINE')
	if nc eq 0 then nc=1
endif
cimage_naxis=fxpar(chead,'NAXIS*')
if naxis(0) ne cimage_naxis(0) then begin
	message,'ERROR: CONTINUUM image isn''t the required X size!',/CON
	return
endif
if naxis(1) ne cimage_naxis(1) then begin
	message,'ERROR: CONTINUUM image isn''t the required Y size!',/CON
	return
endif
if (noparams lt 2) then begin ; try again with astrometry
	extast,chead,astrometry,noparams
	equinox_img=fxpar(chead,'EQUINOX')
	if equinox_img eq 0 then equinox_img=fxpar(chead,'EPOCH')
	equinox_img=float(equinox_img)
	use_coords=1L
;demand that the header contains CDn_m astrometry in order to use WCS info
;it can also contain ROTA + CDELT (AIPS-type) astrometry, as long as the
;CD matrix is filled in a consistent manner.
	if noparams lt 2 then use_coords=0L
	crota=fltarr(2)
	if use_coords then begin
		getrot,chead,rot,cdelt
		crota(0)=0.
		crota(1)=rot
		print,'Found WCS info.'
		if equinox_img lt 1.e-6 then begin
			print,'EQUINOX keyword is missing! Assuming J2000 coords.'
			equinox_img=2000.
		endif
		astrometry_big=astrometry
		astrometry_big.cd=astrometry_big.cd/2.
		astrometry_big.crpix=astrometry_big.crpix*2.
	endif
endif
print,''
endif

;in the case that only the LINE and CONTINUUM images are available
;we will prompt for an appropriate scaling factor in order to
;create a pseudo LINE+CONTINUUM image
if use_cimage and not(use_himage) then begin
	pseudo_linecontinuum,hsimage,cimage,root,h_filename,$
		scaleh,scalec,scalehc,use_coords,astrometry,equinox_img
	use_himage=1L
endif
	
if use_himage then begin
print,'Checking LINE+CONTINUUM image...'
fxread,h_filename,himage,hhead,0,1,0,1
index_bad=where(finite(himage) eq 0,count_bad)
if count_bad gt 0 then himage(index_bad)=0.
if interact then begin
	ghc=fxpar(hhead,'GAIN')
	if ghc eq 0. then ghc=3.
	nhc=fxpar(hhead,'NCOMBINE')
	if nhc eq 0 then nhc=1
endif
himage_naxis=fxpar(hhead,'NAXIS*')
if naxis(0) ne himage_naxis(0) then begin
	message,'ERROR: LINE+CONTINUUM image isn''t the required X size!',/CON
	return
endif
if naxis(1) ne himage_naxis(1) then begin
	message,'ERROR: LINE+CONTINUUM image isn''t the required Y size!',/CON
	return
endif
print,''
endif

if use_bimage then begin
fxread,b_filename,bimage_float,bimage_header
index_bad=where(finite(bimage_float) eq 0,count_bad)
if count_bad gt 0 then bimage(index_bad)=0.
bimage_naxis=fxpar(bimage_header,'NAXIS*')
if naxis(0) ne bimage_naxis(0) then begin
	message,'ERROR: BLANKING image isn''t the right x size!',/CON
	return
endif
if naxis(1) ne bimage_naxis(1) then begin
	message,'ERROR: BLANKING image isn''t the right y size!',/CON
	return
endif
print,''
endif

if use_bimage then bimage=fix(ceil(bimage_float))

;;;;;;;;;;;;;;;;;;;;;; SELECT IMAGE REGIONS (FOR ANALYSIS AND SKY)

get_regions, display, interact, root, hsimage, hshead, naxis, h_filename, $
	use_bimage, bimage, blc, trc, blc_sky, trc_sky, sky_bimage, docongrid
if b_filename eq 'NULL' then bimage_float=float(bimage)

;;;;;;;;;;;;;;;;;;;;;; PAD IF NEEDED (ALONG BOTH AXES, TO SPEED UP FFTs)

if not keyword_set_batch then scr_erase
print,'Evaluating padding along both axes...'

;modified naxis, accounting for analysis section only
naxis(0)=trc(0)-blc(0)+1
naxis(1)=trc(1)-blc(1)+1

case 1 of
	((naxis(0) gt 16) and (naxis(0) lt 32)):	pad_x=32-naxis(0)
	((naxis(0) gt 32) and (naxis(0) lt 64)):	pad_x=64-naxis(0)
	((naxis(0) gt 64) and (naxis(0) lt 128)):	pad_x=128-naxis(0)
	((naxis(0) gt 128) and (naxis(0) lt 256)):	pad_x=256-naxis(0)
	((naxis(0) gt 256) and (naxis(0) lt 512)):	pad_x=512-naxis(0)
	((naxis(0) gt 512) and (naxis(0) lt 1024)):	pad_x=1024-naxis(0)
	((naxis(0) gt 1024) and (naxis(0) lt 2048)):	pad_x=2048-naxis(0)
else: 	begin
		pad_x=0L
		if naxis(0) le 2048 then begin
			print, '  No x-axis padding required.'
		endif else begin
			print, '  Padding would make x dimension too big.'
			print, '  Enjoy a book while waiting for slow FFTs.'
		endelse
	end
endcase
if pad_x gt 0 then print,'  Padding x-axis to ',naxis(0)+pad_x,' pixels.'
case 1 of
	((naxis(1) gt 16) and (naxis(1) lt 32)):	pad_y=32-naxis(1)
	((naxis(1) gt 32) and (naxis(1) lt 64)):	pad_y=64-naxis(1)
	((naxis(1) gt 64) and (naxis(1) lt 128)):	pad_y=128-naxis(1)
	((naxis(1) gt 128) and (naxis(1) lt 256)):	pad_y=256-naxis(1)
	((naxis(1) gt 256) and (naxis(1) lt 512)):	pad_y=512-naxis(1)
	((naxis(1) gt 512) and (naxis(1) lt 1024)):	pad_y=1024-naxis(1)
	((naxis(1) gt 1024) and (naxis(1) lt 2048)):	pad_y=2048-naxis(1)
else: 	begin
		pad_y=0L
		if naxis(1) le 2048 then begin
			print, '  No y-axis padding required.'
		endif else begin
			print, '  Padding would make y dimension too big.'
			print, '  Enjoy a book while waiting for slow FFTs.'
		endelse
	end
endcase
if pad_y gt 0 then print,'  Padding y-axis to ',naxis(1)+pad_y,' pixels.'

;modified naxis, accounting for analysis section and padding
naxis(0)=naxis(0)+pad_x
naxis(1)=naxis(1)+pad_y

wait,1
if not keyword_set_batch then scr_erase

;;;;;;;;;;;;;;;;;;;;;;;; GET MAIN INPUTS AND FIGURE OUT CONSTANTS

get_find_inputs, keyword_set_batch, interact,$
	 dist_mpc, use_coords, scale, cdelt, pc_pix,$
	scalehc, scalec, scaleh, scaleh_em, hs_filename,h_filename,c_filename,$
	cscale, ghc, gc, nhc, nc, sigma_find, PSFfwhm_pix, size_max_pc,$
	termgrad_arr, multiple_cutoffs, correction, NIIpercent

;;;;;;;;;;;;;;;;;;;;;;;; FIND HII REGIONS (Multi-resolution approach)

;choose the smoothed resolutions to be used
;initially increase by 10% each iteration, until this
;increment exceeds half the PSF FWHM.  Then step by 0.5*FWHM.
;only allow the resolution to grow up to 25% of the image size.

if size_max_pc lt 0. then begin
	RESfwhm_values=[PSFfwhm_pix]
	rho_crit=-1.
endif else begin
	RESfwhm_values=(PSFfwhm_pix*1.1^findgen(250))<size_max_pc/pc_pix
	RESfwhm_values=RESfwhm_values<min(naxis)/4
	RESfwhm_values=RESfwhm_values(uniq(RESfwhm_values))
	index_lastgood=where(shift(RESfwhm_values,-1)-RESfwhm_values gt $
		0.5*PSFfwhm_pix,count_lastgood)
	if count_lastgood gt 0 then begin
		RESfwhm_values2=RESfwhm_values(min(index_lastgood))+$
			(PSFfwhm_pix*findgen(250)/2.)
		RESfwhm_values=[RESfwhm_values(0:min(index_lastgood)),$
			RESfwhm_values2]
	endif
	RESfwhm_values=RESfwhm_values<size_max_pc/pc_pix
	RESfwhm_values=RESfwhm_values<(min(naxis-[pad_x,pad_y])/4)
	RESfwhm_values=RESfwhm_values(uniq(RESfwhm_values))
endelse


;set flag so that during 1st FIND iteration we make arrays,
;rather than concatenate
concatenate=0L

;cut sky section from cont-sub image and then forward FFT
hsimage_sky=hsimage(blc_sky(0):trc_sky(0),blc_sky(1):trc_sky(1)) ;supplied unit
fft_hsimage_sky=fft(hsimage_sky,-1) ;in supplied units

;grab all pixels for first try at sigma calculation, exclude blanked pixels
index=where((hsimage_sky gt -1.e12) and (sky_bimage eq 1))

;iterative determination of sky stddev
skysig=stdev(hsimage_sky(index),mean) ;in supplied units, not necess ADU
for i=0,4,1 do begin
        index=where((hsimage_sky ge (mean-3.*skysig)) and $
                 (hsimage_sky le (mean+3.*skysig)) and $
		(sky_bimage eq 1),count)
        skysig=stdev(hsimage_sky(index),mean) ;in supplied units
endfor
skymedian=median(hsimage_sky(index))	;in supplied units, not necess ADU

;use a fudgefactor if the median sky value is negative (in cont-sub img)
if skymedian lt 0. then begin
	print, 'Warning: Median sky value is negative in LINE image!'
	print, '         Adding a constant to the image in order to'
	print, '         correct for the purpose of HIIphot growing.'
	fudge=abs(skymedian)			;in supplied units
	hsimage=temporary(hsimage)+fudge	;in supplied units
	skymedian=skymedian+fudge		;in supplied units
	;cut out sky again (and FFT again) if median was negative
	hsimage_sky=hsimage(blc_sky(0):trc_sky(0),blc_sky(1):trc_sky(1))
	fft_hsimage_sky=fft(hsimage_sky,-1)	;in supplied units
endif else begin
	fudge=0.
endelse

hspadval=skymedian
hsimage=interpolate(hsimage,findgen(naxis(0))+blc(0),$
	findgen(naxis(1))+blc(1),/GRID,MISSING=hspadval);in supplied units,
if pad_x gt 0 then hsimage(trc(0)-blc(0)+1:*,*)=hspadval;not necessarily ADU
if pad_y gt 0 then hsimage(*,trc(1)-blc(1)+1:*)=hspadval;or even EM
fft_hsimage=fft(hsimage,-1)				;in supplied units

cimage_sky=cimage(blc_sky(0):trc_sky(0),blc_sky(1):trc_sky(1)) ;supplied units
fft_cimage_sky=fft(cimage_sky,-1)			       ;supplied units

cpadval=median(cimage_sky(index))
cimage=interpolate(cimage,findgen(naxis(0))+blc(0),$
	findgen(naxis(1))+blc(1),/GRID,MISSING=cpadval)	;in supplied units,
if pad_x gt 0 then cimage(trc(0)-blc(0)+1:*,*)=cpadval	;not necessarily ADU
if pad_y gt 0 then cimage(*,trc(1)-blc(1)+1:*)=cpadval

if use_bimage then bimage_float=interpolate(bimage_float,$
	findgen(naxis(0))+blc(0),$
	findgen(naxis(1))+blc(1),/GRID,MISSING=0.)
if use_bimage and (pad_x gt 0) then bimage_float(trc(0)-blc(0)+1:*,*)=0.
if use_bimage and (pad_y gt 0) then bimage_float(*,trc(1)-blc(1)+1:*)=0.
index_midrange=where(bimage_float gt 0. and bimage_float lt 1.,count_midrange)
if b_filename ne 'NULL' and count_midrange gt 0 then begin
;;;this dies for very large images, rather goes on and on forever...
;uniq_bval_float=bimage_float(uniq(bimage_float(sort(bimage_float))))
;num_uniq_bval_float=n_elements(uniq_bval_float)
read,'Number of unique exposure map values? ',num_uniq_bval_float
if num_uniq_bval_float gt 2 then begin
	print,''
	print,'Blanking image is being used as an exposure map'
	print,strcompress('with '+string(num_uniq_bval_float)+$
		' unique values.')
	print,''
endif
endif

if use_bimage then bimage=interpolate(bimage,$
	findgen(naxis(0))+blc(0),$
	findgen(naxis(1))+blc(1),/GRID,MISSING=0.)
if use_bimage and (pad_x gt 0) then bimage(trc(0)-blc(0)+1:*,*)=0.
if use_bimage and (pad_y gt 0) then bimage(*,trc(1)-blc(1)+1:*)=0.

openw,9,strcompress(root+'.noise_vs_res.dat')

; begin multi-resolution looping
toofewsky=0L
for RESfwhm_index=0,n_elements(RESfwhm_values)-1,1 do begin

RESfwhm_pix=RESfwhm_values(RESfwhm_index)
convl_mult=sqrt((RESfwhm_pix/PSFfwhm_pix)^2 - 1.)
FINDfwhm_pix=sqrt(2.*RESfwhm_pix^2-PSFfwhm_pix^2)

if (convl_mult gt 0.01) then begin

	convl_fwhm=convl_mult*PSFfwhm_pix
	unit_conv=1. ; so that images supplied to FIND have constant <flux/pix>

;could do these convolutions more efficiently with UV tapering of
;the existing FFT(hsimage) and FFT(hsimage_sky), then single reverse FFT each

;convolve the analysis section

	center=size(hsimage)/2
	numpix=N_elements(hsimage)
	fft_psf=fft(psf_Gaussian(npixel=[naxis(0),naxis(1)],$
		FWHM=[convl_fwhm,convl_fwhm],NDIMEN=2,/NORMALIZE),-1)
	product_fft=fft_psf*fft_hsimage*unit_conv*numpix
	;in supplied units, not necess ADU or EM
	hsimage=shift(float(fft(product_fft,1)),center(1),center(2))

;convolve the sky section
	center=size(hsimage_sky)/2
	numpix=N_elements(hsimage_sky)
	fft_psf=fft(psf_Gaussian(npixel=[trc_sky(0)-blc_sky(0)+1,$
		trc_sky(1)-blc_sky(1)+1],$
		FWHM=[convl_fwhm,convl_fwhm],NDIMEN=2,/NORMALIZE),-1)
	product_fft=fft_psf*fft_hsimage_sky*unit_conv*numpix
	;in supplied units, not necessarily ADU or EM
	hsimage_sky=shift(float(fft(product_fft,1)),center(1),center(2))

	product_fft=fft_psf*fft_cimage_sky*unit_conv*numpix
        ;in supplied units, not necessarily ADU or EM
        cimage_sky=shift(float(fft(product_fft,1)),center(1),center(2))
	

endif	

;grab all pixels for first try at sigma calculation
index=where((hsimage_sky gt -1.e12) and (sky_bimage eq 1),indexcount)

;first time through the multi-resolution loop, skysig won't change
;subsquent iterations will cause it to (presumably) drop.
skysig=stdev(hsimage_sky(index),mean)	;in supplied units
for i=0,4,1 do begin
        index=where((hsimage_sky ge (mean-3.*skysig)) and $
                 (hsimage_sky le (mean+3.*skysig)) and $
		(sky_bimage eq 1),count)
        skysig=stdev(hsimage_sky(index),mean) ;in supplied units
	;warn if there are too few independent pixels in the sky area
	if round(count/(RESfwhm_pix^2*1.13)) lt 10. then begin
		toofewsky=1L
	endif
endfor
if (convl_mult le 0.01) then skysig_orig=skysig	;in supplied units
                        ;not necessarily ADU or EM
                        ;determined for the line-only image

;first time through the loop, we determine cmean_orig, cskysig_orig,
;hcmean_orig, hcskysig_orig, offsetc, offsethc

index=where((cimage_sky gt -1.e12) and (sky_bimage eq 1)) 
cskysig=stdev(cimage_sky(index),cmean) ; in supplied units 
for i=0,4,1 do begin 
       	index=where((cimage_sky ge (cmean-3.*cskysig)) and $ 
        	(cimage_sky le (cmean+3.*cskysig)) and $ 
               	(sky_bimage eq 1),count) 
       	cskysig=stdev(cimage_sky(index),cmean) ; in supplied units
endfor 
cmean=cmean/scalec         ;in ADU
cskysig=cskysig/scalec     ;in ADU 
if (convl_mult le 0.01) then begin
	cmean_orig=cmean	;in ADU
	cskysig_orig=cskysig	;in ADU
endif

;recreate the line+continuum sky section in supplied units
hcimage_sky=((hsimage_sky-fudge)/scaleh+cscale*cimage_sky/scalec-$
	correction)*scalehc
index=where((hcimage_sky gt -1.e12) and (sky_bimage eq 1)) 
hcskysig=stdev(hcimage_sky(index),hcmean) ;in supplied units
for i=0,4,1 do begin 
       	index=where((hcimage_sky ge (hcmean-3.*hcskysig)) and $ 
                (hcimage_sky le (hcmean+3.*hcskysig)) and $ 
                (sky_bimage eq 1),count) 
       	hcskysig=stdev(hcimage_sky(index),hcmean)	;in supplied units
endfor 
hcmean=hcmean/scalehc         ;in ADU
hcskysig=hcskysig/scalehc     ;in ADU
if (convl_mult le 0.01) then begin
        hcmean_orig=hcmean        ;in ADU
        hcskysig_orig=hcskysig    ;in ADU

;zero-point offsets are calculated assuming sky-noise limited images
;cmean_orig-offsetc=estimated original bkgd in cimage_sky [ADU]
;hcmean_orig-offsethc=estimated original bkgd in hcimage_sky [ADU]
offsetc=cmean_orig-gc*float(nc)*cskysig_orig^2		;in ADU
offsethc=hcmean_orig-ghc*float(nhc)*hcskysig_orig^2	;in ADU

endif

;determine negative (residual) locations on the orig cont-sub image
;blank the current version of the image in these locations
;working in supplied units
if (convl_mult le 0.01) then index_toolow=where((hsimage lt $
	skymedian-6.*skysig) and (abs(hsimage) ge 1.e-12) ,count_toolow)
if count_toolow gt 0 then hsimage(index_toolow)=0.

;update the user on what is happening
print,'-------------------------------------------------------------------------------'
print, 'Locating regions with FWHM (pix) ~ ',FINDfwhm_pix,' on ',$
	RESfwhm_pix,' image.'
print, '                      FWHM (pc)  ~ ',FINDfwhm_pix*pc_pix,$
	' on ',RESfwhm_pix*pc_pix,' image.'
print,''
print,strcompress('Before rejection, sky region contains ~ '+$
	string(round(indexcount/(RESfwhm_pix^2*1.13)))+' independent points.')
;print warning to the screen if sky region is too small
if toofewsky then begin
	print,'WARNING! Sky region contains < 10 independent points!!!'
	print,'WARNING! Sky noise is likely misjudged this iteration!!!'
	wait,1.
endif
print,''

;initialize FIND variables
rlimit=[-1.e6,1.e6]
slimit=[-1.e6,1.e6]
xpos=0.
ypos=0.
fluxbm=0.
sharp=0.
round=0.
rank=0L
nreg=long(0)

sigma_tolerance_find=0.			 ;this variable can be used to
					 ;account for the fact that
					 ;the s2n in the convolution can
					 ;potentially wind up smaller than
					 ;when evaluated pixel-by-pixel...

					 ;we want a complete sample based
					 ;on the pixel-by-pixel s2n estimate,
					 ;so we might consider including
					 ;a bit more than desired
					 ;at this point, perhaps 0.25?

hmin=(sigma_find-sigma_tolerance_find)*skysig   ;set FIND threshold,
						;supplied units

printf,9,RESfwhm_pix,round(indexcount/(RESfwhm_pix^2*1.13)),skysig,cskysig,$
	(cmean_orig-offsetc)*(RESfwhm_pix^2*1.13),hcskysig,$
	(hcmean_orig-offsethc)*(RESfwhm_pix^2*1.13),$
	format='(E10.3,1X,I6,1X,5(E10.3,1X))'

if (convl_mult le 0.75) and (size_max_pc gt 0.) then begin
;do centroid check during FIND for unconvolved image and convolved images
;up to a 25% increase in beamsize [ note: 0.75 is not derived as 1.0-0.25,
;but instead as sqrt((1.25/1.)^2 - 1.) ]
	find_mod,hsimage,xpos,ypos,fluxbm,sharp,round,hmin,$
		FINDfwhm_pix,rlimit,slimit,/SILENT,/VERBOSE,POINTSOURCE=1
;;;	find_mod,-1.*hsimage,xpos2,ypos2,fluxbm2,sharp2,round2,hmin,$
;;;		FINDfwhm_pix,rlimit,slimit,/SILENT,/VERBOSE,POINTSOURCE=1
endif else begin
;don't do centroid check during FIND when using strongly convolved images,
;the "derivative-type" centroider often fails due to oversampling
	find_mod,hsimage,xpos,ypos,fluxbm,sharp,round,hmin,$
		FINDfwhm_pix,rlimit,slimit,/SILENT,/VERBOSE,POINTSOURCE=0
;;;	find_mod,-1.*hsimage,xpos2,ypos2,fluxbm2,sharp2,round2,hmin,$
;;;		FINDfwhm_pix,rlimit,slimit,/SILENT,/VERBOSE,POINTSOURCE=1
endelse

;;;tmpimage=-1.*hsimage
;;;tmpimage(xpos,ypos)=min(tmpimage)
;;;tmpimage(xpos2,ypos2)=max(tmpimage)

;back to original data values (unsmoothed) for correlation coeff stuff
fxread, hs_filename, hsimage, hshead  ; read again, to save memory!
				      ; in supplied units, not necess ADU or EM
hsimage=temporary(hsimage)+fudge      ; "  "  "  "  "  "
hsimage=interpolate(hsimage,findgen(naxis(0))+blc(0),$
	findgen(naxis(1))+blc(1),/GRID,MISSING=hspadval);in supplied units
index_bad=where(finite(hsimage) eq 0,count_bad)
if count_bad gt 0 then hsimage(index_bad)=hspadval
if pad_x gt 0 then hsimage(trc(0)-blc(0)+1:*,*)=hspadval
if pad_y gt 0 then hsimage(*,trc(1)-blc(1)+1:*)=hspadval

;this detects if no sources were tabulated this iteration
nreg=long(n_elements(fluxbm)) ; the number of HII regions located
index=where(fluxbm lt 1.e-6,count)
if (count eq nreg) or (max(fluxbm) eq 0) then begin
	nreg=long(0)
	goto, NODETECTIONS
endif

s2n_convo=fluxbm/skysig

;sort the HII regions in order of decreasing S/N
rank=reverse(sort(s2n_convo))
s2n_convo=s2n_convo(rank)
RESfwhm_pix_best=replicate(RESfwhm_pix,nreg)
FINDfwhm_pix_best=replicate(FINDfwhm_pix,nreg)
if convl_mult le 0.01 then begin
	xpos=xpos(rank)
	ypos=ypos(rank)
endif else begin ; correct for shift induced by convolution with Gaussian
		 ;                       EVEN CASE                ODD CASE
	if not(naxis(0) mod 2) then xpos=xpos(rank)+0.5 else xpos=xpos(rank)+1.
	if not(naxis(1) mod 2) then ypos=ypos(rank)+0.5 else ypos=ypos(rank)+1.
endelse
;check potential effects of up-rightward shifting
index_temp=where(xpos le trc(0)-blc(0) and ypos le trc(1)-blc(1),nreg)
if nreg eq 0 then goto, NODETECTIONS
s2n_convo=s2n_convo(index_temp)
RESfwhm_pix_best=replicate(RESfwhm_pix,nreg)
FINDfwhm_pix_best=replicate(FINDfwhm_pix,nreg)
xpos=xpos(index_temp)
ypos=ypos(index_temp)

;identify and count those sources inside the edge-effect guard band
index1=where((round(xpos) le trc(0)-blc(0)-ceil(FINDfwhm_pix/2.)) and $
	(round(xpos) ge 0+ceil(FINDfwhm_pix/2.)) and $
	(round(ypos) le trc(1)-blc(1)-ceil(FINDfwhm_pix/2.)) and $
	(round(ypos) ge 0+ceil(FINDfwhm_pix/2.)),count1)
if count1 gt 0 then begin	;if any regions are in-bounds
	xpos=xpos(index1)
	ypos=ypos(index1)
	s2n_convo=s2n_convo(index1)
	RESfwhm_pix_best=RESfwhm_pix_best(index1)
	FINDfwhm_pix_best=FINDfwhm_pix_best(index1)
endif
if long(count1) lt nreg then begin	;if at least 1 source is out-of-bounds
	if long(count1) ne 0 then begin	;...and >= 1 source is in-bounds
		print,' No. of sources rejected as OUT-OF-BOUNDS     ',nreg-$
			long(count1)
	endif else begin		;otherwise all are out-of-bounds
		print,' No. of sources rejected as OUT-OF-BOUNDS     ',nreg
		nreg=long(0)
		if use_bimage then $
	print,' No. of sources rejected by MASK+ROI  criteria       -'
	print,' No. of sources rejected by RESIDUAL  criteria       -'
		goto, NODETECTIONS	
	endelse
endif
nreg=long(n_elements(s2n_convo))	;revised number of regions, after guard band

; throw away detections in masked areas, whether supplied as an image
; or determined using the analysis ROI interface
if use_bimage then begin	;identify and count regions in unmasked areas
	index=where(bimage(round(xpos),round(ypos)) ge 1.,count)
	if count gt 0 then begin	;if any survive MASK+ROI cut
		xpos=xpos(index)
		ypos=ypos(index)
		s2n_convo=s2n_convo(index)
		RESfwhm_pix_best=RESfwhm_pix_best(index)
		FINDfwhm_pix_best=FINDfwhm_pix_best(index)
		print,' No. of sources rejected by MASK+ROI  criteria',$
				nreg-long(count)
		nreg=long(n_elements(s2n_convo))	;revised number of regions
						;after MASK+ROI
	endif else begin	; all detections masked
		print,' No. of sources rejected by MASK+ROI  criteria',nreg
		print,' No. of sources rejected by RESIDUAL  criteria       -'
		nreg=long(0)
		goto, NODETECTIONS
	endelse
endif

;throw away detections within one FINDfwhm of any "toolow" pixel
reject_residuals=1L
if reject_residuals then begin
	;do the distance checking and flag xpos
	for i=long(0),long(n_elements(index_toolow))-long(1) do begin
		cur_pos=index2coord(index_toolow(i),naxis(0))
		d2=(cur_pos(0)-xpos)^2+(cur_pos(1)-ypos)^2
		index=where(d2 le FINDfwhm_pix^2,count)
		if count gt 0 then xpos(index)=-1.
	endfor

	;report the number thrown away due to residual proximity
	index=where(xpos lt 0,count)
	print,' No. of sources rejected by RESIDUAL  criteria',count

	;now update source list vectors appropriately
	index=where(xpos ge 0,count)
	if count gt 0 then begin
		xpos=xpos(index)
		ypos=ypos(index)
		s2n_convo=s2n_convo(index)
		RESfwhm_pix_best=RESfwhm_pix_best(index)
		FINDfwhm_pix_best=FINDfwhm_pix_best(index)
	endif else begin	;if all remaining regions were discarded
		goto, NODETECTIONS
	endelse
endif

;each source making it this far will be compared to stretched/rotated models

nreg=long(n_elements(s2n_convo)) ; the number of HII regions surviving until
			      ; the correlation coeff check. 
			      ; source can still be discarded for low rho
			      ; or any of the footprint criteria
print,strcompress(string(nreg)+' sources survived the FIND procedure.')

;-------------------- distorted models + correlation coeff stuff
;we adopt base morphologies described by rings and gaussians.
;rings are given various ratios of intrinsic radius to ring thickness (r0/sr).
;all models are scaled in size to agree with RESOLUTION used during find,
;gauged via testing to maximize flux/beam.  Don't change the constants
;below unless you really know what you're doing!!!

r0_values=[0.,(RESfwhm_pix/3.1)/1.536,RESfwhm_pix/3.1,$
	(RESfwhm_pix/3.1)/0.769,(RESfwhm_pix/3.1)/0.698,$
	(RESfwhm_pix/3.1)/0.677]
sr_values=[RESfwhm_pix/2.35,r0_values(1)*2.,r0_values(2)*1.,$
	r0_values(3)*0.5,r0_values(4)*0.25,r0_values(5)*0.125]

;keep track of "best-matching" modeltype, axial ratio, and orientation for
;each source.  also, tabulate alpha, beta scaling parameters and correlation
;coefficient, rho, for this match.
modeltype_best=replicate(-1,nreg)
axrat_best=replicate(-1.,nreg)
rotdeg_best=replicate(-1.,nreg)
alpha_best=replicate(-1.,nreg)
beta_best=replicate(-1.,nreg)
rho_best=replicate(-2.,nreg)

;compare the data with each morphology (when possible)
for modeltype=0,5,1 do begin

	;skip model if undersampled
	if sr_values(modeltype) ge 2.0/2.35 then begin

	;generate a symmetrical model on an odd-sized grid about 3x the
	;characteristic source size.
	subsize=2*ceil(3.*FINDfwhm_pix/2.)+1
	subhalf=subsize/2
	sym_model,subsize,r0_values(modeltype),sr_values(modeltype),$
		undistorted

	for axrat=1.0,2.01,0.25 do begin
		print,'checking r0 ='+string(r0_values(modeltype))+', sr ='+$
		string(sr_values(modeltype))+', axrat ='+string(axrat)
		if axrat gt 1. then rotmax=165. else rotmax=0.
		for rotdeg=0.0,rotmax+0.1,15.0 do begin

			;stretch and rotate symmetrical model
			distort,undistorted,distorted,axrat,rotdeg

			;prepare to compute correlation coeff, rho, using
			;only model pix above 1% isophote
			index_signif=where(distorted ge 0.01, count_signif)
			model_pix=distorted(index_signif)
                       	num_pix=float(n_elements(model_pix))

	             	;do the first calculation for this r0,sr,ax,theta
			;combination the slow way
                       	i=0L

			;pull out line-only image data, in supplied units
                       	hdata_subsection=interpolate(hsimage,$
				(xpos(i)+findgen(subsize)-subhalf),$
				(ypos(i)+findgen(subsize)-subhalf),$
				/GRID,MISSING=0.)
                       	hdata_pix=hdata_subsection(index_signif)
			;pull out continuum image data, in supplied units
                       	cdata_subsection=interpolate(cimage,$
                       	        (xpos(i)+findgen(subsize)-subhalf),$
                       	        (ypos(i)+findgen(subsize)-subhalf),$
                       	        /GRID,MISSING=0.)
                       	cdata_pix=cdata_subsection(index_signif)

                       	bimage_float_subsection=interpolate(bimage_float,$
                       	        (xpos(i)+findgen(subsize)-subhalf),$
                       	        (ypos(i)+findgen(subsize)-subhalf),$
                       	        /GRID,MISSING=0.)
                       	bimage_float_pix=bimage_float_subsection(index_signif)
			
			;estimate the formal variance in the line-only image.
			;this is not the std-dev, but an estimate of the
			;true measurement error (based on known errors in
			;line+cont and cont images).
			hcdata_pix=((hdata_pix-fudge)/scaleh+cscale*$
				cdata_pix/scalec-correction)*scalehc ;supplied
								     ;units
			hvar_pix=(scaleh/ghc*$
				sqrt(hcdata_pix/scalehc*ghc/bimage_float_pix+$
				cscale^2*cdata_pix/scalec*$
				gc/bimage_float_pix))^2

                       	Xsum=0.
                       	XXsum=0.
                       	Xbar=0.
			;for this 1st (slow) calculation, compute model sums
			;work in supplied units, so alpha,beta are meaningful
                       	corr_statistics_2d,model_pix,$
                       	        hdata_pix,hvar_pix,$
                       	        num_pix,0.,$
                       	        alpha_cur,beta_cur,rho_cur,$
                       	        Xsum,XXsum,Xbar,/SLOW
			;save model-related values
                       	Msum=Xsum
                       	MMsum=XXsum
                       	Mbar=Xbar
			
			;if this one beats all previous model/data comparisons
			;then save some info
		  	if rho_cur gt rho_best(0) then begin
					rho_best(0)=rho_cur
					alpha_best(0)=alpha_cur
					beta_best(0)=beta_cur
					modeltype_best(0)=modeltype
					axrat_best(0)=axrat
					rotdeg_best(0)=rotdeg
			endif
                       ; reuse the model sums for subsequent detections
                       for i=1L,long(nreg)-1L,1L do begin
				;pull out the line-only image data,
				;work in supplied units
				hdata_subsection=interpolate(hsimage,$
					(xpos(i)+findgen(subsize)-subhalf),$
					(ypos(i)+findgen(subsize)-subhalf),$
					/GRID,MISSING=0.)
				hdata_pix=hdata_subsection(index_signif)
				;pull out the continuum image data,
				;work in supplied units
				cdata_subsection=interpolate(cimage,$
					(xpos(i)+findgen(subsize)-subhalf),$
					(ypos(i)+findgen(subsize)-subhalf),$
					/GRID,MISSING=0.)
				cdata_pix=cdata_subsection(index_signif)

				bimage_float_subsection=interpolate($
					bimage_float,$
					(xpos(i)+findgen(subsize)-subhalf),$
					(ypos(i)+findgen(subsize)-subhalf),$
					/GRID,MISSING=0.)
				bimage_float_pix=bimage_float_subsection($
					index_signif)

				;estimate the variance in the line-only data
			hcdata_pix=((hdata_pix-fudge)/scaleh+cscale*$
				cdata_pix/scalec-correction)*scalehc ;supplied
								     ;units

			hvar_pix=(scaleh/ghc*$
				sqrt(hcdata_pix/scalehc*ghc/bimage_float_pix+$
				cscale^2*cdata_pix/scalec*$
				gc/bimage_float_pix))^2


				;quicker than before, since we have the model
				;sums already done, work in supplied units
				corr_statistics_2d,model_pix,hdata_pix,$
					hvar_pix,n_elements(model_pix),0.,$
					alpha_cur,beta_cur,rho_cur,$
					Msum,MMsum,Mbar

				;if this one beats all previous model/data
				;comparisons then save some info
				if rho_cur gt rho_best(i) then begin
					rho_best(i)=rho_cur
					alpha_best(i)=alpha_cur
					beta_best(i)=beta_cur
					modeltype_best(i)=modeltype
					axrat_best(i)=axrat
					rotdeg_best(i)=rotdeg
				endif			
			endfor
		endfor
	endfor

	endif

endfor

;--------------------create the running list of detections
if not(concatenate) then begin
	concatenate=1L
	xpos_all=xpos
	ypos_all=ypos
	s2n_convo_all=s2n_convo
	RESfwhm_pix_best_all=RESfwhm_pix_best
	FINDfwhm_pix_best_all=FINDfwhm_pix_best
	rho_best_all=rho_best
	alpha_best_all=alpha_best
	beta_best_all=beta_best
	modeltype_best_all=modeltype_best
	axrat_best_all=axrat_best
	rotdeg_best_all=rotdeg_best
endif else begin
	xpos_all=[xpos_all,xpos]
	ypos_all=[ypos_all,ypos]
	s2n_convo_all=[s2n_convo_all,s2n_convo]
	RESfwhm_pix_best_all=[RESfwhm_pix_best_all,RESfwhm_pix_best]
	FINDfwhm_pix_best_all=[FINDfwhm_pix_best_all,FINDfwhm_pix_best]
	rho_best_all=[rho_best_all,rho_best]
	alpha_best_all=[alpha_best_all,alpha_best]
	beta_best_all=[beta_best_all,beta_best]
	modeltype_best_all=[modeltype_best_all,modeltype_best]
	axrat_best_all=[axrat_best_all,axrat_best]
	rotdeg_best_all=[rotdeg_best_all,rotdeg_best]
endelse

NODETECTIONS:

endfor ; end of convolving loop, possibly try another RESfwhm

close,9

skysig=skysig_orig	;in supplied units, not necessarily ADU or EM

;when absolutely no regions survive
if (n_elements(s2n_convo_all) eq 0) then begin
	print,'-------------------------------------------------------------------------------'
	print,'NO REGIONS FOUND! Try adjusting the FIND parameters...'
	print,''
	return
endif	 

;sort the composite database of HII regions in order of decreasing rho
rank_all=reverse(sort(rho_best_all))
s2n_convo_all=s2n_convo_all(rank_all)
xpos_all=xpos_all(rank_all)
ypos_all=ypos_all(rank_all)
RESfwhm_pix_best_all=RESfwhm_pix_best_all(rank_all)
FINDfwhm_pix_best_all=FINDfwhm_pix_best_all(rank_all)
rho_best_all=rho_best_all(rank_all)
alpha_best_all=alpha_best_all(rank_all)
beta_best_all=beta_best_all(rank_all)
modeltype_best_all=modeltype_best_all(rank_all)
axrat_best_all=axrat_best_all(rank_all)
rotdeg_best_all=rotdeg_best_all(rank_all)
nreg=long(size(s2n_convo_all)) & nreg=nreg(1) ; number in composite list
					      ; before footprint/seed criteria
					      ; and simultaneous rho check
r0_best_all=fltarr(nreg)
sr_best_all=fltarr(nreg)

;with all FFTs done, take away padding
naxis(0)=naxis(0)-pad_x
naxis(1)=naxis(1)-pad_y

fxread, hs_filename, hsimage, hshead  ; read again, to save memory!
hsimage=temporary(hsimage)+fudge      ; in supplied units
hsimage=interpolate(hsimage,findgen(naxis(0))+blc(0),$
	findgen(naxis(1))+blc(1),/GRID,MISSING=0.)	;in supplied units
index_bad=where(finite(hsimage) eq 0,count_bad)
if count_bad gt 0 then hsimage(index_bad)=0.

cimage=cimage(0:naxis(0)-1,0:naxis(1)-1)		;in supplied units

if use_bimage then bimage_float=bimage_float(0:naxis(0)-1,0:naxis(1)-1)
if use_bimage then bimage=bimage(0:naxis(0)-1,0:naxis(1)-1)
;this sets things up for outlining BLANKED areas in drawborders code
if use_bimage then index_bimage=where(bimage eq 0,count_bimage) $
	else index_bimage=-1

TRY_FOOTPRINTS_AGAIN:

;intialize footprint and seed integer maps
mapf=lonarr(naxis(0),naxis(1))
maps=lonarr(naxis(0),naxis(1))

;zero counters
n_rhocheck=0L
n_centroidclaimed=0L
n_badmodel=0L
n_f_remain=0L
n_s2n=0L
n_sbmedian=0L
n_mapped=0L

;initialize discard flag array
discardflag_all=s2n_convo_all-s2n_convo_all

;initialize other arrays
f_remain_all=discardflag_all
s2n_initial_all=discardflag_all
sbpeak_all=discardflag_all
sbavg_all=discardflag_all
sbmedian_all=discardflag_all
sb3quartile_all=discardflag_all
sbmin_all=discardflag_all
hbkgd_all=discardflag_all
s2n_linei_all=discardflag_all

;modified version of scaleh_em, in order to go to units of H-alpha EM
;that is, without any contribution from [NII]6548,6584
scaleh_em_mod=scaleh_em*(1.d0/(1.d0+NIIpercent/100.d0))

;place footprints/seeds into the integer maps, starting with the highest rho
;several conditions must be met: unclaimed centroid pixel, S/N, flux fraction 

for cpos=0L,long(nreg)-1L,1L do begin

;setup for model calculations appropriate to char. size of current source
;once again, don't change these constants unless you know what you're doing!!!

	;the source must exceed our critical value for rho, correlation coeff
	if rho_best_all(cpos) lt rho_crit then begin
		n_rhocheck=n_rhocheck+1
		discardflag_all(cpos)=-1.e6
		;;;print,'Rho discard'
		;;;stop
		goto,NOT_IMPORTANT
	endif		 

	r0_values=[0.,(RESfwhm_pix_best_all(cpos)/3.1)/1.536,$
		RESfwhm_pix_best_all(cpos)/3.1,$
		(RESfwhm_pix_best_all(cpos)/3.1)/0.769,$
		(RESfwhm_pix_best_all(cpos)/3.1)/0.698,$
		(RESfwhm_pix_best_all(cpos)/3.1)/0.677]
	sr_values=[RESfwhm_pix_best_all(cpos)/2.35,r0_values(1)*2.,$
		r0_values(2)*1.,$
		r0_values(3)*0.5,r0_values(4)*0.25,r0_values(5)*0.125]

	;keep track of the matching r0 and sr
	r0_best_all(cpos)=r0_values(modeltype_best_all(cpos))
	sr_best_all(cpos)=sr_values(modeltype_best_all(cpos))

	;centroid pixel can't already be claimed
	if maps(round(xpos_all(cpos)),round(ypos_all(cpos))) ne 0 then begin
		n_centroidclaimed=n_centroidclaimed+1
		discardflag_all(cpos)=-1.e6
		;;;print,'Centroid discard'
		;;;stop
		goto,NOT_IMPORTANT
	endif		 

	;recreate the best-matched model that was identified during rho calc.
	subsize=2*ceil(3.*FINDfwhm_pix_best_all(cpos)/2.)+1
	subhalf=subsize/2
	sym_model,subsize,r0_values(modeltype_best_all(cpos)),$
		sr_values(modeltype_best_all(cpos)),undistorted
	distort,undistorted,distorted,axrat_best_all(cpos),$
		rotdeg_best_all(cpos)

	;this time shift the distorted model slightly
	xshift=xpos_all(cpos)-float(round(xpos_all(cpos)))
	yshift=ypos_all(cpos)-float(round(ypos_all(cpos)))
	distorted=interpolate(temporary(distorted),findgen(subsize)-xshift,$
		findgen(subsize)-yshift,/GRID,MISSING=0.)
	
	;now we try to find out the bkgd corrected counts for significant
	;pixels in this model (at the current position), ignore any neighboring
	;detections at this point

	;pull out data section from line-only image, in supplied units
	hdata_subsection=interpolate(hsimage,$
		round(xpos_all(cpos)+findgen(subsize)-subhalf),$
		round(ypos_all(cpos)+findgen(subsize)-subhalf),$
		/GRID,MISSING=0.)

	;pull out seed map section, set missing areas to -1 not 0
	maps_subsection=interpolate(maps,$
		round(xpos_all(cpos)+findgen(subsize)-subhalf),$
		round(ypos_all(cpos)+findgen(subsize)-subhalf),$
		/GRID,MISSING=-1)
	;make a temporary map just for this source
	mapt=maps_subsection
	;find the "significant" pixels ignoring neighboring regions
	index_signif=where(distorted ge footprint_lev, count_signif)
	;flag these pixels in the temporary map
	mapt(index_signif)=-1000
	;get the values of pixels in annulus of width 2, encircling the source
	growth_set,-1001,mapt,hdata_subsection,[subsize,subsize],$
		index_annulus,values_annulus,$
		index_perim,values_perim,D_CRIT=2. ;values in supplied units
	;which pixels in this annulus aren't "missing" within hdata_subsection?
	index_present=where(abs(hdata_subsection(index_annulus)) ge 1.e-6,$
		count_present)
	;bkgd_hdata_subsection is in the supplied units of the line-only input
	if count_present gt 0 then begin
		;find median of pixels in bkgd annulus,excluding missing pix
		guess1=median(values_annulus(index_present))
		;then we try to get the 1st quartile value, as an estimate
		;of the true bkgd for this source (excluding other seeds)
		;the line-only "sky value" is a lower limit
		index_ltmedian=where(values_annulus(index_present) le guess1)
		bkgd_hdata_subsection=median(values_annulus(index_present($
			index_ltmedian)))>hspadval
	endif else begin ;else assume sky bkgd
		bkgd_hdata_subsection=hspadval
	endelse
	;which model-significant pixels in hdata_subsection aren't missing?
	index_present2=where(abs(hdata_subsection) ge 1.e-6 and $
		distorted ge footprint_lev,count_present2)
	;subtract away (at least) this upper limit on bkgd level
	hdata_subsection=hdata_subsection-bkgd_hdata_subsection ;supplied units

	;identify claimed or out-of-image pixels based on maps_subsection
	index_claimed=where((maps_subsection ne 0), count_claimed)
	;set them to zero in distorted model
	if count_claimed gt 0 then distorted(index_claimed)=0.
	;identify available significant pixels for the flux fraction test 
	index_signif=where((distorted ge footprint_lev) and $
		abs(hdata_subsection+bkgd_hdata_subsection) ge 1.e-6,$
		count_signif)	; in supplied units of line-only image
	;we must have ~ 2/3 resolution element worth of pixels avail...
	if count_signif lt ceil(2./3.*(PSFfwhm_pix^2*1.13)) then begin
		n_badmodel=n_badmodel+1
		discardflag_all(cpos)=-1.e6
		;;;print,'Available pixels discard'
		;;;stop
		goto,NOT_IMPORTANT
	endif

	;total bkgd-corrected counts in data for significant model pixels
	;ignoring other regions, working in supplied units of line-only input
	tot_orig=total(hdata_subsection(index_present2))

	;total bkgd-corrected counts in available significant pixels
	tot_signif=total(hdata_subsection(index_signif)) ;in supplied units
	;fraction of flux remaining unclaimed
	f_remain=tot_signif/tot_orig
	f_remain_all(cpos)=f_remain

	;if fraction is less than our threshold, then forget it
	if f_remain lt f_remain_crit then begin
		n_f_remain=n_f_remain+1
		discardflag_all(cpos)=-1.e6
		;;;print,'f_remain discard'
		;;;stop
		goto,NOT_IMPORTANT
	endif
	
	;recut the data sections without previous modifications
	hdata_subsection=interpolate(hsimage,$
		round(xpos_all(cpos)+findgen(subsize)-subhalf),$
		round(ypos_all(cpos)+findgen(subsize)-subhalf),$
		/GRID,MISSING=0.)
	hdata_pix=hdata_subsection(index_signif)	; in supplied units

	cdata_subsection=interpolate(cimage,$
		round(xpos_all(cpos)+findgen(subsize)-subhalf),$
		round(ypos_all(cpos)+findgen(subsize)-subhalf),$
		/GRID,MISSING=0.)
	cdata_pix=cdata_subsection(index_signif) ; in supplied units

	bimage_float_subsection=interpolate(bimage_float,$
		round(xpos_all(cpos)+findgen(subsize)-subhalf),$
		round(ypos_all(cpos)+findgen(subsize)-subhalf),$
		/GRID,MISSING=0.)
	bimage_float_pix=bimage_float_subsection(index_signif)

	;get the values of pixels in annulus of width 2, encircling the source
	growth_set,-1001,mapt,cdata_subsection,[subsize,subsize],$
		index_annulus,values_annulus,$
		index_perim,values_perim,D_CRIT=2. ;values in supplied units
	;which pixels in this annulus aren't "missing" within cdata_subsection?
	index_present=where(abs(cdata_subsection(index_annulus)) ge 1.e-6,$
		count_present)
	;bkgd_cdata_subsection is in the supplied units of the line-only input
	if count_present gt 0 then begin
		;find median of pixels in bkgd annulus,excluding missing pix
		guess1=median(values_annulus(index_present))
		;then we try to get the 1st quartile value, as an estimate
		;of the true bkgd for this source (excluding other seeds)
		;the continuum sky level is a lower limit
		index_ltmedian=where(values_annulus(index_present) le guess1)
		bkgd_cdata_subsection=median(values_annulus(index_present($
			index_ltmedian)))>cpadval
	endif else begin ;else assume sky bkgd
		bkgd_cdata_subsection=cpadval
	endelse

	;neither of the estimated bkgd values can drop below the median sky lev
	hbkgd_pix=bkgd_hdata_subsection	;in supplied units
	cbkgd_pix=bkgd_cdata_subsection	;in supplied units

	hcdata_pix=((hdata_pix-fudge)/scaleh+cscale*cdata_pix/scalec-$
		correction)*scalehc				;supplied units
	hcbkgd_pix=((hbkgd_pix-fudge)/scaleh+cscale*cbkgd_pix/scalec-$
		correction)*scalehc				;supplied units
	nhc_pix=nhc*bimage_float_pix
	nc_pix=nc*bimage_float_pix

	;supply line+cont and cont pixel values in ADU, uncorrected for offset
	;return df, the error estimate, in ADU (to be modified below)
	estimate_error_beg,hcdata_pix/scalehc,cdata_pix/scalec,cscale,ghc,gc,$
		nhc_pix,nc_pix,scalehc,scalec,offsethc,offsetc,$
		hcbkgd_pix/scalehc,cbkgd_pix/scalec,df,hcmean,cmean

	;convert df into supplied units of the line-only image
	df=df(0)*scaleh

	;this leads to a S/N value
	s2n=total(hdata_pix-hbkgd_pix)/df
	s2n_initial_all(cpos)=s2n

	;estimate the line-only initial signal-to-noise ratio
	df_line=sqrt(total((sqrt((hdata_pix>(hcmean_orig-offsethc))/$
		(hcmean_orig-offsethc))*skysig)^2))
	s2n_linei_all(cpos)=total(hdata_pix-hbkgd_pix)/df_line

	;ALL OF THESE VARIABLES WILL CHANGE MEANING BELOW
	;WHEN THE FINAL CATALOG IS BEING GENERATED
	;store the peak surface brightness inside the footprint
	sbpeak_all(cpos)=max(hdata_pix)*scaleh_em_mod ; in EM
	;store the mean surface brightness inside the footprint
	sbavg_all(cpos)=avg(hdata_pix)*scaleh_em_mod ; in EM
	;store the median surface brightness inside the footprint
	sbmedian_all(cpos)=median(hdata_pix)*scaleh_em_mod ; in EM
	;store the 3rd quartile surface brightness inside the footprint
	index3=where(hdata_pix ge sbmedian_all(cpos)/scaleh_em_mod,count3)
	if count3 gt 0 then begin
	sb3quartile_all(cpos)=median(hdata_pix(index3))*scaleh_em_mod ; in EM
	endif else begin
	sb3quartile_all(cpos)=sbmedian_all(cpos)
	endelse
	;store the minimum surface brightness inside the footprint
	sbmin_all(cpos)=min(hdata_pix)*scaleh_em_mod	;in EM
	;store the estimated continuum-sub background for the source
	hbkgd_all(cpos)=hbkgd_pix*scaleh_em_mod ; in EM

	;if the unclaimed bkgd-corrected flux is less than some multiple of
	;the error in this measurement, then forget it

	sigma_tolerance_verify=1.	;here we make use of the adage
					;that (for instance) 5 sigma
					;is appropriate for most cases,
					;but 3 sigma is fine when you already
					;know "where to look".

					;that is, sigma_tolerance_verify
					;permits sources originally detected
					;in the FIND steps at or above
					;sigma_find-sigma_tolerance_find,
					;to drop down in 
					;significance when evaluated on
					;a pixel-by-pixel basis.  In the
					;default case, 5 sigma is used for
					;finding sources, but 4 sigma is
					;allowed during verification.

					;recall that we are more sensitive
					;to large structures when using
					;the low resolution data during find.

	if s2n lt (sigma_find-sigma_tolerance_verify) then begin
		n_s2n=n_s2n+1
		discardflag_all(cpos)=-1.e6
		;;;print,'S/N discard'
		;;;stop
		goto,NOT_IMPORTANT
	endif

	;added after v5.0 to remove large waffle-like low-res detections...
	;generally not real, and associated with some bright image defect
	;factor 1.25 determined empirically... possible consequences for
	;well-sampled low-S.B. point sources?  
	if sbmedian_all(cpos) lt 1.25*skysig+skymedian then begin 
		n_sbmedian=n_sbmedian+1
		discardflag_all(cpos)=-1.e6
		;;;print,'Median S.B. discard'
		;;;stop
		goto,NOT_IMPORTANT
	endif

	;figure out where the unclaimed significant pixels for this source
	;lie within the footprint map
	mapindex_signif=replicate(0L,count_signif)
	for nupdate=0,count_signif-1,1 do begin
		subcoord=index2coord(index_signif(nupdate),subsize)
		mapcoord=coord2coord(subcoord,[round(xpos_all(cpos))-subhalf,$
			round(ypos_all(cpos))-subhalf])
		mapindex_signif(nupdate)=coord2index(mapcoord,naxis(0))
	endfor

	;place footprint into map
	mapf(mapindex_signif)=n_mapped+1

	;pull out continuum-sub values in this region, one last time
	valuesf=hsimage(mapindex_signif)	;supplied units
	medf=median(valuesf)			;supplied units
	minf=hbkgd_pix				;supplied units
	;identify those pixels in the footprint above the seed cutoff
	index_seedpix=where(valuesf ge minf+seed_lev*(medf-minf),count_seedpix)

	;place seed into map
	if count_seedpix gt 0 then maps(mapindex_signif($
		index_seedpix))=n_mapped+1

	;increment counter for number of regions placed
	n_mapped=n_mapped+1

NOT_IMPORTANT:
endfor

print,'-------------------------------------------------------------------------------'
;count the remaining regions
index=where(discardflag_all gt -1.e5,count)
print,nreg-long(count),$
	' multiple detections or insignificant sources eliminated.'

;might eventually want to allow real time adjustment of rho_crit,
;f_remain_crit, sigma_find, footprint_lev, and seed_lev --
;in which case we wouldn't want to actually modify the vectors
;(as done immediately below) until the user had a chance to look
;at the image with foot and seed borders marked, change their mind,
;twiddle parameters, try again, etc...

;;;save,filename=strcompress(root+'.findvar.dat',/REMOVE_ALL)

if count eq 0 then begin
	print,'No significant sources remain...'
	stop
endif

;save some info about the final list
xpos_all=xpos_all(index)
ypos_all=ypos_all(index)
s2n_convo_all=s2n_convo_all(index)
RESfwhm_pix_best_all=RESfwhm_pix_best_all(index)
FINDfwhm_pix_best_all=FINDfwhm_pix_best_all(index)
rho_best_all=rho_best_all(index)
alpha_best_all=alpha_best_all(index)
beta_best_all=beta_best_all(index)
modeltype_best_all=modeltype_best_all(index)
axrat_best_all=axrat_best_all(index)
rotdeg_best_all=rotdeg_best_all(index)
f_remain_all=f_remain_all(index)
s2n_initial_all=s2n_initial_all(index)
sbpeak_all=sbpeak_all(index)
sbavg_all=sbavg_all(index)
sbmedian_all=sbmedian_all(index)
sb3quartile_all=sb3quartile_all(index)
sbmin_all=sbmin_all(index)
hbkgd_all=hbkgd_all(index)
s2n_linei_all=s2n_linei_all(index)
r0_best_all=r0_best_all(index)
sr_best_all=sr_best_all(index)


nreg=long(size(s2n_convo_all)) & nreg=nreg(1) ; the number of HII regions
					   ; in the output catalog
print,nreg,' independent regions remain.'
print,'-------------------------------------------------------------------------------'

;recreate index_toolow (after naxis change above) and skysig
index_toolow=where((hsimage lt skymedian-6.*skysig) and $
	(abs(hsimage) ge 1.e-12),count_toolow)
if count_toolow gt 0 then hsimage(index_toolow)=0.

;output a little info on the characteristic sizes of all detections
for RESfwhm_index=0,n_elements(RESfwhm_values)-1,1 do begin
	RESfwhm_pix=RESfwhm_values(RESfwhm_index)
	tempindex=where(abs(RESfwhm_pix_best_all-RESfwhm_pix) le 1.e-3,$
		tempcount)
	print,tempcount,strcompress(' regions best matched by '+$
		'characteristic FWHM ~ '+string(round(RESfwhm_pix*pc_pix))+$
		' pc')
endfor

print,'-------------------------------------------------------------------------------'

print, 'Creating EM image with footprints marked...' ; in EM units
bigimage=interpolate(hsimage*scaleh_em_mod,findgen(naxis(0)*2+1)/2.-0.5,$
	findgen(naxis(1)*2+1)/2.-0.5,/GRID,MISSING=(skymedian-3.*skysig)*$
	scaleh_em_mod)
drawborders,(skymedian-3.*skysig)*scaleh_em_mod,index_bimage,mapf,naxis,bigimage
if display then imgroam,bigimage ;in EM units
fxhmake,head,float(bigimage)
if use_coords then putast,head,astrometry_big,EQUINOX=equinox_img
fxwrite,strcompress(root+'.FOOT+BORDER.fits',/REMOVE_ALL),head,float(bigimage)
fxhmake,head,mapf
if use_coords then putast,head,astrometry,EQUINOX=equinox_img
fxwrite,strcompress(root+'.FOOT.MAP.fits',/REMOVE_ALL),head,mapf

print, 'Creating EM image with seeds marked...'	;in EM units
bigimage=interpolate(hsimage*scaleh_em_mod,findgen(naxis(0)*2+1)/2.-0.5,$
	findgen(naxis(1)*2+1)/2.-0.5,/GRID,MISSING=(skymedian-3.*skysig)*$
	scaleh_em_mod)
drawborders,(skymedian-3.*skysig)*scaleh_em_mod,index_bimage,maps,naxis,bigimage
if display then imgroam,bigimage ;in EM units
fxhmake,head,float(bigimage)
if use_coords then putast,head,astrometry_big,EQUINOX=equinox_img
fxwrite,strcompress(root+'.SEED+BORDER.fits',/REMOVE_ALL),head,float(bigimage)
fxhmake,head,maps
if use_coords then putast,head,astrometry,EQUINOX=equinox_img
fxwrite,strcompress(root+'.SEED.MAP.fits',/REMOVE_ALL),head,maps

;;;;;;;;;;;;;;;;;;;;;;;;; SETUP FOR THE GROWING PROCESS

;this should keep "too low" residuals from being grown into...
if count_toolow gt 0 then maps(index_toolow)=2^15
;do a similar thing with the bimage...
;;;if n_elements(index_bimage) gt 0 then maps(index_bimage)=2^15

start_lev=fltarr(nreg)
niter=lonarr(nreg)
maxiter=300 ;assume no single region will ever need more than 300 iterations
ncount=lonarr(nreg,maxiter)
deltan=lonarr(nreg,maxiter)
growval_hist=fltarr(nreg,maxiter)
done=lonarr(nreg)-1
termgrad_all=replicate(-1.0,nreg)
alldone=0L
itnum=1L

growth=1L

;rewrite the seed map with "too low" and blanked pixels flagged as 2^15
fxhmake,head,maps
if use_coords then putast,head,astrometry,EQUINOX=equinox_img
fxwrite,strcompress(root+'.SEED.MAP.fits',/REMOVE_ALL),head,maps

;remove pre-existing output catalogs
spawncmd='\rm '+strcompress(root+'.catalog.dat',/REMOVE_ALL)
spawn,spawncmd,spawnout
spawncmd='\rm '+strcompress(root+'.?.catalog.dat',/REMOVE_ALL)
spawn,spawncmd,spawnout

;save the session at this point, so we can loop back EXACTLY in order
;to use multiple growth cutoffs
save,file=strcompress(root+'.pre_growth.dat',/REMOVE_ALL)
CUTOFF_LOOP:
;restore session from disk
restore,file=strcompress(root+'.pre_growth.dat',/REMOVE_ALL)

;figure out which termgrad to use
termgrad=termgrad_arr(0)
root_suffix='.a'
if exist(strcompress(root+root_suffix+'.catalog.dat',/REMOVE_ALL)) then begin
	termgrad=termgrad_arr(1)
	root_suffix='.b'
endif
if exist(strcompress(root+root_suffix+'.catalog.dat',/REMOVE_ALL)) then begin
	termgrad=termgrad_arr(2)
	root_suffix='.c'
endif
if exist(strcompress(root+root_suffix+'.catalog.dat',/REMOVE_ALL)) then begin
	termgrad=termgrad_arr(3)
	root_suffix='.d'
endif
if exist(strcompress(root+root_suffix+'.catalog.dat',/REMOVE_ALL)) then begin
	termgrad=termgrad_arr(4)
	root_suffix='.e'
endif

;decide if we are finally done (with all desired growth cutoffs)
if (termgrad le 1.e-6) or (exist(strcompress($
	root+'.e.catalog.dat',/REMOVE_ALL))) then goto, END_OF_CODE

;;;;;;;;;;;;;;;;;;;;;;;;;;; GROWTH OF REGIONS

;decide which sources will be able to grow, and when they will start...
for i=long(0),nreg-long(1),long(1) do begin  
; get an index pointing to pixels that can potentially be added to
; the current region, also grab the corresponding image values
; we take the median of these values as the "starting level" associated
; with the region, it will not be allowed to start growth until the
; level of pixels being added to other regions drops to this level
	growth_set,i,maps,hsimage*scaleh_em_mod,naxis,index_grow,values_grow,$
		index_perim,values_perim
	size_grow=size(index_grow)
	any_avail=size_grow(0)
	if any_avail then begin
		start_lev(i)=median(values_grow,/EVEN)
	endif else begin ;else use the sky level
		start_lev(i)=hspadval*scaleh_em_mod
		;???might also consider flagging "placeholders" here...
		;by placeholders, I mean semi-significant detections
	endelse
endfor 

;this variable will act as a slowly declining threshold, below which
;no pixels are allowed to augment an HII region (until the threshold drops)
level=max(start_lev,maxindex)

;"turn on" growth for the region with the highest median for neighboring pixels
done(maxindex)=0L

levelnum=0L

print,'-------------------------------------------------------------------------------'
print,'   Iteration     Growing     Started        Done   Level (EM)'

;at the beginning of the "growing" loop, done equals:
;		-1 if the region has not yet been allowed to grow
;		 0 if the region is growing
;		 1 if the region has stopped growing due to threshold

repeat begin

;check to see if any regions are allowed to start growth this iteration
;(they need to have median(values_grow) ge level for this to happen)
index_startgrowth=where((done eq -1) and $
	(start_lev ge level),count_start)

;turn on growth for regions deemed OK to start
if (count_start gt 0) then done(index_startgrowth)=0L

;for display purposes, keep track of the # growing and the # finished
index_growing=where(done eq 0,count_growing)
index_finished=where(done ge 1,count_finished)

;consider adding to any region that is currently growing
;loop through these active regions in order of decreasing rho
for i=long(0),long(nreg)-long(1),long(1) do begin

; only consider a region if it is currently growing
	if (abs(done(i)) ge 1) then goto, NOVECTOR
	
; get an index pointing to pixels that can potentially be added to
; the current region, also grab the corresponding image values.
; the GROWING flag asks for index_perim and values_perim to be filled.
	growth_set,i,maps,hsimage*scaleh_em_mod,naxis,index_grow,values_grow,$
		index_perim,values_perim,/GROWING

;the decision to continue (or stop) growing is based on comparison
;between the median value of pixels just inside and just outside the region
;boundary.  The "perim" vectors specify those pixels which do belong to
;the region but fall adjacent or diagonal to an unclaimed pixel (and
;are therefore on the edge of the region).  The "grow" vectors contain
;information regarding the unclaimed pixels adjacent or diagonal to
;these perimeter pixels.

;check whether or not any pixels are even available to be added
;if none are available, then permanently stop growth of region
	size_grow=size(index_grow)
	any_avail=size_grow(0)
	if not(any_avail) then begin
		done(i)=1L
		goto,NOVECTOR
	endif
;next we try to stop growth via "runaway" fingers...
;;;	numgrow_crit=3
;check whether or not at least numgrow_crit pixels are available to add
;if not this many are available, then permanently stop growth of region
;;;	num_grow=size_grow(1)
;;;	if num_grow lt numgrow_crit then begin
;;;		done(i)=1L
;;;		goto,NOVECTOR
;;;	endif

;only pixels above the global threshold can ever be added...
;just because count_add is zero during an iteration doesn't mean
;that particular region is done growing
	index_add=where(values_grow ge level,count_add)

;add pixels in bunches that consitute >= 1/2 of the "growing" pool
	if float(count_add)/float(n_elements(values_grow)) gt 0.5 then begin

;get the characteristic values for the perimeter and "growth" regions
	growval1=median(values_perim,/EVEN)
	growval2=median(values_grow,/EVEN)

;stop growth if the surface brightness slope is shallower
;than our gradient threshold (in EM/pc)
	if (growval1-growval2)/pc_pix lt termgrad then begin
		done(i)=1L				
		goto,NOVECTOR
	endif
	
	;if we add pixels this iteration...
	;update the terminal gradient value
		termgrad_all(i)=(growval1-growval2)/pc_pix
	;update the integer map
		maps(index_grow(index_add))=i+1
	;keep track of the number of pixels belonging to the region
		index_cur=where(maps eq i+1,count_cur)
		ncount(i,niter(i))=count_cur
	;keep track of the change in this number of pixels
		if niter(i) gt 0 then deltan(i,niter(i)-1)=$
			count_cur-ncount(i,niter(i)-1)
	;keep track of the avg of all pixels added in this iteration
		growval_hist(i,niter(i))=avg(values_grow(index_add))
	;increment the count of iterations in which pixels were added
		niter(i)=niter(i)+1

	endif

NOVECTOR:
endfor

; this controls the rate at which regions are considered for growth.
; the idea is to dig evenly into the faint parts of the image and
; allow neighboring detections to compete fairly with one another.
	level_prev=level
	level_delta=max([skysig*scaleh_em_mod*1.e-2,level-level*(10.^(-0.02))])
	level=level_prev-level_delta

; print an update to the console
	print,itnum,count_growing,count_start,count_finished,level

	itnum=itnum+1

;stop looping if all regions are "done", too many iterations happened,
;or level drops below the sky
	if ((min(done) gt 0) or ((itnum gt maxiter) or $
		(level lt skymedian*scaleh_em_mod))) then alldone=1L

endrep until alldone

index_startgrowth=where((done eq -1) and $
	(start_lev+0.*skysig ge level),count_start)
index_growing=where(done eq 0,count_growing)
;any detection which never finished growing gets termgrad=0.
if count_growing gt 0 then termgrad_all(index_growing)=0.
index_finished=where(done ge 1,count_finished)

; print an update to the console
print,'        ----',count_growing,count_start,count_finished,level
if itnum gt maxiter then print, 'Too many iterations, growth has stagnated!'

print,'-------------------------------------------------------------------------------'

;;;save,ncount,deltan,growval_hist,filename=strcompress(root+'.growvar.dat',$
;;;	/REMOVE_ALL)

;;;;;;;;;;;;;;;;;;;;;;;; WRITE "GROW" DATAPRODUCTS 
 
print,'Creating EM image with region boundaries marked...' 
if count_toolow gt 0 then maps(index_toolow)=0L
bigimage=interpolate(hsimage*scaleh_em_mod,findgen(naxis(0)*2+1)/2.-0.5,$ 
        findgen(naxis(1)*2+1)/2.-0.5,/GRID,MISSING=(skymedian-3.*skysig)*$
	scaleh_em_mod) 
drawborders,(skymedian-3.*skysig)*scaleh_em_mod,index_bimage,maps,naxis,bigimage
if count_toolow gt 0 then maps(index_toolow)=2^15 

if display then imgroam,bigimage
 
fxhmake,head,float(bigimage) 
if use_coords then putast,head,astrometry_big,EQUINOX=equinox_img
fxwrite,strcompress(root+'.GROW+BORDER.fits',/REMOVE_ALL),head,float(bigimage)
fxhmake,head,maps 
if use_coords then putast,head,astrometry,EQUINOX=equinox_img
if count_toolow gt 0 then maps(index_toolow)=0L
fxwrite,strcompress(root+'.GROW.MAP.fits',/REMOVE_ALL),head,maps 
if count_toolow gt 0 then maps(index_toolow)=2^15 

;if count_toolow gt 0 then maps(index_toolow)=0L
;bigimage=interpolate(cimage,findgen(naxis(0)*2+1)/2.-0.5,$ 
;        findgen(naxis(1)*2+1)/2.-0.5,/GRID,MISSING=(skymedian-3.*skysig)*$
;	scaleh_em_mod) 
;drawborders,(skymedian-3.*skysig),index_bimage,maps,naxis,bigimage
;if count_toolow gt 0 then maps(index_toolow)=2^15 
;fxhmake,head,float(bigimage) 
;if use_coords then putast,head,astrometry_big,EQUINOX=equinox_img
;fxwrite,strcompress(root+'.CONT.GROW+BORDER.fits',/REMOVE_ALL),head,$
;	float(bigimage)

;;;;;;;;;;;;;;;;;;;;;;;; SURFACE FITTING OF DIFFUSE BKGD

if (n_elements(bkgd_annulus_pc) ne 1) then begin
	print,''
	print,'Bkgd contributions to the HII region flux will be estimated'
	print,'using median pixel values taken from annuli surrounding all'
	print,'regions.'
	print,''
	print,'The minimum allowable width of such annuli is 5 pix.'
	print,''
	GET_BKGD_ANNULUS_PC:   
	ans = ''
read,'Annulus width for bkgd estimation in pc? (DEFAULT = 250.0 pc): ',ans 
	if ans EQ '' then bkgd_annulus_pc=250. else begin
		bkgd_annulus_pc = getopt(ans,'F')
		if N_elements(bkgd_annulus_pc) NE 1 then begin  
              		message, 'ERROR - Expecting floating-point scalar',/CON
              		goto, GET_BKGD_ANNULUS_PC     
		endif
	endelse
endif 

bkgd_rad=max([5.,bkgd_annulus_pc/pc_pix])
bkgd_rad_pc=bkgd_rad*pc_pix
print,'-------------------------------------------------------------------------------'
if bkgd_annulus_pc/pc_pix lt 5. then begin
	print,'Warning: The minimum allowable width of such annuli is 5 pix.'
	print,'         ',strcompress('This corresponds to '+$
		string(bkgd_rad_pc)+' pc.  Try again...')
	print,''
	goto,GET_BKGD_ANNULUS_PC
endif

print,'Creating diffuse-bkgd EM image with region boundaries marked...' 

num_crit=75
bkgd_spacing=round((32.>(max(naxis)/32.))<(min(naxis)/8.))>4
num_critbkgd=0.75*(acos(-1.d0)*(float(bkgd_rad)/2.)^2)
num_critgrid=0.75*(acos(-1.d0)*(float(bkgd_spacing)/2.)^2)

concatenate=0L
for i=0L,long(nreg)-1L,1L do begin
	;get the indices and values of all pixels within bkgd_rad pix
	;of region i -- hsimage is in user-supplied units
	growth_set,i,maps,hsimage,naxis,index_diffuse_reg,$
		values_diffuse_reg,index_perim,values_perim,D_CRIT=bkgd_rad
	
	;if region i has at least num_crit associated bkgd pixels
	;then it gets used for surface fitting, otherwise rely on
	;"grid" points.  this is needed if we are to avoid winding
	;up with bright "spikes" in the surface near detections without
	;a sufficient number of bkgd pixels to get multiple (independent)
	;control points 
	if n_elements(index_diffuse_reg) ge num_crit then begin

		;use at least 3 control points for each qualified region.
		;if lots of bkgd pixels are available, use about 1 control
		;for each num_crit/1.5 bkgd pixels -- 400 -> 8 ,500 -> 10,etc. 
		npoints=max([3,floor(n_elements(index_diffuse_reg)/$
			(float(num_crit)/1.5))])

		;compute the position angle and distance toward each bkgd pixel
		;based on the region center as origin
		bkgd_pa=fltarr(n_elements(index_diffuse_reg))
		bkgd_dist=fltarr(n_elements(index_diffuse_reg))
		for j=0L,n_elements(index_diffuse_reg)-1L do begin
			coord_cur=index2coord(index_diffuse_reg(j),naxis(0))
			coord_diff=coord_cur-[xpos_all(i),ypos_all(i)]
			coord_dist=sqrt(total(coord_diff^2))
			if asin(coord_diff(0)/coord_dist)*$
				180./acos(-1.) ge 0. then begin
				bkgd_pa(j)=acos(coord_diff(1)/coord_dist)*$
					180./acos(-1.)
			endif else begin
				bkgd_pa(j)=360.-acos(coord_diff(1)/$
					coord_dist)*180./acos(-1.)
			endelse
			bkgd_dist(j)=coord_dist
		endfor

		;sort the bkgd pixels into ascending position angle
		sort_bkgd_pa=sort(bkgd_pa)
		bkgd_pa=bkgd_pa(sort_bkgd_pa)
		bkgd_dist=bkgd_dist(sort_bkgd_pa)
		index_diffuse_reg=index_diffuse_reg(sort_bkgd_pa)

		;start deciding which position angles to try selecting
		;control points from... use 10 deg intervals
		pa_distrib=histogram(bkgd_pa,min=0,max=350,bin=10)
		pa_bincenter=findgen(n_elements(pa_distrib))*10.+5.
		;this vector will hold the bincenter values for each
		;target position angle
		pa_targets=fltarr(npoints)
		;sort the bins into order of decreasing # of bkgd points
		index_pa_sort=reverse(sort(pa_distrib))
		pa_distrib=pa_distrib(index_pa_sort)
		pa_bincenter=pa_bincenter(index_pa_sort)
		;we don't want to include any bins without bkgd points
		index_temp=where(pa_distrib gt 0)
		pa_distrib=pa_distrib(index_temp)
		pa_bincenter=pa_bincenter(index_temp)

		;the 1st target position angle is easy, just use the most
		;populated bin...
		pa_targets(0)=pa_bincenter(0)

		ntargets=1L
		index_cur=1L

		;try to select other target angles as evenly as possible
		repeat begin

			;begin determining the mininum separation between
			;the current bincenter being considered and all the
			;target position angles already accepted
			dist_interbin=pa_targets(0:ntargets-1)-$
				pa_bincenter(index_cur)
			index_temp=where(abs(dist_interbin) le 180.,count_temp)
			if count_temp gt 0 then begin
				dist1=min(abs(dist_interbin(index_temp)))
			endif else begin
				dist1=360.
			endelse

			dist_interbin=360.+pa_bincenter(index_cur)-$
				pa_targets(0:ntargets-1)
			index_temp=where(abs(dist_interbin) le 180.,count_temp)
			if count_temp gt 0 then begin
				dist2=min(abs(dist_interbin(index_temp)))
			endif else begin
				dist2=360.
			endelse
			
			;this is the minimum separation referred to above
			dist_min=min(dist1,dist2)

			;if this bincenter is removed sufficiently (in angle)
			;from all target position angles already accepted,
			;then we include it as well.
			if dist_min ge (360./npoints)/1.5 then begin
				pa_targets(ntargets)=pa_bincenter(index_cur)
				ntargets=ntargets+1
			endif 

			index_cur=index_cur+1

		endrep until (index_cur+1 gt n_elements(pa_distrib)) or $
			(ntargets eq npoints)

		;if we didn't get all the control points we hoped for...
		pa_targets=pa_targets(0:ntargets-1)
		npoints=ntargets

		index_diffuse_sub=lonarr(ntargets)

		;now choose the bkgd point having a typical (median)
		;distance to the source position wheen compared to
		;all the bkgd pixels in the slice of the annulus at
		;any particular target position angle.
		;this effectively tries to place the control points
		;midway out in the annulus...
		for j=0L,long(ntargets)-1L do begin
			index_temp=where((bkgd_pa ge pa_targets(j)-5.) and $
				(bkgd_pa lt pa_targets(j)+5.),count_temp)
			if count_temp gt 0 then begin
				median_distval=median(bkgd_dist(index_temp))
				index_choice=where(abs(bkgd_dist(index_temp)-$
					median_distval) le 1.e-6,count_choice)
				if count_choice eq 0 then stop
			endif else begin
				stop
			endelse
			index_choice=index_choice(0)
			index_diffuse_sub(j)=index_diffuse_reg(index_temp($
				index_choice))
		endfor

		;create/modify the running list of control points and
		;all bkgd points...
		if not(concatenate) then begin
			concatenate=1L
			index_diffuse_fit=index_diffuse_sub
			index_diffuse=index_diffuse_reg
		endif else begin
			index_diffuse_fit=[index_diffuse_fit,index_diffuse_sub]
			index_diffuse=[index_diffuse,index_diffuse_reg]
		endelse
	
	endif

endfor

;remove any duplicate control points
unsorted_index=index_diffuse_fit
sort_order=sort(unsorted_index)
sorted_index=unsorted_index(sort_order)
index_diffuse_fit=sorted_index(uniq(sorted_index))

;remove any duplicate bkgd points
unsorted_index=index_diffuse
sort_order=sort(unsorted_index)
sorted_index=unsorted_index(sort_order)
index_diffuse=sorted_index(uniq(sorted_index))

cur_pos=lonarr(2)

;get the x,y position for each control point
count_diffuse_fit=n_elements(index_diffuse_fit)
x_diffuse_fit=lonarr(count_diffuse_fit)
y_diffuse_fit=lonarr(count_diffuse_fit)
for i=0L,long(count_diffuse_fit)-1L,1L do begin
	cur_pos=index2coord(index_diffuse_fit(i),naxis(0))
	x_diffuse_fit(i)=cur_pos(0)
	y_diffuse_fit(i)=cur_pos(1)
endfor

;get the x,y position and image value at all bkgd points
count_diffuse=n_elements(index_diffuse)
x_diffuse=lonarr(count_diffuse)
y_diffuse=lonarr(count_diffuse)
z_diffuse=fltarr(count_diffuse)
for i=0L,long(count_diffuse)-1L,1L do begin
	cur_pos=index2coord(index_diffuse(i),naxis(0))
	x_diffuse(i)=cur_pos(0)
	y_diffuse(i)=cur_pos(1)
	z_diffuse(i)=(hsimage(x_diffuse(i),y_diffuse(i))*$
		scaleh_em_mod)>((hspadval-1.*skysig)*scaleh_em_mod)
endfor

print,count_diffuse,' available pixels were found in bkgd annuli.'
print,count_diffuse_fit,' control points chosen from these bkgd pixels.'

z_diffuse_fit=fltarr(count_diffuse_fit)
dz_diffuse_fit=fltarr(count_diffuse_fit)
for i=0L,long(count_diffuse_fit)-1L,1L do begin
	if z_diffuse_fit(i) gt -1.e5 then begin
		d2=(float(x_diffuse_fit(i)-x_diffuse))^2+$
			(float(y_diffuse_fit(i)-y_diffuse))^2
		index_vicinity=where(d2 le bkgd_rad^2,count_vicinity)

		;if not the last one, then check to see if any neighboring
		;control points are too close
		if i lt long(count_diffuse_fit)-1L then begin
			d2aux=(float(x_diffuse_fit(i)-$
				x_diffuse_fit(i+1:*)))^2+$
				(float(y_diffuse_fit(i)-$
				y_diffuse_fit(i+1:*)))^2
			;this throws away control points that are too close to
			;each other... one of each set is kept.
			index_vicinityaux=where(d2aux le (0.5*bkgd_rad)^2,$
				count_vicinityaux)
			if count_vicinityaux gt 0 then begin
				z_diffuse_fit(index_vicinityaux)=-1.e6
			endif
		endif

		if (count_vicinity ge num_critbkgd) then begin
			z_diffuse_fit(i)=median(z_diffuse(index_vicinity))
			result_moment=moment(z_diffuse(index_vicinity),$
				SDEV=result_sdev)
			dz_diffuse_fit(i)=result_sdev
		endif else begin
			z_diffuse_fit(i)=-1.e6
		endelse
	endif
endfor

index_goodmedian=where(z_diffuse_fit gt -1.e5,count_goodmedian)

if count_goodmedian gt 0 then begin
	print,count_goodmedian,' control points will be used for surface fitting.'
	x_diffuse_fit=x_diffuse_fit(index_goodmedian)
	y_diffuse_fit=y_diffuse_fit(index_goodmedian)
	z_diffuse_fit=z_diffuse_fit(index_goodmedian)
	dz_diffuse_fit=dz_diffuse_fit(index_goodmedian)
endif else begin
	print,'No suitable control points found for bkgd surface fitting!'
	stop
endelse


x_diffuse_grid=bkgd_spacing*lindgen(floor(naxis(0)/float(bkgd_spacing))+2)
y_diffuse_grid=bkgd_spacing*lindgen(floor(naxis(1)/float(bkgd_spacing))+2)
x_diffuse_grid=x_diffuse_grid<(naxis(0)-1)
y_diffuse_grid=y_diffuse_grid<(naxis(1)-1)
x_diffuse_grid=x_diffuse_grid(uniq(x_diffuse_grid))
y_diffuse_grid=y_diffuse_grid(uniq(y_diffuse_grid))

x_diffuse_grid_fit=lonarr(n_elements(x_diffuse_grid)*$
	n_elements(y_diffuse_grid))
y_diffuse_grid_fit=lonarr(n_elements(x_diffuse_grid)*$
	n_elements(y_diffuse_grid))
z_diffuse_grid_fit=fltarr(n_elements(x_diffuse_grid)*$
	n_elements(y_diffuse_grid))
dz_diffuse_grid_fit=fltarr(n_elements(x_diffuse_grid)*$
	n_elements(y_diffuse_grid))

;this temporary change to maps will insure that no grid points
;are allowed inside bkgd annuli
for i=0L,n_elements(x_diffuse)-1L,1L do begin
	maps(x_diffuse(i),y_diffuse(i))=2^13
endfor

k=0L
for i=0L,n_elements(x_diffuse_grid)-1L,1L do begin
	for j=0L,n_elements(y_diffuse_grid)-1L,1L do begin
		x_diffuse_grid_fit(k)=x_diffuse_grid(i)
		y_diffuse_grid_fit(k)=y_diffuse_grid(j)
		if (i eq 0) or (j eq 0) or $
			(i eq n_elements(x_diffuse_grid)-1L) or $
			(j eq n_elements(y_diffuse_grid)-1L) then begin
			factor=0.2
		endif else begin
			factor=1.0
		endelse
		if maps(x_diffuse_grid_fit(k),$
			y_diffuse_grid_fit(k)) eq 0 then begin

			maps(x_diffuse_grid_fit(k),$
				y_diffuse_grid_fit(k))=nreg+1L
			growth_set,nreg,maps,hsimage*scaleh_em_mod,naxis,$
				index_vicinity,values_vicinity,$
				index_perim,values_perim,D_CRIT=float($
				bkgd_spacing)/1.5 ; use a bit more than
						  ; a circle of radius
						  ; bkgd_spacing, so that
						  ; bkgd regions don't have
						  ; wavy appearance - just
						  ; encourage dependency...
			maps(x_diffuse_grid_fit(k),$
				y_diffuse_grid_fit(k))=0L

			values_vicinity=values_vicinity>((hspadval-3.*skysig)*scaleh_em_mod)
			count_vicinity=n_elements(index_vicinity)

			;grid control points must have at least num_critgrid
			;unclaimed pixels within radius of 1/2*bkgd_spacing,
			;unless they are on the edge then they only need
			;20% of this number
			if count_vicinity ge factor*num_critgrid then begin
				z_diffuse_grid_fit(k)=median(values_vicinity)
				result_moment=moment(values_vicinity,$
					SDEV=result_sdev)
				dz_diffuse_grid_fit(k)=result_sdev
			endif else begin
				z_diffuse_grid_fit(k)=-1.e6
			endelse

		endif else begin
			z_diffuse_grid_fit(k)=-1.e6
		endelse
		k=k+1L
	endfor
endfor

;now we revert maps back to its earlier state
for i=0L,n_elements(x_diffuse)-1L,1L do begin
	maps(x_diffuse(i),y_diffuse(i))=0
endfor

index_goodmedian=where(z_diffuse_grid_fit gt -1.e5,count_goodmedian)

if count_goodmedian gt 0 then begin
	print,count_goodmedian,' grid points will be used for surface fitting.'
	x_diffuse_grid_fit=x_diffuse_grid_fit(index_goodmedian)
	y_diffuse_grid_fit=y_diffuse_grid_fit(index_goodmedian)
	z_diffuse_grid_fit=z_diffuse_grid_fit(index_goodmedian)
	dz_diffuse_grid_fit=dz_diffuse_grid_fit(index_goodmedian)
endif else begin
	print,'No suitable grid points found for bkgd surface fitting!'
	stop
endelse

x_diffuse_fit=[x_diffuse_fit,x_diffuse_grid_fit]
y_diffuse_fit=[y_diffuse_fit,y_diffuse_grid_fit]
z_diffuse_fit=[z_diffuse_fit,z_diffuse_grid_fit]
dz_diffuse_fit=[dz_diffuse_fit,dz_diffuse_grid_fit]

gs=[1,1]
limits=[0,0,naxis(0)-1,naxis(1)-1]

;setup for the surface computation
triangulate,x_diffuse_fit,y_diffuse_fit,tr,b

;this is where the diffuse bkgd and diffuse bkgd error surfaces are computed
hsimage_diffuse=trigrid(x_diffuse_fit,y_diffuse_fit,z_diffuse_fit,$
	tr,gs,limits,MISSING=hspadval*scaleh_em_mod)
index_origimg=where(maps le 0 or maps gt long(nreg),count_origimg)
hsimage_diffuse(index_origimg)=hsimage(index_origimg)*scaleh_em_mod
hsimage_diffuse_err=trigrid(x_diffuse_fit,y_diffuse_fit,dz_diffuse_fit,$
	tr,gs,limits,MISSING=skysig*scaleh_em_mod)

bigimage=interpolate(hsimage_diffuse,$
	findgen(naxis(0)*2+1)/2.-0.5,findgen(naxis(1)*2+1)/2.-0.5,$
	/GRID,MISSING=(skymedian-3.*skysig)*scaleh_em_mod) 
if count_toolow gt 0 then maps(index_toolow)=0L
drawborders,(skymedian-3.*skysig)*scaleh_em_mod,index_bimage,maps,naxis,bigimage
if count_toolow gt 0 then maps(index_toolow)=2^15 

if display then imgroam,bigimage
fxhmake,head,float(bigimage) 
if use_coords then putast,head,astrometry_big,EQUINOX=equinox_img
fxwrite,strcompress(root+'.BKGD+BORDER.fits',/REMOVE_ALL),head,float(bigimage)

;write the bkgd image without borders marked 
fxhmake,head,float(hsimage_diffuse) 
if use_coords then putast,head,astrometry,EQUINOX=equinox_img
fxwrite,strcompress(root+'.BKGD.fits',/REMOVE_ALL),head,float(hsimage_diffuse)

print,'-------------------------------------------------------------------------------'

;;;;;;;;;;;;;;;;;;;;;;;; DEFINE OUTPUT CATALOG ARRAYS

id=lindgen(nreg)+1
comment=strarr(nreg)
alpha=fltarr(nreg)
delta=fltarr(nreg)
n_pix=lonarr(nreg)

fwhm_eff=2.*r0_best_all+2.35*sr_best_all
fwhm_maj=fwhm_eff*sqrt(axrat_best_all)
fwhm_min=fwhm_maj/axrat_best_all

pa_deg=rotdeg_best_all
;here we account for intrinsic rotation of the image
pa_deg=pa_deg-crota(1)
;force all PA measurements into the range [0,180)
index_temp=where(pa_deg ge 180.,count_temp)
if count_temp gt 0 then pa_deg(index_temp)=pa_deg(index_temp)-180.
index_temp=where(pa_deg lt 0.,count_temp)
if count_temp gt 0 then pa_deg(index_temp)=pa_deg(index_temp)+180.

flux_uncorr=fltarr(nreg)
flux_bkgd=fltarr(nreg)
flux_corr=fltarr(nreg)
lum_corr=dblarr(nreg)
delta_flux_corr=fltarr(nreg)
s2n_linef_all=fltarr(nreg)
flux_corrppp=fltarr(nreg)

em2flux=(2.d-18)/1.d0 ; (erg cm^-2 s^-1 arcsec^-2) / (cm^-6 pc) at 10^4 K

for i=0L,long(nreg)-1L,1L do begin

	case modeltype_best_all(i) of
		0:	comment(i)='G '
		1:	comment(i)='S1'
		2:	comment(i)='S2'
		3:	comment(i)='S3'
		4:	comment(i)='S4'
		5:	comment(i)='S5'
	endcase

	index_cur=where(maps eq i+1,count_cur)
	if count_cur gt 0 then begin

		n_pix(i)=count_cur

		hdata_pix=hsimage(index_cur)			;supplied units
		cdata_pix=cimage(index_cur)			;supplied units
		bimage_float_pix=bimage_float(index_cur)
		hcdata_pix=((hdata_pix-fudge)/scaleh+cscale*$
			cdata_pix/scalec-correction)*scalehc	;supplied units
		nhc_pix=nhc*bimage_float_pix
		nc_pix=nc*bimage_float_pix

        	;ALL OF THESE VARIABLES WERE ORIGINALLY DEFINED ABOVE
        	;IN THE CONTEXT OF THE FOOTPRINT FOR EACH SOURCE
        	;store the peak surface brightness inside the region
        	sbpeak_all(i)=max(hdata_pix)*scaleh_em_mod ; in EM
        	;store the mean surface brightness inside the region
        	sbavg_all(i)=avg(hdata_pix)*scaleh_em_mod ; in EM
        	;store the median surface brightness inside the region
        	sbmedian_all(i)=median(hdata_pix)*scaleh_em_mod ; in EM
        	;store the 3rd quartile surface brightness inside the region
        	index3=where(hdata_pix ge sbmedian_all(i)/scaleh_em_mod,count3)
		if count3 gt 0 then begin
        	sb3quartile_all(i)=median(hdata_pix(index3))*scaleh_em_mod ; in EM
		endif else begin
		sb3quartile_all(i)=sbmedian_all(i)
		endelse
		;store the minimum surface brightness inside the region
		sbmin_all(i)=min(hdata_pix)*scaleh_em_mod	; in EM
        	;store the mean continuum-sub background for the source
        	hbkgd_all(i)=avg(hsimage_diffuse(index_cur)) ; in EM

		;get index for brightest pixels (used for PPP method)
		pppval=0.7
		index_ppp=where(hdata_pix*scaleh_em_mod ge hbkgd_all(i)+$
			pppval*(sbpeak_all(i)-hbkgd_all(i)),count_ppp)
		if count_ppp lt 1 then begin ; just in case...
 			print,'No qualified PPP peak pixel.'
			stop
		end

		;to get flux- first convert to EM, then to erg/cm^2/s/arcsec^2,
		;then finally to erg/cm^2/s/pix, and then sum
		flux_uncorr(i)=total(hdata_pix*scaleh_em_mod*em2flux*scale^2)

		;hsimage_diffuse is already in EM units
		flux_bkgd(i)=total(hsimage_diffuse(index_cur)*em2flux*scale^2)

		;compute the bkgd-corrected flux (erg/cm^2/s)
                flux_corr(i)=flux_uncorr(i)-flux_bkgd(i)

		;compute the PPP-style bkgd-corrected flux (erg/cm^2/s)
		flux_corrppp(i)=total(hdata_pix(index_ppp)*scaleh_em_mod*$
			em2flux*scale^2)-total(hsimage_diffuse($
			index_cur(index_ppp))*em2flux*scale^2)
		;correct for flux not included in PPP pixels, see AJ 112, 109
		if not((resFWHM_pix_best_all(i) le min(resFWHM_values)) and $
			(modeltype_best_all(i) eq 0)) then begin
		flux_corrppp(i)=flux_corrppp(i)/(1.-pppval^3.)	;resolved case
		endif else begin
			ppppsf1=psf_Gaussian(NPIXEL=[15,15],$
				FWHM=[PSFfwhm_pix,$
				axrat_best_all(i)*PSFfwhm_pix])
			index_ppppsf1=where(ppppsf1 ge pppval)
			ppppsf2=psf_Gaussian(NPIXEL=[14,14],$
				FWHM=[PSFfwhm_pix,$
				axrat_best_all(i)*PSFfwhm_pix])
			index_ppppsf2=where(ppppsf2 ge pppval)
			pppfix=avg([total(ppppsf1(index_ppppsf1))/$
				total(ppppsf1),$
				total(ppppsf2(index_ppppsf2))/$
				total(ppppsf2)])
			flux_corrppp(i)=flux_corrppp(i)/pppfix ;unresolved case
		endelse

		;compute the luminosity at our assumed distance (erg/s)
		cm_pc=3.084d18
		flux2lum=4.d0*acos(-1.d0)*(double(dist_mpc)*1.d6*cm_pc)^2.d0
		lum_corr(i)=double(flux_corr(i))*flux2lum

		;furnish the L+C and C data in ADU,
		;bkgd uncertainty in electrons (2nd to last argument)
		;the result is returned in ADU
		estimate_error_end,hcdata_pix/scalehc,$
				cdata_pix/scalec,cscale,$
				ghc,gc,nhc_pix,nc_pix,$
				scalehc,scalec,offsethc,offsetc,$
				(hsimage_diffuse_err(index_cur)/$
				(scaleh_em_mod*scaleh))*ghc,$
				flux_corr_err,hcmean,cmean
		delta_flux_corr(i)=flux_corr_err*scaleh*scaleh_em_mod*$
				em2flux*scale^2

		;estimate the line-only final signal-to-noise ratio
		delta_flux_corr_line=sqrt(total((sqrt((hdata_pix>$
			(hcmean_orig-offsethc))/$
			(hcmean_orig-offsethc))*skysig)^2))*$
			scaleh_em_mod*em2flux*scale^2
		s2n_linef_all(i)=flux_uncorr(i)/delta_flux_corr_line
		;this really should be changed to reflect correction
		;for known bkgd emission...

	endif else begin
		print,'Region ',i,' no longer has pixels!'
		stop
	endelse
endfor


;compute the terminal signal-to-noise ratio
;this can be lower than sigma_find, since regions have grown into
;noiser areas after definition of the seeds
s2n=flux_corr/delta_flux_corr


;if we've determined that the images have WCS info, reorder the
;the catalog by decreasing RA (after computing RA and DEC for
;the center of the best-matching model from our FIND step)
if use_coords then begin
	for i=long(0),long(nreg)-long(1),long(1) do begin
		xy2ad,xpos_all+blc(0),ypos_all+blc(1),$
			astrometry,alpha,delta
	endfor
; figure out how to order the catalog by decreasing RA
	sort_order=reverse(sort(alpha))
endif else begin
; figure out how to order the catalog by increasing xpos_all
	sort_order=sort(xpos_all)
endelse

if interact then begin
	print,''
	print,'Dump a composite postage stamp image? (Y/n)'
	repeat begin
		ans=read_key(1)
	endrep until ((ans eq retnum) or (ans eq cnnum) or (ans eq cynum) or $
		(ans eq nnum) or (ans eq ynum))
	if not ((ans eq cynum) or (ans eq ynum) or (ans eq retnum)) then $
		dump_stamps=0L else dump_stamps=1L
endif

if (dump_stamps eq 1) then begin

;;;;;;;;;;;;;;;;;;;;;;;; CREATE LINE-ONLY POSTAGE STAMP IMAGE

;read the line-only image with boundaries marked
fxread,strcompress(root+'.GROW+BORDER.fits'),bigimage

npanels=ceil(sqrt(float(nreg)))
composite=fltarr(65*npanels,65*npanels)
compositename=strcompress(root+'.STAMP.fits',/REMOVE_ALL)
cblc_i=[0,65*(round(npanels)-1)]
ctrc_i=cblc_i+64
cblc=cblc_i
ctrc=cblc+64
borderval=max(bigimage)
;create 'postage stamp' collage, each source gets a 65x65 subsection
for i=long(0),long(nreg)-long(1),long(1) do begin
	i_sorted=sort_order(i)
	sblc=lonarr(2)
	strc=lonarr(2)
	sblc(0)=max([0,2*round(xpos_all(i_sorted))+1-32])
	sblc(1)=max([0,2*round(ypos_all(i_sorted))+1-32])
	sizebig=size(bigimage)
	strc(0)=min([sizebig(1)-1,2*round(xpos_all(i_sorted))+1+32])
	strc(1)=min([sizebig(2)-1,2*round(ypos_all(i_sorted))+1+32])
	subim=fltarr(65,65)
	subim(max([0,-(2*round(xpos_all(i_sorted))+1-32)]): $
		64-max([0,2*round(xpos_all(i_sorted))+1+32-(sizebig(1)-1)]), $
		max([0,-(2*round(ypos_all(i_sorted))+1-32)]): $
		64-max([0,2*round(ypos_all(i_sorted))+1+32-(sizebig(2)-1)]))= $
		bigimage(sblc(0):strc(0),sblc(1):strc(1))
	subim(0,*)=borderval
	subim(*,0)=borderval
	subim(64,*)=borderval
	subim(*,64)=borderval
	composite(cblc(0):ctrc(0),cblc(1):ctrc(1))=subim
;figure out placement of the next postage stamp image (within the collage)
	if cblc(0) ge cblc_i(1) then begin
		cblc(0)=0
		cblc(1)=cblc(1)-65
		ctrc=cblc+64
	endif else begin
		cblc(0)=cblc(0)+65
		ctrc=cblc+64
	endelse
endfor
fxhmake,head,float(composite)
fxwrite,compositename,head,float(composite)

;if use_cimage then begin
;
;;;;;;;;;;;;;;;;;;;;;;;;; CREATE CONTINUUM POSTAGE STAMP IMAGE
;
;fxread,strcompress(root+'.CONT.GROW+BORDER.fits'),bigimage
;
;npanels=ceil(sqrt(float(nreg)))
;composite=fltarr(65*npanels,65*npanels)
;compositename=strcompress(root+'.CONT.STAMP.fits',/REMOVE_ALL)
;cblc_i=[0,65*(round(npanels)-1)]
;ctrc_i=cblc_i+64
;cblc=cblc_i
;ctrc=cblc+64
;borderval=max(bigimage)
;;create 'postage stamp' collage
;for i=long(0),long(nreg)-long(1),long(1) do begin
;	i_sorted=sort_order(i)
;	sblc=lonarr(2)
;	strc=lonarr(2)
;	sblc(0)=max([0,2*round(xpos_all(i_sorted))+1-32])
;	sblc(1)=max([0,2*round(ypos_all(i_sorted))+1-32])
;	sizebig=size(bigimage)
;	strc(0)=min([sizebig(1)-1,2*round(xpos_all(i_sorted))+1+32])
;	strc(1)=min([sizebig(2)-1,2*round(ypos_all(i_sorted))+1+32])
;	subim=fltarr(65,65)
;	subim(max([0,-(2*round(xpos_all(i_sorted))+1-32)]): $
;		64-max([0,2*round(xpos_all(i_sorted))+1+32-(sizebig(1)-1)]), $
;		max([0,-(2*round(ypos_all(i_sorted))+1-32)]): $
;		64-max([0,2*round(ypos_all(i_sorted))+1+32-(sizebig(2)-1)]))= $
;		bigimage(sblc(0):strc(0),sblc(1):strc(1))
;	subim(0,*)=borderval
;	subim(*,0)=borderval
;	subim(64,*)=borderval
;	subim(*,64)=borderval
;	composite(cblc(0):ctrc(0),cblc(1):ctrc(1))=subim
;;figure out placement of the next postage stamp image within the collage
;	if cblc(0) ge cblc_i(1) then begin
;		cblc(0)=0
;		cblc(1)=cblc(1)-65
;		ctrc=cblc+64
;	endif else begin
;		cblc(0)=cblc(0)+65
;		ctrc=cblc+64
;	endelse
;endfor
;fxhmake,head,float(composite)
;fxwrite,compositename,head,float(composite)

;endif

endif

;;;;;;;;;;;;;;;;;;;;;;;; WRITE OUTPUT HII REGION CATALOG

openw,1,strcompress(root+'.catalog.dat',/REMOVE_ALL)
	printf,1,'source          = ',code_version,format='(a18,a)'
if interact then begin
	printf,1,'hs_filename     = ',hs_filename,format='(a18,a)'
	printf,1,'c_filename      = ',c_filename,format='(a18,a)'
	printf,1,'h_filename      = ',h_filename,format='(a18,a)'
	printf,1,'b_filename      = ',b_filename,format='(a18,a)'
	printf,1,'gc              = ',gc,format='(a18,e9.3)'
	printf,1,'nc              = ',nc,format='(a18,e9.3)'
	printf,1,'sc              = ',scalec,format='(a18,e9.3)'
	printf,1,'ghc             = ',ghc,format='(a18,e9.3)'
	printf,1,'nhc             = ',nhc,format='(a18,e9.3)'
	printf,1,'shc             = ',scalehc,format='(a18,e9.3)'
	printf,1,'sh              = ',scaleh,format='(a18,e9.3)'
	printf,1,'sh_em           = ',scaleh_em,format='(a18,e9.3)'	
	printf,1,'[NII]percent    = ',NIIpercent,format='(a18,e9.3)'
	printf,1,'dist_mpc        = ',dist_mpc,format='(a18,e9.3)'
	printf,1,'scale           = ',scale,format='(a18,e9.3)'
	printf,1,'sigma_find      = ',sigma_find,format='(a18,e9.3)'
	printf,1,'PSFfwhm_pix     = ',PSFfwhm_pix,format='(a18,e9.3)'
	printf,1,'size_max_pc     = ',size_max_pc,format='(a18,e9.3)'
	printf,1,'bkgd_annulus_pc = ',bkgd_annulus_pc,format='(a18,e9.3)'
	printf,1,'do_congrid      = ',docongrid,format='(a18,i1)'
	printf,1,'root            = ',root,format='(a18,a)'
	printf,1,'dump_stamps     = ',dump_stamps,format='(a18,i1)'
	printf,1,'footprint_lev   = ',footprint_lev,format='(a18,e9.3)'
	printf,1,'seed_lev        = ',seed_lev,format='(a18,e9.3)'
	printf,1,'f_remain_crit   = ',f_remain_crit,format='(a18,e9.3)'
	printf,1,'rho_crit        = ',rho_crit,format='(a18,e9.3)'
	printf,1,'termgrad        = ',termgrad,format='(a18,e9.3)'
endif else begin
	for i=0L,n_elements(batchinfo)-2L do begin
		printf,1,batchinfo(i)
	endfor
	printf,1,'termgrad        = ',termgrad,format='(a18,e9.3)'
endelse
	printf,1,'RUN_DATE        = ',!stime,format='(a18,a)'
        printf,1,'RUN_COMMENT     = ',format='(a18)'
	printf,1,'NUM_REGIONS     = ',nreg,format='(a18,i6)'
if use_coords then begin
	printf,1,'=============================================================================================================================================================================================================================================================================================================================================='
	printf,1,'  id   h  m    s    d  m    s          x          y   n   fwhm_eff   fwhm_maj   fwhm_min    maj/min         pa   f_uncorr     f_bkgd     f_corr    df_corr  f_corrPPP     l_corr        rho      alpha       beta    peak-EM  median-EM     min-EM    bkgd-EM  term-grad   conv-S/N  LINEi-S/N  LINEf-S/N   init-S/N  final-S/N                   comment'
	printf,1,'=============================================================================================================================================================================================================================================================================================================================================='
endif else begin
	printf,1,'====================================================================================================================================================================================================================================================================================================================='
	printf,1,'  id          x          y   n   fwhm_eff   fwhm_maj   fwhm_min    maj/min         pa   f_uncorr     f_bkgd     f_corr    df_corr  f_corrPPP     l_corr        rho      alpha       beta    peak-EM  median-EM     min-EM    bkgd-EM  term-grad   conv-S/N  LINEi-S/N  LINEf-S/N   init-S/N  final-S/N                   comment'
	printf,1,'====================================================================================================================================================================================================================================================================================================================='
endelse


for i=long(0),long(nreg)-long(1),long(1) do begin
	is=sort_order(i) 
	if use_coords then begin
		printf,1,id(i),adstring(alpha(is),delta(is)),$
			xpos_all(is)+blc(0)+1,ypos_all(is)+blc(1)+1,$
			n_pix(is),fwhm_eff(is),fwhm_maj(is),fwhm_min(is),$
			axrat_best_all(is),pa_deg(is),flux_uncorr(is),$
			flux_bkgd(is),flux_corr(is),delta_flux_corr(is),$
			flux_corrppp(is),lum_corr(is),$
			rho_best_all(is),alpha_best_all(is),beta_best_all(is),$
			sbpeak_all(is),sbmedian_all(is),sbmin_all(is),$
			hbkgd_all(is),termgrad_all(is),$
			s2n_convo_all(is),s2n_linei_all(is),s2n_linef_all(is),$
			s2n_initial_all(is),$
			s2n(is),comment(is),$
			FORMAT='(I4,1X,A24,1X,E10.3,1X,E10.3,1X,I3,24(1X,E10.3),1X,A25)'
	endif else begin
		printf,1,id(i),$
			xpos_all(is)+blc(0)+1,ypos_all(is)+blc(1)+1,$
			n_pix(is),fwhm_eff(is),fwhm_maj(is),fwhm_min(is),$
			axrat_best_all(is),pa_deg(is),flux_uncorr(is),$
			flux_bkgd(is),flux_corr(is),delta_flux_corr(is),$
			flux_corrppp(is),lum_corr(is),$
			rho_best_all(is),alpha_best_all(is),beta_best_all(is),$
			sbpeak_all(is),sbmedian_all(is),sbmin_all(is),$
			hbkgd_all(is),termgrad_all(is),$
			s2n_convo_all(is),s2n_linei_all(is),s2n_linef_all(is),$
			s2n_initial_all(is),$
			s2n(is),comment(is),$
			FORMAT='(I4,1X,E10.3,1X,E10.3,1X,I3,24(1X,E10.3),1X,A25)'
	endelse
endfor
close,1

update_map_files,root,root_suffix,nreg,sort_order

if multiple_cutoffs then change_filenames,root,root_suffix

if multiple_cutoffs then goto, CUTOFF_LOOP

END_OF_CODE:

math_status=CHECK_MATH()

end

;------------------------------------------------------
;------------------------------------------------------
;------------------------------------------------------

pro HIIphot_batch

;this pro allows the user to run through all the input-required points
;of the code and answer any required questions in advance.  The end-result
;is a batch capability.  HIIphot_batch produces *.HIIphot_batch.scr,
;*.HIIphot_batch, *.ROI_analysis.dat, and *.ROI_sky.dat files which
;can be used for hands-off operation of the package -- just
;type something like "______.HIIphot_batch.scr > ______.HIIphot_batch.log &"
;on the UNIX command line to run in the background while creating a
;log file of the batch-session.

code_version='HIIphot_batch v5.1'

;setup constants to interpret (read_key) results
ynum=121 & cynum=89 & nnum=110 & cnnum=78
spcnum=32 & retnum=10 & qnum=113 & cqnum=81

working_dir=getenv('PWD')
astrolib	; load Astronomy User's Library variables
ans=(read_key(0) and exist(working_dir))
scr_erase

;tell the user what is happening
print,''
print,'Running HIIphot_batch:'
print,''
print,'You will be prompted for a number of parameters to use during'
print,'execution of HIIphot.  If in doubt, accept the DEFAULTS.  They'
print,'have been chosen conservatively or using the image header.'
print,''
print,'After this code is finished, you may use the resulting files'
print,'in one of two ways:'
print,''
print,'1 - within IDL'
print,''
print,'    scorpio[91]% idl'
print,'    IDL> .comp $HIIphot/HIIphot'
print,'    IDL> HIIphot, batch=''______.HIIphot_batch'''
print,''
print,'2 - from the command line'
print,''
print,'    scorpio[92]% ______.HIIphot_batch.scr > ______.HIIphot_batch.log &'
print,''
print,'Note: ______ specifies the root name you will soon provide.'
print,''
print,''
print,'Press any key to begin...'
ans=read_key(1)
scr_erase

;setup some flags
batch=1L
interact=1L
display=1L
docongrid=0L

;pick out the images and properly set flags indicating which images were
;actually selected
select_images, display, interact, $
	hs_filename, c_filename, h_filename, b_filename, $
	use_cimage, use_himage, use_bimage

;read in the line image, extract any astrometry
scr_erase
print,'Reading LINE image...'
fxread, hs_filename, hsimage, hshead
index_bad=where(finite(hsimage) eq 0,count_bad)
if count_bad gt 0 then hsimage(index_bad)=0.
ghc=fxpar(hshead,'GAIN')
if ghc eq 0. then ghc=3.
nhc=fxpar(hshead,'NCOMBINE')
if nhc eq 0 then nhc=1
naxis=fxpar(hshead,'NAXIS*')
naxis_orig=naxis
extast,hshead,astrometry,noparams
equinox_img=fxpar(hshead,'EQUINOX')
if equinox_img eq 0 then equinox_img=fxpar(hshead,'EPOCH')
equinox_img=float(equinox_img)
use_coords=1L
;demand that the header contains CDn_m astrometry in order to use WCS info
;it can also contain ROTA + CDELT (AIPS-type) astrometry, as long as the
;CD matrix is filled in a consistent manner.
if noparams lt 2 then use_coords=0L
crota=fltarr(2)
if use_coords then begin
	getrot,hshead,rot,cdelt
	crota(0)=0.
	crota(1)=rot
	print,'Found WCS info.'
	if equinox_img lt 1.e-6 then begin
		print,'EQUINOX keyword is missing! Assuming J2000 coords.'
		equinox_img=2000.
	endif
	astrometry_big=astrometry
	astrometry_big.cd=astrometry_big.cd/2.
	astrometry_big.crpix=astrometry_big.crpix*2.
endif
print,''

;compare the continuum image, extract astrometry if the line image lacked it
if use_cimage then begin
print,'Reading CONTINUUM image...'
fxread,c_filename,cimage,chead 
index_bad=where(finite(cimage) eq 0,count_bad)
if count_bad gt 0 then cimage(index_bad)=0.
gc=fxpar(chead,'GAIN')
if gc eq 0. then gc=3.
nc=fxpar(chead,'NCOMBINE')
if nc eq 0 then nc=1
cimage_naxis=fxpar(chead,'NAXIS*')
if naxis(0) ne cimage_naxis(0) then begin
	message,'ERROR: CONTINUUM image isn''t the required X size!',/CON
	return
endif
if naxis(1) ne cimage_naxis(1) then begin
	message,'ERROR: CONTINUUM image isn''t the required Y size!',/CON
	return
endif
if (noparams lt 2) then begin ; try again with astrometry
	extast,chead,astrometry,noparams
	equinox_img=fxpar(chead,'EQUINOX')
	if equinox_img eq 0 then equinox_img=fxpar(chead,'EPOCH')
	equinox_img=float(equinox_img)
	if equinox_img lt 1.e-6 then begin
		print,'EQUINOX keyword is missing! Assuming J2000 coords.'
		equinox_img=2000.
	endif
	use_coords=1L
;demand that the header contains CDn_m astrometry in order to use WCS info
;it can also contain ROTA + CDELT (AIPS-type) astrometry, as long as the
;CD matrix is filled in a consistent manner.
	if noparams lt 2 then use_coords=0L
	crota=fltarr(2)
	if use_coords then begin
		getrot,chead,rot,cdelt
		crota(0)=0.
		crota(1)=rot
		print,'Found WCS info.'
		astrometry_big=astrometry
		astrometry_big.cd=astrometry_big.cd/2.
		astrometry_big.crpix=astrometry_big.crpix*2.
	endif
endif
print,''
endif

;in the case that only the LINE and CONTINUUM images are available
;we will prompt for an appropriate scaling factor in order to
;create a pseudo LINE+CONTINUUM image
if use_cimage and not(use_himage) then begin
	root='TEMP'
	spawn,'\rm TEMP.PSEUDOhc.fits'
	pseudo_linecontinuum,hsimage,cimage,root,h_filename,$
		scaleh,scalec,scalehc,use_coords,astrometry,equinox_img
	use_himage=1L
endif

;check the line+continuum image
if use_himage then begin
print,'Checking LINE+CONTINUUM image...'
fxread,h_filename,himage,hhead,0,1,0,1 
index_bad=where(finite(himage) eq 0,count_bad)
if count_bad gt 0 then himage(index_bad)=0.
ghc=fxpar(hhead,'GAIN')
if ghc eq 0. then ghc=3.
nhc=fxpar(hhead,'NCOMBINE')
if nhc eq 0 then nhc=1
himage_naxis=fxpar(hhead,'NAXIS*')
if naxis(0) ne himage_naxis(0) then begin
	message,'ERROR: LINE+CONTINUUM image isn''t the required X size!',/CON
	return
endif
if naxis(1) ne himage_naxis(1) then begin
	message,'ERROR: LINE+CONTINUUM image isn''t the required Y size!',/CON
	return
endif
print,''
endif

;check the blanking image
if use_bimage then begin
fxread,b_filename,bimage_float,bimage_header
index_bad=where(finite(bimage_float) eq 0,count_bad)
if count_bad gt 0 then bimage(index_bad)=0.
bimage_naxis=fxpar(bimage_header,'NAXIS*')
if naxis(0) ne bimage_naxis(0) then begin
	message,'ERROR: BLANKING image isn''t the right x size!',/CON
	return
endif
if naxis(1) ne bimage_naxis(1) then begin
	message,'ERROR: BLANKING image isn''t the right y size!',/CON
	return
endif
print,''
endif

if use_bimage then bimage=fix(ceil(bimage_float))


;this is where the user gets to pick what all the output begins with...
;that is, the _____ in _____.HIIphot_batch.scr
root=''
read,'Root name for all output? (DEFAULT=''HIIphot''): ',root
if ((root ne 'TEMP') and exist('TEMP.PSEUDOhc.fits')) then begin
	mvcmd=strcompress('mv TEMP.PSEUDOhc.fits '+strcompress(root+'.PSEUDOhc.fits',/REMOVE_ALL))
	spawn,mvcmd
	h_filename=strcompress(root+'.PSEUDOhc.fits',/REMOVE_ALL)
endif
if root eq '' then root='HIIphot'

;now the user can pick out analysis and sky regions
get_regions, display, interact, root, hsimage, hshead, naxis, h_filename, $
	use_bimage, bimage, blc, trc, blc_sky, trc_sky, sky_bimage, docongrid
if b_filename eq 'NULL' then bimage_float=float(bimage)

;here we prompt for all the FIND related parameters (not yet termgrad_arr
;and multiple cutoffs... they are skipped in the prompting when running
;HIIphot_batch).  the growing paramters will be asked for below...
get_find_inputs, batch, interact, dist_mpc, use_coords, scale, cdelt, pc_pix,$
	scalehc, scalec, scaleh, scaleh_em, hs_filename,h_filename,c_filename,$
	cscale, ghc, gc, nhc, nc, sigma_find, PSFfwhm_pix, size_max_pc,$
	termgrad_arr, multiple_cutoffs, correction, NIIpercent

;this is the spot at which the user selects a width for bkgd annuli
print,''
print,'Bkgd contributions to the HII region flux will be estimated'
print,'using median pixel values taken from annuli surrounding all'
print,'regions.'
print,''
print,'The minimum allowable width of such annuli is 5 pix.'
print,''
GET_BKGD_ANNULUS_PC:   
ans = ''
read, 'Annulus width for bkgd estimation in pc? (DEFAULT = 250.0 pc): ',ans 
if ans EQ '' then bkgd_annulus_pc=250. else begin
	bkgd_annulus_pc = getopt(ans,'F')
	if N_elements(bkgd_annulus_pc) NE 1 then begin  
       		message, 'ERROR - Expecting floating-point scalar',/CON
      		goto, GET_BKGD_ANNULUS_PC     
	endif
endelse

;if the user-selected annulus width was too small then bump it up to 5 pix
bkgd_rad=max([5.,bkgd_annulus_pc/pc_pix])
bkgd_rad_pc=bkgd_rad*pc_pix
print,'-------------------------------------------------------------------------------'
if bkgd_annulus_pc/pc_pix lt 5. then begin
	print,'Warning: The minimum allowable width of such annuli is 5 pix.'
	print,'         ',strcompress('This corresponds to '+$
		string(bkgd_rad_pc)+' pc.  Try again...')
	print,''
	goto,GET_BKGD_ANNULUS_PC
endif

;select if composite postage stamp images are to be created
print,''
print,'Dump a composite postage stamp image? (y/N)'
repeat begin
	ans=read_key(1)
endrep until ((ans eq retnum) or (ans eq cnnum) or (ans eq cynum) or $
	(ans eq nnum) or (ans eq ynum))
if not ((ans eq cnnum) or (ans eq nnum) or (ans eq retnum)) then $
	dump_stamps=1L else dump_stamps=0L
wait,0.5
print,''

;the footprint and seed parameters usually should be modified, but...
print,'Use the DEFAULT values all of parameters regarding initial'
print,'footprints, seeds, and region growth?  Strongly recommended!! (Y/n)'
repeat begin
	ans=read_key(1)
endrep until ((ans eq retnum) or (ans eq cnnum) or (ans eq cynum) or $
	(ans eq nnum) or (ans eq ynum))

if not ((ans eq cnnum) or (ans eq nnum)) then begin ; use the DEFAULTS

	footprint_lev=0.2
	seed_lev=0.50
	f_remain_crit=0.9
	rho_crit=0.25
	termgrad_arr=[1.5,0.0,0.0,0.0,0.0]

endif else begin ;potentially change some of the defaults

	print,''
	GET_FOOTPRINT_LEV:   
	ans = ''
	read,  'Footprint level? (DEFAULT = 0.2): ',ans   
	if ans EQ '' then footprint_lev=0.2 else begin
		footprint_lev = getopt(ans,'F')
		if N_elements(footprint_lev) NE 1 then begin  
              		message, 'ERROR - Expecting floating-point scalar',/CON
              		goto, GET_FOOTPRINT_LEV     
		endif
		if (footprint_lev ge 1.) then begin
			message, 'ERROR - Expecting scalar < 1.',/CON
			goto, GET_FOOTPRINT_LEV
		endif
	endelse

	print,''
	GET_SEED_LEV:   
	ans = ''
	read,  'Seed level? (DEFAULT = 0.50): ',ans   
	if ans EQ '' then seed_lev=0.50 else begin
		seed_lev = getopt(ans,'F')
		if N_elements(seed_lev) NE 1 then begin  
              		message, 'ERROR - Expecting floating-point scalar',/CON
              		goto, GET_SEED_LEV     
		endif
		if (seed_lev ge 1.) then begin
			message, 'ERROR - Expecting scalar < 1.',/CON
			goto, GET_SEED_LEV
		endif
	endelse

	print,''
	GET_F_REMAIN_CRIT:   
	ans = ''
	read,  'Minimum fraction of flux remaining unclaimed? (DEFAULT = 0.9): ',ans   
	if ans EQ '' then f_remain_crit=0.9 else begin
		f_remain_crit = getopt(ans,'F')
		if N_elements(f_remain_crit) NE 1 then begin  
              		message, 'ERROR - Expecting floating-point scalar',/CON
              		goto, GET_F_REMAIN_CRIT     
		endif
		if (f_remain_crit ge 1.) then begin
			message, 'ERROR - Expecting scalar < 1.',/CON
			goto, GET_F_REMAIN_CRIT
		endif
	endelse

	print,''
	GET_RHO_CRIT:   
	ans = ''
	read,  'Minimum correlation coefficient? (DEFAULT = 0.25): ',ans   
	if ans EQ '' then rho_crit=0.25 else begin
		rho_crit = getopt(ans,'F')
		if N_elements(rho_crit) NE 1 then begin  
              		message, 'ERROR - Expecting floating-point scalar',/CON
              		goto, GET_RHO_CRIT     
		endif
		if (rho_crit ge 1.) then begin
			message, 'ERROR - Expecting scalar < 1.',/CON
			goto, GET_RHO_CRIT
		endif
	endelse

        print,'' 
        GET_TERMGRAD:    
        ans = '' 
        read,  'Terminal gradient (EM/pc)? (DEFAULT = 1.5,0,0,0,0): ',ans  
        if ans EQ '' then termgrad_arr=[1.5,0.0,0.0,0.0,0.0] else begin 
                termgrad_arr = getopt(ans,'F') 
                if N_elements(termgrad_arr) NE 5 then begin   
                        message, 'ERROR - Expecting 5 floating-point scalars',/CON 
                        goto, GET_TERMGRAD      
                endif 
                if (min(termgrad_arr) lt 0.) then begin 
                        message, 'ERROR - Expecting all scalars >= 0.',/CON 
                        goto, GET_TERMGRAD 
                endif 
        endelse 


endelse

;now write the *.HIIphot_batch file to disk
openw,1,strcompress(root+'.HIIphot_batch')
printf,1,'source          = ',code_version,format='(a18,a)'
printf,1,'hs_filename     = ',hs_filename,format='(a18,a)'
printf,1,'c_filename      = ',c_filename,format='(a18,a)'
printf,1,'h_filename      = ',h_filename,format='(a18,a)'
printf,1,'b_filename      = ',b_filename,format='(a18,a)'
printf,1,'gc              = ',gc,format='(a18,e9.3)'
printf,1,'nc              = ',nc,format='(a18,e9.3)'
printf,1,'sc              = ',scalec,format='(a18,e9.3)'
printf,1,'ghc             = ',ghc,format='(a18,e9.3)'
printf,1,'nhc             = ',nhc,format='(a18,e9.3)'
printf,1,'shc             = ',scalehc,format='(a18,e9.3)'
printf,1,'sh              = ',scaleh,format='(a18,e9.3)'
printf,1,'sh_em           = ',scaleh_em,format='(a18,e9.3)'
printf,1,'[NII]percent    = ',NIIpercent,format='(a18,e9.3)'
printf,1,'dist_mpc        = ',dist_mpc,format='(a18,e9.3)'
printf,1,'scale           = ',scale,format='(a18,e9.3)'
printf,1,'sigma_find      = ',sigma_find,format='(a18,e9.3)'
printf,1,'PSFfwhm_pix     = ',PSFfwhm_pix,format='(a18,e9.3)'
printf,1,'size_max_pc     = ',size_max_pc,format='(a18,e9.3)'
printf,1,'bkgd_annulus_pc = ',bkgd_annulus_pc,format='(a18,e9.3)'
printf,1,'do_congrid      = ',docongrid,format='(a18,i1)'
printf,1,'root            = ',root,format='(a18,a)'
printf,1,'dump_stamps     = ',dump_stamps,format='(a18,i1)'
printf,1,'footprint_lev   = ',footprint_lev,format='(a18,e9.3)'
printf,1,'seed_lev        = ',seed_lev,format='(a18,e9.3)'
printf,1,'f_remain_crit   = ',f_remain_crit,format='(a18,e9.3)'
printf,1,'rho_crit        = ',rho_crit,format='(a18,e9.3)'
printf,1,'termgrad_arr    = ',termgrad_arr(0),termgrad_arr(1),$
	termgrad_arr(2),termgrad_arr(3),termgrad_arr(4),$
	format='(a18,e9.3,1x,e9.3,1x,e9.3,1x,e9.3,1x,e9.3,1x)'
close,1

;now write the *.HIIphot_batch.scr file to disk
openw,1,strcompress(root+'.HIIphot_batch.scr')
printf,1,'#!/bin/csh'
printf,1,'#'
printf,1,'idl << EOF'
printf,1,'.comp $HIIphot/HIIphot'
printf,1,strcompress('HIIphot,BATCH='''+strcompress(root+'.HIIphot_batch')+'''')
printf,1,'exit'
printf,1,'EOF'
close,1
spawncmd=strcompress('chmod 777 '+root+'.HIIphot_batch.scr')
spawn,spawncmd,spawnout

end

pro select_images, display, interact,$
	hs_filename, c_filename, h_filename, b_filename,$
	use_cimage, use_himage, use_bimage

;setup constants to interpret (read_key) results
ynum=121 & cynum=89 & nnum=110 & cnnum=78
spcnum=32 & retnum=10 & qnum=113 & cqnum=81

;set flags appropriately for a batch session
if not(interact) then begin
	use_cimage=1L
	use_himage=1L
	use_bimage=1L
	if ((c_filename eq 'NULL') or (c_filename eq '')) then use_cimage=0L
	if ((h_filename eq 'NULL') or (h_filename eq '')) then use_himage=0L
	if ((b_filename eq 'NULL') or (b_filename eq '')) then use_bimage=0L
	print,'Image list accepted.'
	return
endif

ASK_AGAIN: ;jump here if the user wasn't happy with the echoed list

scr_erase

working_dir=getenv('PWD')

;determine filenames for all images to be used...

if display then begin

print,'Select the continuum-subtracted LINE image to be analyzed.'
print,''
; get hs_filename interactively
repeat begin
	hs_filename=pickfile($
	TITLE='Select a continuum-subtracted fits file to analyze',$
	FILTER='*')
endrep until ((strlen(hs_filename) ne 0) and $
	(strmid(hs_filename,strlen(hs_filename)-1,132) ne '/') and $
		exist(hs_filename))
print,strcompress('Using '+hs_filename)

;this will help to trim the list of displayed files, to only show images...
extension=strcompress('*'+strmid(hs_filename,rstrpos(hs_filename,'.'),132))

wait,1
scr_erase
print,'Select the CONTINUUM image (if available), else CANCEL...'
print,''
; get c_filename interactively if available
repeat begin
	c_filename=pickfile($
	TITLE='Select the CONTINUUM image (if available)',$
		FILTER=extension)
endrep until ((strlen(c_filename) eq 0) or $
	((strmid(c_filename,strlen(c_filename)-1,132) ne '/') and $
	exist(c_filename)))
if (strlen(c_filename) eq 0) then begin
	use_cimage=0L
	c_filename='NULL'
	print,'CONTINUUM image not selected.'
	print,''
	print,'Currently there is no provision for using only'
	print,'the continuum-subtracted line image, please choose'
	print,'files again.'
	print,''
	print,'Press any key to continue...'
	junk=read_key(1)
	goto,ASK_AGAIN
endif else begin
	use_cimage=1L
	print,strcompress('Using '+c_filename)
endelse

extension=strcompress('*'+strmid(c_filename,rstrpos(c_filename,'.'),132))

wait,1
scr_erase
print,'Select the LINE+CONTINUUM image (if available), else CANCEL...'
print,''
; get h_filename interactively if available
repeat begin
	h_filename=pickfile($
	TITLE='Select the LINE+CONTINUUM image (if available)',$
	FILTER=extension)
endrep until ((strlen(h_filename) eq 0) or $
	((strmid(h_filename,strlen(h_filename)-1,132) ne '/') and $
	exist(h_filename)))
if (strlen(h_filename) eq 0) then begin
	use_himage=0L
	h_filename='NULL'
	print,'LINE+CONTINUUM image not selected.'
endif else begin
	use_himage=1L
	print,strcompress('Using '+h_filename)
endelse

extension=strcompress('*'+strmid(h_filename,rstrpos(h_filename,'.'),132))

; get b_filename interactively if available
wait,1
scr_erase
print,'Select the BLANKING image (if desired), else CANCEL...'
print,''
repeat begin
	b_filename=pickfile($
	TITLE='Select the BLANKING image, else CANCEL...',$
	FILTER=extension)
endrep until ((strlen(b_filename) eq 0) or $
	((strmid(b_filename,strlen(b_filename)-1,132) ne '/') and $
	exist(b_filename)))
if (strlen(b_filename) eq 0) then begin
	use_bimage=0L
	b_filename='NULL'
	print,'BLANKING image not selected.'
endif else begin
	use_bimage=1L
	print,strcompress('Using '+b_filename)
endelse

extension=strcompress('*'+strmid(hs_filename,rstrpos(hs_filename,'.'),132))

wait,1
scr_erase

endif else begin ; no display is available, get the filenames manually

hs_filename=''
c_filename=''
h_filename=''
b_filename=''

read,'Filename of the continuum-subtracted LINE image? ',hs_filename
read,'Filename of the CONTINUUM image? [return to skip] ',c_filename
read,'Filename of the LINE+CONTINUUM image? [return to skip] ',h_filename
read,'Filename of the BLANKING image? [return to skip] ',b_filename
if (strlen(c_filename) eq 0) then begin
	use_cimage=0L
	c_filename='NULL'
	print,''
	print,'CONTINUUM image not selected.'
	print,''
	print,'Currently there is no provision for using only'
	print,'the continuum-subtracted line image, please choose'
	print,'files again.'
	print,''
	print,'Press any key to continue...'
	junk=read_key(1)
	goto,ASK_AGAIN	
endif else begin
	use_cimage=1L
endelse
if (strlen(h_filename) eq 0) then begin
	use_himage=0L
	h_filename='NULL'
endif else begin
	use_himage=1L
endelse
if (strlen(b_filename) eq 0) then begin
	use_bimage=0L
	b_filename='NULL'
endif else begin
	use_bimage=1L
endelse
print,''

;add path information to the names (if needed)
if (strpos(hs_filename,'/') eq -1) then hs_filename=working_dir+'/'+hs_filename
if (strpos(c_filename,'/') eq -1) and (use_cimage) then $
	c_filename=working_dir+'/'+c_filename
if (strpos(h_filename,'/') eq -1) and (use_himage) then $
	h_filename=working_dir+'/'+h_filename
if (strpos(b_filename,'/') eq -1) and (use_bimage) then $
	b_filename=working_dir+'/'+b_filename

endelse

;check with the user to see if they selected the right files
print,'LINE image           : '+hs_filename
print,'CONTINUUM image      : '+c_filename
print,'LINE+CONTINUUM image : '+h_filename
print,'BLANKING image       : '+b_filename
print,''
print,'Is this correct? (y/n)'
repeat begin
	ans=read_key(1)
endrep until ((ans eq retnum) or (ans eq cnnum) or (ans eq cynum) or $
	(ans eq nnum) or (ans eq ynum))

;if the user wasn't happy with the list, then ask again by looping up
if ((ans eq cnnum) or (ans eq nnum)) then goto, ASK_AGAIN

return
end

pro pseudo_linecontinuum,hsimage,cimage,root,h_filename,$
	scaleh,scalec,scalehc,use_coords,astrometry,equinox_img

print,''
print,'The supplied LINE and CONTINUUM images will be used to'
print,'compute a pseudo LINE+CONTINUUM image for use by HIIphot.'
print,''
print,'In order to do this, the continuum scaling factor used'
print,'during original data reduction is required.  That is,'
print,'by what value was the continuum image multiplied'
print,'prior to continuum-subtraction?'
print,''
GET_CSCALE:   
ans = ''
read,  'Continuum scaling factor? (DEFAULT = 1.0): ',ans   
if ans EQ '' then cscale_inp=1.0 else begin
	cscale_inp = getopt(ans,'F')
	if N_elements(cscale_inp) NE 1 then begin  
        	message, 'ERROR - Expecting floating-point scalar',/CON
              		goto, GET_CSCALE     
	endif
endelse

if (n_elements(scaleh) ne 1) then begin
	print,''
	GET_SCALEH:   
	ans = ''
	print, 'Scale factor applied to the continuum-sub (line-only) image'
	print, 'during flux calibration to convert from ADU'
	read,  'to the current units? (DEFAULT = 1.0 -> none applied): ',ans   
	if ans EQ '' then scaleh=1. else begin
		scaleh = getopt(ans,'F')
		if N_elements(scaleh) NE 1 then begin  
              		message, 'ERROR - Expecting floating-point scalar',/CON
              		goto, GET_SCALEH     
		endif
	endelse
endif

if (n_elements(scalec) ne 1) then begin
	print,''
	GET_SCALEC:   
	ans = ''
	print, 'Scale factor applied to the continuum image'
	print, 'during flux calibration to convert from ADU'
	read,  'to the current units? (DEFAULT = 1.0 -> none applied): ',ans   
	if ans EQ '' then scalec=1. else begin
		scalec = getopt(ans,'F')
		if N_elements(scalec) NE 1 then begin  
              		message, 'ERROR - Expecting floating-point scalar',/CON
              		goto, GET_SCALEC     
		endif
	endelse
endif 

hcimg=hsimage/scaleh+cscale_inp*cimage/scalec

fxhmake,headhc,hcimg
if use_coords then putast,headhc,astrometry,EQUINOX=equinox_img
h_filename=strcompress(root+'.PSEUDOhc.fits',/REMOVE_ALL)
fxwrite,h_filename,headhc,hcimg

scalehc=1.

print, ''
print, 'A pseudo-LINE+CONT image has been created.  It has an estimated'
print, 'sky background, which may differ from what was actually observed.'
print, 'This estimate is typically fine, although it obviously cannot'
print, 'account for any genuine sky brightness changes between CONT and'
print, 'LINE+CONT observations.  I would recommend finding the actual'
print, 'LINE+CONT image if at all possible.'
print, ''
print,'Press any key to continue...'
junk=read_key(1)

return
end


pro get_regions, display, interact, root, hsimage, hshead, naxis, h_filename, $
	use_bimage, bimage, blc, trc, blc_sky, trc_sky, sky_bimage, docongrid

;this routine permits the user to select two regions within the data.
;one of these regions delineates the section in which HII regions are
;to be analyzed.  the other is used to characterize the sky properties.
;returned from this routine are the blanking images (either created
;or modified from an existing mask) and the blc,trc variables.

;setup constants to interpret (read_key) results
ynum=121 & cynum=89 & nnum=110 & cnnum=78
spcnum=32 & retnum=10 & qnum=113 & cqnum=81

;if we are running in batch and the user didn't use CONGRID when
;creating the batch file, then skip some of the lines below
if not(interact) and (docongrid le 0) then begin
	aspect_ratio=1.
	ROI_analysis=1L
	ROI_sky=1L
	goto,noninteractive_jump
endif
;if we are running in batch and the user did use CONGRID when
;creating the batch file, then determine aspect_ratio and congrid sizes
if not(interact) and (docongrid gt 0) then begin
	aspect_ratio=float(naxis(1))/float(naxis(0))
	if (naxis(0) gt naxis(1)) then begin
		congrid_xsize=512
		congrid_ysize=512*aspect_ratio
	endif else begin
		congrid_xsize=512*1./aspect_ratio
		congrid_ysize=512
	endelse
	ROI_analysis=1L
	ROI_sky=1L
endif
;batch sessions always use the ROI interface...
	

GET_REGIONS_AGAIN:	;jump here from below when the ROI is incorrect

scr_erase

;if batch mode is not being executed
if interact then begin

	;initialize the aspect ratio
	aspect_ratio=1.

	;when the user asks for real-time display,
	;provide opportunity to inspect image, decide on ROI, blc, trc, etc.
	if display then begin

		print,'-------------------------------------------------------------------------------'
		print,'Please determine BLC and TRC coordinates for a box'
		print,'containing all the HII regions to be analyzed.
		print,''
		print,'Alternatively, use the ''Define ROI'' interface to'
		print,'draw a border for the area of interest.
		print,''
		print,'-------------------------------------------------------------------------------'
		print,''

		;if either dimension of the image is larger than the
		;IMGroam display, allow the possibility of re-gridding
		;the image to fit within 512x512 pixels
		if naxis(0) gt 512 or naxis(1) gt 512 then begin
	print,'With the ''Define ROI'' interface it is only possible to select'
	print,'a 512x512 section of the image if the data is displayed at full'
	print,'resolution.'
	print,''
	print,'Would you like to view a regridded version of the entire image'
	print,'to enable interactive selection of big areas? (y/N)'
	
			;wait until Y,y,N,n, or return is pressed
			repeat begin
				ans=read_key(1)
			endrep until ((ans eq retnum) or (ans eq cnnum) or $
				(ans eq cynum) or (ans eq nnum) or $
				(ans eq ynum))

			;if N,n, or return was pressed, use the original image
			;else regrid to fit within the IMGroam window
			if ((ans eq retnum) or (ans eq cnnum) or $
				(ans eq nnum)) then docongrid=0L $
				else docongrid=1L

		;there is no need to consider regridding if the image
		;already fits within the IMGroam display
		endif else begin
	
			print,'Press any key to begin...'
			ans=read_key(1)
	
			docongrid=0L
	
		endelse
		;now docongrid has been set appropriately for interactive use


		TRY_ANALYSIS_AGAIN:	;jump here if the user wants another
					;try at drawing the analysis ROI

		;remove existing raw ROI files
		spawn,'\rm ROI???.dat',spawnout
		
		;if the interactive user wanted to see original data
		if not(docongrid) then begin
			imgroam,hsimage,hshead
		;else determine congrid parameters and then show regridded data
		endif else begin
			aspect_ratio=float(naxis(1))/float(naxis(0))
			if (naxis(0) gt naxis(1)) then begin
				congrid_xsize=512
				congrid_ysize=512*aspect_ratio
			endif else begin
				congrid_xsize=512*1./aspect_ratio
				congrid_ysize=512
			endelse
			imgroam,congrid(hsimage,congrid_xsize,congrid_ysize)
		endelse

		;check to see if the user made an analysis ROI file
		spawnout=''
		ROI_analysis=0L
		spawn,'ls -1 ROI000.dat',spawnout
		;if they did make one, then accept it and do a name change
		if spawnout(n_elements(spawnout)-1) eq 'ROI000.dat' then begin
			print,'Analysis region accepted.'
			print,''
			spawncmd='\mv ROI000.dat '+strcompress(root+$
				'.ROI_analysis.dat',/REMOVE_ALL)
			spawn,spawncmd,spawnout
			ROI_analysis=1L
		endif

		;if an analysis ROI file was made, ask if it should be used
		if ROI_analysis then begin
			print,'Do you wish to use it? (Y/n)'		
			repeat begin
				ans=read_key(1)
			endrep until ((ans eq retnum) or (ans eq cnnum) or $
				(ans eq cynum) or (ans eq nnum) or $
				(ans eq ynum))
			if ((ans eq cnnum) or (ans eq nnum)) then begin
				;remove the rejected ROI file
				spawncmd='\rm '+strcompress(root+$
				'.ROI_analysis.dat',/REMOVE_ALL)
				spawn,spawncmd,spawnout
				;allow the option of another chance	
				print,'Do you want to try again? (Y/n)'
				repeat begin
					ans=read_key(1)
				endrep until ((ans eq retnum) or $
					(ans eq cnnum) or $
					(ans eq cynum) or (ans eq nnum) or $
					(ans eq ynum))
				if ((ans eq cnnum) or (ans eq nnum)) then begin
					ROI_analysis=0L
				endif else begin
					goto,TRY_ANALYSIS_AGAIN
				endelse
			endif
			scr_erase 
		endif
		;now we know if the analysis ROI was to be used

		;read line+continuum image from disk
		fxread,h_filename,hcimage,hchead
		index_bad=where(finite(hcimage) eq 0,count_bad)
		if count_bad gt 0 then hcimage(index_bad)=0.

		print,'-------------------------------------------------------------------------------'
		print,'Please determine BLC and TRC coordinates for a box'
		print,'that can be used to characterize the sky, eventually'
		print,'leading to a noise model.'
		print,''
		print,'Alternatively, use the ''Define ROI'' interface to'
		print,'draw a border for the sky region.'
		print,''
		print,'Remember to include enough sky pixels so that we still'
		print,'have several independent data points after smoothing.'
		print,''
		print,'-------------------------------------------------------------------------------'
		print,''

		print,'Press any key to begin...'
		ans=read_key(1)
	
		TRY_SKY_AGAIN:	;jump here if the user wants another
				;try at drawing the sky ROI

		;remove existing raw ROI files
		spawn,'\rm ROI???.dat',spawnout
		
		;if the interactive user wanted to see original data
		if not(docongrid) then begin
			imgroam,hcimage,hchead
		;else determine congrid parameters and then show regridded data
		endif else begin
			aspect_ratio=float(naxis(1))/float(naxis(0))
			if (naxis(0) gt naxis(1)) then begin
				congrid_xsize=512
				congrid_ysize=512*aspect_ratio
			endif else begin
				congrid_xsize=512*1./aspect_ratio
				congrid_ysize=512
			endelse
			imgroam,congrid(hcimage,congrid_xsize,congrid_ysize)
		endelse

		;check to see if the user made a sky ROI file
		spawnout=''
		ROI_sky=0L
		spawn,'ls -1 ROI000.dat',spawnout
		;if they did make one, then accept it and do a name change
		if spawnout(n_elements(spawnout)-1) eq 'ROI000.dat' then begin
			print,'Sky region accepted.'
			print,''
			spawncmd='\mv ROI000.dat '+strcompress(root+$
				'.ROI_sky.dat',/REMOVE_ALL)
			spawn,spawncmd,spawnout
			ROI_sky=1L
		endif

		;if a sky ROI file was made, ask if it should be used
		if ROI_sky then begin
			print,'Do you wish to use it? (Y/n)'		
			repeat begin
				ans=read_key(1)
			endrep until ((ans eq retnum) or (ans eq cnnum) or $
				(ans eq cynum) or (ans eq nnum) or $
				(ans eq ynum))
			if ((ans eq cnnum) or (ans eq nnum)) then begin
				;remove the rejected ROI file
				spawncmd='\rm '+strcompress(root+$
				'.ROI_sky.dat',/REMOVE_ALL)
				spawn,spawncmd,spawnout
				;allow the option of another chance	
				print,'Do you want to try again? (Y/n)'
				repeat begin
					ans=read_key(1)
				endrep until ((ans eq retnum) or $
					(ans eq cnnum) or $
					(ans eq cynum) or (ans eq nnum) or $
					(ans eq ynum))
				if ((ans eq cnnum) or (ans eq nnum)) then begin
					ROI_sky=0L
				endif else begin
					goto,TRY_SKY_AGAIN
				endelse
			endif
			scr_erase 
		endif
		;now we know if the sky ROI was to be used
		
	endif else begin ; clearly if the user doesn't want a real-time
			 ; display then we can't use any ROI files...
		ROI_analysis=0L
		ROI_sky=0L
		docongrid=0L
	endelse

endif

;interactive and batch sessions both pick-up here...

noninteractive_jump:
if not(use_bimage) then bimage=fltarr(naxis(0),naxis(1))+1

;if we are to use an analysis ROI file, read it and interpret
if ROI_analysis then begin
	openr,15,strcompress(root+'.ROI_analysis.dat',/REMOVE_ALL)
	junk=strarr(2)
	readf,15,junk
	junk=''
	readf,15,junk
	nverts=round(float(strmid(junk,10,132)))
	junk=strarr(3)
	readf,15,junk
	list_verts=fltarr(2,nverts)
	readf,15,list_verts
	close,15
; need to scale the analysis list_verts for CONGRID taking place above
	if (docongrid) then begin
		list_verts(0,*)=list_verts(0,*)*float(naxis(0))/$
			float(congrid_xsize)
		list_verts(1,*)=list_verts(1,*)*float(naxis(1))/$
			float(congrid_ysize)
	endif

	;list_verts defines a polygon on (1...naxis) grid, not (0...naxis-1)
	;polyfillv takes care of the conversion to (0...naxis-1) coords
	;so if we want to include from x_blc,y_blc to x_trc,y_trc on the
	;(1...naxis) grid -- then list_verts must look like:
	;	 x_blc-delta y_blc-delta
	;	 x_trc+delta y_blc-delta
	;	 x_trc+delta y_trc+delta
	;	 x_blc-delta y_trc+delta , where delta is say 0.01
	index_poly=polyfillv(list_verts(0,*),list_verts(1,*),$
		naxis(0),naxis(1))
	index_checkedge=where(index_poly ge 0 and $
		index_poly le (naxis(0)*naxis(1))-1,count_checkedge)
	if count_checkedge le 0 then begin
		print,'ROI for analysis is apparently outside the image!'
		stop
	endif else begin
		index_poly=index_poly(index_checkedge)
	endelse

; modify existing (analysis) mask or create a new one...
	if use_bimage then begin
		bimage(index_poly)=bimage(index_poly)*1e3
		index_temp=where(bimage lt 1e3,count_temp)
		if count_temp ge 1 then bimage(index_temp)=0L
		index_temp=where(bimage eq 1e3,count_temp)
		if count_temp ge 1 then bimage(index_temp)=1L
		index_temp=where((bimage ne 0) and (bimage ne 1), count_temp)
		if count_temp ge 1 then begin
			print,'Problem incorporating ROI for analysis!' 
			stop
		endif
	endif else begin
		use_bimage=1L
		bimage=lonarr(naxis(0),naxis(1))
		bimage(index_poly)=1L
		fxhmake,bimage_header,bimage
		bimage_naxis=fxpar(bimage_header,'NAXIS*')
	endelse

;choose blc and trc appropriately
	blc=lonarr(2)
	blc(0)=max([floor(min(list_verts(0,*))),0])
	blc(1)=max([floor(min(list_verts(1,*))),0])
	trc=lonarr(2)
	trc(0)=min([ceil(max(list_verts(0,*)))-2,naxis(0)-1])
	trc(1)=min([ceil(max(list_verts(1,*)))-2,naxis(1)-1])
endif
;now the analysis ROI file has been fully interpreted

;if the user wants to use a sky ROI file, read it and interpret
if ROI_sky then begin
	openr,15,strcompress(root+'.ROI_sky.dat',/REMOVE_ALL)
	junk=strarr(2)
	readf,15,junk
	junk=''
	readf,15,junk
	nverts=round(float(strmid(junk,10,132)))
	junk=strarr(3)
	readf,15,junk
	list_verts=fltarr(2,nverts)
	readf,15,list_verts
	close,15
; need to scale list_verts for CONGRID taking place above
	if (docongrid) then begin
		list_verts(0,*)=list_verts(0,*)*float(naxis(0))/$
			float(congrid_xsize)
		list_verts(1,*)=list_verts(1,*)*float(naxis(1))/$
			float(congrid_ysize)
	endif
	index_poly=polyfillv(list_verts(0,*),list_verts(1,*),$
		naxis(0),naxis(1))

	index_checkedge=where(index_poly ge 0 and $
		index_poly le (naxis(0)*naxis(1))-1,count_checkedge)
	if count_checkedge le 0 then begin
		print,'ROI for sky is apparently outside the image!'
		stop
	endif else begin
		index_poly=index_poly(index_checkedge)
	endelse

;choose blc_sky and trc_sky appropriately
	blc_sky=lonarr(2)
	blc_sky(0)=max([floor(min(list_verts(0,*))),0])
	blc_sky(1)=max([floor(min(list_verts(1,*))),0])
	trc_sky=lonarr(2)
	trc_sky(0)=min([ceil(max(list_verts(0,*)))-2,naxis(0)-1])
	trc_sky(1)=min([ceil(max(list_verts(1,*)))-2,naxis(1)-1])
	sky_bimage=lonarr(naxis(0),naxis(1))
	sky_bimage(index_poly)=1L
	sky_bimage=sky_bimage(blc_sky(0):trc_sky(0),blc_sky(1):trc_sky(1))

endif
;now the ROI for the user's sky section has been interpreted

;if ROI files were used above, then we have the masks created/modified
;by this point, also we have blc,trc values set appropriately

;otherwise we need to manually get blc,trc pairs from the user
blctrc_analysis=0L
if ((n_elements(blc) ne 2) or (n_elements(trc) ne 2)) then begin
	print,''
	print,'BLC and TRC for analysis region, accept defaults for whole image:'
	GET_BLC:   
	ans = ''
	read, 'BLC_analysis (DEFAULT = 0,0): ', ans   
	if ans EQ '' then blc = [0,0] else begin
		blc = getopt(ans,'I')
		if N_elements(blc) NE 2 then begin  
              		message, 'ERROR - Expecting 2 scalar values',/CON
              		goto, GET_BLC     
		endif
	endelse
	GET_TRC:   
	ans = ''
	read, 'TRC_analysis (DEFAULT = naxis(0)-1,naxis(1)-1): ', ans   
	if ans EQ '' then trc = [0,0] else begin
		trc = getopt(ans,'I')
		if N_elements(trc) NE 2 then begin  
			message, 'ERROR - Expecting 2 scalar values',/CON
			goto, GET_TRC     
		endif
	endelse
	blctrc_analysis=1L
endif 

;same as above, but now for the sky section
blctrc_sky=0L
if ((n_elements(blc_sky) ne 2) or (n_elements(trc_sky) ne 2)) then begin
	print,''
	print,'DEFAULT values for the sky may cause noise to be overestimated.'
	print,'BLC and TRC for sky region, [RETURN] for defaults:'
	GET_BLC_SKY:   
	ans = ''
	read, 'BLC_sky (DEFAULT = 0,0): ', ans   
	if ans EQ '' then blc_sky = [0,0] else begin
		blc_sky = getopt(ans,'I')
		if N_elements(blc_sky) NE 2 then begin  
	              message, 'ERROR - Expecting 2 scalar values',/CON
	              goto, GET_BLC_SKY     
		endif
	endelse
	GET_TRC_SKY:   
	ans = ''
	read, 'TRC_sky (DEFAULT = naxis(0)-1,naxis(1)-1): ', ans   
	if ans EQ '' then trc_sky = [0,0] else begin
		trc_sky = getopt(ans,'I')
		if N_elements(trc_sky) NE 2 then begin  
	              message, 'ERROR - Expecting 2 scalar values',/CON
	              goto, GET_TRC_SKY     
		endif
	endelse
	print,''
	blctrc_sky=1L
endif


;check the blc, trc, blc_sky, and trc_sky values with respect to the image
if ((blc(0) eq 0) and (blc(1) eq 0) and $
	(trc(0) eq 0) and (trc(1) eq 0)) then begin
	blc=[0,0]
	trc=[naxis(0)-1,naxis(1)-1]
endif
if ((trc(0) eq 0) and (trc(1) eq 0)) then begin
	trc=[naxis(0)-1,naxis(1)-1]
endif
if ((blc_sky(0) eq 0) and (blc_sky(1) eq 0) and $
	(trc_sky(0) eq 0) and (trc_sky(1) eq 0)) then begin
	blc_sky=[0,0]
	trc_sky=[naxis(0)-1,naxis(1)-1]
endif
if ((trc_sky(0) eq 0) and (trc_sky(1) eq 0)) then begin
	trc_sky=[naxis(0)-1,naxis(1)-1]
endif

;if sky_bimage doesn't exist, then make it now (according to blc trc coords)
if (n_elements(sky_bimage) eq 0) then sky_bimage=$
	lonarr(trc_sky(0)-blc_sky(0)+1,trc_sky(1)-blc_sky(1)+1)+1


abort=0L
if (blc(0) lt 0) then begin
	message, 'ERROR - blc(0) lt 0',/CON
	abort=1L
endif
if (blc(1) lt 0) then begin
	message, 'ERROR - blc(1) lt 0',/CON
	abort=1L
endif
if (trc(0) gt naxis(0)-1) then begin
	message, 'ERROR - trc(0) gt naxis(0)-1',/CON
	abort=1L
endif
if (trc(1) gt naxis(1)-1) then begin
	message, 'ERROR - trc(1) gt naxis(1)-1',/CON
	abort=1L
endif
if (blc_sky(0) lt 0) then begin
	message, 'ERROR - blc_sky(0) lt 0',/CON
	abort=1L
endif
if (blc_sky(1) lt 0) then begin
	message, 'ERROR - blc_sky(1) lt 0',/CON
	abort=1L
endif
if (trc_sky(0) gt naxis(0)-1) then begin
	message, 'ERROR - trc_sky(0) gt naxis(0)-1',/CON
	abort=1L
endif
if (trc_sky(1) gt naxis(1)-1) then begin
	message, 'ERROR - trc_sky(1) gt naxis(1)-1',/CON
	abort=1L
endif

;if any of the blc,trc pairs didn't pass the tests above, GET_REGIONS_AGAIN
if (abort gt 0) then begin
	blc=0L
	trc=0L
	blc_sky=0L
	trc_sky=0L
	if not(interact) then return
	goto,GET_REGIONS_AGAIN
endif

delta=0.01
;if the analysis region was defined using blc,trc then
;make a fake ROI file to accomplish the same choice
if blctrc_analysis then begin
	openw,1,strcompress(root+'.ROI_analysis.dat',/REMOVE_ALL)
	printf,1,'; Region of Interest (ROI) Dump File'
	printf,1,'; Polygon Information:'
	printf,1,'NVERTS  = 4'
	printf,1,'NPIXELS = '
	printf,1,'COUNTS  = '
	printf,1,'; List of Vertices follows:'
	printf,1,blc(0)+1-delta,blc(1)+1-delta,format='(2(e12.5))'
	printf,1,trc(0)+1+delta,blc(1)+1-delta,format='(2(e12.5))'
	printf,1,trc(0)+1+delta,trc(1)+1+delta,format='(2(e12.5))'
	printf,1,blc(0)+1-delta,trc(1)+1+delta,format='(2(e12.5))'
	close,1
endif
;likewise for the sky region
if blctrc_sky then begin
	openw,1,strcompress(root+'.ROI_sky.dat',/REMOVE_ALL)
	printf,1,'; Region of Interest (ROI) Dump File'
	printf,1,'; Polygon Information:'
	printf,1,'NVERTS  = 4'
	printf,1,'NPIXELS = '
	printf,1,'COUNTS  = '
	printf,1,'; List of Vertices follows:'
	printf,1,blc_sky(0)+1-delta,blc_sky(1)+1-delta,format='(2(e12.5))'
	printf,1,trc_sky(0)+1+delta,blc_sky(1)+1-delta,format='(2(e12.5))'
	printf,1,trc_sky(0)+1+delta,trc_sky(1)+1+delta,format='(2(e12.5))'
	printf,1,blc_sky(0)+1-delta,trc_sky(1)+1+delta,format='(2(e12.5))'
	close,1
endif

return
end

;------------------

pro get_find_inputs, batch,interact,dist_mpc,use_coords, scale, cdelt, pc_pix,$
	scalehc, scalec, scaleh, scaleh_em, hs_filename,h_filename,c_filename,$
	cscale, ghc, gc, nhc, nc, sigma_find, PSFfwhm_pix, size_max_pc,$
	termgrad_arr, multiple_cutoffs, correction, NIIpercent

GET_FIND_INPUTS_AGAIN:

;ask for the galaxy distance, if not already known
if (n_elements(dist_mpc) ne 1) then begin
	print,''
	GET_DIST_MPC:   
	ans = ''
	read, 'Distance to the galaxy? (Mpc, NO DEFAULT): ',Ans
	if ans EQ '' then goto, GET_DIST_MPC else begin
		dist_mpc = getopt(ans,'F')
		if N_elements(dist_mpc) NE 1 then begin  
              		message, 'ERROR - Expecting floating-point scalar',/CON
              		goto, GET_DIST_MPC     
		endif
	endelse
endif 
print,''

;figure out the plate scale or ask for it from the user
if (n_elements(scale) ne 1) then begin
	if use_coords then begin
		print, 'Using the plate scale of the image ',$
			'header coordinate solution.'
		print, 'scale = ',avg(abs(cdelt))*3600.,' "/pix'
	endif else begin
	GET_SCALE:   
	ans = ''
	read, 'Plate scale of the image? ("/pix, DEFAULT=0.7): ',Ans
	if ans EQ '' then scale=0.7 else begin
		scale = getopt(ans,'F')
		if N_elements(scale) NE 1 then begin  
              		message, 'ERROR - Expecting floating-point scalar',/CON
              		goto, GET_SCALE     
		endif
	endelse
	endelse
	print,''
endif
if (n_elements(scale) eq 1) and use_coords then begin
	print,'Ignoring supplied scale in favor of the image ',$
		'header coordinate solution.'
	print, 'scale = ',avg(abs(cdelt))*3600.,' "/pix'
	print,''
endif
if use_coords then scale=avg(abs(cdelt))*3600.

;compute the number of pc per pixel
pc_pix=((dist_mpc*1.e6)*(1./3600.)*(acos(-1.)/180.))*scale

;by default, assume that the input images are in ADU,
;but ask for confirmation

if (n_elements(scaleh_em) ne 1) then begin
	print,''
	GET_SCALEH_EM:   
	ans = ''
	print, 'Scale factor required to convert the continuum-sub'
	print, '(line-only) image from the current units to EM?'
	read,  '(DEFAULT = 1.0 -> none required): ',ans   
	if ans EQ '' then scaleh_em=1. else begin
		scaleh_em = getopt(ans,'F')
		if N_elements(scaleh_em) NE 1 then begin  
              		message, 'ERROR - Expecting floating-point scalar',/CON
              		goto, GET_SCALEH_EM     
		endif
	endelse
endif 

if (n_elements(NIIpercent) ne 1) then begin
	print,''
	GET_NIIPERCENT:   
	ans = ''
	print, 'Percent contribution from [NII]?'
	read,  '(DEFAULT = 0.0 -> not in bandpass): ',ans   
	if ans EQ '' then NIIpercent=0. else begin
		NIIpercent = getopt(ans,'F')
		if N_elements(NIIpercent) NE 1 then begin  
              		message, 'ERROR - Expecting floating-point scalar',/CON
              		goto, GET_NIIPERCENT
		endif
		;maximum allowed contrib to observed counts is 0.5 of H-alpha
                if (NIIpercent lt 0. or NIIpercent gt 50.) then begin 
                        message, 'ERROR - Expecting (0. <= scalar <= 50.)',/CON 
                        goto, GET_NIIPERCENT 
                endif 
	endelse
endif 

if (n_elements(scaleh) ne 1) then begin
	print,''
	GET_SCALEH:   
	ans = ''
	print, 'Scale factor applied to the continuum-sub (line-only) image'
	print, 'during flux calibration to convert from ADU'
	read,  'to the current units? (DEFAULT = 1.0 -> none applied): ',ans   
	if ans EQ '' then scaleh=1. else begin
		scaleh = getopt(ans,'F')
		if N_elements(scaleh) NE 1 then begin  
              		message, 'ERROR - Expecting floating-point scalar',/CON
              		goto, GET_SCALEH     
		endif
	endelse
endif
 
if (n_elements(scalehc) ne 1) then begin
	print,''
	GET_SCALEHC:   
	ans = ''
	print, 'Scale factor applied to the line+continuum image'
	print, 'during flux calibration to convert from ADU'
	read,  'to the current units? (DEFAULT = 1.0 -> none applied): ',ans   
	if ans EQ '' then scalehc=1. else begin
		scalehc = getopt(ans,'F')
		if N_elements(scalehc) NE 1 then begin  
              		message, 'ERROR - Expecting floating-point scalar',/CON
              		goto, GET_SCALEHC     
		endif
	endelse
endif 

if (n_elements(scalec) ne 1) then begin
	print,''
	GET_SCALEC:   
	ans = ''
	print, 'Scale factor applied to the continuum image'
	print, 'during flux calibration to convert from ADU'
	read,  'to the current units? (DEFAULT = 1.0 -> none applied): ',ans   
	if ans EQ '' then scalec=1. else begin
		scalec = getopt(ans,'F')
		if N_elements(scalec) NE 1 then begin  
              		message, 'ERROR - Expecting floating-point scalar',/CON
              		goto, GET_SCALEC     
		endif
	endelse
endif 

; READ SUBSAMPLED IMAGES TO DETERMINE CONTINUUM SCALING
; FACTOR THAT WAS USED TO PRODUCE hs_filename
; ALSO SERVES AS TEST THAT ALL FILES EXIST
openr,13,hs_filename
fxhread,13, hs_head
hs_naxis=fxpar(hs_head,'NAXIS*')
close,13
sample=floor(float(min(hs_naxis))/32.)
fxread,hs_filename,hs_img,hs_head,0,hs_naxis(0)-1,0,hs_naxis(1)-1,sample
index_use=where(hs_img gt 0. and finite(hs_img) eq 1,count_use)
if count_use eq 0 then begin
	print,'No qualifying pixels for estimation of the cont-sub factor!'
	stop
endif
fxread,h_filename,h_img,h_head,0,hs_naxis(0)-1,0,hs_naxis(1)-1,sample
fxread,c_filename,c_img,c_head,0,hs_naxis(0)-1,0,hs_naxis(1)-1,sample
corr_diff=1.e6
for cscale=0.,2.,0.1 do begin
	corr_diff_cur=abs(1.-correlate(h_img(index_use)/scalehc-$
		cscale*c_img(index_use)/scalec,hs_img(index_use)/scaleh))
	if (corr_diff_cur le corr_diff)  then begin
		corr_diff=corr_diff_cur
		cscale_best=cscale
	endif
endfor
cscale_target=cscale_best
corr_diff=1.e6 
for cscale=cscale_target-0.1,cscale_target+0.1,0.01 do begin
	corr_diff_cur=abs(1.-correlate(h_img(index_use)/scalehc-$
		cscale*c_img(index_use)/scalec,hs_img(index_use)/scaleh))
        if (corr_diff_cur le corr_diff)  then begin 
                corr_diff=corr_diff_cur 
                cscale_best=cscale 
        endif 
endfor 
cscale_target=cscale_best 
corr_diff=1.e6 
for cscale=cscale_target-0.01,cscale_target+0.01,0.001 do begin 
	corr_diff_cur=abs(1.-correlate(h_img(index_use)/scalehc-$
		cscale*c_img(index_use)/scalec,hs_img(index_use)/scaleh))
        if (corr_diff_cur le corr_diff)  then begin 
                corr_diff=corr_diff_cur 
                cscale_best=cscale 
        endif 
endfor 
cscale_target=cscale_best 
corr_diff=1.e6 
for cscale=cscale_target-1.e-4,cscale_target+1.e-4,1.e-5 do begin 
	corr_diff_cur=abs(1.-correlate(h_img(index_use)/scalehc-$
		cscale*c_img(index_use)/scalec,hs_img(index_use)/scaleh))
        if (corr_diff_cur le corr_diff)  then begin 
                corr_diff=corr_diff_cur 
                cscale_best=cscale 
        endif 
endfor 
cscale_target=cscale_best 
corr_diff=1.e6 
for cscale=cscale_target-1.e-5,cscale_target+1.e-5,1.e-6 do begin 
	corr_diff_cur=abs(1.-correlate(h_img(index_use)/scalehc-$
		cscale*c_img(index_use)/scalec,hs_img(index_use)/scaleh))
        if (corr_diff_cur le corr_diff)  then begin 
                corr_diff=corr_diff_cur 
                cscale_best=cscale 
        endif 
endfor 

cscale=cscale_best

maxdiff=max(hs_img(index_use)/scaleh-(h_img(index_use)/scalehc-$
	cscale*c_img(index_use)/scalec))
mindiff=min(hs_img(index_use)/scaleh-(h_img(index_use)/scalehc-$
	cscale*c_img(index_use)/scalec))
correction1=(maxdiff+mindiff)/2.
maxdiff=max(hs_img(index_use)/scaleh-(h_img(index_use)/scalehc-$
	cscale*c_img(index_use)/scalec+correction1))
mindiff=min(hs_img(index_use)/scaleh-(h_img(index_use)/scalehc-$
	cscale*c_img(index_use)/scalec+correction1))
correction2=(maxdiff+mindiff)/2.
correction=correction1+correction2

;having deterimined cscale and correction, it is now possible to
;accurately reconstruct the line+continuum image from only the
;continuum-sub and continuum images (to save memory)

;in order to get the line+continuum data in ADU...
;hs_img/scaleh+cscale*c_img/scalec-correction ---> equivalent to h_img/scalehc

;there may be slight discrepancies at the 1.e-3 of ~100. level
;due to slight error in the estimated continuum scaling

;one might also want to recover the line+continuum data in supplied units
;(hs_img/scaleh+cscale*c_img/scalec-correction)*scalehc --> equiv to h_img

;by default, assume that the header contains an accurate gain value
;and number of images combined, but ask for confirmation
if interact then begin

	print,''
	GET_GHC:   
	ans = ''
	prompttxt=strcompress('Gain for line+continuum image? (DEFAULT = '+$
		string(ghc)+'): ')
	read, ans, PROMPT=prompttxt
	if ans NE '' then begin
		ghc = getopt(ans,'F')
		if N_elements(ghc) NE 1 then begin  
              		message, 'ERROR - Expecting floating-point scalar',/CON
              		goto, GET_GHC     
		endif
	endif

	print,''
	GET_GC:   
	ans = ''
	prompttxt=strcompress('Gain for continuum image? (DEFAULT = '+$
		string(gc)+'): ')
	read, ans, PROMPT=prompttxt
	if ans NE '' then begin
		gc = getopt(ans,'F')
		if N_elements(gc) NE 1 then begin  
              		message, 'ERROR - Expecting floating-point scalar',/CON
              		goto, GET_GC     
		endif
	endif

	print,''
	GET_NHC:   
	ans = ''
	prompttxt=strcompress('Number of images combined to create '+$
		 'line+continuum image? (DEFAULT = '+string(nhc)+'): ')
	read, ans, PROMPT=prompttxt
	if ans NE '' then begin
		nhc = getopt(ans,'F')
		if N_elements(nhc) NE 1 then begin  
              		message, 'ERROR - Expecting floating-point scalar',/CON
              		goto, GET_NHC     
		endif
	endif

	print,''
	GET_NC:   
	ans = ''
	prompttxt=strcompress('Number of images combined to create '+$
		 'continuum image? (DEFAULT = '+string(nc)+'): ')
	read, ans, PROMPT=prompttxt
	if ans NE '' then begin
		nc = getopt(ans,'F')
		if N_elements(nc) NE 1 then begin  
              		message, 'ERROR - Expecting floating-point scalar',/CON
              		goto, GET_NC     
		endif
	endif
endif 

;what is the lowest acceptable S/N (this limit is imposed twice
;during execution of the code: 1-during convolution based FIND
;at original and reduced resolutions, 2-after FIND peaks are tabulated
;and during the explicit pixel-by-pixel S/N computation with
;original data).
if (n_elements(sigma_find) ne 1) then begin
	print,''
	GET_SIGMA_FIND:   
	ans = ''
	read, 'Minimum acceptable S/N? (DEFAULT = 5.0): ', ans   
	if ans EQ '' then sigma_find=5. else begin
		sigma_find = getopt(ans,'F')
		if N_elements(sigma_find) NE 1 then begin  
              		message, 'ERROR - Expecting floating-point scalar',/CON
              		goto, GET_SIGMA_FIND     
		endif
	endelse
endif 

;ask for the fwhm of the PSF in pixels
if (n_elements(PSFfwhm_pix) ne 1) then begin
	print,''
	GET_PSFFWHM_PIX:   
	ans = ''
	read, 'FWHM of image PSF? (DEFAULT = 2.0 pixels): ', ans   
	if ans EQ '' then PSFfwhm_pix=2. else begin
		PSFfwhm_pix = getopt(ans,'F')
		if N_elements(PSFfwhm_pix) NE 1 then begin  
              		message, 'ERROR - Expecting floating-point scalar',/CON
              		goto, GET_PSFFWHM_PIX     
		endif
	endelse
endif 

;what is the largest characteristic size of structures to be searched for?
if (n_elements(size_max_pc) ne 1) then begin
	print,''
	GET_SIZE_MAX_PC:   
	ans = ''
	read, 'FWHM of largest HII region (-1 for single kernel)? (DEFAULT = 500.0 pc): ', ans   
	if ans EQ '' then size_max_pc=500. else begin
		size_max_pc = getopt(ans,'F')
		if N_elements(size_max_pc) NE 1 then begin  
              		message, 'ERROR - Expecting floating-point scalar',/CON
              		goto, GET_SIZE_MAX_PC     
		endif
	endelse
endif 
print,''

;if we are running an interactive session of HIIphot (not HIIphot_batch)
;then ask for the desired terminal gradients of HII region surface brightness
;profiles
if interact and not(batch) then begin 
        GET_TERMGRAD:    
        ans = '' 
        read,  'Terminal gradient (EM/pc)? (DEFAULT = 1.5,0,0,0,0) : ',ans   
        if ans EQ '' then termgrad_arr=[1.5,0.0,0.0,0.0,0.0] else begin 
                termgrad_arr = getopt(ans,'F') 
                if N_elements(termgrad_arr) NE 5 then begin   
                        message, 'ERROR - Expecting 5 floating-point scalars',/CON 
                        goto, GET_TERMGRAD      
                endif 
                if (min(termgrad_arr) lt 0.) then begin 
                        message, 'ERROR - Expecting all scalars >= 0.',/CON 
                        goto, GET_TERMGRAD 
                endif 
        endelse 
	if max(termgrad_arr(1:4)) gt 1.e-6 then multiple_cutoffs=1L else $
		multiple_cutoffs=0L
	print,''
endif

return
end

pro corr_statistics_2d,Array_X,Array_Y,Array_v,num_pix,noise_level,$
        betahat0,betahat1,rhoval,Xsum,XXsum,Xbar,SLOW=slow
 
;this pro computes the noise-corrected linear correlation coefficient
;between model and data. this becomes wildly uncertain if the observed
;variance in the data is low enough to become comparable to the anticipated
;variance (based on our noise model). in such a case, just use the
;canonical (uncorrected for noise) correlation coefficient.
;in any case, return -1. for the coefficient if it is computed to
;be any sort of infinite value (+Inf,-Inf,NaN,etc).

; compute a few sums used below
if keyword_set(SLOW) then begin
        Xsum=total(Array_X)
        XXsum=total(Array_X^2)
        Xbar=Xsum/num_pix
endif
Ysum=total(Array_Y)
YXsum=total(Array_Y*Array_X)
YYsum=total(Array_Y^2)
vsum=total(Array_v)
 
; now calculate Ybar, betahat0, betahat1
Ybar=Ysum/num_pix
betahat1=(YXsum-(Ysum*Xsum)/num_pix)/(XXsum-(Xsum^2)/num_pix)
betahat0=Ybar-betahat1*Xbar  
 
if ((YYsum-num_pix*Ybar^2)/vsum ge 4.0) then begin
        rhoval=(YXsum-num_pix*Xbar*Ybar)/$
                sqrt((XXsum-num_pix*Xbar^2)*(YYsum-num_pix*Ybar^2-vsum))
endif else begin
        rhoval=correlate(Array_X,Array_Y)
endelse
 
if (finite(rhoval) ne 1) then rhoval=-1.

return
end

;--------------------

pro sym_model,n,r0,s,a,NORMALIZE=normalize

;generate a symmetrical model according to the prescription described
;in Thilker, Braun, Walterbos (1999)
 
if keyword_set(NORMALIZE) then normalize=1L else normalize=0L
 
if (s lt 2.0/2.35) then begin
	print,'Warning: undersampled model generated!'
endif

a=fltarr(ceil(n),ceil(n))
naxis=size(a)
c=n/2
 
p=findgen(n)-c
 
for i=0,n-1 do begin
for j=0,n-1 do begin
        d=sqrt(p(i)^2+p(j)^2)-r0
        a(i,j)=exp(-0.5*(d/s)^2)
endfor
endfor
 
if normalize then a=a/total(a)
 
end

;------------------------

pro distort,a,b,r,t

;stretch and rotate an input image (a) according to r and t
 
r=float(r)
naxis=size(a)
n1=naxis(1)
n2=naxis(2)
c1=n1/2
c2=n2/2
m1=sqrt(r)
m2=sqrt(1/r)
 
b=interpolate(a,(findgen(n1)-c1)*m1+c1,$
        (findgen(n2)-c2)*m2+c2,/GRID,MISSING=0.)
 
b=rot(temporary(b),-t,1.,/INTERP,MISSING=0.)
 
end

;-------------------------------

pro update_map_files,root,root_suffix,nreg,sort_order

;this code fixes up integer maps in the case where regions
;have been resorted into a different order than the maps were
;originally written

if root_suffix eq '.a' then begin

print,'Updating sorting order of detections in: ',$
	strcompress(root+'.FOOT.MAP.fits',/REMOVE_ALL)
updatefile=strcompress(root+'.FOOT.MAP.fits',/REMOVE_ALL)

fxread,updatefile,map
map2=map

for i=long(0),long(nreg)-long(1),long(1) do begin
	i_sorted=sort_order(i)
	index_replace=where(map eq i_sorted+1,count)
	if count gt 0 then begin
		map2(index_replace)=i+1
	endif
endfor

if max(map) le 32767 then begin
	map2=fix(map2)
	fxhmake,head,map2
	fxwrite,updatefile,head,map2
endif else begin
	print,'Using long integer FITS format...'
	fxhmake,head,map2
	fxwrite,updatefile,head,map2
endelse

print,'Updating sorting order of detections in: ',$
	strcompress(root+'.SEED.MAP.fits',/REMOVE_ALL)
updatefile=strcompress(root+'.SEED.MAP.fits',/REMOVE_ALL)

fxread,updatefile,map
map2=map

for i=long(0),long(nreg)-long(1),long(1) do begin
	i_sorted=sort_order(i)
	index_replace=where(map eq i_sorted+1,count)
	if count gt 0 then begin
		map2(index_replace)=i+1
	endif
endfor

if max(map) le 32767 then begin 
        map2=fix(map2) 
        fxhmake,head,map2 
        fxwrite,updatefile,head,map2 
endif else begin 
        print,'Using long integer FITS format...' 
        fxhmake,head,map2 
        fxwrite,updatefile,head,map2 
endelse 

endif

print,'Updating sorting order of detections in: ',$
	strcompress(root+'.GROW.MAP.fits',/REMOVE_ALL)
updatefile=strcompress(root+'.GROW.MAP.fits',/REMOVE_ALL)

fxread,updatefile,map
map2=map

for i=long(0),long(nreg)-long(1),long(1) do begin
	i_sorted=sort_order(i)
	index_replace=where(map eq i_sorted+1,count)
	if count gt 0 then begin
		map2(index_replace)=i+1
	endif
endfor

if max(map) le 32767 then begin 
        map2=fix(map2) 
        fxhmake,head,map2 
        fxwrite,updatefile,head,map2 
endif else begin 
        print,'Using long integer FITS format...' 
        fxhmake,head,map2 
        fxwrite,updatefile,head,map2 
endelse 


end

;-----------------------
pro change_filenames,root,root_suffix

;this pro renames a suite of output files in the case of multiple_cutoffs

if exist(strcompress(root+'.catalog.dat',/REMOVE_ALL)) then begin
	spawncmd='\mv '+strcompress(root+'.catalog.dat',/REMOVE_ALL)+' '+strcompress(root+root_suffix+'.catalog.dat',/REMOVE_ALL)
	spawn,spawncmd,spawnout
endif
if exist(strcompress(root+'.GROW+BORDER.fits',/REMOVE_ALL)) then begin
	spawncmd='\mv '+strcompress(root+'.GROW+BORDER.fits',/REMOVE_ALL)+' '+strcompress(root+root_suffix+'.GROW+BORDER.fits',/REMOVE_ALL)
	spawn,spawncmd,spawnout
endif
if exist(strcompress(root+'.CONT.GROW+BORDER.fits',/REMOVE_ALL)) then begin
	spawncmd='\mv '+strcompress(root+'.CONT.GROW+BORDER.fits',/REMOVE_ALL)+' '+strcompress(root+root_suffix+'.CONT.GROW+BORDER.fits',/REMOVE_ALL)
	spawn,spawncmd,spawnout
endif
if exist(strcompress(root+'.BKGD+BORDER.fits',/REMOVE_ALL)) then begin
	spawncmd='\mv '+strcompress(root+'.BKGD+BORDER.fits',/REMOVE_ALL)+' '+strcompress(root+root_suffix+'.BKGD+BORDER.fits',/REMOVE_ALL)
	spawn,spawncmd,spawnout
endif
if exist(strcompress(root+'.BKGD.fits',/REMOVE_ALL)) then begin
	spawncmd='\mv '+strcompress(root+'.BKGD.fits',/REMOVE_ALL)+' '+strcompress(root+root_suffix+'.BKGD.fits',/REMOVE_ALL)
	spawn,spawncmd,spawnout
endif
if exist(strcompress(root+'.GROW.MAP.fits',/REMOVE_ALL)) then begin
	spawncmd='\mv '+strcompress(root+'.GROW.MAP.fits',/REMOVE_ALL)+' '+strcompress(root+root_suffix+'.GROW.MAP.fits',/REMOVE_ALL)
	spawn,spawncmd,spawnout
endif
if exist(strcompress(root+'.STAMP.fits',/REMOVE_ALL)) then begin
	spawncmd='\mv '+strcompress(root+'.STAMP.fits',/REMOVE_ALL)+' '+strcompress(root+root_suffix+'.STAMP.fits',/REMOVE_ALL)
	spawn,spawncmd,spawnout
endif
if exist(strcompress(root+'.CONT.STAMP.fits',/REMOVE_ALL)) then begin
	spawncmd='\mv '+strcompress(root+'.CONT.STAMP.fits',/REMOVE_ALL)+' '+strcompress(root+root_suffix+'.CONT.STAMP.fits',/REMOVE_ALL)
	spawn,spawncmd,spawnout
endif
if exist(strcompress(root+'.growvar.dat',/REMOVE_ALL)) then begin
	spawncmd='\mv '+strcompress(root+'.growvar.dat',/REMOVE_ALL)+' '+strcompress(root+root_suffix+'.growvar.dat',/REMOVE_ALL)
	spawn,spawncmd,spawnout
endif
return
end
;-----------------------

pro estimate_error_beg,hc,c,s,ghc,gc,nhc,nc,scalehc,scalec,offhc,offc,$
	hcbkgd,cbkgd,df,hcmean,cmean

;assume that hc,c,offhc,offc,hcbkgd,cbkgd are all expressed in ADU
;df will be returned in units of ADU

if n_elements(hc) ne n_elements(c) then begin
	print,'ERROR - ESTIMATE_ERROR - hc must the same size as c!!!'
	stop
endif

hc_positive=(hc-offhc)>(hcmean-offhc)		;ADU
c_positive=(c-offc)>(cmean-offc)		;ADU
hcbkgd_positive=(hcbkgd-offhc)>(hcmean-offhc)	;ADU
cbkgd_positive=(cbkgd-offc)>(cmean-offc)	;ADU

db=sqrt(ghc/nhc*hcbkgd_positive+s^2*gc/nc*cbkgd_positive) ;electrons

df=sqrt(total(ghc/nhc*hc_positive)+$
	s^2*total(gc/nc*c_positive)+$
	n_elements(hc)*db^2)/ghc	;ADU

return
end

pro estimate_error_end,hc,c,s,ghc,gc,nhc,nc,scalehc,scalec,offhc,offc,db,df,$
	hcmean,cmean

;assume that hc,c,offhc,offc are all expressed in ADU, but db in electrons
;df will be returned in units of ADU

if n_elements(hc) ne n_elements(c) then begin
	print,'ERROR - ESTIMATE_ERROR - hc must the same size as c!!!'
	stop
endif

hc_positive=(hc-offhc)>(hcmean-offhc)	;ADU
c_positive=(c-offc)>(cmean-offc)	;ADU

df=sqrt(total(ghc/nhc*hc_positive)+$
	s^2*total(gc/nc*c_positive)+$
	total(db^2))/ghc		;ADU

return
end

;----------------------

pro drawborders,bv,index_bimage,map,naxis,bigimage

bv=bv+1.

; if use_bimage (intentional BLANKing) is being employed, then draw a border
; around BLANKed regions in addition to HII regions
if max(index_bimage) ge 0 then begin
	for i=long(0),long(n_elements(index_bimage))-long(1) do begin
		if (map(index_bimage(i)) eq 0) then $
			map(index_bimage(i))=2^14
	endfor
endif

for i=0,naxis(0)-1 do begin
index=where(map(i,1:naxis(1)-1) ne map(i,0:naxis(1)-2), count)
if (count gt 0) then begin
	bigimage(2*i+1,2*index+2)=bv-1.
endif
endfor

for i=0,naxis(1)-1 do begin
index=where(map(1:naxis(0)-1,i) ne map(0:naxis(0)-2,i), count)
if (count gt 0) then begin
	bigimage(2*index+2,2*i+1)=bv-1.
endif
endfor

for j=2,2*naxis(1)-1,2 do begin
for i=2,2*naxis(0)-1,2 do begin
	if ((bigimage(i-1,j) lt bv-0.5) and $
		(bigimage(i+1,j) lt bv-0.5)) then bigimage(i,j)=bv-1.
endfor
endfor

for i=2,2*naxis(0)-1,2 do begin
for j=2,2*naxis(1)-1,2 do begin
	if ((bigimage(i,j-1) lt bv-0.5) and $
		(bigimage(i,j+1) lt bv-0.5)) then bigimage(i,j)=bv-1.
endfor
endfor

for i=2,2*naxis(0)-1,2 do begin
for j=2,2*naxis(1)-1,2 do begin
	k=0
	if (bigimage(i-1,j) lt bv-0.5) then k=k+1
	if (bigimage(i+1,j) lt bv-0.5) then k=k+1
	if (bigimage(i,j-1) lt bv-0.5) then k=k+1
	if (bigimage(i,j+1) lt bv-0.5) then k=k+1
	if k ge 2 then bigimage(i,j)=bv-1.
endfor
endfor

;note-- boundaries falling on the edge of the image will
;       not be drawn at present

; if use_bimage (intentional BLANKing) is being used, we need to undo
; modifications made above to map
if max(index_bimage) ge 0 then begin
	for i=long(0),long(n_elements(index_bimage))-long(1) do begin
		if (map(index_bimage(i)) eq 2^14) then $
			map(index_bimage(i))=0L
	endfor
endif

bv=bv-1.

end

;--------------------------
pro growth_set,i,map,image,naxis,index_out,values_out,index_p,values_p,$
	D_CRIT=d_crit,GROWING=growing

;select a set of pixels adjacent or diagonal to region i if D_CRIT is
;not set, otherwise select those pixels within d_crit pixels

if keyword_set(D_CRIT) then d_crit=d_crit else d_crit=sqrt(2.)+0.001
if keyword_set(GROWING) then doperim=1L else doperim=0L

d2_crit=d_crit^2

; get a subsection containing the region of interest, plus extra...
index_temp=where(map eq i+1,count)
if (count lt 1) then  goto, NO_PIXELS_IN_REG 

max_x=-1
max_y=-1
min_x=naxis(0)
min_y=naxis(1)
for counter=0L,long(count)-1L,1L do begin
	cur_pos=index2coord(index_temp(counter),naxis(0))
	if cur_pos(0) gt max_x then max_x=cur_pos(0)
	if cur_pos(1) gt max_y then max_y=cur_pos(1)
	if cur_pos(0) lt min_x then min_x=cur_pos(0)
	if cur_pos(1) lt min_y then min_y=cur_pos(1)
endfor

padding=ceil(d_crit)
min_x=min_x-padding
min_y=min_y-padding
max_x=max_x+padding
max_y=max_y+padding
x_vector=findgen(max_x-min_x+1)+min_x
y_vector=findgen(max_y-min_y+1)+min_y

sub_image=interpolate(image,x_vector,y_vector,/GRID,MISSING=0.)
sub_map=interpolate(map,x_vector,y_vector,/GRID,MISSING=-1)

num_perim=0L
index_p=lonarr(n_elements(sub_image))
values_p=fltarr(n_elements(sub_image))

index_unclaimed=where(sub_map eq 0,count_unclaimed)
index_region=where(sub_map eq i+1,count_region)
region_pos=lonarr(2,count_region)
for counter=0L,long(count_region)-1L,1L do begin
	cur_pos=index2coord(index_region(counter),max_x-min_x+1)
	region_pos(*,counter)=cur_pos
	if doperim then begin
		index_check=index_region(counter)+[-(max_x-min_x+1)+[-1,0,1],$
			[-1,1],(max_x-min_x+1)+[-1,0,1]]
		index_onimage=where(index_check ge 0 and $
			index_check le ((max_x-min_x+1)*(max_y-min_y+1)-1))
		index_check=index_check(index_onimage)
		values_neighboring=sub_map(index_check)
		index_zero=where(values_neighboring eq 0,count_zero)
		if count_zero gt 0 then begin
			index_p(num_perim)=index_region(counter)
			values_p(num_perim)=sub_image(index_region($
				counter))
			num_perim=num_perim+1
		endif
	endif
endfor
if doperim then begin
	if num_perim gt 0 then begin
		index_p=index_p(0:num_perim-1)
		values_p=values_p(0:num_perim-1)
	endif else begin
		index_p=0L
		values_p=0L
	endelse
endif

index_out=lonarr(n_elements(sub_image))
values_out=fltarr(n_elements(sub_image))
n_set=0L
for counter=0L,long(count_unclaimed)-1L,1L do begin
	cur_pos=index2coord(index_unclaimed(counter),max_x-min_x+1)	
	d2=(cur_pos(0)-region_pos(0,*))^2+(cur_pos(1)-region_pos(1,*))^2
	if min(d2) le d2_crit then begin
		coord_out=coord2coord(cur_pos,[min_x,min_y])
		index_out(n_set)=coord2index(coord_out,naxis(0))
		values_out(n_set)=sub_image(index_unclaimed(counter))
		n_set=n_set+1
	endif
endfor

if n_set gt 0 then begin
	index_out=index_out(0:n_set-1)
	values_out=values_out(0:n_set-1)
endif else begin
	index_out=0L
	values_out=0L
endelse

if doperim and num_perim gt 0 then begin
	for counter=0L,long(num_perim)-1L,1L do begin
		cur_pos=index2coord(index_p(counter),max_x-min_x+1)
		coord_p=coord2coord(cur_pos,[min_x,min_y])
		index_p(counter)=coord2index(coord_p,naxis(0))
		values_p(counter)=image(index_p(counter))
	endfor
endif

return

NO_PIXELS_IN_REG:
print,'ERROR (GROWTH_SET): Region '+string(i)+' has no pixels in map!'
stop

end

;-----------------------------------
pro echo_help

openw,7,'HIIphot.helpfile'
printf,7,'-------------------------------------------------------------------------------'
printf,7,'HELP and EXPLAIN-file for HIIphot.pro'
printf,7,'Still under construction...'
printf,7,'-------------------------------------------------------------------------------'
close,7

spawn,'more HIIphot.helpfile'
spawn,'\rm HIIphot.helpfile'
end

;--------------------

pro find_mod, image, x, y, flux, sharp, roundness,hmin,fwhm,roundlim,sharplim,$
                      PRINT = print, SILENT=silent, VERBOSE=verbose, $
                        POINTSOURCE=pointsource
;+
; NAME:
;       FIND
; PURPOSE:
;       Find positive brightness perturbations (i.e stars) in an image 
; EXPLANATION:
;       Also returns centroids and shape parameters (roundness & sharpness).
;       Adapted from 1986 STSDAS version of DAOPHOT.
;
; CALLING SEQUENCE:
;       FIND, image, [ x, y, flux, sharp, round, hmin, fwhm, roundlim, sharplim 
;               PRINT= , /SILENT , /VERBOSE]
;
; INPUTS:
;       image - 2 dimensional image array (integer or real) for which one
;               wishes to identify the stars present
;
; OPTIONAL INPUTS:
;       FIND will prompt for these parameters if not supplied
;
;       hmin -  Threshold intensity for a point source - should generally 
;               be 3 or 4 sigma above background
;       fwhm  - FWHM to be used in the convolve filter
;       sharplim - 2 element vector giving low and high cutoff for the
;               sharpness statistic (Default: [0.2,1.0] ).   Change this
;               default only if the stars have siginificantly larger or 
;               or smaller concentration than a Gaussian
;       roundlim - 2 element vector giving low and high cutoff for the
;		roundness statistic (Default: [-1.0,1.0] ).   Change this 
;		default only if the stars are significantly elongated.
;
; OPTIONAL INPUT KEYWORDS:
;	SILENT - Normally, FIND will write out each star that meets all
;		selection criteria.   If the SILENT keyword is set and 
;		non-zero, then this printout is suppressed.
;	PRINT - if set and non-zero then T_FIND will also write its results to
;		a file FIND.PRT.   Also one can specify a different output file 
;		name by setting PRINT = 'filename'.
;
; OPTIONAL OUTPUTS:
;	x - vector containing x position of all stars identified by FIND
;	y-  vector containing y position of all stars identified by FIND
;	flux - vector containing flux of identified stars as determined
;		by a gaussian fit.  Fluxes are NOT converted to magnitudes.
;	sharp - vector containing sharpness statistic for identified stars
;	round - vector containing roundness statistic for identified stars
;
; NOTES:
;	The sharpness statistic compares the central pixel to the mean of the
;	surrounding pixels.   If this difference is greater than the originally
;	estimated height of the Gaussian or less than 0.2 the height of the
;	Gaussian (for the default values of SHARPLIM) then the star will be
;	rejected. 
;
; PROCEDURE CALLS:
;	DATATYPE(), GETOPT
; REVISION HISTORY:
;	Written W. Landsman, STX  February, 1987
;	ROUND now an internal function in V3.1   W. Landsman July 1993
;	Change variable name DERIV to DERIVAT    W. Landsman Feb. 1996
;	Use /PRINT keyword instead of TEXTOUT    W. Landsman May  1996
;	Changed loop indices to type LONG       W. Landsman Aug. 1997
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
;
 On_error,2                         ;Return to caller

 npar   = N_params()
 if npar EQ 0 then begin
    print,'Syntax - FIND, image,' + $
          '[ x, y, flux, sharp, round, hmin, fwhm, roundlim, sharplim'
    print,'                      PRINT = , /SILENT ]'
    return
 endif

;;; increased BOX SIZE in this modified version   (D. Thilker)
 maxbox = 251 	;Maximum size of convolution box in pixels 

;Determine if hardcopy output is desired

 type = size(image)
 if ( type[0] NE 2 ) then message, $
     'ERROR - Image array (first parameter) must be 2 dimensional'
 n_x  = type[1] & n_y = type[2]

 if not keyword_set( PRINT )  then print = 0
 if not keyword_set( SILENT ) then silent = 0
 if not keyword_set( VERBOSE ) then verbose = 0
 if ( N_elements(fwhm) NE 1 ) then $
           read, 'Enter approximate FWHM: ', fwhm

 radius = 0.637*FWHM > 2.001             ;Radius is 1.5 sigma
 radsq = radius^2
 nhalf = fix(radius) < (maxbox-1)/2   	;
 nbox = 2*nhalf + 1	;# of pixels in side of convolution box 
 middle = nhalf          ;Index of central pixel

 lastro = n_x - nhalf
 lastcl = n_y - nhalf
 sigsq = ( fwhm/2.35482 )^2
 mask = bytarr( nbox, nbox )   ;Mask identifies valid pixels in convolution box 
 c = fltarr( nbox, nbox )      ;c will contain Gaussian convolution kernel

 dd = indgen(nbox-1) + 0.5 - middle	;Constants need to compute ROUND
 dd2 = dd^2
 w = 1. - 0.5*(abs(dd)-0.5) / (middle-.5)   
 ir = (nhalf-1) > 1

 row2 = (findgen(Nbox)-nhalf)^2

 for i = 0, nhalf do begin
	temp = row2 + i^2
	c[0,nhalf-i] = temp         
        c[0,nhalf+i] = temp                           
 endfor

 mask = fix(c LE radsq)     ;MASK is complementary to SKIP in Stetson's Fortran
 good = where( mask, pixels)  ;Value of c are now equal to distance to center

 c = c*mask               
 c[good] = exp(-0.5*c[good]/sigsq)	;Make c into a gaussian kernel
 sumc = total(c)
 sumcsq = total(c^2) - sumc^2/pixels
 sumc = sumc/pixels
 c[good] = (c[good] - sumc)/sumcsq
 c1 = exp(-.5*row2/sigsq)
 sumc1 = total(c1)/nbox
 sumc1sq = total(c1^2) - sumc1
 c1 = (c1-sumc1)/sumc1sq
 sumc = total(w)                         ;Needed for centroid computation

 if N_elements(hmin) NE 1 then read, $
    'Enter minimum value above background for threshold detection: ',hmin

 if N_elements(sharplim) NE 2 then begin
      print,'Enter low and high cutoffs, press [RETURN] for defaults:'
GETSHARP:   
      ans = ''
      read, 'Image Sharpness Statistic (DEFAULT = 0.2,1.0): ', ans   
      if ans EQ '' then sharplim = [0.2,1.0] else begin
         sharplim = getopt(ans,'F')
          if N_elements(sharplim) NE 2 then begin  
              message, 'ERROR - Expecting 2 scalar values',/CON
              goto, GETSHARP     
          endif
      endelse                                                      

GETROUND: 
  ans = ''
  read, 'Image Roundness Statistic [DEFAULT = -1.0,1.0]: ',ans
  if ans EQ '' then roundlim = [-1.,1.] else begin
      roundlim = getopt( ans, 'F' )
      if N_elements( roundlim ) NE 2 then begin
           message,'ERROR - Expecting 2 scalar values',/CON
           goto, GETROUND   
      endif
 endelse
 endif 

if verbose then message,'Beginning convolution of image', /INF

 h = convol(float(image),c)    ;Convolve image with kernel "c"

    h[0:nhalf-1,*] = 0 & h[n_x-nhalf:n_x-1,*] = 0
    h[*,0:nhalf-1] = 0 & h[*,n_y-nhalf:n_y-1] = 0

if verbose then message,'Finished convolution of image', /INF

 mask[middle,middle] = 0	;From now on we exclude the central pixel
 pixels = pixels -1      ;so the number of valid pixels is reduced by 1
 good = where(mask)      ;"good" identifies position of valid pixels
 xx= (good mod nbox) - middle	;x and y coordinate of valid pixels 
 yy = fix(good/nbox) - middle    ;relative to the center
 offset = yy*n_x + xx
SEARCH: 			    ;Threshold dependent search begins here

 index = where( h GE hmin, nfound)  ;Valid image pixels are greater than hmin
 if nfound EQ 0 then begin          ;Any maxima found?

if verbose then message,'ERROR - No maxima exceed input threshold of ' + $
             string(hmin,'(F9.1)'),/CON
    goto,FINISH    

 endif

 for i= 0L, pixels-1 do begin                             

	stars = where (h[index] GE h[index+offset[i]], nfound)
        if nfound LT 0 then begin  ;Do valid local maxima exist?
             if verbose then message,$
		'ERROR - No maxima exceed input threshold of ' + $
                     string(hmin,'(F9.1)'),/CON
             goto,FINISH  
        endif
	index = index[stars]

 endfor 
 
 ix = index mod n_x              ;X index of local maxima
 iy = index/n_x                  ;Y index of local maxima
 ngood = N_elements(index)       
 if verbose then message,$
	strtrim(ngood,2)+' local maxima located above threshold',/INF

 nstar = 0L       	;NSTAR counts all stars meeting selection criteria
 badround = 0 & badsharp=0  &  badcntrd=0
 if (npar GE 2) or (PRINT) then begin 	;Create output X and Y arrays? 
  	x = fltarr(ngood) & y = x
 endif

 if (npar GE 4) or (PRINT) then begin   ;Create output flux,sharpness arrays?
 	flux = x & sharp = x & roundness = x
 endif

 if PRINT then begin	;Create output file?

         if ( datatype(print) NE 'STR' ) then file = 'find.prt' $
                                         else file = print
         message,'Results will be written to a file ' + file,/INF
         openw,lun,file,/GET_LUN
	printf,lun,' Program: FIND '+ systime()
	printf,lun,format='(/A,F7.1)',' Threshold above background:',hmin
	printf,lun,' Approximate FWHM:',fwhm
	printf,lun,format='(2(A,F6.2))',' Sharpness Limits: Low', $
                sharplim[0], '  High',sharplim[1]
	printf,lun,format='(2(A,F6.2))',' Roundness Limits: Low', $
                roundlim[0],'  High',roundlim[1]
	printf,lun,format='(/A,i6)',' No of sources above threshold',ngood

 endif                      

 if not SILENT then $
  print,format='(/8x,a)','     STAR      X      Y     FLUX     SHARP    ROUND'

;  Loop over star positions; compute statistics

 for i = 0L,ngood-1 do begin   
     temp = float(image[ix[i]-nhalf:ix[i]+nhalf,iy[i]-nhalf:iy[i]+nhalf])
     d = h[ix[i],iy[i]]                  ;"d" is actual pixel intensity        

;  Compute Sharpness statistic

     sharp1 = (temp[middle,middle] - (total(mask*temp))/pixels)/d
     if ( sharp1 LT sharplim[0] ) or ( sharp1 GT sharplim[1] ) then begin
	badsharp = badsharp + 1
	goto, REJECT             ;Does not meet sharpness criteria
     endif

;   Compute Roundness statistic

     dx = total( total(temp,2)*c1)   
     dy = total( total(temp,1)*c1)
;intentionally disable roundness check, in an effort to detect ALL regions
	around=0
;;;     if (dx LE 0) or (dy LE 0) then begin
;;;         badround = badround + 1
;;;	 goto, REJECT           ;Cannot compute roundness
;;;     endif
;;;
;;;     around = 2*(dx-dy) / ( dx + dy )    ;Roundness statistic
     if ( around LT roundlim[0] ) or ( around GT roundlim[1] ) then begin
	badround = badround + 1
	goto,REJECT           ;Does not meet roundness criteria
     endif

	if POINTSOURCE then begin

; Find X centroid
 
     derivat = shift(temp,-1,0) - temp
     derivat = total( derivat[0:nbox-2,middle-ir:middle+ir],2)
     sumd = total(w*derivat)
     sumxd = total(w*dd*derivat)
     sumxsq = total(w*dd2) 
 
     if ( sumxd GE 0. ) then begin
        badcntrd = badcntrd + 1
        goto,REJECT           ;Cannot compute X centroid
     endif
 
     dx =sumxsq*sumd/(sumc*sumxd)
     if abs(dx) GT nhalf then begin
         badcntrd = badcntrd + 1
         goto,REJECT           ;X centroid too far from local X maxima
     endif
 
     xcen = ix[i]-dx               ;Convert back to big image coordinates
 
; Find Y centroid                 
 
     derivat = shift(temp,0,-1) - temp 
     derivat = total( derivat[middle-ir:middle+ir,0:nbox-2], 1 )
     sumd = total( w*derivat )
     sumxd = total( w*dd*derivat )
     sumxsq = total( w*dd2 )
     if (sumxd GE 0) then begin
          badcntrd = badcntrd + 1
          goto, REJECT  
     endif
 
     dy = sumxsq*sumd/(sumc*sumxd)
     if ( abs(dy) GT nhalf ) then begin
        badcntrd = badcntrd + 1
        goto,REJECT
     endif
     
     ycen = iy[i]-dy

	endif else begin ; only do centroiding on demand

	xcen=ix[i]
	ycen=iy[i]

	endelse

;  This star has met all selection criteria.  Print out and save results

   if not SILENT then $
      print,FORM = '(12x,i5,2f7.1,f9.1,2f9.2)', $ 
            nstar, xcen, ycen, d, sharp1, around

   if (npar GE 2) or (PRINT) then begin
              x[nstar] = xcen & y[nstar] = ycen
   endif

   if ( npar GE 4 ) or (PRINT) then begin
	flux[nstar] = d & sharp[nstar] = sharp1 & roundness[nstar] = around
   endif
   
   nstar = nstar+1

REJECT: 
  
 endfor

 nstar = nstar-1		;NSTAR is now the index of last star found

 if PRINT then begin
  printf,lun,' No. of sources rejected by SHARPNESS criteria',badsharp
  printf,lun,' No. of sources rejected by ROUNDNESS criteria',badround
  printf,lun,' No. of sources rejected by CENTROID  criteria',badcntrd
 endif
 
if verbose then print,$
	' No. of sources rejected by SHARPNESS criteria',badsharp
if verbose then print,$
	' No. of sources rejected by ROUNDNESS criteria',badround
if verbose then print,$
	' No. of sources rejected by CENTROID  criteria',badcntrd

  if nstar LT 0 then return               ;Any stars found?

  if (npar GE 2) or (PRINT) then begin
	x=x[0:nstar]  & y = y[0:nstar]
  endif

  if (npar GE 4) or (PRINT) then begin
	flux= flux[0:nstar] & sharp=sharp[0:nstar]  
        roundness = roundness[0:nstar]
  endif

 if PRINT then begin                
   printf,lun, $
      format = '(/8x,a)','     STAR       X       Y     FLUX     SHARP    ROUND'
	for i = 0, nstar do $
	   printf,lun,format='(12x,i5,2f8.2,f9.1,2f9.2)', $
	              i+1, x[i], y[i], flux[i], sharp[i], roundness[i]
        free_lun, lun
 endif

FINISH:

 if SILENT then return

 print,form='(A,F8.1)',' Threshold above background for this pass was',hmin
 ans = ''
 read,'Enter new threshold or [RETURN] to exit: ',ans
 ans = getopt(ans,'F')              
 if ans GT 0. then begin
       hmin = ans
       goto, SEARCH   
 endif

 return                                      
 end


