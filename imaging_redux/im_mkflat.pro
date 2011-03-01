;+
; NAME: 
;  IM_MKFLAT
;
; PURPOSE: 
;  Collect and combine CCD flat frames into a superflat frame
;
; DESCRIPTION:
;
;  The files are assumed to be named as 'root'NNN where 'root' is the
;  root of the file name that you supply and NNN is a three digit number.
;  If your file name has an imbedded '.' then add it to the root.
;
;  The specified range of files are all read in from FITS files.  Then each
;  image has the overscan mean subtracted (if desired), cropped (as indicated),
;  input bias image is subtracted, and then the dark is multipled by each
;  image's exposure time and then subtracted.  These images are then
;  normalized using the mean of the SCALE region and then averaged weighting
;  by the signal-to-noise ratio of each flat image.  The weighted averaging
;  is done with AVGCLIP.PRO which does a robust average of the image
;  stack so that cosmic rays, stars, and other transient image defects are
;  eliminated.  This will do a very good sky flat if enough images are provided
;  of different field pointings.
;
;  When done, the resulting flat image is returned to the caller and the image
;  is saved to a FITS file with the specified output filename.
;
; CATEGORY:
;  CCD data processing
;
; CALLING SEQUENCE:
;  mkflat,root,outsuf,start,nframes,bias,dark,flat
;
; INPUTS:
;  root    - Root of the file name (you must include . in the name).
;  outsuf  - The suffix of the final output file.
;  start   - First frame number to read (integer or long).
;               Start can also be a vector of explicit frame numbers to load.
;               In this case, nframes need not be specified and in fact will
;               be ignored.
;  nframes - Number of frames to average.
;  bias    - Bias frame image name to subtract from each raw flat.
;  dark    - Dark frame image name to subtract from each raw flat.
;
; OPTIONAL INPUT PARAMETERS:
;
; KEYWORD INPUT PARAMETERS:
;
;   CROP     = region of original image to save, default=no cropping.
;                 [x1,x2,y1,y2]
;
;   DDIR    - Path to the raw data, default = ''
;
;   EXCLUDE - Optional vector of image numbers that should be excluded from
;                average.  Default is to include all frames.
;
;   EXPKEY = String - FITS keyword to read to get exposure time, default = EXPTIME
;
;   OVERSCAN = column overscan region to use for frame bias level,
;              default=no overscan subtraction.  jm06nov07nyu -
;              parameter redefined to be the following, as in
;              ARM_FLATCOMBINE: "coordinate array for overscan region
;              (if the input frames are images with NX x NY pixels,
;              then OVERSCAN should be [x1,y1,x2,y2] ([lower-left X
;              and Y coordinates, upper-right X and Y coordinates])"
;
;   SCALE - 4 element vector which, if provide, defines the region of the
;           array dimensions that are used to scale the mean
;           of the arrays before combining (.  If combined in this
;           manner, the arrays are combined weighted by the means.
;                 [x1,x2,y1,y2]
;           These coordinates apply to the pixel locations AFTER cropping.
;           The default is to use the center 50% of the image but not any
;           bigger than 200x200 subsection at the center.
; 
;   FULLFLAT - over-ride the SCALE optional input and use the full
;              image (best if used in combination with BADPIXMASK)
;
; OUTPUTS:
;  flat - Final robust averaged and normalized flat image.
;
; KEYWORD OUTPUT PARAMETERS:
; COMMON BLOCKS:
; SIDE EFFECTS:
; RESTRICTIONS:
; PROCEDURE:
;
; MODIFICATION HISTORY:
;  95/03/10 - Initial crude version written, Marc W. Buie, Lowell Observatory
;  95/06/04, MWB, added EXCLUDE keyword
;  95/06/13, MWB, added SCALE, OVERSCAN, CROP keywords
;  99/06/10, MWB, added EXPTIME keyword and added documentation.
;  2000/02/03, MWB, rewrite to add support for multigroup FITS files.
;  2000/02/28, MWB, added support for frame numbers > 999.
;  2001/02/23, MWB, added option to provide input file list.
;  2001/04/28, MWB, added DDIR keyword
;  2004/05/05, MWB, fixed SCALE keyword default action.  It was not working
;                     right at all and the region selected was much
;                     too large.
;  jm06nov07nyu - added BADPIXMASK optional input and FULLFLAT
;                 keyword; also modified OVERSCAN optional input to
;                 compute the simple median rather than to call
;                 COLBIAS 
;-
pro im_mkflat,root,outsuf,start,nframes,bias,dark,flat, $
   EXCLUDE=exclude,SCALE=in_scale,EXPKEY=expkey, $
   OVERSCAN=in_overscan,CROP=in_crop,DDIR=ddir,badpixmask=in_badpixmask,$
  fullflat=fullflat

   self='MKFLAT: '
   if badpar(root,   7,    0,caller=self+'(root) '   ) then return
   if badpar(outsuf, 7,    0,caller=self+'(outsuf) ' ) then return
   if badpar(start,  [2,3,7],[0,1],caller=self+'(start) ', $
                rank=start_rank,type=start_type) then return

   if badpar(exclude,[0,2,3],[0,1],caller=self+'(exclude) ', $
                                   default=-1) then return

   if badpar(in_overscan,[0,2,3],[1,2],caller=self+'(overscan) ', $
                                       rank=o_rank) then return
   if badpar(in_crop,    [0,2,3],[1,2],caller=self+'(crop) ', $
                                       rank=c_rank) then return
   if badpar(expkey,     [0,7],      0,caller=self+'(EXPKEY) ', $
                                       default='EXPTIME') then return

   if badpar(bias,   [0,7],0,caller=self+'(bias) ', $
                                       default='[[none]]') then return
   if badpar(dark,   [0,7],0,caller=self+'(dark) ', $
                                       default='[[none]]') then return
   if badpar(in_scale,[0,2,3],[0,1],caller=self+'(scale) ', $
                                       default=[-1,-1,-1,-1]) then return
   if badpar(ddir,[0,7],0,caller=self+'(DDIR) ',default='') then return

   if ddir ne '' then ddir=addslash(ddir)

   ; Check to see if it's a sequential list or random list.
   if start_rank eq 0 then begin
      frames=start+indgen(nframes)
      if badpar(nframes,[2,3],   0,caller=self+'(nframes) ') then return
   endif else if start_type ne 7 then begin
      frames=start
      nframes=n_elements(frames)
   endif

   if start_type eq 7 then begin
      fname=start
      nframes=n_elements(fname)
      if exclude[0] eq -1 then $
         bad=intarr(nframes) $
      else $
         bad=exclude
   endif else begin
      ; Setup the file name format string
      digits = fix(ceil(alog10(max(frames)+1))) > 3
      dig = strn(digits)
      fnfmt = '(i'+dig+'.'+dig+')'
      fname=root+string(frames,format=fnfmt)
      bad=intarr(nframes)
      ; Apply the exclusion criteria to the frames list.  Then, find out how
      ;   many frames are being collected and make sure there is something to do.
      for i=0,nframes-1 do begin
         z=where(frames[i] eq exclude,count)
         if count ne 0 then bad[i]=1
      endfor
   endelse

   zg=where(bad eq 0,countg)
   if countg eq 0 then $
      message,'Error ** you have excluded all frames, nothing to do.'

   ; See if we need to append the .fits tag to the file name.
   if exists(ddir+fname[0]+'.fits') then fname=fname+'.fits'

   ; Make the output file name
   outfile = root+outsuf

   ; Check header of image to see if it is a multi-extension image.
   hdr=headfits(ddir+fname[0])
   numext=sxpar(hdr,'NEXTEND')

   ; Setup the overscan/crop control values
   if numext eq 0 then begin
      extend=0
      numext=1
      if o_rank eq 0 then begin
         do_overscan=0
      endif else begin
         do_overscan=1
         overscan = in_overscan
      endelse
      if c_rank eq 0 then begin
         do_crop=0
      endif else begin
         do_crop=1
         crop = in_crop
      endelse
   endif else begin
      extend=1
      if o_rank eq 0 then begin
         do_overscan=0
      endif else if o_rank eq 1 then begin
         do_overscan=1
         overscan = rebin(in_overscan,n_elements(overscan),numext)
      endif else begin
         do_overscan=1
         overscan = in_overscan
      endelse
      if c_rank eq 0 then begin
         do_crop=0
      endif else if o_rank eq 1 then begin
         do_crop=1
         crop = rebin(in_crop,n_elements(crop),numext)
      endif else begin
         do_crop=1
         crop = in_crop
      endelse
   endelse

   ; If it's multi-extension, then the header we've just read is special
   ;   And must be written back out to start the output file.
   if extend then $
      writefits,outfile,0,hdr

   ; Main loop over extension, done just once on "normal" frames
   for ix=1,numext do begin
      ix0=ix-1

      print,'Load ',strn(nframes),' frames.'

      ; Read the bias frame for this image set or extension
      if bias ne '[[none]]' then begin
         if extend then $
            biasim=readfits(bias,exten_no=ix,/silent) $
         else $
            biasim=readfits(bias,/silent)
      endif

      ; Read the dark frame for this image set or extension
      if dark ne '[[none]]' then begin
         if extend then $
            darkim=readfits(dark,exten_no=ix,/silent) $
         else $
            darkim=readfits(dark,/silent)
      endif

      ; Loop over the frames
      j=0
      for i=0,nframes-1 do begin

         if not bad[i] then begin

            ; read in the new image, the dummy line before readfits is to dump
            ;   the storage for the array before its used again.
            image = 0
            if extend then $
               image = float(readfits(ddir+fname[i],hdr,exten_no=ix,/silent)) $
            else $
               image = float(readfits(ddir+fname[i],hdr,/silent))

            if do_overscan and do_crop then begin
; compute a simple average bias level, allowing for multiple columns; jm06nov07nyu
               biasval = median(float(image[overscan[0]:overscan[2],overscan[1]:overscan[3]]))
;              plothist, image[overscan[0]:overscan[2],overscan[1]:overscan[3]], bin=1.0, xsty=3, ysty=3
;              plot, djs_median(image[overscan[0]:overscan[2],overscan[1]:overscan[3]],1), ysty=3, xsty=3, ps=10
;              oplot, !x.crange, biasval*[1,1], line=2, thick=2
;              plot, djs_median(image[overscan[0]:overscan[2],overscan[1]:overscan[3]],2), ysty=3, xsty=3, ps=10
;              oplot, !x.crange, biasval*[1,1], line=2, thick=2
;              print, biasval & cc = get_kbrd(1)
;              biasval = total(image[overscan[0]:overscan[2],overscan[1]:overscan[3]]) / $
;                ((overscan[2]-overscan[0]+1.0)*(overscan[3]-overscan[1]+1.0))
               image = image[crop[0,ix0]:crop[1,ix0],crop[2,ix0]:crop[3,ix0]] - biasval
               if (n_elements(in_badpixmask) ne 0L) then badpixmask = in_badpixmask[crop[0,ix0]:crop[1,ix0],crop[2,ix0]:crop[3,ix0]]
;              image = colbias(image,overscan[0,ix0],overscan[1,ix0], $
;                       crop[0,ix0],crop[1,ix0],crop[2,ix0],crop[3,ix0], $
;                       biasval=biasval)
               if extend then $
                  print,fname[i],' overscan value is ', strtrim(string(biasval),2), $
                                 ' extension ',strn(ix) $
               else $
                  print,fname[i],' overscan value is ', strtrim(string(biasval),2)

               sxaddpar, hdr, 'OVERSCAN', '['+string(overscan[0]+1,format='(I0)')+':'+string(overscan[1]+1,format='(I0)')+','+$
                 string(overscan[2]+1,format='(I0)')+':'+string(overscan[3]+1,format='(I0)')+']', ' overscan region (mean='+$
                 strtrim(string(biasval),2)+')'
               sxaddpar, hdr, 'TRIM', '['+string(crop[0]+1,format='(I0)')+':'+string(crop[1]+1,format='(I0)')+','+$
                 string(crop[2]+1,format='(I0)')+':'+string(crop[3]+1,format='(I0)')+']', ' trim region'
               sxaddhist, "'Overscan region subtracted "+hogg_iso_date()+"'", hdr
               sxaddhist, "'Image trimmed "+hogg_iso_date()+"'", hdr

            endif else if do_overscan then begin
; compute a simple average bias level, allowing for multiple columns; jm06nov07nyu
               biasval = median(image[overscan[0]:overscan[2],overscan[1]:overscan[3]])
;              biasval = total(image[overscan[0]:overscan[2],overscan[1]:overscan[3]]) / $
;                ((overscan[2]-overscan[0]+1.0)*(overscan[3]-overscan[1]+1.0))
               image = image[crop[0,ix0]:crop[1,ix0],crop[2,ix0]:crop[3,ix0]] - biasval
;              image = colbias(image,overscan[0,ix0],overscan[1,ix0], $
;                              biasval=biasval)
               if (n_elements(in_badpixmask) ne 0L) then badpixmask = in_badpixmask
               if extend then $
                  print,fname[i],' overscan value is ', strtrim(string(biasval),2), $
                                 ' extension ',strn(ix) $
               else $
                  print,fname[i],' overscan value is ', strtrim(string(biasval),2)

               sxaddpar, hdr, 'OVERSCAN', '['+string(overscan[0]+1,format='(I0)')+':'+string(overscan[1]+1,format='(I0)')+','+$
                 string(overscan[2]+1,format='(I0)')+':'+string(overscan[3]+1,format='(I0)')+']', ' overscan region (mean='+$
                 strtrim(string(biasval),2)+')'
               sxaddhist, "'Overscan region subtracted "+hogg_iso_date()+"'", hdr

            endif else if do_crop then begin
               image = image[crop[0,ix0]:crop[1,ix0],crop[2,ix0]:crop[3,ix0]]
               if (n_elements(in_badpixmask) ne 0L) then badpixmask = in_badpixmask[crop[0,ix0]:crop[1,ix0],crop[2,ix0]:crop[3,ix0]]
               if extend then $
                  print,fname[i],' extension ',strn(ix) $
               else $
                  print,fname[i]

               sxaddpar, hdr, 'TRIM', '['+string(crop[0]+1,format='(I0)')+':'+string(crop[1]+1,format='(I0)')+','+$
                 string(crop[2]+1,format='(I0)')+':'+string(crop[3]+1,format='(I0)')+']', ' trim region'
               sxaddhist, "'Image trimmed "+hogg_iso_date()+"'", hdr

            endif else begin
               if extend then $
                  print,fname[i],' extension ',strn(ix) $
               else $
                  print,fname[i]

            endelse

            if i eq 0 then begin
               sz=size(image,/dim)
               cube=fltarr(sz[0],sz[1],countg,/nozero)
               if min(in_scale) lt 0 then begin
                  scale=[sz[0]/4,sz[0]/4,sz[1]/4,sz[1]/4]
                  scale=(scale < 100)*[-1,1,-1,1] + $
                           [sz[0]/2,sz[0]/2,sz[1]/2,sz[1]/2]
                  scale = [ scale[0] > 0, $
                            scale[1] < (sz[0]-1), $
                            scale[2] > 0, $
                            scale[3] < (sz[1]-1) ]
               endif else begin
                  scale=in_scale
               endelse
            endif

            if bias ne '[[none]]' then $
               image = temporary(image)-biasim

            if dark ne '[[none]]' then begin
               exptime=float(sxpar(hdr,expkey))
               image = temporary(image) - darkim*exptime
            endif

            cube[*,*,j]=image
            j=j+1

         endif

      endfor 

      if keyword_set(fullflat) then scale = [0L,sz[0]-1L,0L,sz[1]-1L]
      print,'Scaling region: x ',strn(scale[0]),':',strn(scale[1]), $
        '  y ',strn(scale[2]),':',strn(scale[3])

      ; Average the cube down.
      flat=0
      if (n_elements(badpixmask) ne 0L) then $
        im_avgclip,cube,flat,scale=scale,/norm,bad=rebin(badpixmask,[sz,nframes]) else $
          im_avgclip,cube,flat,scale=scale,/norm

;     flat = djs_avsigclip(cube,sigrej=sigrej,maxiter=maxiter,inmask=rebin(badpixmask,[sz,nframes]))

      ; Protect against zeros in the flat.
      z=where(flat le 0.0,count)
      if count ne 0 then begin
         flat[z] = -10000.0
         sflat = median(flat,3)
         flat[z] = sflat[z]
         z=where(flat le 0.0,count)
         if count ne 0 then begin
            flat[z] = max(flat)
         endif
      endif

; set bad pixel values to 0.0 in the output flat

      if (n_elements(badpixmask) ne 0L) then begin
         bad = where(badpixmask eq 1B,nbad)
         if (nbad ne 0L) then flat[bad] = 0.0
      endif
      
      outhdr = hdr
;     mkhdr, outhdr, flat, extend=extend ; jm06nov07nyu

      sxaddpar,outhdr,'NAXIS1',sz[0]
      sxaddpar,outhdr,'NAXIS2',sz[1]
      sxaddpar,outhdr,'BITPIX',-32
      sxdelpar,outhdr,'BSCALE'
      sxdelpar,outhdr,'BZERO'

      if extend then begin
         print,outfile,',  writing extension ',strn(ix)
         writefits,outfile,float(flat),outhdr,/append
      endif else begin
         print,'Saving final flat frame to ',outfile
         writefits,outfile,float(flat),outhdr
      endelse

   endfor

end
