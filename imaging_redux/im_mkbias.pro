;+
; NAME: 
;  mkbias
; PURPOSE: 
;  Collect and combine CCD bias frames into a superbias frame
; DESCRIPTION:
;
;  The files are assumed to be named as 'root'NNN where 'root' is the
;  root of the file name that you supply and NNN is a three digit number.
;  If your file name has an imbedded '.' then add it to the root.
;
;  The specified range of files are all read in from FITS files.  Then each
;  image has the overscan mean subtracted (if desired), cropped (as indicated).
;  These images are then averaged.  The averaging is done with AVGCLIP.PRO
;  which does a robust average of the image stack so that cosmic rays or
;  other transient image defects are eliminated.
;
;  When done, the resulting bias image is returned to the caller and the image
;  is saved to a FITS file with the specified output filename.
;
; CATEGORY:
;  CCD data processing
; CALLING SEQUENCE:
;  mkbias,root,outsuf,start,nframes,bias
; INPUTS:
;  root    - Root of the file name (you must include . in the name).
;  outsuf  - The suffix of the final output file.
;  start   - First frame number to read (integer or long).
;               Start can also be a vector of explicit frame numbers to load.
;               In this case, nframes need not be specified and in fact will
;               be ignored.
;               Additionally, start can also be a string array containing the
;               file names of all files to be read.  In this case, set nframes
;               to 0 or someother innocuous integer.  Exclude is treated
;               differently.  In this case, exclude is a vector of the same
;               length as start and is 0 if the file is to be used, 1 if not.
;  nframes - Number of frames to average.
; OPTIONAL INPUT PARAMETERS:
; KEYWORD INPUT PARAMETERS:
;
;   CROP     = region of original image to save, default=no cropping.
;                 [x1,x2,y1,y2]
;
;   DDIR     = Directory to look for raw data in.  Default = ''
;
;   EXCLUDE - Optional vector of image numbers that should be excluded from
;                average.  Default is to include all frames.
;
;   OVERSCAN = column overscan region to use for frame bias level,
;              default=no overscan subtraction.  jm06nov07nyu -
;              parameter redefined to be the following, as in
;              ARM_FLATCOMBINE: "coordinate array for overscan region
;              (if the input frames are images with NX x NY pixels,
;              then OVERSCAN should be [x1,y1,x2,y2] ([lower-left X
;              and Y coordinates, upper-right X and Y coordinates])"
;
; OUTPUTS:
;  bias - Final robust averaged/cropped bias image.
; KEYWORD OUTPUT PARAMETERS:
; COMMON BLOCKS:
; SIDE EFFECTS:
; RESTRICTIONS:
; PROCEDURE:
; MODIFICATION HISTORY:
;  95/03/09 - Initial crude version written, Marc W. Buie, Lowell Observatory
;  95/06/13, MWB, added OVERSCAN and CROP keywords
;  95/11/22, MWB, add EXCLUDE keyword
;  2000/02/02, MWB, rewrite to add support for multigroup FITS files.
;  2000/02/28, MWB, added support for frame numbers > 999.
;  2001/02/23, MWB, added option to provide input file list.
;  2001/04/28, MWB, added DDIR keyword.
;  jm06nov07nyu - modified OVERSCAN optional input
;-
pro im_mkbias,root,outsuf,start,nframes,bias, $
       DDIR=ddir,OVERSCAN=in_overscan,CROP=in_crop,EXCLUDE=exclude

   if badpar(root,   7,       0,caller='MKBIAS: (root) '   ) then return
   if badpar(outsuf, 7,       0,caller='MKBIAS: (outsuf) ' ) then return
   if badpar(start,  [2,3,7], [0,1],caller='MKBIAS: (start) ', $
                rank=start_rank,type=start_type  ) then return

   if badpar(exclude,[0,2,3], [0,1],caller='MKBIAS: (exclude) ',default=-1) then return

   if badpar(in_overscan,[0,2,3],[1,2],caller='MKBIAS: (overscan) ',rank=o_rank) then return
   if badpar(in_crop,    [0,2,3],[1,2],caller='MKBIAS: (crop) ',rank=c_rank) then return

   if badpar(ddir,[0,7],0,caller='MKBIAS: (DDIR) ',default='') then return
   if ddir ne '' then ddir=addslash(ddir)

   ; Check to see if it's a sequential list or random list.
   if start_rank eq 0 then begin
      frames=start+indgen(nframes)
      if badpar(nframes,[2,3],   0,caller='MKBIAS: (nframes) ') then return
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

      if not exists(ddir+fname[zg[0]]) then begin
         print,'MKBIAS:  ',ddir+fname[zg[0]],' not found, unable to continue.'
         return
      endif

      ; Loop over the frames
      j=0
      for i=0,nframes-1 do begin

         if not bad[i] then begin

            ; read in the new image, the dummy line before readfits is to dump
            ;   the storage for the array before it's used again.
            image = 0
            if i eq 0 then begin
               if extend then $
                  image = float(readfits(ddir+fname[i],hdr,exten_no=ix,/silent)) $
               else $
                  image = float(readfits(ddir+fname[i],hdr,/silent))
            endif else begin
               if extend then $
                  image = float(readfits(ddir+fname[i],exten_no=ix,/silent)) $
               else $
                  image = float(readfits(ddir+fname[i],/silent))
            endelse

            if do_overscan and do_crop then begin
; compute a simple average bias level, allowing for multiple columns; jm06nov07nyu
               biasval = median(image[overscan[0]:overscan[2],overscan[1]:overscan[3]])
;              plothist, image[overscan[0]:overscan[2],overscan[1]:overscan[3]], bin=1.0, xsty=3, ysty=3
;              plot, djs_median(image[overscan[0]:overscan[2],overscan[1]:overscan[3]],1), ysty=3, xsty=3, ps=10
;              oplot, !x.crange, biasval*[1,1], line=2, thick=2
;              plot, djs_median(image[overscan[0]:overscan[2],overscan[1]:overscan[3]],2), ysty=3, xsty=3, ps=10
;              oplot, !x.crange, biasval*[1,1], line=2, thick=2
;              print, biasval & cc = get_kbrd(1)
;              biasval = total(image[overscan[0]:overscan[2],overscan[1]:overscan[3]]) / $
;                ((overscan[2]-overscan[0]+1.0)*(overscan[3]-overscan[1]+1.0))
               image = image[crop[0,ix0]:crop[1,ix0],crop[2,ix0]:crop[3,ix0]] - biasval
;              image = colbias(image,overscan[0,ix0],overscan[1,ix0], $
;                       crop[0,ix0],crop[1,ix0],crop[2,ix0],crop[3,ix0],biasval=biasval)
               if extend then $
                  print,fname[i],' overscan value is ',strtrim(string(biasval),2),' extension ',strn(ix) $
               else $
                  print,fname[i],' overscan value is ',strtrim(string(biasval),2)

               if (i eq 0L) then begin
                  sxaddpar, hdr, 'OVERSCAN', '['+string(overscan[0]+1,format='(I0)')+':'+string(overscan[1]+1,format='(I0)')+','+$
                    string(overscan[2]+1,format='(I0)')+':'+string(overscan[3]+1,format='(I0)')+']', ' overscan region (mean='+$
                    strtrim(string(biasval),2)+')'
                  sxaddpar, hdr, 'TRIM', '['+string(crop[0]+1,format='(I0)')+':'+string(crop[1]+1,format='(I0)')+','+$
                    string(crop[2]+1,format='(I0)')+':'+string(crop[3]+1,format='(I0)')+']', ' trim region'
                  sxaddhist, "'Overscan region subtracted "+hogg_iso_date()+"'", hdr
                  sxaddhist, "'Image trimmed "+hogg_iso_date()+"'", hdr
               endif

            endif else if do_overscan then begin
; compute a simple average bias level, allowing for multiple columns; jm06nov07nyu
               biasval = median(image[overscan[0]:overscan[2],overscan[1]:overscan[3]])
;              biasval = total(image[overscan[0]:overscan[2],overscan[1]:overscan[3]]) / $
;                ((overscan[2]-overscan[0]+1.0)*(overscan[3]-overscan[1]+1.0))
               image = image[crop[0,ix0]:crop[1,ix0],crop[2,ix0]:crop[3,ix0]] - biasval
;              image = colbias(image,overscan[0,ix0],overscan[1,ix0],biasval=biasval)
               if extend then $
                  print,fname[i],' overscan value is ',strtrim(string(biasval),2),' extension ',strn(ix) $
               else $
                  print,fname[i],' overscan value is ',strtrim(string(biasval),2)

               if (i eq 0L) then begin
                  sxaddpar, hdr, 'OVERSCAN', '['+string(overscan[0]+1,format='(I0)')+':'+string(overscan[1]+1,format='(I0)')+','+$
                    string(overscan[2]+1,format='(I0)')+':'+string(overscan[3]+1,format='(I0)')+']', ' overscan region (mean='+$
                    strtrim(string(biasval),2)+')'
                  sxaddhist, "'Overscan region subtracted "+hogg_iso_date()+"'", hdr
               endif

            endif else if do_crop then begin
               image = image[crop[0,ix0]:crop[1,ix0],crop[2,ix0]:crop[3,ix0]]
               if extend then $
                  print,fname[i],' extension ',strn(ix) $
               else $
                  print,fname[i]

               if (i eq 0L) then begin
                  sxaddpar, hdr, 'TRIM', '['+string(crop[0]+1,format='(I0)')+':'+string(crop[1]+1,format='(I0)')+','+$
                    string(crop[2]+1,format='(I0)')+':'+string(crop[3]+1,format='(I0)')+']', ' trim region'
                  sxaddhist, "'Image trimmed "+hogg_iso_date()+"'", hdr
               endif

            endif else begin
               if extend then $
                  print,fname[i],' extension ',strn(ix) $
               else $
                  print,fname[i]

            endelse

            if j eq 0 then begin
               sz=size(image,/dim)
               cube=fltarr(sz[0],sz[1],countg,/nozero)
            endif

            cube[*,*,j]=image
            j=j+1
         endif

      endfor 

      ; Average the cube down.
      bias=0
      im_avgclip,cube,bias

      outhdr = hdr
;     mkhdr, outhdr, bias, extend=extend ; jm06nov07nyu
      
      sxaddpar,outhdr,'NAXIS1',sz[0]
      sxaddpar,outhdr,'NAXIS2',sz[1]
      sxaddpar,outhdr,'BITPIX',-32
      sxdelpar,outhdr,'BSCALE'
      sxdelpar,outhdr,'BZERO'

      if extend then begin
         print,outfile,',  writing extension ',strn(ix)
         writefits,outfile,float(bias),outhdr,/append
      endif else begin
         print,'Saving final bias frame to ',outfile
         writefits,outfile,float(bias),outhdr
      endelse

   endfor 

end
