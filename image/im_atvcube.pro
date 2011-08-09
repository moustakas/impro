;+
; NAME:
;	IM_ATVCUBE
;
; PURPOSE:
;	This routine acts as a wrapper to Aaron Barth's image display
;	program, ATV.  An array of filenames or an image data cube
;	(and optionally an image header cube) can be displayed one at a
;	time with ATV via a simple widget interface.  The routine can
;	be thought of as a generalized BLINK function.
;
; INPUTS:
;	imlist - either an array of filenames or a three-dimensional
;                image data cube whose third dimension is the image
;                index
;
; OPTIONAL INPUTS:
;	headcube - if imlist is passed as an image cube, then the
;                  headers of the images can be passed *as a pointer
;                  array* (see EXAMPLE, below) via this keyword 
;
; OUTPUTS:
;	None.
;
; COMMON BLOCKS:
;	atv_state - see the ATV documentation
;
; EXAMPLE:
;	In this first example ATVCUBE will be passed an image and
;	header data cube.  Imagine that the variable imlist contains a
;	list of 10 FITS files of size (xsize,ysize).  The image and
;	header cube might be created in the following way:
; 
;		IDL> imcube = fltarr(xsize,ysize,10)
;		IDL> headcube = ptrarr(10)
;		IDL> for i = 0L, 9L do begin 
;		IDL>    imcube[*,*,i] = readfits(imlist[i],h)
;		IDL>    headcube[i] = ptr_new(h)
;		IDL> endfor
;
;	Then we would call ATVCUBE with
;
;		IDL> atvcube, imcube, headcube=headcube
;
;	Note that headcube must be a pointer array because in general
;	the header lengths of different files are not the same.
;	Alternatively, we can create a list of FITS files we want to
;	examine one at a time and call ATVCUBE directly:
;
;		IDL> imlist = ['im1.fits','im2.fits','im3.fits']
;		IDL> atvcube, imlist
;
; COMMENTS:
;	The stretch and alignment of subsequent images are carried
;	over from the previous image displayed.
;
;	For the most recent version of ATV go to Aaron's webpage at 
;       http://cfa-www.harvard.edu/~abarth/atv/atv.html.
;
; WARNINGS:
;	If a file list is passed then the first image is read in to
;	create the image cube array.  This procedure assumes that all
;	the subsequent FITS files are the same x- and y-size.
;
; PROCEDURES USED:
;	ATV
;
; MODIFICATION HISTORY:
;	John Moustakas, 2001 February 13, University of Arizona
;
; Copyright (C) 2001, John Moustakas
; 
; This program is free software; you can redistribute it and/or modify 
; it under the terms of the GNU General Public License as published by 
; the Free Software Foundation; either version 2 of the License, or
; (at your option) any later version. 
; 
; This program is distributed in the hope that it will be useful, but 
; WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
; General Public License for more details. 
;
;-

pro atvcube_shutdown, base_id
; destroy all the widgets and free the header pointers

    widget_control, base_id, get_uvalue = info
    for k = 0L, n_elements(info.jarray)-1L do ptr_free, info.headcube[k]
    if xregistered('atvcube') then widget_control, base_id, /destroy
    atv_shutdown
    
return
end

pro atvcube_event, event
; atvcube event handler

    widget_control, event.id, get_uvalue = event_name
    widget_control, event.top, get_uvalue = info, /no_copy

; if atv is dead, then kill the wrapper
    
    valid_atv = widget_info(info.state.base_id,/valid_id)
    if not valid_atv then event_name = 'done'

    case event_name of
       'next': begin
          info.jarray = shift([info.jarray],-1)
          info.j = info.jarray[0]

          atv, info.imcube[*,*,info.j], head=*(info.headcube[info.j]), /align, /stretch

          if (size(info.imlist[0],/type) eq 7L) then $
            widget_control, info.field, set_value = info.imlist[info.j] else $
            widget_control, info.field, set_value = info.j
       end
       'previous': begin
          info.jarray = shift([info.jarray],1)
          info.j = info.jarray[0]

          atv, info.imcube[*,*,info.j], head=*(info.headcube[info.j]), /align, /stretch

          if (size(info.imlist[0],/type) eq 7L) then $
            widget_control, info.field, set_value = info.imlist[info.j] else $
            widget_control, info.field, set_value = info.j
       end
       'help': help = dialog_message(['Press the "Next" and "Previous" buttons to scroll',$
                                      'through the image cube.  Press "Done" to quit ATV.'],$
                                     dialog_parent=event.top,/info,title='ATV Wrapper Help')
       'done': begin
          widget_control, event.top, set_uvalue = info, /no_copy
          atvcube_shutdown, event.top
       end
    endcase    

; if the wrapper widget hasn't been killed then continue

    valid_atvcube = widget_info(event.top,/valid_id)
    if valid_atvcube then widget_control, event.top, set_uvalue = info, /no_copy
    
return
end

pro im_atvcube, imlist, headcube=headcube

    on_error, 2
    resolve_routine, ['atv'], /compile_full_file

    common atv_state, state
    
    imsize = size([imlist])
    nimage = imsize[3]

    if keyword_set(headcube) then if (nimage ne (size(headcube))[3]) then $
      message, 'Image and header cube dimensions do not agree!'
    
    jarray = lindgen(nimage)   ; image cube index array
    j = jarray[0]              ; display the zeroth image

    imtype = size(imlist,/type)
    if imtype eq 7L then begin  ; string file list

       fxread, imlist[0], temp, htemp, 1, 2, 1, 2 ; read a corner to extract image dimensions

       imcube = fltarr(fxpar(htemp,'NAXIS1'),fxpar(htemp,'NAXIS2'),nimage)
       headcube = ptrarr(nimage)

       for k = 0L, nimage-1L do begin ; read in the fits files
          imcube[*,*,k] = readfits(imlist[k],h,/silent)
          headcube[k] = ptr_new(h)
       endfor       

       imid = imlist[j]            ; image ID name
       xsize = max(strlen(imlist))
       
       delvarx, temp, htemp    ; blow away temporary variables

    endif else begin

       imcube = imlist         ; imlist is a fits cube
       imid = j
       xsize = 3

    endelse

; check the headcube keyword

    if not keyword_set(headcube) then headcube = replicate(ptr_new(0L),nimage) else begin ; no header information
       for i = 0L, nimage-1L do if not ptr_valid(headcube[i]) then message, 'Headcube is not a valid pointer!'
    endelse

    atv, imcube[*,*,j], head = *(headcube[j])

; create and realize the wrapper widget

    base = widget_base(title='ATV Wrapper',/row,/align_right,kill_notify='atvcube_shutdown')
    field = cw_field(base,/row,title='Image:',value=imid,xsize=xsize,/noedit)
    button = widget_button(base,value='Next',uvalue='next',xsize=2,units=2)
    button = widget_button(base,value='Previous',uvalue='previous',xsize=2,units=2)
    button = widget_button(base,value='Help',uvalue='help',xsize=2,units=2)
    button = widget_button(base,value='Done',uvalue='done',xsize=2,units=2)

    info = {base_id: base, field: field, imlist: imlist, jarray: jarray, j: j, $
            imcube: imcube, headcube: headcube, state: state}
    
    device, get_screen_size = scrsize
    atvsize = widget_info(state.base_id,/geometry)

    widget_control, base, /realize, event_pro='atvcube_event', set_uvalue=info, /no_copy, $
      xoffset = (atvsize.xsize+atvsize.xoffset)<scrsize[0], yoffset = atvsize.yoffset
    xmanager, 'atvcube', base, event_handler='atvcube_event', /no_block, group_leader=base, $
      cleanup = 'atvcube_shutdown'

return
end





