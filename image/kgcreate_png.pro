pro kgcreate_png, fits_file, png_file, int_range=int_range, thumbnail=thumbnail, $
  scaled_image=scaled_image, ss_scale=ss_scale, log_scale=log_scale, $
  spec_thumbnail=spec_thumbnail, file_ext=file_ext

; if both SS_SCALE and LOG_SCALE are set then SS_SCALE takes
; precedence 

    if keyword_set(ss_scale) and keyword_set(log_scale) then log_scale = 0

    if (keyword_set(file_ext)) then begin

       fits_read, fits_file, image, header, exten_no=file_ext

    endif else begin

       fits_read, fits_file, image, header

    endelse

;tot_time = fxpar(header,'EXPTIME')
;if (tot_time NE 0.0) then image = image/tot_time

    if (keyword_set(thumbnail)) then begin

       size_image = size(image)
       if (size_image[1] GT size_image[2]) then begin
          ysize = 150
          xsize = ysize*size_image[1]/size_image[2]
       endif else begin
          xsize = 150
          ysize = xsize*size_image[2]/size_image[1]
       endelse
       disp_image = congrid(image,xsize,ysize)

    endif else if (keyword_set(spec_thumbnail)) then begin

       size_image = size(image)
       if (size_image[1] GT size_image[2]) then begin
          ysize = 150
          xsize = ysize*size_image[1]/size_image[2]
       endif else begin
          xsize = 150
          ysize = xsize*size_image[2]/size_image[1]
       endelse
       disp_image = congrid(image,xsize,ysize)

    endif else begin

       disp_image = image

    endelse

stop

    if (keyword_set(ss_scale)) then begin ; sqrt(sqrt) scaling

;      image_size = size(disp_image)
;      ave = median(disp_image)
;      std_dev = sqrt(total((disp_image - ave)^2)/(image_size[1]*image_size[2]))
;      if ((ave - 1.0*std_dev) LT 0.0) then begin
;          print,ave,std_dev
;          disp_image = disp_image + 1.0*std_dev
;      endif

       indxs = where(disp_image GT 0,n_indxs)
       if (n_indxs GT 0) then disp_image[indxs] = (disp_image[indxs])^0.25
       min_val = min(disp_image[indxs],max=max_val)

;      image_size = size(disp_image)
;      ave = median(disp_image)
;      std_dev = sqrt(total((disp_image - ave)^2)/(image_size[1]*image_size[2]))
;      int_range = [min_val,ave+1.*std_dev]

       int_range = [min_val,max_val]

    endif

    
    if (n_elements(int_range) eq 0L) then begin

       stats = im_stats(image,sigrej=3.0)
       int_range = stats.median + stats.sigma_rej*[-3.0,3.0]
       
;      image_size = size(image)
;      ave = median(image)
;      std_dev = sqrt(total((image - ave)^2)/(image_size[1]*image_size[2]))
;      int_range = ave + [-3*std_dev,3.*std_dev]

    endif else if (int_range[0] EQ int_range[1]) then begin

       min_val = min(disp_image,max=max_val)
       int_range = [min_val,max_val]
       
    endif

    disp_image = imgscl(disp_image,min=int_range[0],max=int_range[1],log=log_scale)
;   disp_image = bytscl(disp_image,min=int_range[0],max=int_range[1])

    loadct, 3, /silent
    tvlct, r, g, b, /get
    write_png, png_file, disp_image, r, g, b

    scaled_image = disp_image

return
end
