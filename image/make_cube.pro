function make_cube, incube, sigmamap=sigmamap, sky=sky, mask=mask
; jm0feb14uofa
; jm03dec8uofa - added SKY keyword

    nim = (size(incube))[3]
    imsize = size(incube[0].image,/dimension)
    xsize = imsize[0]
    ysize = imsize[1]

    if keyword_set(mask) then $
      outcube = make_array(xsize,ysize,nim,/integer) else $
      outcube = make_array(xsize,ysize,nim,/float)

    if keyword_set(sigmamap) then for k = 0L, nim-1L do outcube[*,*,k] = incube[k].sigmamap
    if keyword_set(sky) then for k = 0L, nim-1L do outcube[*,*,k] = incube[k].sky 
    if keyword_set(mask) then for k = 0L, nim-1L do outcube[*,*,k] = incube[k].mask 

    if (not keyword_set(sigmamap)) and (not keyword_set(sky)) and (not keyword_set(mask)) then $
      for k = 0L, nim-1L do outcube[*,*,k] = incube[k].image

return, outcube
end
