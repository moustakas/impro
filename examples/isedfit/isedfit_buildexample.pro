PRO isedfit_buildexample, _extra=extra
  field = 'cosmos'
  nobjects = 10
  cat = im_mrdfits(getfilename(field))
  ii = where(cat.isgood AND $
             cat.ismass AND $
             cat.zquality EQ 4, nii) ;; slightly more strict
  cat =  cat[ii[0:nobjects-1]]
  filters = getfilters(field)
  match, strtrim(cat[0].filters, 2), filters, nn, mm, count=nfilters
  
  example = {ra:0d, $
             dec:0d, $
             z:0.0, $
             filterlist:strarr(nfilters), $
             maggies:fltarr(nfilters), $
             ivarmaggies:fltarr(nfilters)}
  example = replicate(example, n_elements(cat))
  example.ra =  cat.ra
  example.dec = cat.dec
  example.z = cat.z
  example.filterlist = filters[mm]
  example.maggies = cat.maggies[nn]
  example.ivarmaggies = cat.ivarmaggies[nn]

  im_mwrfits, example, 'cosmos_example.fits', _extra=extra
END
