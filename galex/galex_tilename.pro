function galex_tilename, tile
; jm10apr22ucsd - internal support routine for BUILD_GALEX_CATALOG
; that rebuilds the GALEX tilename from the tile structure (see
; BUILD_GALEX_TILELIST)

    nobj = n_elements(tile)
    if (nobj gt 1L) then begin
       name = strarr(nobj)
       for ii = 0L, nobj-1L do name[ii] = galex_tilename(tile[ii])
       return, name
    endif
    
    galex_dir = getenv('GALEX_DIR')
    if (strtrim(galex_dir,2) eq '') then message, $
      '$GALEX_DIR environment variable required!'
    gr = strtrim(tile.gr,2)
    vsn = string(tile.vsn,format='(I2.2)')
; the AIS tiles are a special case because there are many possible
; "sub-visits" for each tile
    tilename = strtrim(tile.tile_name,2)
    if strmatch(tilename,'*AIS_*') then tilename = $
      strmid(tilename,0,strpos(tilename,'_sg'))
    name = file_search(galex_dir+'/'+gr+'/pipe/'+vsn+'-vsn/'+$
      string(tile.tile_num,format='(I5.5)')+'-'+tilename+$
      '/d/01-main/*/*/*'+strtrim(tile.tile_name,2)+'*mcat*fits*')
    if (n_elements(name) gt 1) or (strtrim(name,2) eq '') then $
      message, 'Problem finding tile!'
return, name
end

