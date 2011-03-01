function galex_read_tileinfo, gr=gr
; jm10apr26ucsd - read the GALEX tile info file on all the available
; tiles
; jm10jul22ucsd - generalized GR release

    if (n_elements(gr) eq 0) then gr = 'gr6'

    infofile = getenv('CATALOGS_DIR')+'/galex/'+gr+'.dat'
    splog, 'Reading '+infofile
    readcol, infofile, survey, tilename, tilenum, ra, dec, $
      exptime, nuv, comment, comment='#', delimiter='|', $
      format='A,A,L,D,D,F,A,A', /silent
    ntile = n_elements(survey)
    info = {survey: '', tile_name: '', tile_num: 0L, $
      tile_ra: 0.0D, tile_dec: 0.0D, exptime: 0.0, $
      nuv_comment: '', comment: ''}
    info = replicate(info,ntile)
    info.survey = survey
    info.tile_name = tilename
    info.tile_num = tilenum
    info.tile_ra = ra
    info.tile_dec = dec
    info.exptime = exptime
    info.nuv_comment = nuv
    info.comment = comment

;   ii = galex_read_info()
;   jj = mrdfits('galex_tiles_gr4_gr5.fits.gz',1)    
;   im_plotconfig, 0, pos, psfile='galex_gr45.ps'
;   djs_plot, ii.tile_ra, ii.tile_dec, ps=6, sym=0.1, xsty=1, ysty=1, xtitle='RA', $
;     ytitle='Dec', xrange=[0,360], yrange=[-90,90], position=pos
;   djs_oplot, jj.tile_ra, jj.tile_dec, ps=6, sym=0.1, color='red'
;   im_plotconfig, /psclose, /gzip, psfile='galex_gr45.ps'
    
return, info
end
