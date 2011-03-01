;+
; NAME:
;   read_galex_tileinfo()
; PURPOSE:
;   Read basic information on all the tiles given a GALEX data
;   release. 
; CALLING SEQUENCE:
;   info = read_galex_tileinfo()
; KEYWORDS:
;   qaplot - optionally make a QAplot showing the sky distribution of
;     all the tiles
; INPUTS: 
;   gr - galex release number (default 'gr6')
; OUTPUTS: 
;   info - basic information structure
; COMMENTS:
;   See http://galex.stsci.edu/GR4/doc/GR6.txt
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Jul 26, UCSD
;-

function read_galex_tileinfo, gr=gr, qaplot=qaplot

    if (n_elements(gr) eq 0) then gr = 'gr6'

    galexpath = getenv('CATALOGS_DIR')+'/galex/'
    infofile = galexpath+gr+'.dat'
    if (file_test(infofile) eq 0) then begin
       splog, 'Tile information file '+infofile+' not found'
       return, -1
    endif
    
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

    if keyword_set(qaplot) then begin
       psfile = galexpath+'galex_alltiles_'+gr+'.ps'
       im_plotconfig, 0, pos, psfile=psfile, height=5.5, xmargin=[1.2,0.3]
       djs_plot, [0], [0], /nodata, position=pos, xsty=1, $
         ysty=1, xtitle='\alpha_{J2000} (degree)', $
         ytitle='\delta_{J2000} (degree)', $
         xrange=[0,360], yrange=[-90,90]
       surv = strtrim(info.survey,2)
       dis = where(surv eq 'DIS')
       mis = where(surv eq 'MIS')
       gii = where(surv eq 'GII')
       cai = where(surv eq 'CAI')
       if (dis[0] ne -1L) then djs_oplot, info[dis].tile_ra, $
         info[dis].tile_dec, psym=6, sym=0.5, color='red'
       if (mis[0] ne -1L) then djs_oplot, info[mis].tile_ra, $
         info[mis].tile_dec, psym=6, sym=0.2, color='blue'
       if (gii[0] ne -1L) then djs_oplot, info[gii].tile_ra, $
         info[gii].tile_dec, psym=6, sym=0.5, color='orange'
       if (cai[0] ne -1L) then djs_oplot, info[cai].tile_ra, $
         info[cai].tile_dec, psym=6, sym=0.9, color='dark green'
       legend, ['DIS','MIS','GII','CAI'], /left, /top, box=0, $
         charsize=1.3, psym=6, margin=0, $
         color=djs_icolor(['red','blue','orange','dark green'])
       im_plotconfig, /psclose, /gzip, psfile=psfile
    endif
    
return, info
end
