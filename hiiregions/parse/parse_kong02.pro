pro parse_kong02, myredcorr=myredcorr
; jm07nov06nyu - parse Kong et al. 2002; abundances are from Shi et
;                al. 2005; this paper is a nightmare; first of all, I
;                can't reproduce Table 3 in the Kong paper from Table
;                2; if I do my own reddening correction, some of the
;                line-fluxes are totally different, especially for
;                4363, which is worrisome [and even the E(B-V) values
;                are totally different!]; second, even when I use the
;                published "intrinsic" 4363 values, I only get the
;                published Te and (O/H) values for about half the
;                sample: OK it turns out that this is because Shi used
;                the Pilyugin (2001) empirical R23-Te relation, which
;                is crap!  note that Shi also cut on EW(4363)>2 A and
;                S/N[4363]>5!

; THROUGH EXTENSIVE EXPERIMENTATION BY COMPARING MY VALUES WITH SHI I
; HAVE DETERMINED THAT THE 4363 VALUES IN THE KONG PAPER (AND HENCE
; THE Te AND O/H VALUES IN THE SHI PAPER) ARE TOTAL CRAP!    
    
    myredcorr = 1L ; NOTE!
    
    outpath = hiiregions_path()
    table1 = im_read_vizier_tsv(outpath+'02kong/2002_kong_table1.dat')
    table2 = im_read_vizier_tsv(outpath+'02kong/2002_kong_table2.dat')
    table = struct_addtags(table1,struct_trimtags(table2,except=['RECNO','NAME','TYPE']))
    
    table3 = im_read_vizier_tsv(outpath+'02kong/2002_kong_table3.dat')
    table3 = im_struct_trimtags(table3,select=tag_names(table3),newtags=repstr(tag_names(table3),'FHB','EBV'))

    rem = where(strmatch(table.name,'*VIIZw631*',/fold),comp=keep) ; this object only appears in Tables 1&2
    table = table[keep]
    table = struct_addtags(table,struct_trimtags(table3,except=['RECNO','NAME','TYPE']))

    table1_shi = rsex(outpath+'02kong/2005_shi_table1.dat')

; match against the Shi abundances; AGN have been removed
    
    match, strlowcase(strtrim(table.name,2)), strlowcase(strtrim(table1_shi.shi_name,2)), aa, bb 
    table = table[aa] & table1_shi = table1_shi[bb]
;   niceprint, table.name, table1_shi.shi_name

    table = struct_addtags(table,struct_trimtags(table1_shi))

    nobject = n_elements(table)
    out = init_hii_region_structure(nobject)

    out.hii_galaxy  = strtrim(table.name,2)
    out.hii_region  = 'H1'

; ZW0855 (also called ZW0855+06) is not recognized by NED; replace
; with UGC04703; also ZW2220=UGC12011 and ZW2335=UGCA441=MRK328 

    zw = where(strmatch(out.hii_galaxy,'*zw0855*',/fold))
    out[zw].hii_region = out[zw].hii_galaxy
    out[zw].hii_galaxy = 'UGC04703'

    zw = where(strmatch(out.hii_galaxy,'*zw2220*',/fold))
    out[zw].hii_region = out[zw].hii_galaxy
    out[zw].hii_galaxy = 'UGC12011'

    zw = where(strmatch(out.hii_galaxy,'*zw2335*',/fold))
    out[zw].hii_region = out[zw].hii_galaxy
    out[zw].hii_galaxy = 'UGCA441'

    good = where((table.ewhb gt 0.0),ngood)
    if (ngood ne 0L) then begin
       out[good].ewhb = table[good].ewhb
    endif

; TEST
;   good = where((table.ewo3c gt 0.0),ngood)
;   if (ngood ne 0L) then begin
;      out[good].ewhb = table[good].ewo3c 
;   endif

; THESE VALUES ARE CRAP -- DO NOT USE    
;   good = where((table.shi_toiii gt 0.0),ngood)
;   if (ngood ne 0L) then begin
;      out[good].t4363 = table[good].shi_toiii/1D4
;      out[good].t4363_err = 0.1*table[good].shi_toiii/1D4 ; assume 10%
;   endif
;
;   good = where((table.shi_log12oh gt 0.0),ngood)
;   if (ngood ne 0L) then begin
;      out[good].log12oh = table[good].shi_log12oh
;      out[good].log12oh_err = 0.05 ; assume 0.5 dex
;   endif

    if keyword_set(myredcorr) then begin
       
; reddening-correct the data myself; I *assume* that the lines have
; already been corrected for Galactic extinction!!

       line = {$
         linename:         ['OII_3727','H_BETA','OIII_4363','OIII_5007','H_ALPHA','NII_6584','SII_6716','SII_6731'], $
         oii_3727_wave:    3727.42, $
         h_beta_wave:      4861.33, $
         oiii_4363_wave:   4363.21, $
         oiii_5007_wave:   5006.84, $
         h_alpha_wave:     6562.80, $
         nii_6584_wave:    6583.46, $
         sii_6716_wave:    6716.14, $
         sii_6731_wave:    6730.81, $
         oii_3727:      [0.0,-2.0], $
         h_beta:        [0.0,-2.0], $
         oiii_4363:     [0.0,-2.0], $
         oiii_5007:     [0.0,-2.0], $
         h_alpha:       [0.0,-2.0], $
         nii_6584:      [0.0,-2.0], $
         sii_6716:      [0.0,-2.0], $
         sii_6731:      [0.0,-2.0]}
       line = replicate(line,nobject)

       k_oii_3727  = k_lambda(line[0].oii_3727_wave,/odonnell)
       k_h_beta    = k_lambda(line[0].h_beta_wave,/odonnell)
       k_oiii_4363 = k_lambda(line[0].oiii_4363_wave,/odonnell)
       k_oiii_5007 = k_lambda(line[0].oiii_5007_wave,/odonnell)
       k_h_alpha   = k_lambda(line[0].h_alpha_wave,/odonnell)
       k_nii_6584  = k_lambda(line[0].nii_6584_wave,/odonnell)
       k_sii_6716  = k_lambda(line[0].sii_6716_wave,/odonnell)
       k_sii_6731  = k_lambda(line[0].sii_6731_wave,/odonnell)
       
       good = where((table.foii gt 0.0) and (table.e_foii gt 0.0),ngood)
       if (ngood ne 0L) then begin
          line[good].oii_3727[0]  = table[good].foii*1D-15
          line[good].oii_3727[1]  = table[good].e_foii*1D-15
       endif

       good = where((table.fo3b gt 0.0) and (table.e_fo3b gt 0.0),ngood)
       if (ngood ne 0L) then begin
          line[good].oiii_5007[0]  = table[good].fo3b*1D-15
          line[good].oiii_5007[1]  = table[good].e_fo3b*1D-15
       endif

       good = where((table.fo3c gt 0.0) and (table.e_fo3c gt 0.0),ngood)
       if (ngood ne 0L) then begin
          line[good].oiii_4363[0]  = table[good].fo3c*1D-15
          line[good].oiii_4363[1]  = table[good].e_fo3c*1D-15
       endif

       good = where((table.fha gt 0.0) and (table.e_fha gt 0.0),ngood)
       if (ngood ne 0L) then begin
          line[good].h_alpha[0]  = table[good].fha*1D-15
          line[good].h_alpha[1]  = table[good].e_fha*1D-15
       endif

       good = where((table.fhb gt 0.0) and (table.e_fhb gt 0.0),ngood)
       if (ngood ne 0L) then begin
          line[good].h_beta[0]  = table[good].fhb*1D-15
          line[good].h_beta[1]  = table[good].e_fhb*1D-15
       endif

       good = where((table.fnii gt 0.0) and (table.e_fnii gt 0.0),ngood)
       if (ngood ne 0L) then begin
          line[good].nii_6584[0]  = table[good].fnii*1D-15
          line[good].nii_6584[1]  = table[good].e_fnii*1D-15
       endif

       good = where((table.fsiia gt 0.0) and (table.e_fsiia gt 0.0),ngood)
       if (ngood ne 0L) then begin
          line[good].sii_6716[0]  = table[good].fsiia*1D-15
          line[good].sii_6716[1]  = table[good].e_fsiia*1D-15
       endif

       good = where((table.fsiib gt 0.0) and (table.e_fsiib gt 0.0),ngood)
       if (ngood ne 0L) then begin
          line[good].sii_6731[0]  = table[good].fsiib*1D-15
          line[good].sii_6731[1]  = table[good].e_fsiib*1D-15
       endif

; correct for reddening; two objects (IIZW70 and IZW123) have
; unphysical Balmer decrements (2.665+/-0.0741 and 2.739+/-0.0725,
; respectively; allow these to enter the sample by increasing CUTSIG
; to 3-sigma

       line_nodust = iunred_linedust(line,snrcut=1.0,/silent,cutsig=3.0)
;      niceprint, table.name, line_nodust.ebv_hahb, line_nodust.ebv_hahb_err, table.ebv

; construct the final line ratios, relative to H-beta

       good = where((line_nodust.oii_3727[1] gt 0.0) and (line_nodust.h_beta[1] gt 0.0),ngood)
       if (ngood ne 0L) then begin
          out[good].oii_3727 = line_nodust[good].oii_3727[0]/line_nodust[good].h_beta[0]
          out[good].oii_3727_err = im_compute_error(line_nodust[good].oii_3727[0],line_nodust[good].oii_3727[1],$
            line_nodust[good].h_beta[0],line_nodust[good].h_beta[1],/quotient)
       endif

       good = where((line_nodust.oiii_4363[1] gt 0.0) and (line_nodust.h_beta[1] gt 0.0),ngood)
       if (ngood ne 0L) then begin
          out[good].oiii_4363 = line_nodust[good].oiii_4363[0]/line_nodust[good].h_beta[0]
          out[good].oiii_4363_err = im_compute_error(line_nodust[good].oiii_4363[0],line_nodust[good].oiii_4363[1],$
            line_nodust[good].h_beta[0],line_nodust[good].h_beta[1],/quotient)
       endif

       good = where((line_nodust.oiii_5007[1] gt 0.0) and (line_nodust.h_beta[1] gt 0.0),ngood)
       if (ngood ne 0L) then begin
          out[good].oiii_5007 = line_nodust[good].oiii_5007[0]/line_nodust[good].h_beta[0]
          out[good].oiii_5007_err = im_compute_error(line_nodust[good].oiii_5007[0],line_nodust[good].oiii_5007[1],$
            line_nodust[good].h_beta[0],line_nodust[good].h_beta[1],/quotient)
       endif
       
       good = where((line_nodust.h_alpha[1] gt 0.0) and (line_nodust.h_beta[1] gt 0.0),ngood)
       if (ngood ne 0L) then begin
          out[good].ha = line_nodust[good].h_alpha[0]/line_nodust[good].h_beta[0]
          out[good].ha_err = im_compute_error(line_nodust[good].h_alpha[0],line_nodust[good].h_alpha[1],$
            line_nodust[good].h_beta[0],line_nodust[good].h_beta[1],/quotient)
       endif
       
       good = where((line_nodust.nii_6584[1] gt 0.0) and (line_nodust.h_beta[1] gt 0.0),ngood)
       if (ngood ne 0L) then begin
          out[good].nii_6584 = line_nodust[good].nii_6584[0]/line_nodust[good].h_beta[0]
          out[good].nii_6584_err = im_compute_error(line_nodust[good].nii_6584[0],line_nodust[good].nii_6584[1],$
            line_nodust[good].h_beta[0],line_nodust[good].h_beta[1],/quotient)
       endif
       
       good = where((line_nodust.sii_6716[1] gt 0.0) and (line_nodust.h_beta[1] gt 0.0),ngood)
       if (ngood ne 0L) then begin
          out[good].sii_6716 = line_nodust[good].sii_6716[0]/line_nodust[good].h_beta[0]
          out[good].sii_6716_err = im_compute_error(line_nodust[good].sii_6716[0],line_nodust[good].sii_6716[1],$
            line_nodust[good].h_beta[0],line_nodust[good].h_beta[1],/quotient)
       endif
       
       good = where((line_nodust.sii_6731[1] gt 0.0) and (line_nodust.h_beta[1] gt 0.0),ngood)
       if (ngood ne 0L) then begin
          out[good].sii_6731 = line_nodust[good].sii_6731[0]/line_nodust[good].h_beta[0]
          out[good].sii_6731_err = im_compute_error(line_nodust[good].sii_6731[0],line_nodust[good].sii_6731[1],$
            line_nodust[good].h_beta[0],line_nodust[good].h_beta[1],/quotient)
       endif

; these two objects have a null "observed" detection of 4363, but a
; non-zero "intrinsic" detection!!!       
       
       w = where((out.oiii_4363 lt -900.0) and (table.io3c gt -900.0),nw)
       if (nw ne 0L) then begin
          niceprint, out[w].hii_galaxy, out[w].oiii_4363, table[w].fo3c, table[w].e_fo3c, $
            table[w].fhb, table[w].fha, table[w].io3c/100.0, line_nodust[w].ebv_hahb
       endif
       
    endif else begin
       
; use the published "intrinsic" (reddening-corrected) line-fluxes, but
; no uncertainties were published, so scale the uncertainties in the
; extinction-corrected fluxes by the S/N of the observed fluxes
; (this doesn't work for all galaxies!!)

       good = where((table.ioii gt 0.0) and (table.foii gt 0.0) and (table.e_foii gt 0.0),ngood)
       if (ngood ne 0L) then begin
          out[good].oii_3727     = table[good].ioii/100.0
          out[good].oii_3727_err = (table[good].ioii/100.0)/(table[good].foii/table[good].e_foii)
       endif

       good = where((table.io3c gt 0.0) and (table.fo3c gt 0.0) and (table.e_fo3c gt 0.0),ngood)
       if (ngood ne 0L) then begin
          out[good].oiii_4363     = table[good].io3c/100.0
;         out[good].oiii_4363_err = (table[good].io3c/100.0)/(table[good].ewo3c/table[good].e_ewo3c)
;         niceprint, (table[good].io3c/100.0), (table[good].ewo3c/table[good].e_ewo3c)
          out[good].oiii_4363_err = (table[good].io3c/100.0)/(table[good].fo3c/table[good].e_fo3c)
       endif
       
       good = where((table.io3b gt 0.0) and (table.fo3b gt 0.0) and (table.e_fo3b gt 0.0),ngood)
       if (ngood ne 0L) then begin
          out[good].oiii_5007     = table[good].io3b/100.0
          out[good].oiii_5007_err = (table[good].io3b/100.0)/(table[good].fo3b/table[good].e_fo3b)
       endif

       good = where((table.iha gt 0.0) and (table.fha gt 0.0) and (table.e_fha gt 0.0),ngood)
       if (ngood ne 0L) then begin
          out[good].ha     = table[good].iha/100.0
          out[good].ha_err = (table[good].iha/100.0)/(table[good].fha/table[good].e_fha)
       endif

       good = where((table.inii gt 0.0) and (table.fnii gt 0.0) and (table.e_fnii gt 0.0),ngood)
       if (ngood ne 0L) then begin
          out[good].nii_6584     = table[good].inii/100.0
          out[good].nii_6584_err = (table[good].inii/100.0)/(table[good].fnii/table[good].e_fnii)
       endif

       good = where((table.isiia gt 0.0) and (table.fsiia gt 0.0) and (table.e_fsiia gt 0.0),ngood)
       if (ngood ne 0L) then begin
          out[good].sii_6716     = table[good].isiia/100.0
          out[good].sii_6716_err = (table[good].isiia/100.0)/(table[good].fsiia/table[good].e_fsiia)
       endif

       good = where((table.isiib gt 0.0) and (table.fsiib gt 0.0) and (table.e_fsiib gt 0.0),ngood)
       if (ngood ne 0L) then begin
          out[good].sii_6731     = table[good].isiib/100.0
          out[good].sii_6731_err = (table[good].isiib/100.0)/(table[good].fsiib/table[good].e_fsiib)
       endif

    endelse

; TOTAL CRAP -- DO NOT USE
    out.oiii_4363     = -999.0
    out.oiii_4363_err = -999.0
    
; write out

    filename = '2002_kong.sex'
    reference = 'Kong & Cheng 2002, A&A, 389, 845'

    openw, lun, outpath+filename, /get_lun
    printf, lun, '## Written by PARSE_KONG02 '+im_today()+'.'
    printf, lun, '## '+reference
    printf, lun, '## '
    printf, lun, '## See PARSE_KONG02 for details.'
    printf, lun, '## '
    tags = tag_names(out)
    for i = 0L, n_tags(out)-1L do printf, lun, '# '+string(i+1L,format='(I2)')+' '+strlowcase(tags[i])
    struct_print, out, lun=lun, /no_head
    free_lun, lun
    
return
end
