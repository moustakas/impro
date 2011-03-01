pro match_sources, catname, refcat, datapath=datapath, rcrit=rcrit, $
         master_catalog=master_catalog, header=header
;+
; NAME:
;	MATCH_SOURCES
;
; PURPOSE:
;
;
; CALLING SEQUENCE:
;
;
; INPUTS:
;	catname - names of the ASCII catalogs to match (assumed
;                 format) 
;	refcat  - reference index number of the catalog to match *to*
;                 (zero-indexed)
;
; OPTIONAL INPUTS:
;	rcrit   - critical source matching radius
;	master_catalog - name of the output catalog (default:
;                        master_catalog.txt)
;	header - two-line header to provide for the text output of the
;                master catalog (first line should identify the
;                catalog, and the second should list the filters) 
;	
; KEYWORD PARAMETERS:
;
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
;
; EXAMPLE:
;
;
; PROCEDURES USED:
;	REMOVE, SRCOR, READFAST
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2001 August 31, U of A
;-

    if n_params() ne 2L then begin
       print, 'Syntax - match_sources, catname, refcat, [datapath=], '+ $
         '[rcrit=], [master_catalog=], [header=]'
       return
    endif

    if size(catname,/type) ne 7L then message, 'Catalog name is not type string!'
    if not keyword_set(datapath) then datapath = './'

    ncat = n_elements(catname)
    if ncat eq 1L then message, 'Matching requires two or more catalogs!'

    refcat = long(refcat)
    if (refcat ge ncat) or (refcat lt 0L) then message, 'Reference catalog number is out of range!'
    
    if not keyword_set(rcrit) then rcrit = 5.0 ; pixels
    
; read the catalogs

    datarray = ptrarr(ncat)
    
    for i = 0L, ncat-1L do begin

       readfast, datapath+catname[i], data, nlines=nsources, skip=1
       datarray[i] = ptr_new(data)
       
    endfor

; match the reference catalog to every other catalog by position
    
    catarray = lindgen(ncat) ; catalog array
    remove, refcat, catarray ; remove the reference catalog from the catalog array

    for j = 0L, ncat-2L do begin

       xyref = (*(datarray[refcat]))[1:2,*]      ; reference xy positions
       xycat = (*(datarray[catarray[j]]))[1:2,*] ; catalog xy positions

       srcor, xyref[0,*], xyref[1,*], xycat[0,*], xycat[1,*], rcrit, match1, match2, /option

; trim the lists

       refmatch = (*(datarray[refcat]))[*,match1]      ; reference catalog
       catmatch = (*(datarray[catarray[j]]))[*,match2] ; matched catalog

       ptr_free, datarray[refcat]      ; manage memory
       ptr_free, datarray[catarray[j]]

       datarray[refcat] = ptr_new(refmatch) 
       datarray[catarray[j]] = ptr_new(catmatch) 

       if j gt 0L then for k = j-1L, 0L, -1L do begin ; recursively trim previous catalogs

          backcat = (*(datarray[catarray[k]]))[*,match1] 
          ptr_free, datarray[catarray[k]]
          datarray[catarray[k]] = ptr_new(backcat) 

       endfor 

    endfor 

    nsources = n_elements(match1)
    print, format='("There are ", I0," sources in the final catalog.")', nsources

; sanity check
    
;   xycheck = fltarr(2*ncat,nsources)
;   for i = 0L, ncat-1L do xycheck[2*i:2*(i+1)-1L,*] = (*datarray[i])[1:2,*]
;   print, xycheck
    
; write the final catalog: x- and y-positions in the reference catalog
; plus flux and flux errors in all the catalogs

    finalcat = fltarr(2+2*ncat,nsources)
    finalcat[0:3,*] = (*(datarray[refcat]))[1:4,*]
    for i = 0L, ncat-2L do finalcat[4+2*i:4+2*(i+1)-1,*] = (*datarray[catarray[i]])[3:4,*]

    if not keyword_set(master_catalog) then master_catalog = 'master_catalog.txt'
    if not keyword_set(header) then header = ['# Master simulated images catalog.','#']

stop
    
    openw, lun, datapath+master_catalog, /get_lun
    printf, lun, header[0]
    printf, lun, header[1]
    for k = 0L, nsources-1L do printf, lun, k+1, finalcat[*,k], format='(I7,2F12.5,'+$
      string(2*ncat)+'G17.5)'
    free_lun, lun

; clean up memory

    if long(total(ptr_valid(datarray)))/ncat then ptr_free, datarray

return
end
