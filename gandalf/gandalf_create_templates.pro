FUNCTION gandalf_CREATE_TEMPLATES, EMISSION_SETUP=emission_setup, PARS=pars, NPIX=npix         ,$
                           LSTEP_GAL=lstep_gal, INT_DISP_PIX=int_disp_pix, LOG10=log10 ,$
                           FOR_ERRORS=for_errors

; Take the emission-setup structure and the input pars parameter array
; to make emission-line single or multi-Gaussian templates.
;
; The input pars array should contiant only the Gaussian paramenter for the lines
; that we are actually fitting, such as Hb or [OIII]5007, containing
; only the V_gas and S_gas parameters, like [V_Hb, S_Hb, V_OIII5007, S_OIII5007, ...]
;
; On the other hand if we wish to fit also the amplitudes with MPFIT,
; then the pars array should be [A_Hb, V_Hb, S_Hb, A_OIII5007, ...]

; This happens when we are evaluating the errors in the Gaussian
; parameters around the best solution.

i_lines   = where(emission_setup.kind eq 'l')
nlines   = n_elements(i_lines)
if not keyword_set(for_errors) then n_pars = nlines*2 else n_pars = nlines*3
if (n_pars ne n_elements(pars)) then message,'Hey, this is not the right emission-line parameter array'
; array that will contain the emission-line templates, including multiplets
gaus = dblarr(npix,nlines)

; a) First create the emission-line templates corresponding to each
; single line (like Hb) or to each main line of a multiplet (like
; [OIII]5007)
for i = 0,nlines-1 do begin
    ; Create the emission-line templates. 
    if not keyword_set(for_errors) then begin
        ; Use the amplitude from the emission-line setup. If that is set to unity, 
        ; than the NNLS weight assigned to each template will actually correspond 
        ; to the emission-line amplitude.
       ampl_i = emission_setup[i_lines[i]].a
        gaus[*,i]=create_gaussn(X=findgen(npix),PARS=[ampl_i,pars[2*i:2*i+1]],INT_DISP_PIX=int_disp_pix)
    endif else begin
        ; Make Gaussian templates amplitudes specified by the input pars array
        gaus[*,i]=create_gaussn(X=findgen(npix),PARS=[pars[3*i:3*i+2]],INT_DISP_PIX=int_disp_pix)
    endelse
endfor

; b) Then find all the satellite lines belonging to multiplets
; (like [OIII]4959), and add them the main emission-line template 
; we just created in a)
i_slines = where(strmid(emission_setup.kind,0,1) eq 'd')
n_slines = n_elements(i_slines)

if i_slines[0] ne -1 then begin
    ; loop over the satellite lines
    for i = 0,n_slines-1 do begin
        ; Current index in the emission-line setup structure for 
        ; the current satellite line (e.g. 1 for [OIII]4959)
        j = i_slines[i]
        ; Find the reference line index, as given in the "kind" tag of the
        ; emission setup structure, which points to the main line in the
        ; present multiplet (e.g. 2 for [OIII]5007 if kind was d2 for [OIII]4959)
        k_mline = fix(strmid(emission_setup[j].kind,1)) 
        ; which correspond to the following position in the emission setup
        ; structure that was passed to this function (e.g. still 2, if
        ; the indices of the lines in the emission setup start at 0 and
        ; increase by 1)
        j_mline = where(emission_setup.i eq k_mline)
        ; Get the wavelengths of both satellite and main lines
        ; and compute the offset (in pix) that we need to apply in order
        ; to correctly place the satellite line.     
        l_sline = emission_setup[j].lambda
        l_mline = emission_setup[j_mline].lambda 
        offset  = alog(l_mline/l_sline)/lstep_gal   
        ; to deal with log10-lambda rebinned data, instead of ln-lambda
        if keyword_set(log10) then offset  = alog10(l_mline/l_sline)/lstep_gal
        ; Get the index in the array of the to fit, corresponding
        ; to the main line of the present multiplet, so that we can add the 
        ; satellite emission-line template in the right place in the
        ; gaussian templates array
        i_mline = where(emission_setup[i_lines].i eq k_mline)
        ; Finally, create the satellite template, and add it to that of
        ; the corresponding main line of this multiplet.
        ; Use the amplitude given in the emission setup structure as the
        ; relative strength of the satellite line w.r.t to the amplitude
        ; of the main lines. 
        ; In turn these are as specified either by the emission-lines
        ; setup or by the input pars array.
        if (keyword_set(for_errors) eq 0) then begin
            a_sline = emission_setup[j].a*emission_setup[j_mline].a
            gaus_sline = create_gaussn(X=findgen(npix), $
                                       PARS=[a_sline,pars[i_mline*2]-offset,pars[i_mline*2+1]], $
                                       INT_DISP_PIX=int_disp_pix)
        endif else begin
            a_sline = emission_setup[j].a*pars[i_mline*3]
            gaus_sline = create_gaussn(X=findgen(npix), $
                                       PARS=[a_sline,pars[i_mline*3+1]-offset,pars[i_mline*3+2]], $
                                       INT_DISP_PIX=int_disp_pix)
        endelse
        gaus[*,i_mline] = gaus[*,i_mline] + gaus_sline 
    endfor
endif

return, gaus
END

;------------------------------------------------------------------------------------
