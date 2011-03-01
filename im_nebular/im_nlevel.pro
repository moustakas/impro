;+
; NAME:
;     IM_NLEVEL
;
; PURPOSE:
;     Calculates the relative populations and emissivities for the lines
;     of a desired ion at a given density and temperature.  
;
; EXPLANATION:
;     The emission line strengths are calculated from an N-level
;     atom at a given density and temperature.  The atomic data
;     for each ion is stored in the procedure 'GETMATRIX.PRO'.
;
; CALLING SEQUENCE:
;      NLEVEL, IONSTR, DENS, TEMP, WAVE, EMIS [, HBETA]
;
; INPUT:
;      IONSTR: A string containing the name of the ion to be calculated.
;           In general the designation is given by [ABBREV]_[ROMAN]; i.e.,
;           doubly-ionized oxygen is 'O_III', neutral silicon is 'SI_I', &c.
;           See the procedure 'GETMATRIX' for specifics.
;      DENS: The electron density for the emission-line region, in 1/cm^3.
;      TEMP: The electron temperature for the emission-line region, in K.
;
; OUTPUT:
;
;      WAVE: The NxN matrix containing the wavelengths of the relevant
;            transitions.  The transition from level U to level L is
;            in the array element WAVE[U,L].
;      EMIS: The array of line emissivities per ion per electron relative to
;            H-beta at the corresponding wavelengths.
;
; OPTIONAL INPUT KEYWORD:
;
;      HBETA: The H-beta emissivity at the given density and temperature.
;             The emissivity is given in units of erg/cm^3/s/[N(p) N(e)]
;
; EXAMPLE:
;
;      The matrix of wavelengths and line emissivities for the O++ ion
;      at electron density 100/cm^3 and electron temperature 10,000 K
;      is given by
;     
;      NLEVEL, 'O_III', 1d2, 1d4, WAVE, EMIS
;
;
; NOTES:
;
;      Note that it is the emissivities per ion per electron that is
;      given.  Thus, for the emissivity EMIS[U,L], the line strength
;      in a given spectrum relative to H-beta 4861 is given by
;      EMIS[U,L] x N(Ion)/N(H+).
;
;      Note also that a few ions require a bit more work -- namely,
;      [Ar III], [Ne III], and [S III].  These ions are problematic
;      when compared to the CLOUDY outputs. The N-level atomic data
;      has been tailored to maximize agreement, requiring the use
;      of branching ratios to recast those results into those expected
;      from the output.  That part of the code will only execute
;      if N < 5 for the ions in question.
;
;
; REVISION HISTORY:
;       1.0	B.D.Moore	Rice Univ.	September 2003
;-
;      J. Moustakas, 2007-Nov-27, NYU: now a function; pack everything
;         into a structure; default 10,000 K and 100 cm^-3 temperature
;         and density assumed
;
;  Get the atomic parameters for the desired ion.
;  N.B.  The energy levels 'eij' are in eV.
;

function im_nlevel, ionstr, dens=dens1, temp=temp1, matrix=matrix

    if (n_elements(ionstr) eq 0L) then begin
       doc_library, 'im_nlevel'
       return, -1L
    endif

    if (n_elements(temp1) eq 0L) then temp = 1D4 else temp = temp1 ; [K]
    if (n_elements(dens1) eq 0L) then dens = 1D2 else dens = dens1 ; [cm^-3]

    ntemp = n_elements(temp)
    ndens = n_elements(dens)

; call this routine recursively if TEMP and/or DENS are arrays; in
; this case LEVEL (and, optionally, MATRIX) are [NTEMP,NDENS] arrays 

    if (ntemp gt 1L) or (ndens gt 1L) then begin
       for idens = 0L, ndens-1L do begin
          for itemp = 0L, ntemp-1L do begin
             level1 = im_nlevel(ionstr,dens=dens[idens],$
               temp=temp[itemp],matrix=matrix1)
             if (n_elements(level) eq 0L) then begin
                level = level1
                if arg_present(matrix) then matrix = matrix1
             endif else begin
                level = [level,level1]
                if arg_present(matrix) then matrix = [matrix,matrix1]
             endelse
          endfor
       endfor
       if arg_present(matrix) then matrix = reform(matrix,ntemp,ndens)       
       return, reform(level,ntemp,ndens)
    endif else begin
       temp = temp[0] ; this prevents crashing on some of the matrix math
       dens = dens[0]
    endelse
    
    matrix = im_getmatrix(ionstr,temp)
    qij = matrix.qij
    eij = matrix.eij
    aij = matrix.aij
    nelem = fix(sqrt(n_elements(qij))) - 1
;
;  Make the vector and matrix to be solved by LU_SOL.
;
    bb = reform(dens*qij[0,1:nelem])
    aa = -dens*qij[1:nelem,1:nelem] - aij[1:nelem,1:nelem]
    for i = 1, nelem do aa[i-1,i-1] = total(dens*qij[i,*] + aij[i,*])
;
;     Solve for x in A#x = b.
    ludc, aa, index
    xx = lusol(aa, index, bb)
;
;     And now calculate the spectrum.
;
    yy = [1,xx]
    yy = yy/total(yy)
;
;
    wave = (eij ne 0)*12398.5/(eij > 1e-30)
    emis = 0.*eij
    hbeta = 1.387d-25*(1d4/temp)^(0.983)*exp(-976.3/temp)
    norm = dens*hbeta
    for i = 1, nelem do begin
       for j = 0, i do begin
          emis[i,j] = yy[i] * aij[i,j] * eij[i,j] * 1.602d-12 / norm
       endfor
    endfor

    if (n_elements(emis) lt 10) then begin
       case strupcase(ionstr) of
;
          'NE_III' : BEGIN
             ff = fltarr(5,5)
             wave = ff
             wave[3,0] = 3869.
             wave[3,1] = 3968.
             wave[4,3] = 3343.
             ff[3,0] = emis[1,0]*3.319/(1+3.319)
             ff[3,1] = emis[1,0]/(1+3.319)
             ff[4,3] = emis[2,1]
             emis = ff
          END
;
          'S_III' : BEGIN 
             ff = fltarr(5,5)
             wave = ff
             wave[3,2] = 9532.
             wave[3,1] = 9069.
             wave[4,3] = 6312.
             ff[3,2] = emis[1,0]*2.486/(1+2.486)
             ff[3,1] = emis[1,0]/(1+2.486)
             ff[4,3] = emis[2,1]
             emis = ff
          END
;
          'AR_III' : BEGIN 
             ff = fltarr(5,5)
             wave = ff
             wave[3,0] = 7135.
             wave[3,1] = 7751.
             wave[4,3] = 5192.
             ff[3,0] = emis[1,0]*4.140/(1+4.140)
             ff[3,1] = emis[1,0]/(1+4.140)
             ff[4,3] = emis[2,1]
             emis = ff
          END
;
          ELSE   : BEGIN 
          END
;
       endcase
    endif

; pack everything into a structure and return

    level = {ion: ionstr, density: dens, temp: temp, linewave: wave, $
      emissivity: emis, hbeta: hbeta}
    
return, level
end
