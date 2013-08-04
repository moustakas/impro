;+
; NAME:
;     IM_IONSPEC
; PURPOSE:
;     Calculates the emission-line spectrum of an ion from an N-level
;     atom, and returns an array containing the wavelengths and
;     line emissivities.
;
; EXPLANATION:
;     This procedure calls the IDL NEBULAR procedure 'NLEVEL', and
;     extracts the non-zero emissivities.  These are written to a
;     2xN array of the wavelengths and emissivities.  Unless otherwise
;     specified, the procedure writes to the screen the transition
;     in question, the wavelength, the emissivity relative to H-beta,
;     and the emissivity in erg/cm^3/s/[N(Ion) N(e)].
;
;
; CALLING SEQUENCE:
;     IM_IONSPEC, IONSTR, DENS, TEMP, RESULT, HBETA, ANNOT [, /SILENT ]
;
; INPUT:
;      IONSTR: A string containing the name of the ion to be calculated.
;           In general the designation is given by [ABBREV]_[ROMAN]; i.e.,
;           doubly-ionized oxygen is 'O_III', neutral silicon is 'SI_I', &c.
;           See the procedure 'GETMATRIX' for specifics.
;      DENS: The electron density for the emission-line region, in 1/cm^3.
;      TEMP: The electron temperature for the emission-line region, in K.
;
;
; OUTPUT:
;
;      RESULT: A 2xN matrix containing the wavelengths and emmisivities
;           of levels with emissivities > 0 *relative* to H-beta
;      HBETA: The emissivity of H-beta at the input electron temperature,
;                in units of erg/cm^3/s/[ N(H+) N(e) ].
;      ANNOT: String array containing the annotation of the transitions
;                in terms of the upper level (UL) and lower level (LL),
;                written in the form 'UL --> LL'.  The annotation uses
;                LL = 1 for the ground state.
;
; OPTIONAL INPUT KEYWORD:
;      SILENT: Setting this keyword suppresses the output to the screen.
;
;
; EXAMPLE:
;
;      The wavelengths and line emissivities for the O++ ion    
;      at electron density 100/cm^3 and electron temperature 10,000 K       
;      is given by
;      
;      IONSPEC, 'O_III', 1d2, 1d4, RES, HBETA, ANNOT, /SILENT
;
;
; NOTES:
;     This procedure is a more user-friendly version of 'NLEVEL',
;     which it invokes to calculate the spectrum.
;
;      J. Moustakas, 2007-Nov-27, NYU: default 10,000 K and 100 cm^-3
;         temperature and density assumed; LEVEL is optionally
;         returned from IM_NLEVEL, which also contains the H-BETA
;         emissivity at the specific Te and ne
;
; REVISION HISTORY:
;       1.0	B.D.Moore	Rice Univ.	September 2003
;
;       jm07dec04nyu - hacked significantly
;-
;

pro im_ionspec, ionstr, result, dens=dens, temp=temp, level=level, $
  matrix=matrix, full_annot=full_annot, silent=SILENT, $
  sortbystrength=sortbystrength

    if (n_elements(ionstr) eq 0L) then begin
       doc_library, 'im_ionspec'
       return
    endif

    if (n_elements(temp) eq 0L) then temp = 1D4 ; [K]
    if (n_elements(dens) eq 0L) then dens = 1D2 ; [cm^-3]

    level = im_nlevel(ionstr,dens=dens,temp=temp,matrix=matrix)

    wave = level.linewave
    emis = level.emissivity
    hbeta = level.hbeta
    aij = matrix.aij
    
    nelem = n_elements(emis[0,*])
    full_annot = strarr(nelem,nelem)
    full_index = strarr(nelem,nelem)
    for i = 0, nelem - 1 do begin
       for j = 0, nelem - 1 do begin
          full_annot[i,j] = strcompress(string(i+1) + ' -->' + string(j+1),/remove)
          full_index[i,j] = '['+strcompress(string(i) + ',' + string(j),/remove)+']'
       endfor
    endfor

    q = where((wave ne 0) and (emis ne 0),nq)
    result = replicate({annot: '', index: '', vacwave: 0.0d, $
      airwave: 0.0d, aij: 0.0, relative_emissivity: 0.0, absolute_emissivity: 0.0},nq)

    vacwave = wave[q]
    airwave = vacwave
    vactoair, airwave

    result.annot = full_annot[q]
    result.index = full_index[q]
    result.vacwave = vacwave
    result.airwave = airwave
    result.relative_emissivity = emis[q]
    result.absolute_emissivity = emis[q]*hbeta
    result.aij = aij[q]

    result = result[sort(result.vacwave)]
    if keyword_set(sortbystrength) then result = result[reverse(sort(result.relative_emissivity))]
    if (not keyword_set(silent)) then struct_print, result

return
end
