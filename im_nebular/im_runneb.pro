;+
; NAME:
;     IM_RUNNEB
; PURPOSE:
;     Derives the ionic abundances from the input emission-line spectra.
;
;
; EXPLANATION:
;     This is the main procedure of the IDL NEBULAR package.  It analyzes
;     the emission lines in the input spectra, using standard diagnostic
;     line ratios to infer the physical connditions in the nebulae.
;     The nebulae are then modeled with a multizone model, where each
;     ion is presumed to exist solely in one of the zones.  The nebular
;     diagnostics for the ions in each zone are used to derive average
;     densities and temperatures for each zone.  The conditions in each
;     zone, combined with the strength of the emission lines in the spectra,
;     are used to infer the abundances of the ions in the nebulae.
;
; CALLING SEQUENCE:
;     IM_RUNNEB, LIST, SPEC, ABUND, DENS, TEMP, DSOL, TSOL, DRATS, 
;               TRATS, [, DDIR=DDIR ]
; INPUT:
;     LIST: A vector of strings designating the files containing the spectra.
;           If the list is empty, it is assumed that the spectra already
;           are contained in the variable SPEC.
;
; OUTPUT:
;     SPEC: An array of structures containing the input spectra.
;     ABUND: An array of structures containing the inferred ionic abundances.
;     DENS:  An array of structures containing the electron densities of
;              each ion for each spectrum.
;     TEMP: An array of structures containing the electron temperatures of
;              each ion for each spectrum.
;     DSOL: An array of structures containing the detailed solutions for
;              the density diagnostics of each ion for each spectrum.
;     TSOL: An array of structures containing the detailed solutions for
;              the temperature diagnostics of each ion for each spectrum.
;     DRATS: The density diagnostic ratios of each ion for each spectrum.
;     TRATS: The density diagnostic ratios of each ion for each spectrum.
;
; OPTIONAL INPUT KEYWORD:
;
;     SPEC: If the array of spectra structures already exist, then one
;           can use them instead of the files by passing the null character
;           in for the variable LIST.
;
;     DDIR: If the database of diagnostics is in a subfolder other than
;           'data', the path to the other folder is passed as a string
;           with "ddir = DDIR", where DDIR is the appropriate string.
;
; EXAMPLE:
;
;     With the database in the default directory, and the list of files
;     to be analyzed in LIST, the IDL NEBULAR package is run with the call
;
;     IM_RUNNEB, LIST, SPEC, ABUND, DENS, TEMP, DSOL, TSOL, DRATS, TRATS
;
;     If one wanted to run the IDL NEBULAR package again on the same spectra,
;     but using an alternate database of diagnostic ratios in the folder
;     'alt_database', the second call would be 
;
;     IM_RUNNEB, '', SPEC, ABUND, DENS, TEMP, DSOL, TSOL, DRATS, TRATS,
;               ddir = 'alt_database'
;
;     assuming that SPEC still held the spectra to be analyzed.
;
; NOTES:
;
;     The ability to run on different databases of diagnostic ratios may
;     be useful to some, hence the option to change directories at runtime.
;     HOWEVER, it should be noted that this implies a different set of
;     atomic data for one or more ions.  Thus, for the analysis to make
;     sense one must also compile a corresponding version of the IDL
;     NEBULAR procedure 'GETMATRIX' before running from the alternate
;     database.  CAVEAT EXECUTOR.
;
; REVISION HISTORY:
;       1.0	B.D.Moore	Rice Univ.	September 2003
;       J. Moustakas, 2005-Jun-01, U of A - renamed IM_RUNNEB; added
;          error checking, generalized I/O; LIST and SPEC were
;-

pro im_runneb, list, abund, dens, temp, dsol, tsol, drats, trats, ddir = DDIR, verbose=verbose

    if (n_elements(ddir) eq 0L) then ddir = filepath('',root=getenv('NEBULAR_DIR'))
;   if not keyword_set(DDIR) then ddir = '/Users/bemo/idldir/nebular/data/'

;     Read in diagnostic ratios of 

    if (file_test(ddir+'ddiag.dat') eq 0L) or (file_test(ddir+'tdiag.dat') eq 0L) then begin
       print, 'Required file ddiag.dat and/or tdiag.dat not found -- run NEBSETUP.'
       return
    endif
    
    ddiag = struc_init('DDIAG')
    openr, fun, ddir + 'ddiag.dat', /get_lun
    readu, fun, ddiag
    free_lun, fun

    tdiag = struc_init('TDIAG')
    openr, fun, ddir + 'tdiag.dat', /get_lun
    readu, fun, tdiag
    free_lun, fun

;    If the variable 'list' is empty, then assume that the variable
;    'spec' already contains the spectral data in the appropriate format.

    spec = im_wparse(list)

    specdiag, spec, ddiag, tdiag, drats, trats, dsol, tsol
    denstemp, dsol, tsol, dens, temp
    ionic, spec, dens, temp, abund, verbose=verbose


return
end
