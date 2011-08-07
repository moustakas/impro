;+
; NAME:
;     IM_WPARSE()
;
; PURPOSE:
;  reads the files with input spectra into the structure used in
;  the IDL NEBULAR package.
; 
; EXPLANATION:
;  Parses input spectra by assigning the fluxes for a particular
;  wavelength into the corresponding tag in a structure.  The
;  resulting structure is returned by the function.
;
;
; CALLING SEQUENCE:
;
;   spec_structure = IM_WPARSE(infiles [, outfiles, LLIST = ])
;
; INPUT:  
;     INFILES  - list of strings containing the names of the input
;                 spectral files.
;
;
; OUTPUT:
;     SPEC - array of structure containing the line fluxes of each
;             matched transition of interest and the name of the
;             associated infile and outfile.
;
; OPTIONAL INPUT KEYWORD:
;     OUTFILES - list of strings containing the names of the output
;                 spectral files.  Any relevant data for a given element
;                 of 'SPEC' can be written to these files in this and 
;                 other procedures in NEBULAR.  If 'OUTFILES' not given,
;                 the string defaults to that set in 'struc_init.pro',
;                 usually a null string.  For elements of 'OUTFILES'
;                 with a length of zero, the data will be written to the
;                 screen.
;
;     LLIST - the file 'linelist' contains the user-defined database
;             of wavelengths corresponding to a given transition.
;             Alternately, the user can use another file with
;             'LLIST = foo_list'.
;        
; KEYWORDS:
;     VERBOSE - turns on the print commands, writing either to the screen
;               or to each file in the variable 'outfiles'.  The printout
;               consists of the wavelength from the spectrum file along
;               with an annotation declaring the ion and transition the
;               program assumes it represents.
;
; EXAMPLE:
;
;   Read in the spectra listed in the variable 'infiles', assign it to the
;   appropriate transition in 'spec' as outlined in 'llist', assign
;   to each spectrum an output file according to the variable 'outfiles',
;   line fluxes, and print the relevant information to the outfiles by
;   setting the keyword /VERBOSE:
;
;   spec_structure = WPARSE(infiles, outfiles, LLIST = llist, /VERBOSE )
;
; NOTES:
;
;
; REVISION HISTORY:
;       1.0	B.D.Moore	Rice Univ.	September 2003
;-
;
; 
;  Get the user-defined wavelengths from the procedure 'translate'.
;

function im_wparse, infiles, outfiles, llist = LLIST, verbose = VERBOSE, nebdir=nebdir

    if (n_elements(nebdir) eq 0L) then nebdir = filepath('',root=getenv('NEBULAR_DIR'))
;   nebdir = '/Users/bemo/idldir/nebular/'

    if (n_elements(llist) eq 0L) then llist = 'linelist'
;   if (not keyword_set(LLIST)) then llist = nebdir + 'data/linelist'

    loud = keyword_set(VERBOSE)
    translate, nebdir+llist, name, wneb, wuser, annot

;  The structure is created to store the input spectra.  While it may
;  seem needlessly large and complex, creating the structure this way
;  improves the readability of the code as well as simplifying the
;  inclusion of other line fluxes in the future.
;
;  The names of the input and output files are stored in the tags 'infile'
;  and 'outfile', respectively.  The flux at an individual wavelength 'W'
;  (or that of a typically unresolved multiplet) is referenced within the
;  structure with the tag 'L' + 'W'.  The flux of a multiplet referenced
;  at wavelength 'MW' is given by 'LL' + 'MW'.  All wavelengths are
;  quoted in Angstroms unless followed by 'm', which denotes the microns.
;  A detailed annotation for each wavelength is given in 'linelist'.
;
;
    ii = n_elements(infiles)
    spec = struc_init('SPEC',ii)
    spec.infile = infiles
    if keyword_set(OUTFILES) then spec.outfile = outfiles
;
;
    taglist = tag_names(spec)
;   
;   The program transfers data from the input files to the structure.
;   Files are read in by the procedure 'readspec.pro', which allows
;   the user to customize it to fit the layout of his/her spectral files.
;
;   The wavelengths and fluxes are read from the file into the variables
;   'wavein' and 'fluxin'.  The wavelengths are then referenced to the
;   corresponding index wavelength.  It does this by comparing the input
;   wavelength to that defined by the user for the given transition.
;   This is done by adding 'X' to the end of the comparison string,
;   preventing the match between the wavelength and any of 'taglist'.
;   In the (unlikely) event that the input wavelength is exactly between
;   two user-defined wavelengths, the lower of the two is used.  A
;   message is printed to the screen indicating the transition assumed
;   for the wavelength or a failure to do so.
;
;   The nebular package is designed to operate for spectra which do
;   not contain the complete suite of lines in its database.  For
;   example, one often has the total flux from a multiplet but not
;   not those of the individual components.  Since the total flux
;   of a multiplet is referenced by a wavelength intermediate to its
;   components, the program must be able to distinguish between spectra
;   that give only the total flux, only the individual fluxes, or both.
;   This is done by setting a maximum tolerance for deviation, within
;   which the closest input wavelength is assumed to be the target
;   wavelength.  If no lines lie within the tolerance, it is assumed that
;   the flux for that line or multiplet is not given in the input spectrum.
;   For simplicity, the tolerance is set here as equal to one half of the
;   MINIMUM FRACTIONAL spread among the user-defined wavelengths.
;
    maxtol = 1.01
    wsorted = 1.*wuser[sort(wuser)]
    left = (abs(wsorted - shift(wsorted,-1)) < maxtol)/wsorted
    right = (abs(wsorted - shift(wsorted,1)) < maxtol)/wsorted
    toler = (left < right)
;
;   The default Hbeta flux is given by 'def_hbeta'.
;   It is currently set to 1.00, changeable at user discretion.
;
    def_hbeta = 1.00
;
    for jj = 0, ii-1 do begin
       if ( (jj mod 20) eq 0 ) then print, jj, format = '("Reading spectrum", I)'
;
       if ( loud AND (strlen(spec[jj].outfile) gt 0) ) $
         then openu, outlun, spec[jj].outfile, /get_lun, /append $
       else outlun = -1
;
       readspec, infiles[jj], wavein, fluxin
       if loud then begin
          printf, outlun, '---------------------------------------------------'
          printf, outlun, 'Nebular: ', infiles[jj]
          printf, outlun, '---------------------------------------------------'
       endif
       nwave = n_elements(wavein)
       if (nwave le 0) then begin
          print,'No spectral data found in ' + infiles[jj]
       endif else begin
          for kk = 0, nwave-1 do begin
             diff = abs(wuser - wavein[kk])/wavein[kk]
             qmin = min(where(diff eq min(diff)))
             wstr = 'L'+ strtrim(string(wneb[qmin]),2) $
               + strmid('X',0,(diff[qmin] gt toler[kk]))
             wtag = total(1 + where(taglist eq wstr) $
               + where(taglist eq ('L'+ wstr)))
             if (wtag lt 0) then begin
                if loud then printf, outlun, wavein[kk], $
                  format = '("Input wavelength ", F7.2, " unknown")'
             endif else begin
                if loud then printf, outlun, wavein[kk], annot[qmin], $
                  format = '("Input wavelength ", F7.2, " ", A25)'
                spec[jj].(wtag) = fluxin[kk]
             endelse
          endfor                ;         loop kk
       endelse
;
;   The spectrum is normalized relative to Hbeta, with the normalization
;   factor stored in 'spec.norm'.   If the input Hbeta flux us zero, then
;   the spectrum is normalized by the default Hbeta flux 'def_hbeta'.
;
       if (spec[jj].l4861 ne 0.0) then spec[jj].norm = spec[jj].l4861  $
       else spec[jj].norm = def_hbeta
       for hh = 0, n_elements(taglist) - 1 do begin
          if (strmid(taglist[hh], 0, 1) eq 'L') then $
            spec[jj].(hh) = spec[jj].(hh)/spec[jj].norm
       endfor                   ;         loop hh
;
       if loud then begin
          printf, outlun, '---------------------------------------------------'
          printf, outlun, '---------------------------------------------------'
       endif
       if (outlun gt 0) then free_lun, outlun
;
    endfor                      ;         loop jj



return, spec
end
