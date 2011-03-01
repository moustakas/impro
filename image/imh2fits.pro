pro imh2fits, file_list
;+
; NAME:
;	IMH2FITS
;
; PURPOSE:
;	Convert a list of IRAF format images to fits files.
;
; INPUTS:
;	file_list : text file with the files, one per line, to be
;		    converted 
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; PROCEDURES USED:
;	RDTXT(), IRAFRD(), WRITEFITS
;
; MODIFICATION HISTORY:
;	John Moustakas, 2000 February 8, UCB
;-


	spawn, ['pwd'], datapath
        imnames = rdtxt(file_list)

        for j = 0L, n_elements(imnames)-1L do begin
            fname = datapath[0]+'/'+imnames[j]
            irafrd, image, header, fname, /silent
            writefits, fname+'fits', image, header
        endfor
            
return
end
