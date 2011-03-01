pro fits2imh, file_list
;+
; NAME:
;	FITS2IMH
;
; PURPOSE:
;	Convert a list of fits files to IRAF format.
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
;	RDTXT(), READFITS(), IRAFWRT
;
; MODIFICATION HISTORY:
;	John Moustakas, 2000 February 8, UCB
;-

	spawn, ['pwd'], datapath
        imnames = rdtxt(file_list)

        for j = 0L, n_elements(imnames)-1L do begin
            fname = datapath+'/'+imnames[j]
            image = readfits(fname[0]+'.fits',header)
            print, 'Writing '+imnames[j]+'.pix'
            irafwrt, image, header, imnames[j], pixdir=datapath[0]
        endfor

return
end
