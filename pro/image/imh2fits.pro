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
;
; Copyright (C) 2000, John Moustakas
; 
; This program is free software; you can redistribute it and/or modify 
; it under the terms of the GNU General Public License as published by 
; the Free Software Foundation; either version 2 of the License, or
; (at your option) any later version. 
; 
; This program is distributed in the hope that it will be useful, but 
; WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
; General Public License for more details. 
;-

pro imh2fits, file_list

	spawn, ['pwd'], datapath
        imnames = rdtxt(file_list)

        for j = 0L, n_elements(imnames)-1L do begin
            fname = datapath[0]+'/'+imnames[j]
            irafrd, image, header, fname, /silent
            writefits, fname+'fits', image, header
        endfor
            
return
end
