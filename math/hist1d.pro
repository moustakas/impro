;+
; NAME:
;        HIST1D
;
; PURPOSE:
;        Returns the weighted histogram of a 1D array.
;
; CATEGORY:
;        Image processing, statistics, probability.
;
; CALLING SEQUENCE:
;        Result = HIST1D( Array [, Weight, OBIN=Obin, ..other HISTOGRAM keywords] )
;
; INPUTS:
;        Array:    1D array which holds the data to be histogrammed.
;
; OPTIONAL INPUTS:
;       Weight:    1D array of the same dimension as Array which holds the
;             weighted values associated with each Array element.
;
; INPUT KEYWORD PARAMETERS:
;      BINSIZE:    The size of the bin to use. If this keyword is not specified,
;             a bin size of 1 is used.
;
;        INPUT:    1D array to be added to the output of HIST1D.
;
;          MAX:    The maximum value to consider in the histogram.  If this
;             keyword is not specified, Array is searched for its largest value.
;
;          MIN:    The minimum value to consider in the histogram.  If this
;             keyword is not specified, and Array is of type byte, 0 is used.
;             If this keyword is not specified and Array is not of byte type,
;             Array is searched for its smallest value.
;
;      BINEDGE:    This keyword specfies the location of the bin values returned
;             by OBIN. The values are:
;
;                  0 : Center of each bin, [DEFAULT].
;                 -1 : Left  or lower edge of each bin.
;                  1 : Right or upper edge of each bin.
;
; OUTPUTS:
;        1D weighted histogram of Array using the Weight array.  If Weight
;        is not specified, then the usual 1D histogram or density function
;        is returned.
;
; OUTPUT KEYWORD PARAMETERS:
;
;      DENSITY:    Density function of Array; i.e. the unweighted 1D histogram.
;
;         OBIN:    An array holding the bin values of the histogram.
;
;         OMAX:    A named variable that, upon exit, contains the maximum data
;             value used in constructing the histogram.
;
;         OMIN:    A named variable that, upon exit, contains the minimum data
;             value used in constructing the histogram.
;
; REVERSE_INDICES: Set this keyword to a named variable in which the list of
;             reverse indices is returned.  This list is returned as a longword
;             vector whose number of elements is the sum of the number of elements
;             in the histogram, N, and the number of array elements included in
;             the histogram.
;
;             The subscripts of the original array elements falling in the
;             ith bin, 0 le i lt N, are given by the expression: R(R(i):R(i+1)-1),
;             where, R is the reverse index list.  If R(i) is equal to R(i+1),
;             no elements are present in the ith bin.
;
; EXAMPLE:
;        Here'a a very simple example of a weighted histogram:
;
;        z = indgen(10)
;        w = z
;        y = HIST1D( z,w,OBIN=x)
;        plot,x,y,psym=10
;
; MODIFICATION HISTORY:
;        Written by:    Han Wen, July, 1994.
;        01-AUG-1994    Forced H to be double, H=double(Density) avoids
;                       the intermittent bug when H returns as a INTEGER
;                       histogram!
;        31-JAN-1995    Added the ANONYMOUS_ keyword to avoid a very subtle/
;                       rare bug with _EXTRA keyword. (Has no functionality
;                       for the USER -> Internal user only.)
;        21-FEB-1995    Bugfix: didn't add Input to histogram for 1D,
;                       UN-weighted histogram. Also made use of TEMPORARY.
;        09-JAN-1995    Bugfix: change for loop index to long instead of fix.
;                       Improperly indexed R(R(i)) pointers.
;        28-FEB-1996    Added the BINEDGE keyword.
;        06-APR-1996    Assume the data type of the Weight parameter if
;                       specified. (previously defaulted to double)
;-
function HIST1D, Array, Weight, BINSIZE=Binsize, INPUT=Input, $
                   _EXTRA=INPUTKeywords, OMAX=Omax, OMIN=Omin, $
                   OBIN=Obin, REVERSE_INDICES=R, DENSITY=Density, $
                   ANONYMOUS_=Dummy_, BINEDGE=Binedge

         ON_ERROR, 2          ; on error, return control to caller

         if keyword_set( BINSIZE ) then dbin = binsize $
         else dbin = 1.0

         Density = HISTOGRAM( Array, BINSIZE=dbin, _EXTRA=INPUTKeywords, $
                              OMAX=Omax, OMIN=Omin, REVERSE_INDICES=R  )

;   Determine value at the center of each bin
         if (N_ELEMENTS( Binedge ) eq 0) then Binedge = 0
         offset    = (Binedge+1)*0.5
         nbin      = N_ELEMENTS( Density )
         Obin      = Omin + (LINDGEN(nbin)+offset)*dbin
         if N_ELEMENTS( Weight ) eq 0 then begin
              if KEYWORD_SET(INPUT) then return, Density + Input $
              else return, Density
         endif

;   Determine which bins are non-zero
         i_bin     = INDGEN(nbin)
         diff      = R(1:nbin) - R(0:nbin-1)     ; bugfix 1/9/95
         i_gt0     = WHERE( diff ne 0, Ngt0 )
         if Ngt0 eq 0 then return, Density

;   Using reverse-indice information add the weighted values for each
;   data point in each bin.
         dim       = SIZE(Weight)
         data_type = dim(dim(0)+1)
         case data_type of
              1    : H = BYTARR(nbin)
              2    : H = INTARR(nbin)
              3    : H = LONARR(nbin)
              4    : H = FLTARR(nbin)
              5    : H = DBLARR(nbin)
              6    : H = COMPLEXARR(nbin)
              9    : H = DCOMPLEXARR(nbin)
              else : message,'Invalid Weight data type.'
         endcase
         for j=0L,Ngt0-1 do begin
              j_bins = R( R(i_gt0(j)) : R(i_gt0(j)+1)-1 )
              H(i_gt0(j))   = TOTAL( Weight(j_bins) )
         endfor

         if keyword_set( INPUT ) then H = TEMPORARY(H) + Input
         return, H
end
