;+
; NAME:
;     IM_GETMATRIX()
;
; PURPOSE:
;     Returns the energy levels, Einstein A-values, and the matrix
;     of collision strengths for the desired N-level atom at the
;     desired temperature.
;
; EXPLANATION:
;     This procedure is invoked by the IDL NEBULAR procedure 'NLEVEL'.
;     It returns the atomic data necessary to calculate the quasi-static
;     N-level atom for the requested ion.  This in turn is used to
;     calculate the emission-line spectrum for the ion in an electron
;     gas of the requested temperature.
;
; CALLING SEQUENCE:
;     GETMATRIX, IONSTR, TT, QVAL, EN, EINAS
;
; INPUT:
;     IONSTR: string designating the ion for which the atomic data 
;               will be used in the matrix.
;     TT: electron temperature, in K.
;
;
; OUTPUT:
;     QVAL: The matrix of rate coefficients between the various N levels.
;     EN: A vector containg the N energy levels of the ion, in eV.
;     EINAS: Einstein A-values for the N energy levels, in 1/s.
;
;
; OPTIONAL INPUT KEYWORD:
;     None.
;
; EXAMPLE:
;    To get the matrix of Q-values, energy levels, and Einstein A-values
;    of the O++ ion at 10,000 K, the command would be
;
;     GETMATRIX, 'O_III', 1d4, QVAL, EN, EINAS
;
; NOTES:
;     The ions are called by the spectroscopic notation for the ion. The
;     underscore is used to separate the atomic symbol from the Roman
;     numeral designation, to avoid any ambiguities.  So, for example,
;     the S+ ion is called with 'S_II', O++ is called with 'O_III', &c.
;     The calling string is not case-sensitive.
;
;
; REVISION HISTORY:
;       1.0	B.D.Moore	Rice Univ.	September 2003
;-
;      J. Moustakas, 2007-Nov-27, NYU: now a function; pack everything
;         into a structure; default 10,000 K temperature assumed
;      jm07dec05nyu - atomic data update: see comments below, Froese
;                     Fischer et al. 2004, and
;                     http://atoms.vuse.vanderbilt.edu 
;      j11sep07ucsd - bug fix on P3C energy level diagram (for N_I)
;
;
;
;      *****     CLOUDY 94.00 "COMPATIBLE"     *****
;
;  The atomic data has been drawn from the same literature noted in the
;  relevant subroutines in CLOUDY 94.00.  In cases where this was not
;  possible the relevant references are included.
;
;  Calculate log10(tt), and set conversion from various energy units to eV.
;

function im_getmatrix, ionstr, tt, cloudy94=cloudy94

    if (n_elements(ionstr) eq 0L) then begin
       doc_library, 'im_getmatrix'
       return, -1L
    endif

    if (n_elements(tt) eq 0L) then tt = 1D4 ; [K]
    
    logt = alog10(tt)
    invcm = 12398.5d-8
    ryd = invcm*109733.5522d    ; from http://hf8.vuse.vanderbilt.edu/
;   ryd =  13.606d0 
    ktoev = 1.0d0/1.16045d4
    yaff = 1.438783933d0
;
    istr = strcompress(strupcase(ionstr),/remove_all)
    case strupcase(istr) of
;
;
       'C_I'  : BEGIN
;
;  Set ion type, energy levels and Einstein A-values.  Energy levels,
;  A-values, and collision strengths from Mendoza 1982 except P?-P?
;  A-values, from Tielens and Hollenbach 1985.
;
          iontype = 'p2'
          ee = invcm*[0, 16.4, 43.4, 10192.6, 21648.0, 33735.2]
;                       ^    ^     ^     ^        ^        ^
;                      3P0  3P1   3P2   1D2      1S0      5S2
;
          einas = dblarr(6,6)
          einas[1,0] = 7.93d-8
          einas[2,0:1] = [2.00d-14, 2.68d-7]
          einas[3,0:2] = [7.77d-8, 8.21d-5, 2.44d-4]
          einas[4,1:3] = [2.71d-3, 2.00d-5, 5.28d-1]
          einas[5,1:3] = [6.94d0, 1.56d1, 2.35d-5]
;
;  Collision strengths, interpolated for temperature 'logt'.  P?-P?
;  collision strengths from Johnson et al. 1985, with those for 15000 K
;  and 20000 K extrapolated from the results for 5000 K and 10000 K.
;
          tgrid = alog10(5.0d3 + 5.0d3*findgen(4))
          p1p0 = spline(tgrid, [2.43, 3.71, 4.66, 6.12], logt)*1d-1
          p2p0 = spline(tgrid, [1.82, 2.46, 3.00, 3.95], logt)*1d-1
          p2p1 = spline(tgrid, [7.14, 10.2, 12.6, 16.7], logt)*1d-1
          d2px = spline(tgrid, [6.03, 11.4, 16.0, 19.6], logt)*1d-1/9.
          s0px = spline(tgrid, [1.49, 2.52, 3.20, 3.65], logt)*1d-1/9.
          s0d2 = spline(tgrid, [1.96, 2.77, 3.40, 3.92], logt)*1d-1
          s2px = spline(tgrid, [4.75, 6.71, 8.22, 9.50], logt)*1d-1/9.
;
       END
;
;
       'C_II' : BEGIN
;
;
;  Set ion type, energy levels and Einstein A-values.
;
;  Energy levels from Edlen 1981 as cited in Lennon et al. 1985, Table 4.
;
          iontype = 'p1'
          ee = invcm*[0,   63.4,  43003.3,  43025.3, 43053.6, $
;                       ^     ^       ^         ^        ^
;                      2P1/2 2P3/2   4P1/2     4P3/2    4P5/2
;
            74930.1,   74932.6]
;                               ^          ^
;                              2D5/2      2D3/2
;
;
;  Einstein A-values from Lennon et al., ApJ, 294, 200, with
;  additional values from Pradhan and Peng, STScI #8.
;  A(2P3/2-2P1/2) transition rate from Froese Fischer 1983.
;
          einas = dblarr(7,7)
          einas[1,0] = 2.0494d-6
          einas[2,0:1] = [7.440d1, 7.780d1]
          einas[3,0:2] = [1.700d0, 1.240d1, 2.390d-7]
          einas[4,1:3] = [5.390d1, 3.490d-14, 3.670d-7]
          einas[5,0:2] = [0.000d0, 2.840d8, 0.000d0]
          einas[6,0:2] = [2.380d8, 4.700d7, 0.000d0]
          
;  Collision strengths, interpolated for temperature 'logt', from
;  Blum and Pradhan 1992.
          
          tgrid = alog10(4.0d3 + 4.0d3*findgen(5))
          pl3pl1 = spline(tgrid, [1.8034,2.0756,2.2021,2.2562,2.2776], $
            logt)*1d0
          pu1pl1 = spline(tgrid, [2.459, 2.413, 2.439, 2.464, 2.475], $
            logt)*1d-1
          pu1pl3 = spline(tgrid, [1.748, 1.749, 1.793, 1.826, 1.842], $
            logt)*1d-1
          pu3pl1 = spline(tgrid, [3.650, 3.595, 3.645, 3.687, 3.707], $
            logt)*1d-1
          pu3pl3 = spline(tgrid, [4.763, 4.728, 4.821, 4.893, 4.928], $
            logt)*1d-1
          pu3pu1 = spline(tgrid, [6.358, 7.564, 8.885, 9.869,10.609], $
            logt)*1d-1
          pu5pl1 = spline(tgrid, [2.304, 2.315, 2.382, 2.429, 2.453], $
            logt)*1d-1
          pu5pl3 = spline(tgrid, [1.0315,1.0170,1.0317,1.0440,1.0500], $
            logt)*1d0
          pu5pu1 = spline(tgrid, [7.078, 8.070, 8.910, 9.423, 9.705], $
            logt)*1d-1
          pu5pu3 = spline(tgrid, [1.5926,1.8499,2.0969,2.2721,2.3880], $
            logt)*1d0
          d5pl1  = spline(tgrid, [1.3094,1.3669,1.4240,1.4749,1.5204], $
            logt)*1d0
          d5pl3  = spline(tgrid, [8.793, 9.143, 9.449, 9.652, 9.778], $
            logt)*1d-1
          d5pu1  = spline(tgrid, [4.404, 4.614, 4.793, 4.936, 5.046], $
            logt)*1d-1
          d5pu3  = spline(tgrid, [6.835, 7.187, 7.478, 7.708, 7.885], $
            logt)*1d-1
          d5pu5  = spline(tgrid, [5.318, 5.678, 5.946, 6.152, 6.310], $
            logt)*1d-1
          d3pl1  = spline(tgrid, [5.145, 5.341, 5.501, 5.585, 5.615], $
            logt)*1d-1
          d3pl3  = spline(tgrid, [2.7685,2.8877,3.0032,3.1016,3.1856], $
            logt)*1d0
          d3pu1  = spline(tgrid, [2.495, 2.669, 2.797, 2.896, 2.971], $
            logt)*1d-1
          d3pu3  = spline(tgrid, [6.963, 7.379, 7.703, 7.955, 8.149], $
            logt)*1d-1
          d3pu5  = spline(tgrid, [1.5379,1.6171,1.6826,1.7343,1.7741], $
            logt)*1d0
          d3d5   = spline(tgrid, [1.9846,1.9086,1.8848,1.8735,1.8663], $
            logt)*1d0
;
       END
;
;
       'C_III': BEGIN
;
;  Set ion type, energy levels and Einstein A-values.
;
;
;  Energies and Einstein A's from Fleming et al., MNRAS, 279, 1289
;  except 3P1-1S0, from Kwong et al., ApJ 411, 341.
;
          iontype = 's2'
          ee = 2.7212d1*[0, 0.238350, 0.238456, 0.238710, 0.466343]
;                          ^     ^         ^         ^         ^
;                         1S0   3P0       3P1       3P2       1P1U
;
          einas = dblarr(5,5)
          einas[1,0] = 0.000d0
          einas[2,0:1] = [1.210d1, 2.270d-7]
          einas[3,0:2] = [5.19d-3, 1.370d-13, 2.330d-6]
          einas[4,0] = 1.790d9
;
;  Collision strengths, interpolated for temperature 'logt'.
;  Berrington et al., ADNDT, 33. 195.  Note that the collision
;  strengths are extrapolated for T < 12500 K.
;
          tgrid = 4.1d0 + findgen(9)/5.0d0
          pxs0 = spline(tgrid, [10.4, 10.3, 10.3, 9.66, $
            9.16, 8.08, 6.92, 5.79, 4.73], logt)*1d-1/9.
          p1p0 = spline(tgrid, [9.62, 10.3, 10.9, 11.1, $
            10.6, 9.53, 8.15, 6.73, 5.41], logt)*1d-1
          p2p0 = spline(tgrid, [7.18, 8.57, 10.1, 11.1, $
            10.9, 9.88, 8.43, 7.03, 5.88], logt)*1d-1
          p2p1 = spline(tgrid, [2.78, 3.17, 3.60, 3.83, $
            3.74, 3.37, 2.88, 2.39, 1.97], logt)*1d0
          pus0 = spline(tgrid, [4.00, 4.11, 4.24, 4.42, $
            4.65, 5.00, 5.50, 6.17, 7.01], logt)*1d0
          pupx = spline(tgrid, [3.70, 3.54, 3.29, 2.93, $
            2.50, 2.07, 1.68, 1.38, 1.20], logt)*1d0/9.
;
       END
;
;
       'N_I'  : BEGIN
;
;  Set ion type, energy levels and Einstein A-values.
;
;     Energy levels and A-values from Butler and Zeippen 1984,
;     assuming that A = A(E2) + A(M1r).  A(2P3/2-2P1/2) from
;     Pradhan and Peng 1992.
;
          iontype = 'p3c'
          ee = ryd*[0, 0.175189, 0.175273, 0.262810, 0.262183]
;                     ^     ^         ^         ^         ^
;                    4S3/2 2D5/2     2D3/2     2P1/2     2P3/2
;
          einas = dblarr(5,5)
          einas[1,0] = 6.124d-6
          einas[2,0:1] = [2.278d-5, 1.239d-8]
          einas[3,0:2] = [2.716d-3, 3.135d-2, 4.804d-2]
          einas[4,0:3] = [6.604d-3, 5.589d-2, 2.523d-2, 5.17d-13]
;
;  Collision strengths, interpolated for temperature 'logt'.  CS from
;  Dopita et al. 1976.
;
          tgrid = alog10(5.0d3 + 5.0d3*findgen(4))
          dxs3 = spline(tgrid, [2.34, 4.80, 6.63, 8.00], logt)*1d-1/10.
          pxs3 = spline(tgrid, [8.30, 17.1, 24.2, 29.8], logt)*1d-2/6.
          d3d5 = spline(tgrid, [1.28, 2.69, 3.79, 4.65], logt)*1d-1
          p1d5 = spline(tgrid, [6.26, 10.9, 15.1, 19.0], logt)*1d-2
          p3d5 = spline(tgrid, [1.62, 2.66, 3.52, 4.38], logt)*1d-1
          p1d3 = spline(tgrid, [6.01, 9.70, 12.8, 15.7], logt)*1d-2
          p3d3 = spline(tgrid, [8.56, 14.7, 20.2, 25.2], logt)*1d-2
          p3p1 = spline(tgrid, [3.29, 7.10, 11.1, 15.3], logt)*1d-2
;
       END
;
;
       'N_II' : BEGIN

; jm07dec05nyu: Energy levels from NIST http://physics.nist.gov) and
; Einstein A values from the Chianti project
; (http://www.arcetri.astro.it/science/chianti/chianti.html), as
; tabulated in Stasinska 2007 (astro-ph/0704.0348); A-values for 3071,
; 3177, and 76 micron do not appear in the Chianti database, and have
; been set to 1d-20; additional info here:
; http://wwwsolar.nrl.navy.mil/database/n/n_2_table.html

          iontype = 'p2'
          ee = invcm*[0.0, 48.67, 130.80, 15316.19, 32688.64, 46784.56]
;                      ^     ^       ^         ^         ^         ^
;                     3P0   3P1     3P2       1D2       1S0       5S2

          einas = dblarr(6,6)
          einas[1,0]   = [2.080d3-6]                  ; [205 mu]
          einas[2,0:1] = [1d-20,7.460d-6]             ; [76 mu,121 mu]
          einas[3,0:2] = [1.928d-6,9.819d-4,3.015d-3] ; [6529,6549,6585]
          einas[4,1:3] = [4.588d-6,1d-20,9.234d-1]    ; [3063,3071,5756]
          einas[5,1:3] = [7.17d1,1.77d2,1d-20]        ; [2139,2143,3177]

;;;  Set ion type, energy levels and Einstein A-values.  jm07dec05nyu
;;;  update: Energy levels and Einstein A values from Froese Fischer et
;;;  al. 2004 [Table 3, C-like, Z=7]
;;
;;          iontype = 'p2'
;;          ee = invcm*[0.0, 48.74, 130.69, 15316.61, 32688.35, 46784.03]
;;;                      ^     ^       ^         ^         ^         ^
;;;                     3P0   3P1     3P2       1D2       1S0       5S2
;;
;;          einas = dblarr(6,6)
;;          einas[1,0]   = 2.08d3-6                       ; [205 mu]
;;          einas[2,0:1] = [1.116d-12, 7.420d-6]          ; [76 mu,122 mu]
;;          einas[3,0:2] = [5.253d-7, 9.842d-4, 2.905d-3] ; [6528,6549,6585]
;;          einas[4,1:3] = [3.185d-2, 1.547d-4, 1.136d0]  ; [3063,3071,5756]
;;          einas[5,1:3] = [5.155d1, 1.266d2, 8.949d-4]   ; [2139,2143,3177]

;;;  A-values from Nussbaumer and Rusca 1979, except A(5S2-3P{2,1}) from
;;;  Brage et al. 1997
;;;
;;            iontype = 'p2'
;;            ee = invcm*[0, 48.7, 130.8, 15316.3, 32688.9, 46857.7]
;;;                       ^    ^     ^      ^        ^        ^
;;;                      3P0  3P1   3P2    1D2      1S0      5S2
;;
;;            einas = dblarr(6,6)
;;            einas[1,0] = 2.08d-6
;;            einas[2,0:1] = [1.16d-12, 7.46d-6]
;;            einas[3,0:2] = [5.35d-7, 9.24d-4, 2.73d-3]
;;            einas[4,1:3] = [3.16d-2, 1.51d-4, 1.17d0]
;;            einas[5,1:3] = [5.36d1, 1.306d2, 4.92d-4]

;  Collision strengths, interpolated for temperature 'logt'.  CS from
;  Stafford et al. 1994.
;
          tgrid = alog10(5.0d3 + 5.0d3*findgen(4))
          p1p0 = spline(tgrid, [4.097, 4.232, 4.371, 4.491], logt)*1d-1
          p2p0 = spline(tgrid, [2.439, 2.657, 2.859, 3.019], logt)*1d-1
          p2p1 = spline(tgrid, [1.061, 1.127, 1.190, 1.241], logt)*1d0
          s2px = spline(tgrid, [1.155, 1.146, 1.152, 1.157], logt)*1d0/9.
          d2px = spline(tgrid, [3.047, 3.007, 2.995, 2.990], logt)*1d0/9.
          s0px = spline(tgrid, [3.812, 3.694, 3.656, 3.642], logt)*1d-1/9.
          s0d2 = spline(tgrid, [5.027, 5.014, 5.071, 5.140], logt)*1d-1
;
       END
;
;
       'N_III': BEGIN
;
;  Set ion type, energy levels and Einstein A-values.
;
;  Energy levels and collision strengths from Blum and Pradhan 1992.
;
;
          iontype = 'p1'
          ee = ryd*[0, 0.001590, 0.521173, 0.521719, 0.522459, $
;                     ^      ^         ^         ^         ^
;                    2P1/2  2P3/2     4P1/2     4P3/2     4P5/2
;
            0.920597, 0.920667]
;                           ^         ^
;                          2D5/2     2D3/2
;
;  Einstein A-values from Brage et al. 1995, except A(2P3/2-2P1/2) from
;  Pradhan and Peng, STScI #8.  Note that CLOUDY uses Stafford et al. 1994,
;  which only treats total level transition rates. jm07nov27nyu update
;  for 6548&6584: Storey & Zeippen (2000) 
;
          einas = dblarr(7,7)
          einas[1,0] = 4.77d-5
          einas[2,0:1] = [3.61d2, 3.72d1]
          einas[3,0:2] = [9.11d0, 6.51d1, 0.00d0]
          einas[4,1:3] = [2.82d2, 0.00d0, 0.00d0]
          einas[5,0:2] = [4.27d8, 8.34d7, 0.00d0]
          einas[6,1:3] = [5.08d8, 0.00d0, 0.00d0]
;
;  Collision strengths, interpolated for temperature 'logt'.  CS from
;  Blum and Pradhan 1992.
;
          tgrid = alog10(4.0d3 + 4.0d3*findgen(5))
          pl3pl1 = spline(tgrid, [1.2976,1.3937,1.4922,1.5729,1.6421], $
            logt)*1d0
          pu1pl1 = spline(tgrid, [1.879, 1.944, 2.009, 2.049, 2.068], $
            logt)*1d-1
          pu1pl3 = spline(tgrid, [1.323, 1.445, 1.558, 1.637, 1.687], $
            logt)*1d-1
          pu3pl1 = spline(tgrid, [2.784, 2.910, 3.034, 3.114, 3.157], $
            logt)*1d-1
          pu3pl3 = spline(tgrid, [3.619, 3.866, 4.100, 4.259, 4.353], $
            logt)*1d-1
          pu3pu1 = spline(tgrid, [0.9904,1.0745,1.1207,1.1446,1.1581], $
            logt)*1d0
          pu5pl1 = spline(tgrid, [1.741, 1.923, 2.091, 2.210, 2.285], $
            logt)*1d-1
          pu5pl3 = spline(tgrid, [7.864, 8.242, 8.610, 8.848, 8.980], $
            logt)*1d-1
          pu5pu1 = spline(tgrid, [5.959, 6.502, 6.804, 6.983, 7.108], $
            logt)*1d-1
          pu5pu3 = spline(tgrid, [1.8323,1.9927,2.0813,2.1301,2.1605], $
            logt)*1d0
          d5pl1  = spline(tgrid, [1.7651,1.8062,1.8263,1.8371,1.8433], $
            logt)*1d0
          d5pl3  = spline(tgrid, [1.0242,1.0194,1.0232,1.0237,1.0198], $
            logt)*1d0
          d5pu1  = spline(tgrid, [4.968, 5.010, 5.023, 5.016, 5.002], $
            logt)*1d-1
          d5pu3  = spline(tgrid, [7.840, 7.882, 7.896, 7.886, 7.867], $
            logt)*1d-1
          d5pu5  = spline(tgrid, [6.520, 6.475, 6.470, 6.465, 6.457], $
            logt)*1d-1
          d3pl1  = spline(tgrid, [5.593, 5.485, 5.483, 5.469, 5.426], $
            logt)*1d-1
          d3pl3  = spline(tgrid, [3.6247,3.6900,3.7260,3.7442,3.7520], $
            logt)*1d0
          d3pu1  = spline(tgrid, [3.085, 3.059, 3.056, 3.053, 3.051], $
            logt)*1d-1
          d3pu3  = spline(tgrid, [8.267, 8.257, 8.261, 8.253, 8.238], $
            logt)*1d-1
          d3pu5  = spline(tgrid, [1.7641,1.7734,1.7766,1.7743,1.7700], $
            logt)*1d0
          d3d5   = spline(tgrid, [1.2191,1.3452,1.4475,1.5183,1.5663], $
            logt)*1d0
;
       END
;
;
       'O_I'  : BEGIN
;
;  Set ion type, energy levels and Einstein A-values.
;
          iontype = 'p4'
          ee = invcm*[0, 158.3, 227.0, 15867.9, 33792.6]
;                       ^   ^      ^       ^        ^
;                      3P2 3P1    3P0     1D2      1S0
;
;  A-values from Mendoza 1982.
;
          einas = dblarr(5,5)
          einas[1,0] = 8.92d-5
          einas[2,0:1] = [1.275d-10, 1.74d-5]
          einas[3,0:2] = [6.34d-3, 2.11d-3, 7.23d-7]
          einas[4,0:3] = [2.88d-4, 7.32d-2, 0.00d0, 1.22d0]
;
;  Collision strengths, interpolated for temperature 'logt'.  CS from
;  Berrington 1988 for 5000K and 1000K, and Mendoza 1982 for 20000 K.
;
          tgrid = alog10([5.0d3, 1.0d4, 2.0d4])
          p1p2 = spline(tgrid, [4.52, 10.6, 20.7], logt)*1d-2
          p0p2 = spline(tgrid, [1.53, 3.21, 5.36], logt)*1d-2
          p0p1 = spline(tgrid, [9.12, 28.3, 69.3], logt)*1d-3
          d2px = spline(tgrid, [1.24, 2.66, 5.01], logt)*1d-1/9.
          s0px = spline(tgrid, [1.53, 3.24, 6.07], logt)*1d-2/9.
          s0d2 = spline(tgrid, [7.32, 10.5, 14.8], logt)*1d-2
;
       END
;
;
       'O_II' : BEGIN

; jm07dec05nyu: Energy levels from NIST http://physics.nist.gov) and
; Einstein A values from the Chianti project
; (http://www.arcetri.astro.it/science/chianti/chianti.html), as
; tabulated in Stasinska 2007 (astro-ph/0704.0348); the A-value for
; 5030 micron does not appear in the Chianti database, and has been
; set to 1d-20; additional info here:
; http://wwwsolar.nrl.navy.mil/database/o/o_2_table.html

          iontype = 'p3a'
          ee = invcm*[0.0, 26810.55, 26830.57, 40468.01, 40470.00]
;                      ^        ^         ^         ^         ^
;                     4S3/2   2D5/2     2D3/2     2P3/2     2P1/2

          einas = dblarr(5,5)
          einas[1,0]   = [3.588d-5]                         ; [3729]
          einas[2,0:1] = [1.810d-4,1.320d-7]                ; [3727,499 mu]
          einas[3,0:2] = [5.800d-2,1.074d-1,5.800d-2]       ; [2471,7322,7332]
          einas[4,0:3] = [2.380d-2,5.630d-2,9.410d-2,1d-20] ; [2470,7320,7331,5030 mu]
          
;;; jm07dec05nyu: Energy levels and Einstein A values from Froese
;;;  Fischer et al. 2004 [Table 4, N-like, Z=8]
;;
;;          iontype = 'p3a'
;;          ee = invcm*[0.0, 26810.73, 26830.45, 40468.36, 40470.96]
;;;                      ^        ^         ^         ^         ^
;;;                     4S3/2   2D5/2     2D3/2     2P3/2     2P1/2
;;
;;          einas = dblarr(5,5)
;;          einas[1,0]   = 3.382d-5                                ; [3729]
;;          einas[2,0:1] = [1.414d-4,1.241d-7]                     ; [3727,507 mu]
;;          einas[3,0:2] = [5.646d-2,1.018d-1,4.313d-2]            ; [2471,7321,7332]
;;          einas[4,0:3] = [2.265d-2,5.824d-2,8.694d-2,3.158d-10]  ; [2470,7320,7331,3843 mu]
          
;;; A-values from NIST, as quoted in CLOUDY (cooloxyg.c)
;;
;;            iontype = 'p3a'
;;            ee = 1.d0*[0, 3.32506, 3.32756, 5.01888, 5.01913]
;;;                      ^     ^        ^        ^        ^
;;;                     4S3/2 2D5/2    2D3/2    2P3/2    2P1/2
;;
;;            einas = dblarr(5,5)
;;            einas[1,0] = 3.058d-5
;;            einas[2,0:1] = [1.776d-4, 1.30d-7]
;;            einas[3,0:2] = [5.22d-2, 9.907d-2, 5.34d-2]
;;            einas[4,0:3] = [2.12d-2, 5.19d-2, 8.672d-2, 1.41d-10]

;  Collision strengths, interpolated for temperature 'logt'.  First set
;  of CS from McLaughlin and Bell 1993.  Second set of CS from Pradhan
;  and Peng 1992, with scale factor for P?-D? suggested by McLaughlin
;  and Bell 1993 as used in CLOUDY (cooloxyg.c).
;
          tgrid = alog10([5.0d3, 1.0d4, 1.6d4, 2.0d4])
          dxs3 = spline(tgrid, [1.361, 1.375, 1.391, 1.401], logt)*1d0/10.
          pxs3 = spline(tgrid, [4.052, 4.147, 4.243, 4.304], logt)*1d-1/6.
;
          tgrid = alog10(5.0d3 + 5.0d3*findgen(4))
          d3d5 = spline(tgrid, [1.22, 1.17, 1.14, 1.11], logt)*1d0
          p1p3 = spline(tgrid, [2.80, 2.87, 2.93, 3.00], logt)*1d-1
          p3d5 = spline(tgrid, 1.207*[7.18, 7.30, 7.41, 7.55], logt)*1d-1
          p1d5 = spline(tgrid, 1.207*[2.90, 2.95, 3.00, 3.05], logt)*1d-1
          p3d3 = spline(tgrid, 1.207*[4.01, 4.08, 4.14, 4.22], logt)*1d-1
          p1d3 = spline(tgrid, 1.207*[2.70, 2.75, 2.81, 2.84], logt)*1d-1
;
       END
;
;
       'O_III' : BEGIN

; jm07dec05nyu update: Energy levels from NIST
; http://physics.nist.gov) and Einstein A values from the Chianti
; project (http://www.arcetri.astro.it/science/chianti/chianti.html),
; as tabulated in Stasinska 2007 (astro-ph/0704.0348), except for
; 4364, whose A-value is given as 0.4760 in the Chianti database, and
; 1.71 in NIST (makes a ~1000 K difference in the inferred Te
; values!!); A-values for 2332 and 2496 do not appear in the Chianti
; database, and have been set to 1d-20; additional info here:
; http://wwwsolar.nrl.navy.mil/database/o/o_3_table.html

          iontype = 'p2'
          ee = invcm*[0.0, 113.18, 306.17, 20273.27, 43185.74, 60324.79]
;                      ^      ^       ^         ^         ^         ^
;                     3P0    3P1     3P2       1D2       1S0       5S2

          einas = dblarr(6,6)
          einas[1,0]   = [2.70d-5]                    ; [88 mu]
          einas[2,0:1] = [3.032d-11,1.010d-4]         ; [32 mu,51 mu]
          einas[3,0:2] = [7.250d-6,6.791d-3,2.046d-2] ; [4932,4960,5008]
          einas[4,1:3] = [2.26d-1,1d-20,1.71d0]       ; [2321,2332,4364]
          einas[5,1:3] = [1.45d2,4.26d2,1d-20]        ; [1660,1666,2496]

;;;  jm07dec05nyu: Energy levels and Einstein A values from Froese
;;;  Fischer et al. 2004  [Table 3, C-like, Z=8]
;;
;;          iontype = 'p2'
;;          ee = invcm*[0.0, 113.03, 305.61, 20273.11, 43186.75, 60324.71]
;;;                      ^      ^       ^         ^         ^         ^
;;;                     3P0    3P1     3P2       1D2       1S0       5S2
;;
;;          einas = dblarr(6,6)
;;          einas[1,0]   = [2.596d-5]                    ; [88 mu]
;;          einas[2,0:1] = [3.032d-11,9.632d-5]          ; [32 mu,52 mu]
;;          einas[3,0:2] = [2.322d-6,6.946d-3,2.025d-2]  ; [4932,4960,5008]
;;          einas[4,1:3] = [2.255d-1,6.998d-4,1.685d0]   ; [2321,2332,4364]
;;          einas[5,1:3] = [2.308d2,5.765d2,5.777d-3]    ; [1660,1666,2496]

;;;  Energy levels from Mendoza 1982; A values from Wiese et al. 1996,
;;;  except 5S2-3P{1,2} from Mendoza 1982.
;;;
;;            iontype = 'p2'
;;            ee = invcm*[0, 113.2, 306.2, 20273.3, 43185.7, 60325.0]
;;;                       ^   ^      ^       ^        ^        ^
;;;                      3P0 3P1    3P2     1D2      1S0      5S2
;;
;;            einas = dblarr(6,6)
;;            einas[1,0] = 2.62d-5
;;            einas[2,0:1] = [3.17d-11, 9.76d-5]
;;            einas[3,0:2] = [2.41d-6, 6.21d-3, 1.81d-2]
;;            einas[4,1:3] = [2.15d-1, 6.34d-4, 1.71d0]
;;            einas[5,1:3] = [2.12d2, 5.22d2, 6.36d-3]

;  Collision strengths, interpolated for temperature 'logt'.  CS from
;  Lennon and Burke 1994, except 1S0-1D2 from Burke et al. 1989.
;
          tgrid = 3.0d0 + findgen(11)/5.0d0
          p1p0 = spline(tgrid, [4.975, 5.066, 5.115, 5.180, 5.296, $
            5.454, 5.590, 5.678, 5.788, 5.918, 5.938], logt)*1d-1
          p2p0 = spline(tgrid, [2.455, 2.493, 2.509, 2.541, 2.609, $
            2.713, 2.832, 2.955, 3.101, 3.254, 3.314], logt)*1d-1
          p2p1 = spline(tgrid, [1.173, 1.193, 1.203, 1.218, 1.248, $
            1.291, 1.335, 1.373, 1.419, 1.468, 1.482], logt)*1d0
          s2px = spline(tgrid, [9.760, 9.673, 9.712, 10.224, 11.196, $
            12.074, 12.574, 12.720, 12.451, 11.704, 10.600], logt)*1d-1/9.
          
          d2px = spline(tgrid, [2.2233, 2.1888, 2.1416, 2.1117, 2.1578, $
            2.2892, 2.4497, 2.5851, 2.6730, 2.7019, 2.6594], logt)*1d0/9.
          s0px = spline(tgrid, [2.754, 2.738, 2.713, 2.693, 2.747, 2.925, $
            3.174, 3.405, 3.563, 3.621, 3.571], logt)*1d-1/9.
          
          tgrid = alog10([6000., 10000., 15000., 20000.])
          s0d2 = spline(tgrid, [6.17, 6.77, 6.80, 6.64], logt)*1d-1
;
          d2px = (3.0211144d0 - 101.57536d0/sqrt(tt) + 817.57670d0*logt/tt)/9.
          s0px = (0.0784d0*tt^0.143 < 0.36d0)/9.
          s0d2 = 0.32412181d0 + 79.051672d0/sqrt(tt) - 4374.7816d0/tt

       END
;
;
       'O_IV' : BEGIN
;
;  Set ion type, energy levels and Einstein A-values.
;
;  Energy levels and collision strengths from Blum and Pradhan 1992.
;  It is assumed that the energies in Table 5 are out of order, and
;  are entered here in increasing order.
;
          iontype = 'p1'
          ee = ryd*[0, 0.003522, 0.648613, 0.649810, 0.651491, $
;                     ^      ^         ^         ^         ^
;                    2P1/2  2P3/2     4P1/2     4P3/2     4P5/2
;
            1.156729, 1.156856]
;                           ^         ^
;                          2D5/2     2D3/2
;
;  Einstein A-values from Brage et al. 1996, except 2P3/2-2P1/2
;  from Pradhan and Peng 1992.       
;
          einas = dblarr(7,7)
          einas[1,0] = 5.18d-4
          einas[2,0:1] = [1.47d3, 1.43d3]
          einas[3,0:2] = [3.84d1, 2.94d2, 0.00d0]
          einas[4,1:3] = [1.17d3, 0.00d0, 1.02d-4]
          einas[5,0:2] = [6.08d8, 1.17d8, 0.00d0]
          einas[6,1:3] = [7.18d8, 0.00d0, 0.00d0]
;
;  Collision strengths, interpolated for temperature 'logt'.  CS from
;  Blum and Pradhan 1992.  CLOUDY uses Zhang et al. 1994, but spot
;  checks show good (few %) agreement.
;
          tgrid = alog10(4.0d3 + 4.0d3*findgen(5))
          pl3pl1 = spline(tgrid, [1.9103,2.3072,2.4862,2.5482,2.5795], $
            logt)*1d0
          pu1pl1 = spline(tgrid, [1.195, 1.284, 1.375, 1.442, 1.489], $
            logt)*1d-1
          pu1pl3 = spline(tgrid, [0.849, 0.958, 1.080, 1.179, 1.252], $
            logt)*1d-1
          pu3pl1 = spline(tgrid, [1.774, 1.924, 2.082, 2.202, 2.288], $
            logt)*1d-1
          pu3pl3 = spline(tgrid, [2.314, 2.560, 2.827, 3.039, 3.195], $
            logt)*1d-1
          pu3pu1 = spline(tgrid, [1.0418,1.0684,1.1067,1.1417,1.1709], $
            logt)*1d0
          pu5pl1 = spline(tgrid, [1.119, 1.275, 1.453, 1.597, 1.706], $
            logt)*1d-1
          pu5pl3 = spline(tgrid, [5.012, 5.451, 5.911, 6.264, 6.519], $
            logt)*1d-1
          pu5pu1 = spline(tgrid, [7.359, 6.883, 6.923, 7.138, 7.397], $
            logt)*1d-1
          pu5pu3 = spline(tgrid, [2.0672,2.0349,2.0813,2.1467,2.2114], $
            logt)*1d0
          d5pl1  = spline(tgrid, [1.6472,1.7811,1.8913,2.0046,2.1065], $
            logt)*1d0
          d5pl3  = spline(tgrid, [0.6397,0.9123,1.1009,1.2349,1.3296], $
            logt)*1d0
          d5pu1  = spline(tgrid, [3.747, 4.078, 4.259, 4.325, 4.333], $
            logt)*1d-1
          d5pu3  = spline(tgrid, [5.841, 6.360, 6.647, 6.757, 6.778], $
            logt)*1d-1
          d5pu5  = spline(tgrid, [4.631, 5.052, 5.293, 5.404, 5.449], $
            logt)*1d-1
          d3pl1  = spline(tgrid, [2.586, 4.634, 6.021, 6.950, 7.568], $
            logt)*1d-1
          d3pl3  = spline(tgrid, [3.1719,3.5768,3.8862,4.1643,4.3971], $
            logt)*1d0
          d3pu1  = spline(tgrid, [2.178, 2.376, 2.490, 2.544, 2.567], $
            logt)*1d-1
          d3pu3  = spline(tgrid, [6.008, 6.548, 6.852, 6.981, 7.021], $
            logt)*1d-1
          d3pu5  = spline(tgrid, [1.3142,1.4311,1.4955,1.5203,1.5250], $
            logt)*1d0
          d3d5   = spline(tgrid, [2.9983,3.5450,3.9576,4.0895,4.0731], $
            logt)*1d0
;
       END
;
;
       'O_V' : BEGIN
;
;  Set ion type, energy levels and Einstein A-values.
;  Energies and Einstein A's from Fleming et al., MNRAS, 279, 1289
;
          iontype = 's2'
          ee = 2.7212d1*[0, 0.373359, 0.373982, 0.375362, 0.723559]
;                          ^     ^         ^         ^         ^
;                         1S0   3P0       3P1       3P2       1P1U
;
          einas = dblarr(5,5)
          einas[1,0] = 0.000d0
          einas[2,0:1] = [2.280d3, 4.620d-5]
          einas[3,0:2] = [2.160d-2, 1.280d-10, 3.800d-4]
          einas[4,0] = 2.800d9
;
;  Collision strengths, interpolated for temperature 'logt'.  CS from
;  Berrington et al. 1985.
;
          tgrid = 4.5d0 + findgen(9)/5.0d0
          pxs0 = spline(tgrid, [5.98, 5.55, 5.02, 4.40, 3.76, $
            3.16, 2.63, 2.16, 1.74], logt)*1d-1/9.
          p1p0 = spline(tgrid, [8.70, 8.42, 7.76, 6.80, 5.71, $
            4.63, 3.68, 2.88, 2.22], logt)*1d-1
          p2p0 = spline(tgrid, [10.5, 10.2, 9.38, 8.11, 6.63, $
            5.23, 4.07, 3.20, 2.58], logt)*1d-1
          p2p1 = spline(tgrid, [3.32, 3.24, 2.98, 2.59, 2.14, $
            1.71, 1.35, 1.06, 0.852], logt)*1d0
          pus0 = spline(tgrid, [2.73, 2.79, 2.85, 2.94, 3.06, $
            3.24, 3.49, 3.82, 4.21], logt)*1d0
          pupx = spline(tgrid, [15.8, 15.1, 13.8, 11.5, 9.54, $
            7.77, 6.26, 4.98, 3.89], logt)*1d-1/9.
;
;           Note that for typical nebular temperatures the collision
;           strengths must be extrapolated rather than interpolated.
;
       END
;
;
       'F_II' : BEGIN
;
;  Set ion type, energy levels and Einstein A-values.  Experimental
;  energy levels cited in Butler and Zeippen 1994.
;
          iontype = 'p4'
          ee = invcm*[0, 341.0, 489.9, 20873.4, 44918.1]
;                       ^   ^      ^       ^        ^
;                      3P2 3P1    3P0     1D2      1S0
;
;  A-values calculated in Galavis et al. 1997.
;
          einas = dblarr(5,5)
          einas[1,0] = 8.912d-4
          einas[2,0:1] = [1.985d-9, 1.784d-4]
          einas[3,0:2] = [3.941d-2, 1.249d-2, 2.618d-6]
          einas[4,0:3] = [1.257d-3, 4.690d-1, 0.00d0, 1.828d0]
;
;  Collision strengths, interpolated for temperature 'logt'.
;
          tgrid = 3.0d0 + findgen(11)/5.0d0
          p1p2 = spline(tgrid, [7.23, 6.90, 6.53, 6.20, 6.02, $
            6.00, 6.05, 6.12, 6.27, 6.61, 7.21], logt)*1d-1
          p0p2 = spline(tgrid, [1.74, 1.66, 1.58, 1.51, 1.50, $
            1.52, 1.56, 1.60, 1.65, 1.75, 1.91], logt)*1d-1
          p0p1 = spline(tgrid, [2.65, 2.53, 2.39, 2.24, 2.13, $
            2.06, 2.02, 2.01, 2.04, 2.14, 2.33], logt)*1d-1
          d2px = spline(tgrid, [1.387, 1.339, 1.281, 1.236, 1.215, $
            1.215, 1.230, 1.258, 1.303, 1.365, 1.442], logt)*1d0/9.
          s0px = spline(tgrid, [1.32, 1.32, 1.32, 1.32, 1.33, $
            1.35, 1.38, 1.42, 1.47, 1.56, 1.68], logt)*1d-1/9.
          s0d2 = spline(tgrid, [3.63, 3.62, 3.61, 3.58, 3.54, $
            3.48, 3.42, 3.38, 3.39, 3.47, 3.58], logt)*1d-1
;
       END
;
;
       'F_IV' : BEGIN
;
;  Set ion type, energy levels and Einstein A-values.  Experimental energy
;  levels cited in Galavis et al. 1997.  The absence of an energy level
;  and A-values for the 5S2 level requires that its energy be set
;  sufficiently high to effectively remove it from the calculation.  
;
          iontype = 'p2'
          ee = invcm*[0, 226, 615, 25236, 53538, 999999.9]
;                       ^   ^    ^     ^      ^       ^
;                      3P0 3P1  3P2   1D2    1S0     5S2
;
;  A-values calculated in Galavis et al. 1997.  CLOUDY does not specify
;  its source for A-values.
;

          einas = dblarr(6,6)
          einas[1,0] = 2.075d-4
          einas[2,0:1] = [4.727d-10, 7.937d-4]
          einas[3,0:2] = [6.124d-6, 3.398d-2, 9.731d-2]
          einas[4,1:3] = [1.092d0, 2.129d-3, 2.114d0]
;           einas[5,1:3] = [8.13d2, 2.03d3, 5.51d-2]
          einas[5,1:3] = [0.00d0, 0.00d0, 0.00d0]

;  Collision strengths, interpolated for temperature 'logt'.  CS from
;  Lennon and Burke 1994.

          tgrid = 3.0d0 + findgen(11)/5.0d0
          p1p0 = spline(tgrid, [3.687, 4.047, 4.518, 4.959, 5.268, $
            5.536, 5.962, 6.494, 6.920, 7.108, 6.977], logt)*1d-1
          p2p0 = spline(tgrid, [1.413, 1.596, 1.895, 2.277, 2.630, $
            2.922, 3.251, 3.677, 4.128, 4.455, 4.508], logt)*1d-1
          p2p1 = spline(tgrid, [7.778, 8.638, 9.901, 11.31, 12.49, $
            13.49, 14.76, 16.38, 17.92, 18.86, 18.79], logt)*1d-1
          d2px = spline(tgrid, [2.5218, 2.5073, 2.5559, 2.6392, 2.7069, $
            2.7338, 2.7201, 2.6763, 2.6139, 2.5344, 2.4159], logt)*1d0/9.
          s0px = spline(tgrid, [2.907, 2.956, 3.000, 3.078, 3.189, $
            3.296, 3.357, 3.357, 3.312, 3.228, 3.094], logt)*1d-1/9.
          s0d2 = spline(tgrid, [2.977, 3.024, 3.054, 3.081, 3.128, $
            3.245, 3.461, 3.782, 4.240, 4.745, 5.050], logt)*1d-1
          s2px = spline(tgrid, [9.493, 9.583, 9.315, 9.193, 9.430, $
            9.786, 10.135, 10.440, 10.517, 10.136, 9.277], logt)*1d0/9.
;
       END
;
;
       'NE_III' : BEGIN

          if keyword_set(cloudy94) then begin ; micromanaged to match CLOUDY 94 results

             iontype = '3lev'
             ee = ktoev*[0, 35830., 78840.]
             einas = dblarr(3,3)
             einas[1,0] = 2.09d-1
             einas[2,0] = 2.052d0
             einas[2,1] = 2.677d0
             G0 = 9.00d0
             G1 = 5.00d0
             G2 = 1.00d0
             O10 = 1.348d0
             O20 = 1.51d-1
             O21 = 2.69d-1
             if (logt gt 4) then O21 = (0.1142d0*tt^0.093) < 3.33d-1
;

          endif else begin
             
;  Set ion type, energy levels and Einstein A-values.  Experimental
;  energy levels cited in Butler and Zeippen 1994.
             
             iontype = 'p4'
             ee = invcm*[0, 642.9, 920.4, 25840.8, 55750.6]
;                        ^   ^      ^       ^        ^
;                       3P2 3P1    3P0     1D2      1S0
             
             einas = dblarr(5,5)
             einas[1,0] = 5.974d-3
             einas[2,0:1] = [2.081d-8, 1.159d-3]
             einas[3,0:2] = [1.730d-1, 5.344d-2, 8.269d-6]
             einas[4,0:3] = [3.985d-3, 2.028d0, 0.00d0, 2.563d0]
             
;  Collision strengths, interpolated for temperature 'logt'.  CS from
;  Butler and Zeippen 1994.
             
             tgrid = 3.0d0 + findgen(11)/5.0d0
             p1p2 = spline(tgrid, [4.81, 5.45, 6.34, 7.08, 7.52, $
               7.74, 7.78, 7.71, 7.64, 7.73, 7.94], logt)*1d-1
             p0p2 = spline(tgrid, [1.28, 1.49, 1.74, 1.94, 2.04, $
               2.08, 2.08, 2.05, 2.03, 2.05, 2.11], logt)*1d-1
             p0p1 = spline(tgrid, [1.54, 1.68, 1.94, 2.18, 2.35, $
               2.44, 2.47, 2.47, 2.46, 2.49, 2.56], logt)*1d-1
             d2px = spline(tgrid, [1.349, 1.377, 1.387, 1.380, 1.368, $
               1.357, 1.348, 1.343, 1.348, 1.374, 1.407], logt)*1d0/9.
             s0px = spline(tgrid, [1.50, 1.50, 1.50, 1.50, 1.51, $
               1.51, 1.51, 1.54, 1.62, 1.71, 1.80], logt)*1d-1/9.
             s0d2 = spline(tgrid, [2.66, 2.66, 2.66, 2.66, 2.67, $
               2.69, 2.77, 2.92, 3.10, 3.25, 3.33], logt)*1d-1
             d2px = 1.348d0/9.
             s0px = 0.151d0/9.
             if (logt lt 4.0) then s0d2 = 0.269d0  $
             else s0d2 = 0.1142d0*tt^0.993 < 0.333d0

          endelse
          
       END
;
;
       'NE_IV' : BEGIN
;
;  Set ion type, energy levels and Einstein A-values.  Experimental
;  energy levels cited in Zeippen 1982.
;
          iontype = 'p3c'
          ee = ryd*[0, 0.375767, 0.376177, 0.568961, 0.569022]
;                     ^      ^         ^         ^         ^
;                    4S3/2  2D5/2     2D3/2     2P1/2     2P3/2
;
;  A-values (=A(E2)+A(M1)) from Zeippen 1982.
;
          einas = dblarr(5,5)
          einas[1,0] = 4.84d-4
          einas[2,0:1] = [5.54d-3, 1.48d-6]
          einas[3,0:2] = [5.21d-1, 1.15d-1, 3.93d-1]
          einas[4,0:3] = [1.27d0, 4.00d-1, 4.37d-1, 2.68d-9]
;
;  Collision strengths, interpolated for temperature 'logt'.  CS from
;  Giles 1981.
;
          tgrid = alog10(6.0d3 + 2.0d3*findgen(8))
          dxs3 = spline(tgrid, [1.405, 1.400, 1.397, 1.393, 1.390, $
            1.385, 1.380, 1.374], logt)*1d0/10.
          pxs3 = spline(tgrid, [4.620, 4.670, 4.690, 4.690, 4.680, $
            4.680, 4.660, 4.640], logt)*1d-1/6.
          d3d5 = spline(tgrid, [1.365, 1.362, 1.362, 1.359, 1.355, $
            1.348, 1.340, 1.331], logt)*1d0
          p1d5 = spline(tgrid, [3.47, 3.61, 3.68, 3.71, 3.73, $
            3.74, 3.74, 3.74], logt)*1d-1
          p3d5 = spline(tgrid, [8.67, 8.89, 9.00, 9.05, 9.08, $
            9.09, 9.10, 9.09], logt)*1d-1
          p1d3 = spline(tgrid, [3.27, 3.34, 3.36, 3.38, 3.39, $
            3.39, 3.39, 3.39], logt)*1d-1
          p3d3 = spline(tgrid, [4.82, 5.00, 5.09, 5.13, 5.15, $
            5.16, 5.16, 5.16], logt)*1d-1
          p3p1 = spline(tgrid, [3.23, 3.35, 3.43, 3.50, 3.55, $
            3.61, 3.66, 3.70], logt)*1d-1
;
       END
;
;
       'NE_V' : BEGIN
;
;  Set ion type, energy levels and Einstein A-values.  Experimental
;  energy levels cited in Galavis et al. 1997, except 5S2 inferred
;  from reference wavelength in Pradhan and Peng 1992.
;
          iontype = 'p2'
          ee = invcm*[0, 413, 1111, 30291, 63915, 88364]
;                       ^   ^     ^     ^      ^      ^
;                      3P0 3P1   3P2   1D2    1S0    5S2
;
;  A-values from Baluja 1985, except 5S2-3P? from Mendoza 1982 and
;  5S2-1S0 from Bhatia and Dorshek 1993.
;
          einas = dblarr(6,6)
          einas[1,0] = 1.19d-3
          einas[2,0:1] = [4.59d-9, 4.33d-3]
          einas[3,0:2] = [2.74d-5, 1.37d-1, 3.84d-1]
          einas[4,1:3] = [4.34d0, 6.41d-3, 2.76d0]
          einas[5,1:3] = [2.37d3, 6.06d3, 1.973d-1]
;
;  Collision strengths, interpolated for temperature 'logt'.  CS from
;  Lennon and Burke 1994.
;
          tgrid = 3.0d0 + findgen(11)/5.0d0
          p1p0 = spline(tgrid, [18.30, 18.40, 18.20, 17.50, 16.08, $
            14.08, 11.92, 10.06, 8.752, 7.879, 7.157], logt)*1d-1
          p2p0 = spline(tgrid, [32.28, 31.95, 29.82, 26.43, 22.35, $
            18.10, 14.16, 10.96, 8.709, 7.205, 6.098], logt)*1d-1
          p2p1 = spline(tgrid, [9.551, 9.488, 8.984, 8.133, 7.038, $
            5.832, 4.675, 3.722, 3.049, 2.597, 2.254], logt)*1d0
          d2px = spline(tgrid, [1.7871, 2.0001, 2.1171, 2.1442, 2.1153, $
            2.0863, 2.1118, 2.1761, 2.2105, 2.1748, 2.0668], logt)*1d0/9.
          s0px = spline(tgrid, [3.155, 2.967, 2.765, 2.597, 2.491, $
            2.457, 2.486, 2.538, 2.564, 2.538, 2.452], logt)*1d-1/9.
          s0d2 = spline(tgrid, [5.289, 6.767, 7.354, 7.009, 6.251, $
            5.774, 6.103, 6.878, 7.299, 7.148, 6.621], logt)*1d-1
          s2px = spline(tgrid, [2.806, 3.435, 5.658, 9.418, 12.803, $
            14.258, 13.933, 12.871, 11.764, 10.638, 9.343], logt)*1d-1/9.
;
       END
;
;
       'NA_IV' : BEGIN
;
;  Set ion type, energy levels and Einstein A-values.  Experimental
;  energy values cited in Galavis et al. 1997.
;
          iontype = 'p4'
          ee = invcm*[0, 1107, 1576, 30841, 66496]
;                       ^    ^     ^     ^      ^
;                      3P2  3P1   3P0   1D2    1S0
;
;  A-values from Galavis et al. 1997.
;
          einas = dblarr(5,5)
          einas[1,0] = 3.048d-2
          einas[2,0:1] = [1.599d-7, 5.562d-3]
          einas[3,0:2] = [6.149d-1, 1.837d-1, 2.202d-5]
          einas[4,0:3] = [1.066d-2, 7.152d0, 0.00d0, 3.315d0]
;
;  Collision strengths, interpolated for temperature 'logt'.  CS from
;  Butler and Zeippen 1994.
;
          tgrid = 3.0d0 + findgen(11)/5.0d0
          p1p2 = spline(tgrid, [6.99, 7.51, 7.82, 7.90, 7.91, $
            8.02, 8.27, 8.51, 8.61, 8.78, 9.19], logt)*1d-1
          p0p2 = spline(tgrid, [1.78, 1.94, 2.02, 2.03, 2.02, $
            2.05, 2.13, 2.21, 2.25, 2.33, 2.50], logt)*1d-1
          p0p1 = spline(tgrid, [2.38, 2.51, 2.61, 2.66, 2.69, $
            2.73, 2.79, 2.83, 2.83, 2.83, 2.86], logt)*1d-1
          d2px = spline(tgrid, [1.096, 1.105, 1.104, 1.086, 1.071, $
            1.098, 1.168, 1.251, 1.320, 1.361, 1.374], logt)*1d0/9.
          s0px = spline(tgrid, [1.82, 1.81, 1.80, 1.79, 1.78, $
            1.78, 1.80, 1.82, 1.83, 1.83, 1.83], logt)*1d-1/9.
          s0d2 = spline(tgrid, [2.08, 2.08, 2.08, 2.08, 2.08, $
            2.11, 2.18, 2.27, 2.39, 2.54, 2.71], logt)*1d-1
;
       END
;
;
       'NA_V' : BEGIN
;
;  Set ion type, energy levels and Einstein A-values.  Experimental
;  energy levels cited in Becker et al. 1989.
;
          iontype = 'p3c'
          ee = ryd*[0, 0.440275, 0.440693, 0.667081, 0.667395]
;                     ^      ^         ^         ^         ^
;                    4S3/2  2D5/2     2D3/2     2P1/2     2P3/2
;
;  A-values from Kaufamn and Sugar 1986, as cited by Pradhan and Peng 1992.
;
          einas = dblarr(5,5)
          einas[1,0] = 1.73d-3
          einas[2,0:1] = [1.78d-2, 7.50d-7]
          einas[3,0:2] = [1.96d0, 1.41d-1, 1.43d0]
          einas[4,0:3] = [4.74d0, 1.40d0, 1.91d0, 4.55d-7]
;
;  Collision strengths, interpolated for temperature 'logt'.  CS from
;  Mendoza 1982.
;
          tgrid = alog10(5.0d3 + 5.0d3*findgen(4))
          dxs3 = spline(tgrid, [9.19, 9.19, 9.19, 9.19], logt)*1d-1/10.
          pxs3 = spline(tgrid, [3.59, 3.59, 3.59, 3.59], logt)*1d-1/6.
          d3d5 = spline(tgrid, [6.96, 6.96, 6.96, 6.96], logt)*1d-1
          p1d5 = spline(tgrid, [2.01, 2.01, 2.01, 2.01], logt)*1d-1
          p3d5 = spline(tgrid, [5.02, 5.02, 5.02, 5.02], logt)*1d-1
          p1d3 = spline(tgrid, [1.90, 1.90, 1.90, 1.90], logt)*1d-1
          p3d3 = spline(tgrid, [2.79, 2.79, 2.79, 2.79], logt)*1d-1
          p3p1 = spline(tgrid, [4.38, 4.38, 4.38, 4.38], logt)*1d-1
;
;           Collision strengths given only for 10000K.
;
       END
;
;
       'NA_VI' : BEGIN
;
;  Set ion type, energy levels and Einstein A-values.  Experimental energy
;  levels cited in Galavis et al. 1997.  The absence of an energy level
;  and A-values for the 5S2 level requires that its energy be set
;  sufficiently high to effectively remove it from the calculation.  
;
;
          iontype = 'p2'
          ee = invcm*[0, 698, 1858, 35506, 74423, 999999.9]
;                       ^   ^     ^     ^      ^       ^
;                      3P0 3P1   3P2   1D2    1S0     5S2
;
;  A-values calculated in Galavis et al. 1997.  CLOUDY does not specify
;  its source for A-values.
;
          einas = dblarr(6,6)
          einas[1,0] = 6.111d-3
          einas[2,0:1] = [3.788d-8, 2.104d-2]
          einas[3,0:2] = [5.007d-5, 4.131d-1, 1.119d0]
          einas[4,1:3] = [1.284d1, 1.589d-2, 3.389d0]
          einas[5,1:3] = [0.00d0, 0.00d0, 0.00d0]
;
;  Collision strengths, interpolated for temperature 'logt'.  CS from
;  Lennon and Burke 1994.
;
          tgrid = 3.0d0 + findgen(11)/5.0d0
          p1p0 = spline(tgrid, [7.506, 6.948, 6.803, 7.052, 7.422, $
            7.697, 7.725, 7.428, 6.922, 6.372, 5.825], logt)*1d-1
          p2p0 = spline(tgrid, [5.294, 4.828, 4.702, 4.899, 5.140, $
            5.214, 5.082, 4.807, 4.504, 4.246, 3.995], logt)*1d-1
          p2p1 = spline(tgrid, [2.130, 1.955, 1.908, 1.984, 2.084, $
            2.133, 2.102, 1.995, 1.853, 1.718, 1.589], logt)*1d0
          d2px = spline(tgrid, [1.7356, 1.7075, 1.6532, 1.5870, 1.5182, $
            1.4505, 1.3947, 1.3602, 1.3466, 1.3417, 1.3232], logt)*1d0/9.
          s0px = spline(tgrid, [1.769, 1.757, 1.748, 1.738, 1.727, $
            1.718, 1.720, 1.746, 1.789, 1.823, 1.815], logt)*1d-1/9.
          s0d2 = spline(tgrid, [9.780, 9.950, 10.21, 10.53, 10.94, $
            11.57, 12.76, 14.98, 18.87, 24.10, 28.76], logt)*1d-2
          s2px = spline(tgrid, [2.976, 2.943, 2.925, 3.009, 3.483, $
            4.538, 5.838, 6.886, 7.449, 7.461, 6.954], logt)*1d-1/9.
;
       END
;
;
       'MG_V' : BEGIN
;
;  Set ion type, energy levels and Einstein A-values.  Experimental
;  energy levels cited in Butler and Zeippen 1994.
;
          iontype = 'p4'
          ee = invcm*[0, 1783.1, 2521.8, 35926.0, 77289.0]
;                       ^   ^       ^        ^       ^
;                      3P2 3P1     3P0      1D2     1S0
;
;  A-values calculated in Galavis et al. 1997.  CLOUDY does not specify
;  its source for A-values.
;
          einas = dblarr(5,5)
          einas[1,0] = 1.273d-1
          einas[2,0:1] = [9.683d-7, 2.166d-2]
          einas[3,0:2] = [1.868d0, 5.348d-1, 5.146d-5]
          einas[4,0:3] = [2.499d-2, 2.157d1, 0.00d0, 4.090d0]
;
;  Collision strengths, interpolated for temperature 'logt'.  CS from
;  Butler and Zeippen 1994.
;
          tgrid = 3.0d0 + findgen(11)/5.0d0
          p1p2 = spline(tgrid, [5.01, 5.10, 5.54, 6.37, 7.48, $
            8.52, 9.08, 9.32, 9.74, 10.29, 10.58], logt)*1d-1
          p0p2 = spline(tgrid, [1.16, 1.21, 1.36, 1.63, 1.97, $
            2.25, 2.40, 2.46, 2.62, 2.84, 2.97], logt)*1d-1
          p0p1 = spline(tgrid, [1.92, 1.91, 1.97, 2.15, 2.45, $
            2.76, 2.96, 3.02, 3.07, 3.11, 3.11], logt)*1d-1
          d2px = spline(tgrid, [1.308, 1.296, 1.296, 1.302, 1.306, $
            1.317, 1.335, 1.335, 1.306, 1.262, 1.217], logt)*1d0/9.
          s0px = spline(tgrid, [1.27, 1.36, 1.50, 1.62, 1.65, $
            1.64, 1.59, 1.53, 1.51, 1.51, 1.51], logt)*1d-1/9.
          s0d2 = spline(tgrid, [1.81, 1.82, 1.82, 1.82, 1.81, $
            1.82, 1.86, 1.99, 2.21, 2.51, 2.78], logt)*1d-1
;
       END
;
;
       'MG_VI' : BEGIN
;
;  Set ion type, energy levels and Einstein A-values.  Energy levels
;  derived from wavelengths cited in Kafatos and Lynch 1980.
;
          iontype = 'p3c'
          ee = 1.d0*[0, 6.86327, 6.86555, 10.40492, 10.41803]
;                      ^     ^        ^         ^         ^
;                     4S3/2 2D5/2    2D3/2     2P1/2     2P3/2
;
;
;  A values from Kafatos and Lynch 1980.
;
          einas = dblarr(6,6)
          einas[1,0] = 5.40d-3
          einas[2,0:1] = [1.20d-1, 1.50d-7]
          einas[3,0:2] = [5.30d0, 1.50d-1, 2.50d0]
          einas[4,0:3] = [1.30d1, 2.40d0, 3.80d0, 1.6d-5]
;
;  Collision strengths, interpolated for temperature 'logt'.  CS from
;  Kafatos and Lynch 1980, specified at a single temperature.
;
;
          tgrid = alog10(5.0d3 + 5.0d3*findgen(4))
          dxs3 = spline(tgrid, [6.52, 6.52, 6.52, 6.52], logt)*1d-1/10.
          pxs3 = spline(tgrid, [2.89, 2.89, 2.89, 2.89], logt)*1d-1/6.
          d3d5 = spline(tgrid, [5.60, 5.60, 5.60, 5.60], logt)*1d-1
          p3p1 = spline(tgrid, [2.03, 2.03, 2.03, 2.03], logt)*1d-1
          p3d5 = spline(tgrid, [4.28, 4.28, 4.28, 4.28], logt)*1d-1
          p1d5 = spline(tgrid, [1.67, 1.67, 1.67, 1.67], logt)*1d-1
          p3d3 = spline(tgrid, [2.26, 2.26, 2.26, 2.26], logt)*1d-1
          p1d3 = spline(tgrid, [1.49, 1.49, 1.49, 1.49], logt)*1d-1
;
       END
;
;
       'MG_VII' : BEGIN
;
;  Set ion type, energy levels and Einstein A-values.  Experimental
;  energy levels cited in Galavis et al. 1997.  Energy for 5S2
;  level unknown, so set to arbitrarily large value for calculation.
;
          iontype = 'p2'
          ee = invcm*[0, 1118, 2933, 40957, 85163, 999999.9]
;                       ^    ^     ^     ^      ^       ^
;                      3P0  3P1   3P2   1D2    1S0     5S2
;
;  A-values from Galavis et al. 1997.  CLOUDY does not cite its source
;  for A-values.
;
          einas = dblarr(6,6)
          einas[1,0] = 2.510d-2
          einas[2,0:1] = [2.319d-7, 8.052d-2]
          einas[3,0:2] = [1.163d-4, 1.192d0, 3.107d0]
          einas[4,1:3] = [3.585d1, 3.611d-2, 3.958d0]
          einas[5,1:3] = [0.00d0, 0.00d0, 0.00d0]
;
;  Collision strengths, interpolated for temperature 'logt'.
;
          tgrid = 3.0d0 + findgen(11)/5.0d0
          p1p0 = spline(tgrid, [2.414, 2.565, 2.597, 2.633, 2.858, $
            3.366, 3.952, 4.333, 4.476, 4.488, 4.396], logt)*1d-1
          p2p0 = spline(tgrid, [1.787, 1.795, 1.690, 1.688, 2.115, $
            3.005, 3.879, 4.303, 4.307, 4.132, 3.873], logt)*1d-1
          p2p1 = spline(tgrid, [7.038, 7.245, 7.049, 7.086, 8.298, $
            10.79, 13.24, 14.46, 14.56, 14.20, 13.57], logt)*1d-1
          d2px = spline(tgrid, [7.362, 7.344, 7.496, 7.776, 8.135, $
            8.567, 9.105, 9.739, 10.330, 10.670, 10.651], logt)*1d-1/9.
          s0px = spline(tgrid, [1.998, 2.150, 2.225, 2.162, 2.006, $
            1.849, 1.753, 1.705, 1.670, 1.627, 1.564], logt)*1d-1/9.
          s0d2 = spline(tgrid, [3.201, 4.201, 5.110, 5.417, 5.099, $
            4.461, 3.902, 3.735, 3.859, 3.965, 3.901], logt)*1d-1
          s2px = spline(tgrid, [8.217, 9.880, 10.484, 9.827, 8.426, $
            6.970, 5.973, 5.560, 5.506, 5.474, 5.211], logt)*1d-1/9.
;
       END
;
;
       'AL_VI' : BEGIN
;
;  Set ion type, energy levels and Einstein A-values.  Experimental
;  energy levels cited in Bulter and Zeippen 1994.
;
          iontype = 'p4'
          ee = invcm*[0, 2732, 3829, 41167, 88213]
;                       ^    ^     ^     ^      ^
;                      3P2  3P1   3P0   1D2    1S0
;
;  A-values from Galavis et al. 1997.  CLOUDY does not cite source
;  for A-values.
;
          einas = dblarr(5,5)
          einas[1,0] = 4.580d-1
          einas[2,0:1] = [4.849d-6, 7.051d-2]
          einas[3,0:2] = [5.030d0, 1.366d0, 1.085d-4]
          einas[4,0:3] = [5.292d-2, 5.761d1, 0.00d0, 4.891d0]
;
;  Collision  strengths, interpolated for temperature 'logt'.  CS from
;  Butler and Zeippen 1994.  CLOUDY does not cite source for CS.
;
          tgrid = 3.0d0 + findgen(11)/5.0d0
          p1p2 = spline(tgrid, [5.713, 5.442, 5.419, 5.475, 5.271, $
            4.673, 3.831, 2.969, 2.242, 1.711, 1.365], logt)*1d0
          p0p2 = spline(tgrid, [2.011, 1.899, 1.856, 1.833, 1.725, $
            1.498, 1.204, 0.914, 0.675, 0.504, 0.395], logt)*1d0
          p0p1 = spline(tgrid, [9.50, 9.36, 9.93, 10.81, 11.11, $
            10.41, 8.97, 7.30, 5.78, 4.62, 3.81], logt)*1d-1
          d2px = spline(tgrid, [11.99, 12.27, 12.16, 11.94, 11.70, $
            11.35, 10.90, 10.44, 10.09, 9.87, 9.66], logt)*1d-1/9.
          s0px = spline(tgrid, [0.890, 0.910, 1.01, 1.17, 1.33, $
            1.41, 1.44, 1.45, 1.44, 1.42, 1.38], logt)*1d-2/9.
          s0d2 = spline(tgrid, [2.58, 3.01, 3.66, 4.18, 4.29, $
            4.17, 4.25, 4.63, 4.87, 4.71, 4.27], logt)*1d-1
;
       END
;
;
       'AL_VIII' : BEGIN
;
;  Set ion type, energy levels and Einstein A-values.  Experimental
;  energy levels cited in Galavis et al. 1997.  Energy for 5S2
;  level unknown, so set to arbitrarily large value for calculation.
;
          iontype = 'p2'
          ee = invcm*[0, 1715, 4419, 46729, 96243, 999999.9]
;                       ^    ^     ^     ^      ^       ^
;                      3P0  3P1   3P2   1D2    1S0     5S2
;
;  A-values from Galavis et al. 1997.  CLOUDY does not cite source
;  for A-values.
;
          einas = dblarr(6,6)
          einas[1,0] = 9.052d-2
          einas[2,0:1] = [1.185d-6, 2.660d-1]
          einas[3,0:2] = [2.495d-4, 3.116d0, 7.746d0]
          einas[4,1:3] = [8.983d1, 7.605d-2, 4.554d0]
          einas[5,1:3] = [0.00d0, 0.00d0, 0.00d0]
;
;  Collision strengths, interpolated for temperature 'logt'.  CS from
;  Lennon and Burke 1994.
;
          tgrid = 3.0d0 + findgen(11)/5.0d0
          p1p0 = spline(tgrid, [1.756, 1.880, 2.063, 2.303, 2.538, $
            2.735, 2.994, 3.365, 3.726, 3.912, 3.887], logt)*1d-1
          p2p0 = spline(tgrid, [0.714, 0.780, 0.869, 0.997, 1.130, $
            1.255, 1.462, 1.839, 2.280, 2.585, 2.679], logt)*1d-1
          p2p1 = spline(tgrid, [3.787, 4.088, 4.503, 5.102, 5.693, $
            6.220, 6.995, 8.269, 9.648, 10.50, 10.62], logt)*1d-1
          d2px = spline(tgrid, [2.3957, 2.3562, 2.2370, 2.0213, 1.7733, $
            1.5566, 1.3854, 1.2519, 1.1437, 1.0532, 0.9751], logt)*1d0/9.
          s0px = spline(tgrid, [0.921, 1.002, 1.043, 1.059, 1.064, $
            1.067, 1.066, 1.070, 1.083, 1.098, 1.100], logt)*1d-1/9.
          s0d2 = spline(tgrid, [9.172, 10.60, 11.13, 10.61, 9.461, $
            8.047, 6.548, 5.142, 4.018, 3.301, 2.923], logt)*1d-1/9.
          s2px = spline(tgrid, [1.552, 1.540, 1.565, 1.773, 2.401, $
            3.367, 4.258, 4.800, 5.025, 5.007, 4.734], logt)*1d-1
;
       END
;
;
       'SI_II' : BEGIN
;
;  Set ion type, energy levels and Einstein A-values.
;
;  Experimental energy levels cited in Dufton et al. 1981 (P-levels) and 
;  inferred from reference wavelengths in Bergeson and Lawler, 1993
;  (D-levels).
;
          iontype = 'p1'
          ee = invcm*[0, 287.0, 42824.0, 42933.0, 43108.0, $
;                       ^     ^     ^        ^        ^
;                      2P1/2 2P3/2 4P1/2    4P3/2    4P5/2
;
            55309.0,    55304.1]
;                             ^           ^
;                            2D5/2       2D3/2
;
; Einstein A-values for P-levels from Calamai et al. 1993.  A-values
; for D-levels inferred from lifetimes in Bergeson and Lawler 1993
; and 2D3/2 branching ratio in Nussbaumer, MNRAS, 58, 291.
;
          einas = dblarr(7,7)
          einas[1,0] = 2.170d-4
          einas[2,0:1] = [5.20d3, 4.41d3]
          einas[3,0:2] = [1.32d1, 1.22d3, 0.00d0]
          einas[4,1:3] = [2.46d3, 0.00d0, 0.00d0]
          einas[5,0:2] = [1.94d6, 3.32d5, 0.00d0]
          einas[6,0:2] = [0.00d0, 2.24d6, 0.00d0]
;
;  Collision strengths, interpolated for temperature 'logt'.  CS from
;  Dufton and Kingston 1991.
;
          tgrid = 3.6d0 + findgen(6)/5.0d0
          pl3pl1 = spline(tgrid, [5.58, 5.61, 5.70, 5.79, $
            5.75, 5.47], logt)*1d0
          pu1pl1 = spline(tgrid, [5.58, 5.42, 5.16, 4.84, $
            4.49, 4.13], logt)*1d-1
          pu1pl3 = spline(tgrid, [4.42, 4.24, 4.02, 3.78, $
            3.52, 3.23], logt)*1d-1
          pu3pl1 = spline(tgrid, [8.44, 8.19, 7.80, 7.31, $
            6.80, 6.26], logt)*1d-1
          pu3pl3 = spline(tgrid, [11.5, 11.1, 10.5, 9.90, $
            9.21, 8.43], logt)*1d-1
          pu3pu1 = spline(tgrid, [5.03, 4.81, 4.51, 4.14, $
            3.73, 3.33], logt)*1d0
          pu5pl1 = spline(tgrid, [5.81, 5.61, 5.34, 5.04, $
            4.71, 4.32], logt)*1d-1
          pu5pl3 = spline(tgrid, [2.35, 2.30, 2.19, 2.06, $
            1.92, 1.76], logt)*1d0
          pu5pu1 = spline(tgrid, [1.68, 1.68, 1.67, 1.63, $
            1.51, 1.32], logt)*1d0
          pu5pu3 = spline(tgrid, [7.46, 7.26, 6.94, 6.53, $
            6.10, 5.69], logt)*1d0
          d5pl1  = spline(tgrid, [2.77, 2.76, 2.74, 2.66, $
            2.49, 2.24], logt)*1d0
          d5pl3  = spline(tgrid, [3.50, 3.49, 3.48, 3.41, $
            3.23, 2.93], logt)*1d0
          d5pu1  = spline(tgrid, [12.0, 12.1, 12.0, 11.4, $
            10.4, 9.27], logt)*1d-1
          d5pu3  = spline(tgrid, [1.84, 1.87, 1.85, 1.76, $
            1.62, 1.44], logt)*1d0
          d5pu5  = spline(tgrid, [1.38, 1.40, 1.40, 1.35, $
            1.27, 1.15], logt)*1d0
          d3pl1  = spline(tgrid, [2.45, 2.45, 2.44, 2.38, $
            2.23, 2.00], logt)*1d0
          d3pl3  = spline(tgrid, [6.90, 6.86, 6.79, 6.55, $
            6.05, 5.35], logt)*1d0
          d3pu1  = spline(tgrid, [6.43, 6.54, 6.53, 6.33, $
            5.93, 5.39], logt)*1d-1
          d3pu3  = spline(tgrid, [1.84, 1.87, 1.86, 1.78, $
            1.65, 1.49], logt)*1d0
          d3pu5  = spline(tgrid, [4.14, 4.20, 4.16, 3.96, $
            3.64, 3.25], logt)*1d0
          d3d5   = spline(tgrid, [6.09, 5.99, 5.92, 5.83, $
            5.67, 5.45], logt)*1d0
;
       END
;
;
       'SI_III' : BEGIN
;
;  Set ion type, energy levels and Einstein A-values.
;
;  Energies and Einstein A's from Nussbaumer 1986 and Pradhan and Peng 1992.
;
          iontype = 's2'
          ee = 1d0*[0, 6.537225, 6.553161, 6.585618, 10.276419]
;                     ^      ^         ^         ^          ^
;                    1S0    3P0       3P1       3P2        1P1U
;
          einas = dblarr(5,5)
          einas[1,0] = 0.000d0
          einas[2,0:1] = [1.800d4, 3.860d-5]
          einas[3,0:2] = [1.300d-2, 3.200d-9, 2.420d-4]
          einas[4,0] = 2.590d9
;
;  Collision strengths, interpolated for temperature 'logt'.  CS from
;  Berrington et al., 1994.  CLOUDY claims to have used Hallaway 1994
;  (same ADNDT, 57, 9 vs. 57, 273) but citation makes no sense.
;
          tgrid = alog10(5.0d3 + 5.0d3*findgen(4))
          pxs0 = spline(tgrid, [6.960, 5.456, 4.816, 4.411], logt)*1d0/9.
          p1p0 = spline(tgrid, [1.778, 1.814, 1.830, 1.832], logt)*1d0
          p2p0 = spline(tgrid, [3.660, 3.624, 3.529, 3.432], logt)*1d0
          p2p1 = spline(tgrid, [10.42, 10.42, 10.24, 10.02], logt)*1d0
          pus0 = spline(tgrid, [5.302, 5.600, 5.927, 6.219], logt)*1d0
          pupx = spline(tgrid, [8.532, 8.295, 8.121, 7.917], logt)*1d0/9.
;
;
       END

       'S_II' : BEGIN

          if keyword_set(cloudy94) then begin ; micromanaged to match CLOUDY 94 results
             
             iontype = '5lev'
             ee = yaff*ktoev*[0, 14851.9, 14883.4, 24851.9, 24572.8]
             einas = dblarr(5,5)
             einas[1,0] = 8.82d-4
             einas[2,0:1] = [2.60d-4, 3.35d-7]
             einas[3,0:2] = [9.06d-2, 1.63d-1, 7.80d-2]
             einas[4,0:3] = [2.25d-1, 1.33d-1, 1.79d-1, 1.03d-6]
             G0 = 4.00d0
             G1 = 4.00d0
             G2 = 6.00d0
             G3 = 2.00d0
             G4 = 4.00d0

;   Collision strengths valid for 3.5 < logt < 5.0

             O20 = ( (8.1458628d0 - 0.5389108d0*logt^1.4486586d0) $
               < 4.77d0) > 2.54d0
             O10 = O20/1.5d0
             O40 = ( ( 2.6106784d0 - 3.2766908d-5*logt^6.5105436d0) $
               < 2.46d0) > 1.45d0
             O30 = O40/2.0d0
             O21 = ( ( -27.497273d0 + 247.27405d0/logt - 429.9142d0/logt/logt) $
               < 8.01d0) > 4.79d0
             O31 = ( ( 8.274085d0 - 2.6223732d0*logt + 0.2502924d0*logt*logt) $
               < 2.14d0) > 1.42d0
             O41 = ( ( 6.690242d0 - 1.061514d0*logt + 0.034535506d0*logt*logt) $
               < 3.38d0) > 2.24d0
             O32 = ( ( 4.2250081d0 - 0.46549935d0*logt $
               - 0.010172139d0*logt*logt) < 2.46d0) > 1.64d0
             O42 = ( (18.335524d0 - 5.1180248d0*logt + 0.44482438d0*logt*logt) $
               < 5.82d0) > 3.87d0
             O43 = ( ( -5.1994665d0 + 49.334586d0/logt - 70.93344d0/logt/logt) $
               < 3.07d0) > 1.85d0

          endif else begin

;  jm07dec05nyu: Energy levels from NIST (http://physics.nist.gov) and
;  Einstein A values from the Chianti project
;  (http://www.arcetri.astro.it/science/chianti/chianti.html), as
;  tabulated in Stasinska 2007 (astro-ph/0704.0348)

           iontype = 'p3d'
           ee = invcm*[0.0, 14852.94, 14884.73, 24524.83, 24571.54] ; NIST
;                       ^        ^         ^         ^         ^
;                     4S3/2    2D3/2     2D5/2     2P1/2     2P3/2

           einas = dblarr(5,5)
           einas[1,0]   = [1.231d-3]                           ; [6732]
           einas[2,0:1] = [3.338d-4,3.452d-7]                  ; [6718,314 mu]
           einas[3,0:2] = [1.076d-1,1.812d-1,7.506d-2]         ; [4077,10339,10373]
           einas[4,0:3] = [2.670d-1,1.644d-1,1.938d-1,9.16d-7] ; [4069,10289,10323,214 mu]

;;;  Experimental energy levels cited in Mendoza and Zeippen 1982;
;;;  A-values from Mendoza and Zeippen 1982 (assuming A = A(E2) +
;;;  A(M1).)
;;
;;           iontype = 'p3d'
;;           ee = ryd*[0, 0.135343, 0.135630, 0.223485, 0.223928]
;;;                     ^      ^         ^         ^         ^
;;;                    4S3/2  2D3/2     2D5/2     2P1/2     2P3/2
;;
;;           einas = dblarr(5,5)
;;           einas[1,0] = 8.82d-4
;;           einas[2,0:1] = [2.60d-4, 3.35d-7]
;;           einas[3,0:2] = [9.06d-2, 1.63d-1, 7.80d-2]
;;           einas[4,0:3] = [2.25d-1, 1.33d-1, 1.79d-1, 1.03d-6]
           
;  Collision strengths, interpolated for temperature 'logt'.  CS from
;  Ramsbottom et al. 1996.

           tgrid = 1.0d0*[3.5, 3.6, 3.7, 3.8, 3.9, 4.0, $
                          4.1, 4.2, 4.4, 4.6, 4.8, 5.0]
           dxs3 = spline(tgrid, [7.95, 7.82, 7.64, 7.42, 7.17, 6.90, $
                   6.63, 6.37, 5.90, 5.45, 4.90, 4.23], logt)*1d0/10.
           pxs3 = spline(tgrid, [3.69, 3.70, 3.69, 3.66, 3.60, 3.52, $
                   3.43, 3.33, 3.13, 2.895, 2.579, 2.17], logt)*1d0/6.
           d5d3 = spline(tgrid, [7.99, 8.01, 7.98, 7.87, 7.69, 7.47, $
                   7.22, 6.95, 6.42, 5.92, 5.39, 4.79], logt)*1d0
           p1d3 = spline(tgrid, [2.14, 2.08, 2.01, 1.94, 1.87, 1.79, $
                   1.72, 1.66, 1.57, 1.51, 1.46, 1.42], logt)*1d0
           p3d3 = spline(tgrid, [3.38, 3.32, 3.25, 3.17, 3.08, 3.00, $
                   2.91, 2.82, 2.68, 2.55, 2.40, 2.24], logt)*1d0
           p1d5 = spline(tgrid, [2.46, 2.42, 2.37, 2.32, 2.26, 2.20, $
                   2.14, 2.08, 1.97, 1.88, 1.76, 1.64], logt)*1d0
           p3d5 = spline(tgrid, [5.82, 5.68, 5.52, 5.35, 5.17, 4.99, $
                   4.81, 4.65, 4.40, 4.21, 4.03, 3.87], logt)*1d0
           p3p1 = spline(tgrid, [3.07, 3.03, 2.98, 2.90, 2.81, 2.71, $
                   2.61, 2.50, 2.31, 2.16, 2.01, 1.85], logt)*1d0

        endelse 

       END

       'S_III' : BEGIN

          if keyword_set(cloudy94) then begin ; micromanaged to match CLOUDY 94 results

             iontype = '3lev'
             ee = ktoev*[0, 1.55d4, 3.83d4]
             einas = dblarr(3,3)
             einas[1,0] = 7.97e-2
             einas[2,0] = 8.07d-1
             einas[2,1] = 2.22d0
             G0 = 9.0d0
             G1 = 5.0d0
             G2 = 1.0d0
             O10 = 7.98d0
             O20 = 1.14d0
             O21 = (8.21d-2*tt^0.3) < 2.05

          endif else begin

; jm07dec09nyu: these A-values seem to be wrong, as observationally
; the 9532/9069 ratio is 2.5-3, but the Chianti project give a ratio
; of 5.98!  so use the Mendoza & Zeippen (1982) values for now
             
;  jm07dec05nyu: Energy levels from NIST (http://physics.nist.gov) and
;  Einstein A values from the Chianti project
;  (http://www.arcetri.astro.it/science/chianti/chianti.html), as
;  tabulated in Stasinska 2007 (astro-ph/0704.0348); the 1713, 1728,
;  2111, and 3798 Einstein A values are from the references cited by
;  B. Moore, commented out below 

;;             iontype = 'p2'
;;             ee = invcm*[0.0, 298.68, 833.06, 11322.70, 27161.00, 58671.92]
;;;                         ^      ^       ^         ^         ^         ^
;;;                        3P0    3P1     3P2       1D2       1S0       5S2
;;             
;;             einas = dblarr(6,6)
;;             einas[1,0]   = [8.429d-4]                   ; [33 mu]
;;             einas[2,0:1] = [5.008d-8,1.877d-3]          ; [12 mu,18 mu]
;;             einas[3,0:2] = [2.570d-4,3.107d-2,1.856d-1] ; [8831,9071,9533]
;;             einas[4,1:3] = [1.016,1.05d-2,3.045]        ; [3722,3798,6313]
;;             einas[5,1:3] = [5.01d3,1.38d4,4.43d0]       ; [1713,1728,2111]
             
;  Experimental energy levels from NIST cited in Kohstall et
;  al. 1998 (5S2 energy level calculated); A-values from Mendoza and
;  Zeippen 1982, except 5S2-?? from Kohstall et al. 1998.
             
             iontype = 'p2'
;            ee = ryd*[0, 0.002708, 0.007586, 0.103157, 0.247532, 0.543985]
             ee = invcm*[0.0, 298.68, 833.06, 11322.70, 27161.00, 58671.92]
;                         ^      ^       ^         ^         ^         ^
;                        3P0    3P1     3P2       1D2       1S0       5S2
             
             einas = dblarr(6,6)
             einas[1,0]   = [4.72d-4]                 ; [33 mu]
             einas[2,0:1] = [4.61d-8,2.07d-3]         ; [12 mu,18 mu]
             einas[3,0:2] = [5.82d-6,2.21d-2,5.76d-2] ; [8831,9071,9533]
             einas[4,1:3] = [7.96d-1,1.05d-2,2.22d0]  ; [3722,3798,6313]
             einas[5,1:3] = [5.01d3,1.38d4,4.43d0]    ; [1713,1728,2111]
             
;  Collision strengths, interpolated for temperature 'logt'.  CS for
;  3P?-3P? and S2-P? from Tayal and Gupta 1999.  CS for 1D2-?? and 1S0-??
;  from Galavis et al. 1995.  CS for 5S2-3P? from Tayal 1997.
             
             tgrid = alog10([5.0d3, 8.0d3, 1.0d4, 2.0d4])
             p1p0 = spline(tgrid, [4.44, 4.17, 3.98, 3.24], logt)*1d0
             p2p0 = spline(tgrid, [1.41, 1.35, 1.31, 1.28], logt)*1d0
             p2p1 = spline(tgrid, [8.72, 8.20, 7.87, 6.88], logt)*1d0
             
             s2px = spline(tgrid, [2.45, 2.74, 2.84, 3.12], logt)*1d0/9.
             
             tgrid = 3.0d0 + findgen(11)/5.0d0
             d2px = spline(tgrid, [6.023, 6.932, 7.689, 8.038, 8.023, $
               7.951, 7.989, 8.006, 7.829, 7.332, 6.478], logt)*1d0/9.
             s0px = spline(tgrid, [1.372, 1.297, 1.214, 1.150, 1.117, $
               1.110, 1.131, 1.193, 1.261, 1.257, 1.146], logt)*1d0/9.
             s0d2 = spline(tgrid, [1.094, 1.037, 0.984, 0.982, 1.093, $
               1.301, 1.540, 1.772, 1.962, 2.054, 2.016], logt)*1d0
             

          endelse
          
       END
;
;
       'CL_II' : BEGIN
;
;  Set ion type, energy levels and Einstein A-values.  Experimental
;  nergy levels cited in Mendoza and Zeippen 1983.
;
          iontype = 'p4'
          ee = ryd*[0, 0.006343, 0.009081, 0.106197, 0.254047]
;                       ^    ^         ^         ^         ^
;                      3P2  3P1       3P0       1D2       1S0
;
;  A-values from Mendoza and Zeippen 1983.  CLOUDY specifies this
;  reference to [CL III] IR lines, presumably a mistake.
;
          einas = dblarr(5,5)
          einas[1,0] = 7.57d-3
          einas[2,0:1] = [4.57d-7, 1.46d-3]
          einas[3,0:2] = [1.04d-1, 2.92d-2, 9.82d-6]
          einas[4,0:3] = [1.97d-2, 1.31d0, 0.00d0, 2.06d0]
;
;  Collision strengths, interpolated for temperature 'logt'.  CS
;  from Mendoza 1982 given for 10000K only.
;
          tgrid = alog10(5.0d3 + 5.0d3*findgen(4))
          p1p2 = spline(tgrid, [2.17, 2.17, 2.17, 2.17], logt)*1d0
          p0p2 = spline(tgrid, [4.43, 4.43, 4.43, 4.43], logt)*1d-1
          p0p1 = spline(tgrid, [9.33, 9.33, 9.33, 9.33], logt)*1d-1
          d2px = spline(tgrid, [3.86, 3.86, 3.86, 3.86], logt)*1d0/9.
          s0px = spline(tgrid, [4.56, 4.56, 4.56, 4.56], logt)*1d-1/9.
          s0d2 = spline(tgrid, [1.15, 1.15, 1.15, 1.15], logt)*1d0
;
;
       END
;
;
       'CL_III' : BEGIN
;
;  Set ion type, energy levels and Einstein A-values.  Energy levels
;  from Mendoza et al. 1982.
;
          iontype = 'p3d'
          ee = invcm*[0.0,  18053.0, 18118.6, 29812.0, 29907.0]
;                        ^      ^        ^        ^        ^
;                       4S3/2  2D3/2    2D5/2   2P1/2     2P3/2
;
;  A-values from Mendoza 1982.
;
          einas = dblarr(5,5)
          einas[1,0] = 4.83d-3
          einas[2,0:1] = [7.04d-4, 3.22d-6]
          einas[3,0:2] = [3.05d-1, 3.03d-1, 1.00d-1]
          einas[4,0:3] = [7.54d-1, 3.23d-1, 3.16d-1, 7.65d-6]
;
;  Collision strengths, interpolated for temperature 'logt'.  CS from
;  Mendoza 1982.
;
          tgrid = 5.0d3 + 5.0d3*findgen(4)
          dxs3 = spline(tgrid, [3.14, 3.14, 3.14, 3.14], logt)*1d0/10.
          pxs3 = spline(tgrid, [1.89, 1.89, 1.89, 1.89], logt)*1d0/6.
          d5d3 = spline(tgrid, [3.19, 3.19, 3.19, 3.19], logt)*1d0
          p1d3 = spline(tgrid, [1.24, 1.24, 1.24, 1.24], logt)*1d0
          p3d3 = spline(tgrid, [1.91, 1.91, 1.91, 1.91], logt)*1d0
          p1d5 = spline(tgrid, [1.38, 1.38, 1.38, 1.38], logt)*1d0
          p3d5 = spline(tgrid, [3.33, 3.33, 3.33, 3.33], logt)*1d0
          p3p1 = spline(tgrid, [1.34, 1.34, 1.34, 1.34], logt)*1d0
;
       END
;
;
       'CL_IV' : BEGIN
;
;  Set ion type, energy levels and Einstein A-values.  Experimental
;  energy values cited in Mendoza 1982.
;
          iontype = 'p2'
          ee = invcm*[0, 492.0 , 1341.9, 13767.6, 32547.8, 999999.9]
;                       ^   ^        ^       ^        ^         ^
;                      3P0 3P1      3P2     1D2      1S0       5S2
;
;  A-values from Mendoza 1982.  CLOUDY does not give source for its
;  values, but appears to agree with Mendoza 1982.
;
          einas = dblarr(6,6)
          einas[1,0] = 2.14d-3
          einas[2,0:1] = [2.70d-7, 8.25d-3]
          einas[3,0:2] = [1.54d-5, 7.23d-2, 1.79d-1]
          einas[4,1:3] = [2.47d0, 2.62d-2, 2.80d0]
          einas[5,1:3] = [0.00d0, 0.00d0, 0.00d0]
;
;  Collision strengths, interpolated for temperature 'logt'.  CS from
;  Galavis et al. 1995.
;
          tgrid = 3.0d0 + findgen(11)/5.0d0
          p1p0 = spline(tgrid, [2.694, 2.672, 2.579, 2.376, 2.099, $
            1.828, 1.645, 1.595, 1.631, 1.658, 1.597], logt)*1d0
          p2p0 = spline(tgrid, [2.040, 2.054, 2.078, 2.046, 1.927, $
            1.753, 1.592, 1.489, 1.436, 1.392, 1.308], logt)*1d0
          p2p1 = spline(tgrid, [7.956, 7.961, 7.899, 7.574, 6.959, $
            6.229, 5.639, 5.344, 5.269, 5.204, 4.938], logt)*1d0
          d2px = spline(tgrid, [9.272, 8.872, 8.272, 7.527, 6.856, $
            6.437, 6.337, 6.489, 6.619, 6.443, 5.867], logt)*1d0/9.
          s0px = spline(tgrid, [1.743, 1.733, 1.718, 1.721, 1.780, $
            1.922, 2.089, 2.160, 2.066, 1.836, 1.534], logt)*1d0/9.
          s0d2 = spline(tgrid, [0.697, 0.717, 0.743, 0.803, 0.960, $
            1.240, 1.578, 1.866, 2.042, 2.090, 2.010], logt)*1d0
          tgrid = alog10([5000., 10000., 15000., 20000.])
          s2px = spline(tgrid, [1.00, 1.00, 1.00, 1.00], logt)*1d-50/9.
;
;           No parameters for 5S2 level given.  Values set arbitrarily
;           to values expected to yield extremely small level population.
;
       END
;
;
       'AR_III' : BEGIN

          if keyword_set(cloudy94) then begin ; micromanaged to match CLOUDY 94 results

             iontype = '3lev'
             ee = ktoev*[0, 19550., 47250.]
             einas = dblarr(3,3)
             einas[1,0] = 3.963d-1
             einas[2,0] = 3.952d0
             einas[2,1] = 2.59d0
             G0 = 9.00d0
             G1 = 5.00d0
             G2 = 1.00d0
             O10 = 4.825d0
             O20 = 8.41d-1
             O21 = (3.296d0/tt^0.108) < 1.30d0

          endif else begin

;  Set ion type, energy levels and Einstein A-values.  Experimental
;  energy levels from Mendoza 1982.
             
             iontype = 'p4'
             ee = invcm*[0, 458.1, 1570.2, 14010.0, 33265.7]
;                        ^    ^      ^       ^        ^
;                       3P2  3P1    3P0     1D2      1S0
             
;  A-values from Mendoza 1982.  
             
             einas = dblarr(5,5)
             einas[1,0] = 3.06d-2
             einas[2,0:1] = [2.37d-6, 5.31d-3]
             einas[3,0:2] = [3.14d-1, 8.23d-2, 2.21d-5]
             einas[4,0:3] = [4.17d-2, 3.952d0, 0.00d0, 2.59d0]
             
;  Collision strengths, interpolated for temperature 'logt'.  CS from
;  Galavis et al. 1995.
             
             tgrid = 3.0d0 + findgen(11)/5.0d0
             p1p2 = spline(tgrid, [3.726, 3.507, 3.342, 3.226, 3.137, $
               3.087, 3.117, 3.209, 3.315, 3.343, 3.164],logt)*1d0
             p0p2 = spline(tgrid, [7.25, 6.96, 6.77, 6.65, 6.60, 6.71, $
               7.15, 7.81, 8.54, 9.06, 8.97],logt)*1d-1
             p0p1 = spline(tgrid, [1.676, 1.552, 1.456, 1.384, 1.322, $
               1.261, 1.207, 1.162, 1.128, 1.089, 1.006],logt)*1d0
             d2px = spline(tgrid, [5.092, 5.055, 5.030, 4.988, 4.916, $
               4.825, 4.738, 4.687, 4.663, 4.534, 4.153],logt)*1d0/9.
             s0px = spline(tgrid, [9.12, 9.27, 9.13, 8.84, 8.58, 8.41, $
               8.28, 8.15, 7.97, 7.57, 6.77],logt)*1d-1/9.
             s0d2 = spline(tgrid, [1.587, 1.506, 1.395, 1.302, 1.248, $
               1.219, 1.188, 1.145, 1.096, 1.038, 0.951],logt)*1d0

          endelse
          
       END
;
;
       'AR_IV' : BEGIN
;
;  Set ion type, energy levels and Einstein A-values.  Experimental
;  energy values cited in Mendoza 1982.
;
          iontype = 'p3d'
          ee = invcm*[0.0, 21090.4, 21219.3, 34855.5, 35032.6]
;                        ^     ^        ^        ^        ^
;                       4S3/2 2D3/2    2D5/2    2P1/2    2P3/2
;
;  A-values from Mendoza and Zeippen 1982.  CLOUDY does not cite source,
;  but uses identical values.
;
          einas = dblarr(5,5)
          einas[1,0] = 2.23d-2
          einas[2,0:1] = [1.77d-3, 2.30d-5]
          einas[3,0:2] = [8.62d-1, 6.03d-1, 1.19d-1]
          einas[4,0:3] = [2.11d0, 7.89d-1, 5.98d-1, 4.94d-5]
;
;
;  Collision strengths, interpolated for temperature 'logt'.  CS from
;  Ramsbottom et al. 1997.
;
          tgrid = alog10([5.0d3, 8.0d3, 1.0d4, 1.2d4, 1.5d4, 1.7d4, 2.0d4])
          dxs3 = spline(tgrid, [1.913, 1.901, 1.902, 1.914, $
            1.918, 1.931, 1.947], logt)*1d0/10.
          pxs3 = spline(tgrid, [12.57, 11.70, 11.78, 12.00, $
            12.39, 12.63, 12.98], logt)*1d-1/6.
          d5d3 = spline(tgrid, [7.10, 7.09, 7.06, 6.99, 6.87, $
            6.79, 6.65], logt)*1d0
          p1d3 = spline(tgrid, [1.55, 1.52, 1.51, 1.51, 1.51, $
            1.52, 1.53], logt)*1d0
          p3d3 = spline(tgrid, [2.12, 2.13, 2.14, 2.15, 2.16, $
            2.17, 2.18], logt)*1d0
          p1d5 = spline(tgrid, [1.51, 1.53, 1.53, 1.54, 1.55, $
            1.55, 1.56], logt)*1d0
          p3d5 = spline(tgrid, [4.00, 3.95, 3.94, 3.94, 3.96, $
            3.98, 4.01], logt)*1d0
          p3p1 = spline(tgrid, [1.46, 1.84, 2.07, 2.25, 2.47, $
            2.58, 2.70], logt)*1d0
;
       END
;
;
       'AR_V' : BEGIN
;
;  Set ion type, energy levels and Einstein A-values.  Experimental
;  energy levels from Mendoza 1982.  No data for 5S2, so value is set
;  arbitrarily high to remove level from solution.
;
          iontype = 'p2'
          ee = invcm*[0.0, 763.9, 2029.2, 16299.4, 37912.5, 999999.9]
;                        ^    ^       ^       ^        ^         ^
;                       3P0  3P1     3P2     1D2      1S0       5S2
;
;  A-values from Mendoza 1982.
;
          einas = dblarr(6,6)
          einas[1,0] = 7.99d-3
          einas[2,0:1] = [1.24d-6, 2.72d-2]
          einas[3,0:2] = [3.50d-5, 2.04d-1, 4.76d-1]
          einas[4,1:3] = [6.55d0, 5.69d-2, 3.29d0]
          einas[5,1:3] = [0.00d0, 0.00d0, 0.00d0]
;
;  Collision strengths, interpolated for temperature 'logt'.
;
          tgrid = 3.0d0 + findgen(11)/5.0d0
          p1p0 = spline(tgrid, [4.162, 3.959, 3.747, 3.519, 3.257, $
            2.941, 2.582, 2.240, 1.981, 1.816, 1.682], logt)*1d0
          p2p0 = spline(tgrid, [1.869, 1.871, 1.931, 1.979, 1.953, $
            1.837, 1.672, 1.514, 1.407, 1.349, 1.285], logt)*1d0
          p2p1 = spline(tgrid, [9.409, 9.159, 9.030, 8.852, 8.465, $
            7.811, 6.989, 6.207, 5.642, 5.305, 4.994], logt)*1d0
          d2px = spline(tgrid, [5.026, 4.142, 3.697, 3.456, 3.281, $
            3.207, 3.319, 3.622, 4.043, 4.401, 4.454], logt)*1d0/9.
          s0px = spline(tgrid, [4.02, 4.66, 5.16, 5.40, 5.50, 5.59, $
            5.85, 6.37, 6.98, 7.33, 7.10], logt)*1d-1/9.
          s0d2 = spline(tgrid, [2.319, 2.030, 1.847, 1.745, 1.685, $
            1.648, 1.649, 1.688, 1.742, 1.775, 1.741], logt)*1d0
          tgrid = alog10([5000., 10000., 15000., 20000.])
          s2px = spline(tgrid, [1.00, 1.00, 1.00, 1.00], logt)*1d-50/9.
;
       END
;
;
       'K_IV' : BEGIN
;
;  Set ion type, energy levels and Einstein A-values.  Experimental
;  energy levels cited in Mendoza 1982.
;
          iontype = 'p4'
          ee = invcm*[0.0, 1671.4, 2321.2, 16384.1, 38546.3]
;                       ^      ^       ^       ^        ^
;                      3P2    3P1     3P0     1D2      1S0
;
;  A-values from Mendoza 1982.
;
          einas = dblarr(5,5)
          einas[1,0] = 1.04d-1
          einas[2,0:1] = [1.01d-5, 1.48d-2]
          einas[3,0:2] = [8.14d-1, 1.98d-1, 4.54d-5]
          einas[4,0:3] = [8.17d-2, 1.00d1, 0.00d0, 3.18d0]
;
;  Collision strengths, interpolated for temperature 'logt'.  CS from
;  Galavis et al. 1997.
;
          tgrid = 3.0d0 + findgen(11)/5.0d0
          p1p2 = spline(tgrid, [3.045, 3.568, 3.976, 4.173, 4.181, $
            4.124, 4.170, 4.333, 4.490, 4.491, 4.200], logt)*1d0
          p0p2 = spline(tgrid, [0.660, 0.818, 0.982, 1.112, 1.185, $
            1.219, 1.258, 1.304, 1.330, 1.308, 1.208], logt)*1d0
          p0p1 = spline(tgrid, [1.248, 1.382, 1.413, 1.338, 1.212, $
            1.105, 1.071, 1.119, 1.199, 1.239, 1.186], logt)*1d0
          d2px = spline(tgrid, [4.719, 4.988, 5.310, 5.652, 5.909, $
            6.013, 6.065, 6.139, 6.123, 5.842, 5.197], logt)*1d0/9.
          s0px = spline(tgrid, [1.718, 2.132, 2.572, 2.792, 2.673, $
            2.311, 1.885, 1.523, 1.254, 1.046, 0.854], logt)*1d0/9.
          s0d2 = spline(tgrid, [0.768, 0.791, 0.790, 0.769, 0.760, $
            0.789, 0.875, 1.032, 1.225, 1.348, 1.314], logt)*1d0
;
       END
;
;
       'CA_V' : BEGIN
;
;  Set ion type, energy levels and Einstein A-values.  Experimental
;  energy levels cited in Mendoza 1982.
;
          iontype = 'p4'
          ee = invcm*[0.0, 2404.7, 3275.6, 18830.3, 43836.5]
;                        ^     ^       ^       ^        ^
;                       3P2   3P1     3P0     1D2      1S0
;
;  A-values from Mendoza 1982.
;
          einas = dblarr(5,5)
          einas[1,0] = 3.10d-1
          einas[2,0:1] = [3.67d-5, 3.54d-2]
          einas[3,0:2] = [1.90d0, 4.26d-1, 8.42d-5]
          einas[4,0:3] = [1.45d-1, 2.31d1, 0.00d0, 3.73d0]
;
;  Collision strengths, interpolated for temperature 'logt'.  CS from
;  Galavis et al. 1997.
;
          tgrid = 3.0d0 + findgen(11)/5.0d0
          p1p2 = spline(tgrid, [2.748, 2.436, 2.251, 2.196, 2.220, $
            2.298, 2.440, 2.681, 3.001, 3.278, 3.348], logt)*1d0
          p0p2 = spline(tgrid, [0.668, 0.603, 0.579, 0.590, 0.615, $
            0.648, 0.704, 0.789, 0.887, 0.962, 0.972], logt)*1d0
          p0p1 = spline(tgrid, [0.996, 0.863, 0.758, 0.695, 0.670, $
            0.671, 0.685, 0.725, 0.804, 0.891, 0.928], logt)*1d0
          d2px = spline(tgrid, [3.135, 2.868, 2.712, 2.741, 2.869, $
            3.067, 3.346, 3.691, 3.987, 4.066, 3.829], logt)*1d0
          s0px = spline(tgrid, [0.147, 0.165, 0.199, 0.270, 0.379, $
            0.522, 0.690, 0.825, 0.870, 0.823, 0.713], logt)*1d0
          s0d2 = spline(tgrid, [1.090, 1.198, 1.288, 1.339, 1.353, $
            1.348, 1.354, 1.379, 1.402, 1.381, 1.275], logt)*1d0
;
       END
;
;
       ELSE : BEGIN
;
;  To test the matrix manipulations.
;
          print, ionstr, ' not in database'
          iontype = 'test'
          ee = findgen(6)
;
          einas = dblarr(6,6)
          einas[1,0] = 1.00d0
          einas[2,0:1] = 1.00d0 + dblarr(2)
          einas[3,0:2] = 1.00d0 + dblarr(3)
          einas[4,0:3] = 1.00d0 + dblarr(4)
          einas[5,0:4] = 1.00d0 + dblarr(5)
;
          tgrid = alog10(5.0d3 + 5.0d3*findgen(4))
          p1p0 = spline(tgrid, [1.00, 1.00, 1.00, 1.00], logt)*1d0
          p2p0 = spline(tgrid, [1.00, 1.00, 1.00, 1.00], logt)*1d0
          p2p1 = spline(tgrid, [1.00, 1.00, 1.00, 1.00], logt)*1d0
          d2px = spline(tgrid, [1.00, 1.00, 1.00, 1.00], logt)*1d0
          s0px = spline(tgrid, [1.00, 1.00, 1.00, 1.00], logt)*1d0
          s0d2 = spline(tgrid, [1.00, 1.00, 1.00, 1.00], logt)*1d0
          s2px = spline(tgrid, [1.00, 1.00, 1.00, 1.00], logt)*1d0
          s2d2 = spline(tgrid, [1.00, 1.00, 1.00, 1.00], logt)*1d0
          s2s0 = spline(tgrid, [1.00, 1.00, 1.00, 1.00], logt)*1d0
;
       END
;
    endcase
;
    case strupcase(iontype) of
;
;
       'S2':     BEGIN          ;     Beryllium-like ions
          ombyw = reform($
            [[ 0.00 , 1*pxs0, 3*pxs0, 5*pxs0,  pus0 ]/1. , $
            [1*pxs0,  0.00 ,  p1p0 ,  p2p0 , 1*pupx]/1. , $
            [3*pxs0,  p1p0 ,  0.00 ,  p2p1 , 3*pupx]/3. , $
            [5*pxs0,  p2p0 ,  p2p1 ,  0.00 , 5*pupx]/5. , $
            [ pus0 , 1*pupx, 3*pupx, 5*pupx,  0.00 ]/3.], 5, 5)
;                  ^       ^       ^       ^       ^     ^     ^
;                 1S0     3P0     3P1     3P2     1P1U   wi
;
          en = reform([ (ee-ee(0) > 0), (ee-ee(1) > 0), (ee-ee(2) > 0), $
            (ee-ee(3) > 0), (ee-ee(4) > 0)], 5, 5)
       END
;
       'P1':     BEGIN          ;     boron-like ions -- CII config
          ombyw = reform($
            [[ 0.00 ,pl3pl1,pu1pl1,pu3pl1,pu5pl1, d5pl1, d3pl1 ]/2. , $
            [pl3pl1, 0.00 ,pu1pl3,pu3pl3,pu5pl3, d5pl3, d3pl3 ]/4. , $
            [pu1pl1,pu1pl3, 0.00 ,pu3pu1,pu5pu1, d5pu1, d3pu1 ]/2. , $
            [pu3pl1,pu3pl3,pu3pu1, 0.00 ,pu5pu3, d5pu3, d3pu3 ]/4. , $
            [pu5pl1,pu5pl3,pu5pu1,pu5pu3, 0.00 , d5pu5, d3pu5 ]/6. , $
            [ d5pl1, d5pl3, d5pu1, d5pu3, d5pu5, 0.00 , d3d5  ]/6. , $
            [ d3pl1, d3pl3, d3pu1, d3pu3, d3pu5, d3d5 , 0.00  ]/4.], 7, 7)
;                  ^      ^      ^      ^      ^      ^      ^     ^
;                2P1/2  2P3/2  4P1/2  4P3/2  4P5/2  2D5/2  2D3/2   wi
;
          en = reform([ (ee-ee(0) > 0), (ee-ee(1) > 0), (ee-ee(2) > 0), $
            (ee-ee(3) > 0), (ee-ee(4) > 0), (ee-ee(5) > 0), $
            (ee-ee(6) > 0)], 7, 7)
       END
;
;
       'P2':     BEGIN
          ombyw = reform($
            [[ 0.00 ,  p1p0 ,  p2p0 , 1*d2px, 1*s0px, 1*s2px]/1. , $
            [ p1p0 ,  0.00 ,  p2p1 , 3*d2px, 3*s0px, 3*s2px]/3. , $
            [ p2p0 ,  p2p1 ,  0.00 , 5*d2px, 5*s0px, 5*s2px]/5. , $
            [1*d2px, 3*d2px, 5*d2px,  0.00 ,  s0d2 ,  0.00 ]/5. , $
            [1*s0px, 3*s0px, 5*s0px,  s0d2 ,  0.00 ,  0.00 ]/1. , $
            [1*s2px, 3*s2px, 5*s2px,  0.00 ,  0.00 ,  0.00 ]/5.], 6, 6)
;                  ^       ^       ^       ^       ^       ^     ^
;                 3P0     3P1     3P2     1D2     1S0     5S2    wi
;
          en = reform([ (ee-ee(0) > 0), (ee-ee(1) > 0), (ee-ee(2) > 0), $
            (ee-ee(3) > 0), (ee-ee(4) > 0), (ee-ee(5) > 0)], 6, 6)
       END
;
       'P3A':    BEGIN          ;     Order of energy levels is OII-like
          ombyw = reform($
            [[ 0.00 , 6*dxs3, 4*dxs3, 4*pxs3, 2*pxs3]/4. ,$
            [6*dxs3,  0.00 ,  d3d5 ,  p3d5 ,  p1d5 ]/6. ,$
            [4*dxs3,  d3d5 ,  0.00 ,  p3d3 ,  p1d3 ]/4. ,$
            [4*pxs3,  p3d5 ,  p3d3 ,  0.00 ,  p1p3 ]/4. ,$
            [2*pxs3,  p1d5 ,  p1d3 ,  p1p3 ,  0.00 ]/2.], 5, 5)
;                    ^       ^       ^       ^       ^     ^
;                  4S3/2   2D5/2   2D3/2   2P3/2   2P1/2   wi
          en = reform([ (ee-ee(0) > 0), (ee-ee(1) > 0), (ee-ee(2) > 0), $
            (ee-ee(3) > 0), (ee-ee(4) > 0)], 5, 5)
       END
;
       'P3B':    BEGIN          ;     Order of energy levels is
          ombyw = reform($
            [[ 0.00 , 4*dxs3, 6*dxs3, 4*pxs3, 2*pxs3]/4. ,$
            [4*dxs3,  0.00 ,  d5d3 ,  p3d3 ,  p1d3 ]/4. ,$
            [6*dxs3,  d5d3 ,  0.00 ,  p3d5 ,  p1d5 ]/6. ,$
            [4*pxs3,  p3d3 ,  p3d5 ,  0.00 ,  p1p3 ]/4. ,$
            [2*pxs3,  p1d3 ,  p1d5 ,  p1p3 ,  0.00 ]/2.], 5, 5)
;                    ^       ^       ^       ^       ^     ^
;                  4S3/2   2D3/2   2D5/2   2P3/2   2P1/2   wi
          en = reform([ (ee-ee(0) > 0), (ee-ee(1) > 0), (ee-ee(2) > 0), $
            (ee-ee(3) > 0), (ee-ee(4) > 0)], 5, 5)
       END
;
       'P3C':    BEGIN          ;     Order of energy levels is NI-like
          ombyw = reform($
            [[ 0.00 , 6*dxs3, 4*dxs3, 2*pxs3, 4*pxs3]/4. ,$
            [6*dxs3,  0.00 ,  d3d5 ,  p1d5 ,  p3d5 ]/6. ,$
            [4*dxs3,  d3d5 ,  0.00 ,  p1d3 ,  p3d3 ]/4. ,$
            [2*pxs3,  p1d5 ,  p1d3 ,  0.00 ,  p3p1 ]/2. ,$
            [4*pxs3,  p3d5 ,  p3d3 ,  p3p1 ,  0.00 ]/4.], 5, 5)
;                    ^       ^       ^       ^       ^     ^
;                  4S3/2   2D5/2   2D3/2   2P1/2   2P3/2   wi
          en = reform([ (ee-ee(0) > 0), (ee-ee(1) > 0), (ee-ee(2) > 0), $
            (ee-ee(3) > 0), (ee-ee(4) > 0)], 5, 5)
       END
;
       'P3D':    BEGIN          ;     Order of energy levels is SII-like
          ombyw = reform($
            [[ 0.00 , 4*dxs3, 6*dxs3, 2*pxs3, 4*pxs3]/4. ,$
            [4*dxs3,  0.00 ,  d5d3 ,  p1d3 ,  p3d3 ]/4. ,$
            [6*dxs3,  d5d3 ,  0.00 ,  p1d5 ,  p3d5 ]/6. ,$
            [2*pxs3,  p1d3 ,  p1d5 ,  0.00 ,  p3p1 ]/2. ,$
            [4*pxs3,  p3d3 ,  p3d5 ,  p3p1 ,  0.00 ]/4.], 5, 5)
;                    ^       ^       ^       ^       ^     ^
;                  4S3/2   2D3/2   2D5/2   2P1/2   2P3/2   wi
          en = reform([ (ee-ee[0] > 0), (ee-ee[1] > 0), (ee-ee[2] > 0), $
            (ee-ee[3] > 0), (ee-ee[4] > 0)], 5, 5)

       END
;
       'P4':     BEGIN
          ombyw = reform($
            [[ 0.00 ,  p1p2 ,  p0p2 , 5*d2px, 5*s0px]/5. , $
            [ p1p2 ,  0.00 ,  p0p1 , 3*d2px, 3*s0px]/3. , $
            [ p0p2 ,  p0p1 ,  0.00 , 1*d2px, 1*s0px]/1. , $
            [1*d2px, 3*d2px, 5*d2px,  0.00 ,  s0d2 ]/5. , $
            [1*s0px, 3*s0px, 5*s0px,  s0d2 ,  0.00 ]/1.], 5, 5)
;                  ^       ^       ^       ^       ^     ^
;                 3P2     3P1     3P0     1D2     1S0    wi
;
          en = reform([ (ee-ee(0) > 0), (ee-ee(1) > 0), (ee-ee(2) > 0), $
            (ee-ee(3) > 0), (ee-ee(4) > 0)], 5, 5)
       END
;
       '3LEV':   BEGIN
          ombyw = reform($
            [[ 0.0, O10, O20]/G0 , $
            [ O10, 0.0, O21]/G1 , $
            [ O20, O21, 0.0]/G2], 3, 3)
;
          en = reform([ (ee-ee(0) > 0), (ee-ee(1) > 0), (ee-ee(2) > 0)], 3, 3)
;
       END
;
       '5LEV':   BEGIN
          ombyw = reform($
            [[ 0.0, O10, O20, O30, O40 ]/G0 , $
            [ O10, 0.0, O21, O31, O41 ]/G1 , $
            [ O20, O21, 0.0, O32, O42 ]/G2 , $
            [ O30, O31, O32, 0.0, O43 ]/G3 , $
            [ O40, O41, O42, O43, 0.0  ]/G4], 5, 5)
;
          en = reform([ (ee-ee(0) > 0), (ee-ee(1) > 0), (ee-ee(2) > 0), $
            (ee-ee(3) > 0), (ee-ee(4) > 0)], 5, 5)
;
       END
;
       'TEST':   BEGIN
          ombyw = reform($
            [[ 0.00 ,  p1p0 ,  p2p0 ,  d2px ,  s0px ,  s2px ] , $
            [ p1p0 ,  0.00 ,  p2p1 ,  d2px ,  s0px ,  s2px ] , $
            [ p2p0 ,  p2p1 ,  0.00 ,  d2px ,  s0px ,  s2px ] , $
            [ d2px ,  d2px ,  d2px ,  0.00 ,  s0d2 ,  s2d2 ] , $
            [ s0px ,  s0px ,  s0px ,  s0d2 ,  0.00 ,  s2s0 ] , $
            [ s2px ,  s2px ,  s2px ,  s2d2 ,  s2s0 ,  0.00 ]], 6, 6)
;                  ^       ^       ^       ^       ^       ^     ^
;                 3P0     3P1     3P2     1D2     1S0     5S2    wi
;
          en = reform([ (ee-ee(0) > 0), (ee-ee(1) > 0), (ee-ee(2) > 0), $
            (ee-ee(3) > 0), (ee-ee(4) > 0), (ee-ee(5) > 0)], 6, 6)
       END
;
       ELSE :    BEGIN
          ombyw = identity(5)
          en = ombyw
          einas = en
       END

    endcase

    tev = tt/1.16045d4
    qval = 8.629d-6*transpose(ombyw*exp(-en/tev))/sqrt(tt)


; pack everything into a structure and return

    matrix = {ion: ionstr, temp: tt, qij: qval, eij: en, aij: einas}
    
return, matrix
end
