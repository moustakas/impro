pro jd_spherematch
; written by JD Smith; Match stars in one list to another, within some
; tolerance.  Pre-bin into a 2D histogram, and use DUAL HISTOGRAM
; matching to select

    n1=1000000                  ;number of stars
    x1=randomu(sd,n1)           ;points to find matches near
    y1=randomu(sd,n1)

    n2=n1
    x2=randomu(sd,n2)           ;points to search in
    y2=randomu(sd,n2)

    t1=systime(1)

    max_r=.0005                 ;maximum allowed radius for a match
    bs=2*max_r                  ;this is the smallest binsize allowed
    h=hist_nd([1#x2,1#y2],bs,MIN=0,MAX=1,REVERSE_INDICES=ri)
    bs=bs[0]
    d=size(h,/DIMENSIONS)

;; Bin location of X1,Y1 in the X2,Y2 grid
    xoff=x1/bs       &  yoff=y1/bs
    xbin=floor(xoff) &  ybin=floor(yoff)
    bin=(xbin + d[0]*ybin)<(d[0]*d[1]-1L) ;The bin it's in

;; We must search 4 bins worth for closest match, depending on
;; location within bin (towards any of the 4 quadrants).
    xoff=1-2*((xoff-xbin) lt 0.5) ;add bin left or right
    yoff=1-2*((yoff-ybin) lt 0.5) ;add bin down or up

    min_pos=make_array(n1,VALUE=-1L)
    min_dist=fltarr(n1,/NOZERO)

    for i=0,1 do begin ;; Loop over 4 bins in the correct quadrant direction
       for j=0,1 do begin
          b=0L>(bin+i*xoff+j*yoff*d[0])<(d[0]*d[1]-1) ;current bins (offset)

          ;; Dual HISTOGRAM method, loop by repeat count in bins
          h2=histogram(h[b],MIN=1,REVERSE_INDICES=ri2)

          ;; Process all bins with the same number of repeats >= 1
          for k=0L,n_elements(h2)-1 do begin
             if h2[k] eq 0 then continue
             these_bins=ri2[ri2[k]:ri2[k+1]-1] ;the points with k+1 repeats in bin

             if k eq 0 then begin ; single point (n)
                these_points=ri[ri[b[these_bins]]]
             endif else begin   ; range over k+1 points, (n x k+1)
                these_points=ri[ri[rebin(b[these_bins],h2[k],k+1,/SAMPLE)]+ $
                  rebin(lindgen(1,k+1),h2[k],k+1,/SAMPLE)]
                these_bins=rebin(temporary(these_bins),h2[k],k+1,/SAMPLE)
             endelse

             dist=(x2[these_points]-x1[these_bins])^2 + $
               (y2[these_points]-y1[these_bins])^2

             if k gt 0 then begin ;multiple point in bin: find closest
                dist=min(dist,DIMENSION=2,p)
                these_points=these_points[p] ;index of closest point in bin
                these_bins=ri2[ri2[k]:ri2[k+1]-1] ;original bin list
             endif

             ;; See if a minimum is already set
             set=where(min_pos[these_bins] ge 0, nset, $
               COMPLEMENT=unset,NCOMPLEMENT=nunset)

             if nset gt 0 then begin
                ;; Only update those where the new distance is less
                closer=where(dist[set] lt min_dist[these_bins[set]], cnt)
                if cnt gt 0 then begin
                   set=set[closer]
                   min_pos[these_bins[set]]=these_points[set]
                   min_dist[these_bins[set]]=dist[set]
                endif
             endif

             if nunset gt 0 then begin ;; Nothing set, closest by default
                min_pos[these_bins[unset]]=these_points[unset]
                min_dist[these_bins[unset]]=dist[unset]
             endif
          endfor
       endfor
    endfor


    print,systime(1)-t1

return
end
    
