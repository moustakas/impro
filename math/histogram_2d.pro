pro histogram_2d,x,y,hist,xrange,yrange,nxbins,nybins,silent=silent

;+
; NAME:
;       histogram_2d
;
;PURPOSE:
;       takes two dimentional data x,y and makes a histogram, 
;       hist,with hist_2d.  hist is a structure that also contains
;       the scaling parameters so one can display
;       it with the right x,y, scale etc.
;
;CALLING SEQUENCE:
;
;       -syntax histogram_2d,x,y,hist,xrange,yrange,nxbins,nybins,silent=silent
;
;INPUTS:
;       x ,y : the two one-dimentional arrays that you want histogrammed
;       
;OUTPUTS:
;       hist : the histogram structure contains 3 fields
;       .map (the 2D array) , .xrange, .yrange (the two ranges)
;       these ranges permit mapping from the data to the histogram map
;       and vice versa.  also contains the mean x and mean y (xbins,ybins)
;
;KEYWORDS:
;       xrange,yrange : the range to include data from
;       these are output and saved in the hist structure
;       the default is min, max
;       if these are flipped like [12,3]
;       it will use 3 as min and 12 as max and then flip the histogram
;       
;       nxbins, nybins : the number of bins for each dimention
;
;       silent : speak no evil  
;
;EXTERNAL CALLS:
;       none
;
;METHOD:
;       uses hist2d rather than the built in
;       IDL routine hist_2d because it has to work with floating
;       point numbers as well
;
;EXAMPLE:
;  IDL> histogram2d,radius,mag,hist,[0,6],[24,14],250,250
;  IDL> tvim2,hist.map,range=[0,100],xrange=hist.xrange,yrange=hist.yrange
;
;NOTES
;
;HISTORY:  written by David Johnston -University of Chicago
;       Physics Department 1999
;       Erin Sheldon, UofChicago. Added calculation of mean in x and y direction
;             30-Jul-2003
;       Returns bin ranges and bin middles as well as mean in each bin from
;       data.  July 2006 Erin Sheldon , NYU
;-
;
;
;
;  Copyright (C) 2006  Dave Johnston
;
;    This program is free software; you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation; either version 2 of the License, or
;    (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.
;
;    You should have received a copy of the GNU General Public License
;    along with this program; if not, write to the Free Software
;    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
;
;

if n_params() eq 0 then begin
        print,'-syntax histogram_2d,x,y,hist,xrange,yrange,nxbins,nybins,silent=silent'
        return
endif

if n_elements(xrange) eq 0 then xrange=[min(x),max(x)]
if n_elements(yrange) eq 0 then begin
        yrange=[min(y),max(y)]
        if not keyword_set(silent) then begin
                print,'xrange= ',xrange
                print,'yrange= ',yrange
        endif
endif

if n_elements(nxbins) eq 0 then nxbins=100
if n_elements(nybins) eq 0 then nybins=100
;the default size

min1=xrange(0) < xrange(1)
min2=yrange(0) < yrange(1)
max1=xrange(0) > xrange(1)
max2=yrange(0) > yrange(1)

;this is because xrange may be something like [12.0, 3.0]
;in which case the caller wants 3,0 to be the min and
;12.0 to be the max BUT the output array should be flipped

w=where(x gt min1 and x lt max1 and y gt min2 and y lt max2,wif)
;the relevent data

if not keyword_set(silent) then print,wif,' in histogram'

xx=x(w)-min1
yy=y(w)-min2

rangex=float(abs(xrange(1)-xrange(0)))
rangey=float(abs(yrange(1)-yrange(0)))

hist2d,xx,yy,h,[0,rangex],[0,rangey],nxbins,nybins

h=float(h)
xrange=float(xrange)
yrange=float(yrange)

if xrange(0) gt xrange(1) then h=rotate(h,5)
;flip x
if yrange(0) gt yrange(1) then h=rotate(h,7)
;flip y

;; Calculate the middles of the bins
xbins_low = fltarr(nxbins)
xbins_high = xbins_low
ybins_low = fltarr(nybins)
ybins_high = ybins_low

xstep = (max1-min1)/nxbins
xbins_low[0] = min1
xbins_high[0] = xbins_low[0] + xstep
for i=1l, nxbins-1 do begin 
    xbins_low[i] = xbins_low[i-1] + xstep
    xbins_high[i] = xbins_high[i-1] + xstep
endfor 
xbins = (xbins_low + xbins_high)/2.0

ystep = (max2-min2)/nybins
ybins_low[0] = min2
ybins_high[0] = ybins_low[0] + ystep
for i=1l, nybins-1 do begin 
    ybins_low[i] = ybins_low[i-1] + ystep
    ybins_high[i] = ybins_high[i-1] + ystep
endfor 
ybins = (ybins_low + ybins_high)/2.0


;; calculate the mean x and mean y
xmeans = dblarr(nxbins)
ymeans = dblarr(nybins)

xmin = min(xrange,max=xmax)
ymin = min(yrange,max=ymax)



hx = histogram(x, min=xmin, max=xmax, nbins=nxbins, rev=xrev)
hy = histogram(y, min=ymin, max=ymax, nbins=nybins, rev=yrev)

for i=0l, nxbins-1 do begin 
    if xrev[i] ne xrev[i+1] then begin 
        w2=xrev[ xrev[i]:xrev[i+1]-1  ]
        nw2 = n_elements(w2)
        
        xmeans[i] = total(x[w2])/nw2
        
    endif 
endfor 

for i=0l, nybins-1 do begin 
    if yrev[i] ne yrev[i+1] then begin 
        w2=yrev[ yrev[i]:yrev[i+1]-1  ]
        nw2 = n_elements(w2)

        ymeans[i] = total(y[w2])/nw2

    endif 
endfor 


if xrange[0] gt xrange[1] then begin 
    xbins_low=reverse(xbins_low)
    xbins_high=reverse(xbins_high)
    xbins=reverse(xbins)
    xmeans = reverse(xmeans)
endif 
if yrange[0] gt yrange[1] then begin 
    ybins_low=reverse(ybins_low)
    ybins_high=reverse(ybins_high)
    ybins=reverse(ybins)
    ymeans = reverse(ymeans)
endif 


hist={map:h,$
      xrange:xrange,$
      yrange:yrange,$
      $
      xbins_low: xbins_low, $
      xbins_high: xbins_high, $
      xbins: xbins,$
      $
      ybins_low: ybins_low, $
      ybins_high: ybins_high, $
      ybins: ybins,$      
      $
      xmeans:xmeans,$
      ymeans:ymeans}
;make the hist structure

return
end
