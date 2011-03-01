;+
; NAME:
;
;  GANG_PLOT_POS
;
; PURPOSE:
;
;  Simplify abutting multiple plots.
;
; CATEGORY:
;
;  Plotting.
;
; CALLING SEQUENCE:
;
;  struct=gang_plot_pos(nrow,ncol,col,row, [XTICKFORMAT=,YTICKFORMAT=,
;                       OFFSET=,SIZE=, /POSITION_ONLY]
;
; INPUTS:
;
;  nrow: The number of rows of plots to abut.
;
;  ncol: The number of columns of plots to abut.
;
;  col_ind: The column index starting at 0, or, if row isn't passed,
;    the 1D index in the grid [0 .. ncol*nrow-1].
;
; OPTIONAL INPUTS:
;
;  row: The row index starting at 0.
;
; KEYWORD PARAMETERS:
;
;  (X|Y)TICKFORMAT: Format of tick labeling used (if different from
;    IDL's default).  Note that tick labels will be placed on the
;    bottom and left of the plot only.
;
;  OFFSET: The normalized offset of the lower left plot ([0.1,0.1] by default)
;
;  SIZE: The normalized x and y size of the region in which the plots
;    should be positioned. (default [0.8,0.85] for default offset).
;
;  POSITION_ONLY: If set, the position will be returned rather than a
;    full _EXTRA structure.
;
; OUTPUTS:
;
;  struct: A structure in the form of EXTRA structure for keyword
;    inheritence, which can be passed directly to plot and will set
;    the tick labelling correctly (see example).
;
; RESTRICTIONS: 
; 
;  With the _EXTRA type format, the plot is not erased, so may need
;  to issue ERASE beforehand to clean the screen.
;
; EXAMPLE:
;   
;  Gang two plots vertically:
;  plot,findgen(10),  _EXTRA=gang_plot_pos(2,1,0)
;  plot,findgen(10)^2,_EXTRA=gang_plot_pos(2,1,1)
;
;  Gang 3x3 plots:
;
;   for c=0,2 do for r=0,2 do $
;     plot,findgen(10)^(float(c+3*r)/5),YRANGE=[0,9.^(((r+1)*3-1)/5.)],
;          _EXTRA=gang_plot_pos(3,3,c,r)
;
; MODIFICATION HISTORY:
;
;  Feb 7, 2006, J.D. Smith: Written
;
;-
;##############################################################################
;
; LICENSE
;
;  Copyright (C) 2006 J.D. Smith
;
;  This file is free software; you can redistribute it and/or modify
;  it under the terms of the GNU General Public License as published
;  by the Free Software Foundation; either version 2, or (at your
;  option) any later version.
;
;  This file is distributed in the hope that it will be useful, but
;  WITHOUT ANY WARRANTY; without even the implied warranty of
;  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
;  General Public License for more details.
;
;  You should have received a copy of the GNU General Public License
;  along with this file; see the file COPYING.  If not, write to the
;  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
;  Boston, MA 02110-1301, USA.
;
;##############################################################################

function gang_plot_pos,nrow,ncol,col_ind,row, XTICKFORMAT=xt, $
                       YTICKFORMAT=yt,OFFSET=off,SIZE=sz, POSITION_ONLY=po
  
  if n_params() eq 3 then begin ;; 1D index, wrap
     col_ind mod= nrow*ncol
     row=col_ind/ncol
     col=col_ind mod ncol
  endif else col=col_ind
  
  if n_elements(off) eq 0 then off=[.1,.1] ; normalized offset of plot position
  if n_elements(sz) eq 0 then $
     sz=((1.-2.*off)+[0.,.05])<(1.-off-.05) ; size of plot win, normalized
  
  pane_sz=sz/[ncol,nrow]
  
  left=off[0]+col*pane_sz[0]
  right=off[0]+(col+1)*pane_sz[0]
  top=off[1]+(nrow-row)*pane_sz[1]
  bottom=off[1]+(nrow-1-row)*pane_sz[1]
  
  pos=[left,bottom,right,top]
  
  if keyword_set(po) then return,pos
  
  at_top=row eq 0
  at_bottom=row eq (nrow-1)
  at_left=col eq 0
  at_right=col eq (ncol-1)
  
  if n_elements(xt) eq 0 then xt=''
  if n_elements(yt) eq 0 then yt=''
  
  extra={NOERASE:1, $
         POSITION: pos, $
         YTICKFORMAT: at_left?yt:'(A1)', $
         XTICKFORMAT: at_bottom?xt:'(A1)'}

  return,extra
end
