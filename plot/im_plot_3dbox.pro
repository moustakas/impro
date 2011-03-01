;  $Id: plot_3dbox.pro,v 1.10 2003/02/03 18:13:22 scottm Exp $
;
; Copyright (c) 1994-2003, Research Systems, Inc.  All rights reserved.
;       Unauthorized reproduction prohibited.
;+
; NAME:
;	Plot_3dbox
;
; PURPOSE:
;	This procedure plots data in a 3-dimensional box, with options
;	to have the data displayed on the walls surrounding the plot area.
;
; CATEGORY:
;	Plotting, Three-dimensional
;
; CALLING SEQUENCE:
;	Plot_3dbox, X, Y, Z
; 
; INPUTS:
;	X:	A one dimensional array that contains the X coordinats
;
;	Y:	A one dimensional array that contains the Y coordinates
;
;	Z:	A one dimensional array that contains the Z coordinates
;
; OPTIONAL INPUTS:
;	None.
;	
; KEYWORD PARAMETERS:
;	COLOR:		The color for the Grid and Lines or the Color
;		 	for the box walls when the keyword SOLID_WALLS
;			is set.
;
;	BACKGROUND:	The background color of the plot or the color
;			of the Grid and Plot data when the SOLID_WALLS
;			keyword is set.
;
;	XY_PLANE:	Setting this keyword will cause the X and Y values
;			of the data to be plotted on the Z=0 axis plane.
;
;	XZ_PLANE:	Setting this keyword will cause the X and Z values
;			of the data to be plotted on the Y=Ymax axis plane.
;
;	YZ_PLANE:	Setting this keyword will cause the Y and Z values
;			of the data to be plotted on the X=Xmax axis plane.
;
;	SOLID_WALLS:	Setting this keyword causes the axis "walls" of 
;			the plot box to be filled with the value of COLOR.
;
;	PSYM:		The plotting symbol that the data is draw with.
;
;	GRIDSTYLE:	Set this keyword to the linestyle that will be 
;			used in drawing the gridlines.
;
;	TITLE:		Set this keyword to the Main plot title
;	
;	XTITLE:		Set this keyword to the X axis title.
;
;	YTITLE:		Set this keyword to the Y axis title.
;
;	ZTITLE:		Set this keyword to the Z axis title.
;
; 	SUBTITLE:	Set this keyword to the Sub-Title 
;
;	LINESTYLE:	The linestyle used to plot the data.
;
;	XYSTYLE:	The linesytle used to draw the plot in the XY plane.
;			If this keyword is not set, the value of LINESTYLE
;			is used.
;
;	XZSTYLE:	The linesytle used to draw the plot in the XZ plane.
;			If this keyword is not set, the value of LINESTYLE
;			is used.
;
;	YZSTYLE:	The linesytle used to draw the plot in the YZ plane.
;			If this keyword is not set, the value of LINESTYLE
;			is used.
;
;	Surface         All other keywords available to SURFACE are also
;	Keywords:	used by this procedure.
;
; OUTPUTS:
;	None.
;
; OPTIONAL OUTPUTS:
;	None.
;
; COMMON BLOCKS:
;	None.
;
; SIDE EFFECTS:
;	Plotting on the current device is performed.
;
; RESTRICTIONS:
;	None.
;
; PROCEDURE:
;	Straightforward.  Unrecognized keywords are passed to the SURFACE
;	procedure.  
;
; EXAMPLE:
;       Create some data that can be passed to Plot_3dbox
;
;       x = Replicate(5., 10)
;       x1 = cos(findgen(36)*10.*!dtor)*2.+5.     
;       x=[x,x1,x]     
;       y = findgen(56)     
;       z = Replicate(5., 10)
;       z1 =sin(findgen(36)*10.*!dtor)*2.+5.     
;       z=[z,z1,z]     
;
;     ; Plot this data in a "plot box" 
;
;       Plot_3dbox, X, Y, Z, /XY_PLANE, /YZ_PLANE, /XZ_PLANE, $
;                 /SOLID_WALLS, GRIDSTYLE=1, XYSTYLE=3, XZSTYLE=4, $
;                 YZSTYLE=5, AZ=40, TITLE="Example Plot Box",      $
;                 Xtitle="X Coodinate", Ytitle="Y Coodinate",      $
;                 Ztitle="Z Coodinate", SubTitle="Sub Title",      $
;                 /YSTYLE, ZRANGE=[0,10], XRANGE=[0,10],Charsize=1.6
;
;     ; Then to plot symbols on the locations of the above plot
;
;       plots, X, Y, Z, /T3D, PSYM=4, COLOR=!p.background
;
; MODIFICATION HISTORY:
;       6/94   KDB, RSI   - Initial Coding and Testing
;
;-

PRO im_Plot_3dbox, X, Y, Z, COLOR=COLOR, BACKGROUND=BACKGROUND,        $  
                    XY_PLANE=XY_PLANE, YZ_PLANE=YZ_PLANE,        $  
                    XZ_PLANE=XZ_PLANE, SOLID_WALLS=SOLID_WALLS,  $  
                    PSYM=PSYM, GRIDSTYLE=GRIDSTYLE, TITLE=TITLE, $
		    XTITLE=XTITLE, YTITLE=YTITLE, ZTITLE=ZTITLE, $
		    SUBTITLE=SUBTITLE, LINESTYLE=LINESTYLE,      $
                    XYSTYLE=XYSTYLE, YZSTYLE=YZSTYLE, XZSTYLE=XZSTYLE, $
                    _EXTRA=e, nogrid=nogrid

;  Set up simple error handling
       
   On_ERROR, 2

;  Lets make sure that all arrays are the same size     
       
   Xcnt = N_Elements(X)      
   Ycnt = N_Elements(Y)     
   Zcnt = N_Elements(Z)     
     
   if (Xcnt ne Ycnt) or (Xcnt ne Zcnt) then $
          Message, "X, Y and Z arrays must have same number of elements "      
     
;  Check the values of the keywords  
     
   if(N_Elements(PSYM)  eq 0)then      PSYM = 0  
   if(N_Elements(COLOR) eq 0)then      COLOR = !P.COLOR  
   if(N_Elements(BACKGROUND) eq 0)then BACKGROUND = !P.BACKGROUND  
   if(N_Elements(GRIDSTYLE) eq 0)then  GRIDSTYLE = 1  
   if(N_Elements(LINESTYLE) eq 0)then  LINESTYLE = 0  
   if(N_Elements(XYSTYLE) eq 0)then    XYSTYLE = LINESTYLE  
   if(N_Elements(XZSTYLE) eq 0)then    XZSTYLE = LINESTYLE
   if(N_Elements(YZSTYLE) eq 0)then    YZSTYLE = LINESTYLE

   if(not KeyWord_Set(TITLE))then      TITLE=''
   if(not KeyWord_Set(SUBTITLE))then   SUBTITLE=''
   if(not KeyWord_Set(XTITLE))then     XTITLE=''
   if(not KeyWord_Set(YTITLE))then     YTITLE=''
   if(not KeyWord_Set(ZTITLE))then     ZTITLE=''
   if(not KeyWord_Set(ZRANGE))then     ZRANGE=[Min(Z, MAX=zmax), zmax]

   if keyword_set(nogrid) then ticklen = 0L else ticklen = 1L
   
;  Use SURFACE to set up the coordinates system, handle titles ect...     
     
   Surface, FltArr(Xcnt,Xcnt), X, Y, /NODATA, /SAVE, TICKLEN=ticklen,        $ 
      COLOR=COLOR, XGRIDSTYLE=GRIDSTYLE, YGRIDSTYLE=GRIDSTYLE,         $ 
      ZGRIDSTYLE=GRIDSTYLE, BACKGROUND=BACKGROUND, XTICK_GET=xt,       $ 
      YTICK_GET=yt, ZTICK_GET=zt,ZRANGE=ZRANGE, SUBTITLE=SUBTITLE,     $
      TITLE=TITLE, YTITLE=YTITLE, ZTITLE=ZTITLE, XTITLE=XTITLE, _EXTRA=e    

   name = replicate(' ' ,30)   ; Make up null tick names  
     
;  See if the user wants to have "Solid" box walls  
  
   if(KeyWord_Set(SOLID_WALLS))then BEGIN  
  
   ;  Using the values of Crange, fill in the box walls with the value  
   ;  of color  
  
      PolyFill, [!X.Crange[0], !X.Crange, !X.Crange[1]],  $;bottom  
                [!Y.Crange, !Y.Crange[1], !Y.Crange[0]],  $  
                Replicate(!Z.Crange[0], 4), /T3D, COLOR=COLOR  
      PolyFill, [!X.Crange[0], !X.Crange, !X.Crange[1]],  $ ; Back  
                Replicate(!Y.Crange[1],4),                $  
                [!Z.Crange, !Z.Crange[1], !Z.Crange[0]],  $  
                /T3D, COLOR=COLOR  
      PolyFill, Replicate(!X.Crange[1],4), /T3D, COLOR=COLOR, $ ;side  
                [!Y.Crange, !Y.Crange[1], !Y.Crange[0]],  $  
                [!Z.Crange[0], !Z.Crange, !Z.Crange[1]]  

   ;  Now Replot the surface data  
  
      COLOR=BACKGROUND  ; reverse colors 

      Surface, FltArr(Xcnt,Xcnt), X, Y, /NODATA, /SAVE, TICKLEN=ticklen,       $ 
              /NOERASE, COLOR=COLOR, XTICKNAME=name, YTICKNAME=name,     $  
              ZTICKNAME=name, XGRIDSTYLE=GRIDSTYLE, YGRIDSTYLE=GRIDSTYLE,$  
              ZGRIDSTYLE=GRIDSTYLE, ZRANGE=ZRANGE, _EXTRA=e     

   ENDif  ;if solid walls set.  
   
;  Now complete the drawing of the axis box     
     
   Axis, !X.Crange[1], !Y.Crange[1], !Z.Crange[1], /YAXIS, /T3D,    $ 
         YTICKNAME = name, COLOR=COLOR, YTICKLEN=0     
   Axis, !X.Crange[1], !Y.Crange[0], !Z.Crange[0], ZAXIS=-1, /T3D,  $ 
         ZTICKNAME = name, COLOR=COLOR, ZTICKLEN=0     
   Axis, !X.Crange[0], !Y.Crange[1], !Z.Crange[0], /XAXIS, /T3D,    $ 
         XTICKNAME = name, COLOR=COLOR, XTICKLEN=0     
   Axis, !X.Crange[0], !Y.Crange[1], !Z.Crange[1], /XAXIS, /T3D,    $  
         XTICKNAME = name, COLOR=COLOR, XTICKLEN=0     
   Axis, !X.Crange[1], !Y.Crange[1], !Z.Crange[0], ZAXIS=-1, /T3D,  $ 
         ZTICKNAME = name, COLOR=COLOR, ZTICKLEN=0     
   Axis, !X.Crange[1], !Y.Crange[0], !Z.Crange[1], YAXIS=-1, /T3D,  $ 
         YTICKNAME = name, COLOR=COLOR, YTICKLEN=0     
   Axis, !X.Crange[1], !Y.Crange[0], !Z.Crange[0], /YAXIS, /T3D,    $ 
         YTICKNAME = name, COLOR=COLOR, YTICKLEN=0     

; now plot the data     
     
 Plots, X, Y, Z, /T3D , PSYM=PSYM, COLOR=COLOR ,LINESTYLE=LINESTYLE, $
                        _EXTRA=e

; And now plot the data along the walls of the box if requested

  if(KeyWord_Set(XY_PLANE))then $     
      Plots, X, Y, Replicate(!Z.Crange[0], Zcnt), /T3D, COLOR=COLOR, $
               LINESTYLE=XYSTYLE, psym=psym
  if(KeyWord_Set(YZ_PLANE))then $     
      Plots, Replicate(!X.Crange[1], Xcnt), Y, Z, /T3D, COLOR=COLOR, $
               LINESTYLE=YZSTYLE, psym=psym
  if(KeyWord_Set(XZ_PLANE))then $     
      Plots, X, Replicate(!Y.Crange[1], Ycnt), Z, /T3D, COLOR=COLOR, $
               LINESTYLE=XZSTYLE, psym=psym

; now draw the grid lines on the plot that were not drawn by surface     

  if not keyword_set(nogrid) then begin
     
     for i=0, N_Elements(yt)-1 do $     
       Plots, [!X.Crange[1],!X.Crange[1]], [yt[i],yt[i]], !Z.Crange, /T3D, $  
         LINESTYLE=GRIDSTYLE, COLOR=COLOR  
     for i=0, N_Elements(zt)-1 do $     
       Plots, [!X.Crange[1],!X.Crange[1]], !Y.Crange, [zt[i],zt[i]], /T3D, $  
         LINESTYLE=GRIDSTYLE, COLOR=COLOR  
     for i=0, N_Elements(xt)-1 do $     
       Plots, [xt[i],xt[i]], [!Y.Crange[1],!Y.Crange[1]], !Z.Crange, /T3D, $  
         LINESTYLE=gridstyle, COLOR=COLOR  

  endif
 
END


