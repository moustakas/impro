; readplot.pro (beta version), A.R.Marble

; a final version of this routine will be available in the future


function yesorno

    read:
    input='' & read,input
    if input eq 'n' or input eq 'N' then return,'n'
    if input eq 'y' or input eq 'Y' then return,'y'
    print & print,"Please select either 'y' or 'n'."
    goto,read
    
end

pro readplot,infile,outfile

    wsize=600
    tvlct,getcolor('white'),0
    tvlct,getcolor('green'),1
    tvlct,getcolor('blue'),2
    tvlct,getcolor('yellow'),3
    tvlct,getcolor('red'),4
    tvlct,getcolor('black'),10
    npts=49
    xsym = cos(findgen(npts)/(npts-1)*360.*(!pi/180.))
    ysym = sin(findgen(npts)/(npts-1)*360.*(!pi/180.))
    USERSYM, xsym, ysym, /fill

;   CHECK FOR VALID FITS FILE
    check=findfile(infile)
    if check[0] eq '' then begin
       print & print,'The file '+infile+' does not exist.'
       return
    endif else begin
       parts=strsplit(infile,'.',/extract)
       if n_elements(parts) ne 2 then begin
          print & print,'The file name '+infile+' is invalid.'
          print,'Valid names are "file.fits" or "file.FITS".' & print
          return
       endif else if parts[1] ne 'FITS' and parts[1] ne 'fits' then begin
          print & print,'The file name '+infile+' is invalid.'
          print,'Valid names are "file.fits" or "file.FITS".' & print
          return
       endif
    endelse 

    if not keyword_set(outfile) then outfile=parts[0]+'.dat'
    if infile eq outfile then begin
       print & print,'The input and output file names must be different.'
       print
       return
    endif
    check=findfile(outfile)
    if check[0] ne '' then begin
       print & print,'The output file '+outfile+' already exists.'
       print,'Do you wish to write over it? (y/n)'
       input=yesorno()
       if input eq 'n' then return
    endif

    im=readfits(infile)
    nx=n_elements(im[*,0])
    ny=n_elements(im[0,*])
    im=im[0:nx-1<800,0:ny-1<800]
    nx=n_elements(im[*,0])
    ny=n_elements(im[0,*])

    im=1.*(im-min(im))
    im=round(im/max(im)*245L+10L)


    if nx gt ny then begin
       nx2=wsize & ny2=wsize*ny/nx
    endif else begin
       nx2=wsize*nx/ny & ny2=wsize
    endelse
;    im=congrid(im,nx2,ny2)
nx2=nx & ny2=ny
;;;;;;;;;;;;;;;
;    window,xsize=nx2,ysize=ny2,color=2    
;    tv,im

;    print & print,'Do you wish to crop the image? (y/n)'
;    input=yesorno()
input='n'
    if input eq 'y' then begin
       zoom:
       print & print,'Select lower-left and upper-right'
       print,'corners of desired region with mouse...'
       cursor,xx1,yy1,3,/device
       xx1=xx1>0
       yy1=yy1>0
       cursor,xx2,yy2,3,/device
       xx2=xx2<nx2-1
       yyw=yy2<ny2-1
       copy=im
       copy[xx1:xx1,yy1:yy2]=4
       copy[xx1:xx2,yy2:yy2]=4
       copy[xx2:xx2,yy1:yy2]=4
       copy[xx1:xx2,yy1:yy1]=4
       tv,copy

       print & print,'Click inside the red box to crop,'
       print,'or click outside to select a different region...'
       cursor,xx3,yy3,3,/device
       
       if xx3 gt xx1 and xx3 lt xx2 and yy3 gt yy1 and yy3 lt yy2 then begin
;          im2=congrid(im[xx1:xx2,yy1:yy2],nx2,ny2-25)
          im2=im[xx1:xx2,yy1:yy2]
          nx2=xx2-xx1+1
          ny2=yy2-yy1+1
       endif else begin
          tv,im
          goto,zoom
       endelse
    endif else im2=im;congrid(im,nx2,ny2-25)
    
    im=fltarr(nx2,ny2+25)+0
    im[*,0:ny2-1]=im2
    window,0,xsize=nx2,ysize=ny2+25
    tv,im

    rotate:
    print & print,'Do you wish to rotate the image? (y/n)'
    input=yesorno()
    if input eq 'y' then begin
       print,'Select two points on what SHOULD be the horizontal axis...'
       cursor,x1,y1,3,/device
       cursor,x2,y2,3,/device
       angle=atan((y2-y1)*1./(x2-x1))*180./!pi
print,x1,x2
print,y1,y2
print,angle
       im2=rot(im2,angle)
       im=fltarr(nx2,ny2+25)+0
       im[*,0:ny2-1]=im2
       tv,im
       goto,rotate
    endif

    im=fltarr(nx2,ny2+25)+3
    im[*,0:ny2-1]=im2
    im[*,ny2:ny2+1]=10
    im[*,0:1]=10
    im[*,ny2+23:ny2+24]=10
    im[0:1,*]=10
    im[nx2-2:nx2-1,*]=10
    tv,im

    xyouts,nx2/2.,ny2+7,align=.5,'FOLLOW INSTRUCTIONS TO ESTABLISH SCALING....',charsize=1.5,color=10,/device

    print & print,'Select a point on x axis with mouse...'
    cursor,xpt1,dummy,3,/device
    print & print,'Enter x coordinate of selected point...'
    read,xnum1
    print & print,'Select another point on x axis with mouse...'
    cursor,xpt2,dummy,3,/device
    print & print,'Enter x coordinate of selected point...'
    read,xnum2
    print & print,'Select a point on y axis with mouse...'
    cursor,dummy,ypt1,3,/device
    print & print,'Enter y coordinate of selected point...'
    read,ynum1
    print & print,'Select another point on y axis with mouse...'
    cursor,dummy,ypt2,3,/device
    print & print,'Enter y coordinate of selected point...'
    read,ynum2


    im[2:nx2-3,ny2+2:ny2+22]=4    
    tv,im
    xyouts,nx2/2.,ny2+7,align=.5,'CLICK INSIDE RED RECTANGLE WHEN FINISHED SELECTING DATA POINTS...',charsize=1.5,color=10,/device

    xnum1=xnum1*1.
    xnum2=xnum2*1.
    ynum1=ynum1*1.
    ynum2=ynum2*1.
    xpt1=xpt1*1.
    xpt2=xpt2*1.
    ypt1=ypt1*1.
    ypt2=ypt2*1.

    print & print,'Is the x-axis log scale? (y/n)'
    xinput=yesorno()
    print & print,'Is the y-axis log scale? (y/n)'
    yinput=yesorno()

    if xinput eq 'y' then begin
       xnum1=alog10(xnum1)
       xnum2=alog10(xnum2)
    endif
    if yinput eq 'y' then begin
       ynum1=alog10(ynum1)
       ynum2=alog10(ynum2)
    endif

ny2=ny2+25


    xxvala=xnum1+(0-xpt1)/(xpt2-xpt1)*(xnum2-xnum1)
    xxvalb=xnum1+(nx2-1-xpt1)/(xpt2-xpt1)*(xnum2-xnum1)
    yyvala=ynum1+(0-ypt1)/(ypt2-ypt1)*(ynum2-ynum1)
    yyvalb=ynum1+(ny2-1-ypt1)/(ypt2-ypt1)*(ynum2-ynum1)

    if xinput eq 'y' then begin
       xxvala=10.^xxvala
       xxvalb=10.^xxvalb
    endif 
    if yinput eq 'y' then begin
       yyvala=10.^yyvala
       yyvalb=10.^yyvalb
    endif

    if xinput eq 'y' and yinput eq 'y' then begin
       plot,[0,0],/nodata,xr=[xxvala,xxvalb],xst=5,yr=[yyvala,yyvalb],yst=5, $
         position=[0,0,nx2-1,ny2-1],/noerase,color=2,/xlog,/ylog
    endif else if xinput eq 'y' and yinput ne 'y' then begin
       plot,[0,0],/nodata,xr=[xxvala,xxvalb],xst=5,yr=[yyvala,yyvalb],yst=5, $
         position=[0,0,nx2-1,ny2-1],/noerase,color=2,/xlog
    endif else if xinput ne 'y' and yinput ne 'y' then begin
       plot,[0,0],/nodata,xr=[xxvala,xxvalb],xst=5,yr=[yyvala,yyvalb],yst=5, $
         position=[0,0,nx2-1,ny2-1],/noerase,color=2
    endif else if xinput ne 'y' and yinput eq 'y' then begin
       plot,[0,0],/nodata,xr=[xxvala,xxvalb],xst=5,yr=[yyvala,yyvalb],yst=5, $
         position=[0,0,nx2-1,ny2-1],/noerase,color=2,/ylog
    endif


    print & print,"Okay, we're ready to select data points now."
    print,'Use the mouse to select a point (a blue dot'
    print,'will appear).  To de-select a point, simply'
    print,'click on the blue dot.  When you are done,'
    print,'click anywhere on the red bar at the top.'

    xdata=[-999]
    ydata=[-999]
    xcord=[-999]
    ycord=[-999]
    next:
    copy=im
    cursor,xx,yy,3,/device

    xxval=xnum1+(xx-xpt1)/(xpt2-xpt1)*(xnum2-xnum1)
    yyval=ynum1+(yy-ypt1)/(ypt2-ypt1)*(ynum2-ynum1)

    if xinput eq 'y' then xxval=10.^xxval
    if yinput eq 'y' then yyval=10.^yyval
    
    if yy ge ny2-26 then goto,done
    r=sqrt((xcord-xx)^2.+(ycord-yy)^2.)
    wh=where(r le 1)
    if wh[0] eq -1 then begin
       xcord=[xcord,xx]
       xdata=[xdata,xxval]
       ycord=[ycord,yy]
       ydata=[ydata,yyval]
       oplot,[xxval,xxval],[yyval,yyval],psym=8,color=1,symsize=1.
    endif else begin
       wh=(where(r eq min(r)))[0]
       oplot,[xdata[wh],xdata[wh]],[ydata[wh],ydata[wh]],psym=8,color=4,symsize=1.
       if wh[0] eq n_elements(xdata)-1 then begin
          xdata=xdata[0:wh-1]
          xcord=xcord[0:wh-1]
          ydata=ydata[0:wh-1]
          ycord=ycord[0:wh-1]
       endif else begin
          xdata=[xdata[0:wh-1],xdata[wh+1:n_elements(xdata)-1]]
          xcord=[xcord[0:wh-1],xcord[wh+1:n_elements(xcord)-1]]
          ydata=[ydata[0:wh-1],ydata[wh+1:n_elements(ydata)-1]]
          ycord=[ycord[0:wh-1],ycord[wh+1:n_elements(ycord)-1]]
       endelse
    endelse

    goto,next
    
    done:
    
    xdata=xdata[1:n_elements(xdata)-1]
    ydata=ydata[1:n_elements(ydata)-1]

;    ord=sort(xdata)
;    xdata=xdata[ord]
;    ydata=ydata[ord]
    
    oplot,xdata,ydata,color=1

    close,2
    openw,2,outfile
    for i=0,n_elements(xdata)-1 do printf,2,f='(e14.7,3x,e14.7)',xdata[i],ydata[i]
    close,2

stop

    return
    
 end

