load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
;
; this plots LES output data using command line arguments to control plots
; arguments include:
;    xmin, xmax, ymin, ymax
;    ffld : field name
;    tbeg : start time for averaging of ps data
;    tend : end time for averaging of ps data
;    foname : name of output file
;    fname : name of input file
;    ncols : number of columns
;    npages: number of pages
;
; example call statements can be found in plotfld.csh
;
begin

  colors = (/"black","red","dodgerblue","cyan","magena","orange"/)
  dashes = (/0,1,2,3,4,5/)
  if (.not.isvar("foname")) then
    foname = "plotfld"
  end if
  
  f1 = addfile(fname(0),"r")
  nplots = dimsizes(ffld)

  compare = False
  if (dimsizes(fname).gt.1) then
    compare = True
  end if

  wks = gsn_open_wks("pdf",foname)
  plots = new ( (/nplots/), "graphic")

  time = f1->time/3600.
  time@long_name = "Time [hrs]"
;
; set up default plot style
;
  xyres = True
  xyres@gsnDraw = False
  xyres@gsnFrame = False
  xyres@xyLineThicknessF = 2.
;
; remove upper and right borders
;
  xyres@tmYRBorderOn = False
  xyres@tmYROn = False
  xyres@tmXTBorderOn = False
  xyres@tmXTOn = False
;
; define plots
;
  do k=0,nplots-1
    dmy = f1->$ffld(k)$

    if (isvar("ymax")) then
      xyres@trYMaxF = ymax
    end if
    if (isvar("xmax")) then
      xyres@trXMaxF = xmax
    end if
    if (isvar("ymin")) then
      xyres@trYMinF = ymin
    end if
    if (isvar("xmin")) then
      xyres@trXMinF = xmin
    end if

    xyres@gsnRightString = f1->$ffld(k)$@longname
    xyres@tiYAxisString = f1->$ffld(k)$@units

    if ( dimsizes(dimsizes(dmy)) .gt. 1) then
      if (.not.isvar("tbeg")) then
        tbeg = f1->time(dimsizes(f1->time)-1)
      end if
      if (.not.isvar("tend")) then
        tend = f1->time(dimsizes(f1->time)-1)
      end if

      zvar = f1->$ffld(k)$!1
      xfld = dim_avg(dmy({$zvar$|0:},{time|tbeg:tend}))
      yfld = f1->$zvar$({0:})
      plots(k) = gsn_csm_xy(wks,xfld,yfld,xyres)
      do n=1,dimsizes(fname)-1
        fn = addfile(fname(n),"r")
        delete (dmy)
        delete (xfld)
        delete (yfld)
        delete (zvar)
        dmy = fn->$ffld(k)$
        zvar = fn->$ffld(k)$!1
        xfld = dim_avg(dmy({$zvar$|0:},{time|tbeg:tend}))
        yfld = fn->$zvar$({0:})
        xyres@xyLineColor = colors(n)
        xyres@xyDashPattern = dashes(n)
        ovrly = gsn_csm_xy(wks,xfld,yfld,xyres)
        xyres@xyLineColor = "black"
        xyres@xyDashPattern = dashes(0)
        overlay (plots(k),ovrly)
      end do
    else
      xyres@vpHeightF = 0.15
      xfld = f1->time
      yfld = dmy
      plots(k) = gsn_csm_xy(wks,xfld,yfld,xyres)
      do n=1,dimsizes(fname)-1
        fn = addfile(fname(n),"r")
        delete (dmy)
        delete (xfld)
        delete (yfld)
        dmy  = fn->$ffld(k)$
        xfld = fn->time
        yfld = dmy
        xyres@xyLineColor = colors(n)
        xyres@xyDashPattern = dashes(n)
        ovrly = gsn_csm_xy(wks,xfld,yfld,xyres)
        xyres@xyLineColor = "black"
        xyres@xyDashPattern = dashes(0)
        overlay (plots(k),ovrly)
      end do
    end if
    delete (xfld)
    delete (yfld)
    delete (dmy)
  end do
;
; lay out plots
;
  pres = True
  pres@gsnPanelYWhiteSpacePercent = 5.
  pres@gsnFrame = False

  tres = True
  tres@txJust = "BottomRight"
  tres@txFontHeightF = 0.0075
  tres@txFuncCode = "#"
  tstring = "Created with plotfld.ncl on "+systemfunc("date")+" from " \
  +fname(0)

  if (.not.isvar("npages")) then
    npages = 1
  end if
  do k=1,npages
    if (.not.isvar("ncols")) then
      ncols=1
    end if
    nrows = nplots/ncols/npages
    n1 = (k-1)*nrows*ncols
    n2 = min((/nplots-1,k*nrows*ncols/))
    if (k.eq.npages .and. n2 .lt. nplots-1) then
      n2 = nplots-1
      nrows = nrows+1
    end if
    gsn_panel(wks,plots(n1:n2),(/nrows,ncols/),pres)
    gsn_text_ndc(wks,tstring,1.,0.,tres)
    frame(wks)
  end do
end





