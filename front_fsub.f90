MODULE front_fsub

CONTAINS
!*********************************************************
SUBROUTINE mk_a1diagrat_sphere(a1lon, a1lat, nx, ny, a1diagrat)
implicit none
!--------------------------
! make ratio of area relative to a grid located on the equator
!------------------------
integer                                 nx,ny
!-- in ------------------
real,dimension(nx)                   :: a1lon
!f2py intent(in)                        a1lon

real,dimension(ny)                   :: a1lat
!f2py intent(in)                        a1lat
!-- out -----------------
real,dimension(ny)                   :: a1diagrat
!f2py intent(out)                       a1diagrat
!-- calc ----------------
integer                                 iy
real                                    lats, latn, lonw, lone
real                                    diag_eq
!!------------------------
lonw = a1lon(1)
lone = a1lon(2)
diag_eq  = hubeny_real(a1lat(int(ny*0.5)), lonw, a1lat(int(ny*0.5)+1), lone)
do iy =2,ny-1
  latn = (a1lat(iy+1) + a1lat(iy))*0.5
  lats = (a1lat(iy-1) + a1lat(iy))*0.5
  a1diagrat(iy) = hubeny_real(lats, lonw, latn, lone) / diag_eq
end do
do iy =1,ny,ny-1
  latn = (a1lat(2) + a1lat(1))*0.5   !<- same for iy=1 and iy=ny
  lats = -90.0                       !<- save for iy=1 and iy=ny
  a1diagrat(iy) = hubeny_real(lats, lonw, latn, lone) / diag_eq
end do
!------------------------
return
END SUBROUTINE mk_a1diagrat_sphere

!*********************************************************
SUBROUTINE mk_a1arearat_sphere(a1lat, ny, a1arearat)
implicit none
!--------------------------
! make ratio of area relative to a grid located on the equator
!------------------------
integer                                 ny
!-- in ------------------
real,dimension(ny)                   :: a1lat
!f2py intent(in)                        a1lat
!-- out -----------------
real,dimension(ny)                   :: a1arearat
!f2py intent(out)                       a1arearat
!-- calc ----------------
integer                                 iy
real                                    lats, latn, lonw, lone
real                                    area_eq
!!------------------------
area_eq  = cal_area_sphere(a1lat(int(ny/2)), a1lat(int(ny/2)+1), 0.0, 1.0)
do iy =2,ny-1
  latn = (a1lat(iy+1) + a1lat(iy))*0.5
  lats = (a1lat(iy-1) + a1lat(iy))*0.5
  lonw = 0.0
  lone = 1.0
  a1arearat(iy) = cal_area_sphere(lats, latn, lonw, lone) / area_eq
end do
do iy =1,ny,ny-1
  latn = (a1lat(2) + a1lat(1) )*0.5  !<- same for iy=1 and iy=ny
  lats = -90.0                       !<- save for iy=1 and iy=ny
  lonw = 0.0
  lone = 1.0
  a1arearat(iy) = cal_area_sphere(lats, latn, lonw, lone) / area_eq
end do
!------------------------
return
END SUBROUTINE mk_a1arearat_sphere
!!*********************************************************
FUNCTION cal_area_sphere(lats, latn, lonw, lone)
!----------------------
! estimate area (km2) assuming that the eath is a sphere
! S = r2 * dlon * pi / 180 * (sin(lat2) - sin(lat1))
! lat, lon are in degree
!----------------------
implicit none
!--- in -----------
real                                    lats, latn, lonw, lone  ! (deg.)
!f2py intent(in)                        lats, latn, lonw, lone

!--- out ----------
real                                    cal_area_sphere
!f2py intent(out)                       cal_area_sphere
!--- para ---------
real,parameter                       :: r = 6371.012  ! (km)
real,parameter                       :: pi = 3.1416
!--- calc ---------
real                                    dlon
!----------------------
dlon = min( abs(lone-lonw), abs(360.0-(lone-lonw)))

cal_area_sphere = r**2.0 *dlon*pi/180.0 *(sin(latn/180.0*pi)-sin(lats/180.0*pi))
!----------------------
return
END FUNCTION cal_area_sphere



!*********************************************************
SUBROUTINE mk_8gridsxy(nx, ny, ix, iy, a1x, a1y)
implicit none
!--- in ---------
integer                  nx, ny
!f2py intent(in)         nx, ny

integer                  ix, iy
!f2py intent(in)         ix, iy
!--- out --------
integer,dimension(8)  :: a1x, a1y
!f2py intent(in)         a1x, a1y
!--- calc -------
integer                  ixn, ixs, ixe, ixw
integer                  iyn, iys, iye, iyw
!----------------
call ixy2iixy(nx, ny, ix, iy+1, ixn, iyn)
call ixy2iixy(nx, ny, ix, iy-1, ixs, iys)
call ixy2iixy(nx, ny, ix-1, iy, ixw, iyw)
call ixy2iixy(nx, ny, ix+1, iy, ixe, iye)
!------
a1x(1) = ixn
a1x(2) = ixe
a1x(3) = ixe
a1x(4) = ixe
a1x(5) = ixw
a1x(6) = ixw
a1x(7) = ixw
a1x(8) = ixs     
!------
a1y(1) = iyn
a1y(2) = iyn
a1y(3) = iy
a1y(4) = iys
a1y(5) = iyn
a1y(6) = iy
a1y(7) = iys
a1y(8) = iys     
!------
return
END SUBROUTINE mk_8gridsxy


!*********************************************************
SUBROUTINE del_front_lesseq_ngrids_wgt(a2in, a1lat, miss, thnum, nx, ny, a2out)
implicit none
!--- in ---------
integer                     nx, ny
integer                     thnum
!f2py intent(in)            thnum
real,dimension(nx,ny)    :: a2in
!f2py intent(in)            a2in
real,dimension(ny)       :: a1lat
!f2py intent(in)            a1lat
real                        miss
!f2py intent(in)            miss
!--- out --------
real,dimension(nx,ny)    :: a2out
!f2py intent(out)           a2out
!--- calc -------
integer                     ix,iy,ik
integer                     icount_surr
integer                     iix,iiy
integer                     id, id_min, ik_tmp
integer                     stopflag
integer                     idMAX, idMIN, MAXmin, MINmin
integer,dimension(8)     :: a1x, a1y
integer,dimension(nx,ny) :: a2id
integer,dimension(nx*ny) :: a1connect
real,dimension(nx*ny)    :: a1num
real,dimension(nx,ny)    :: a2in_tmp
real,dimension(ny)       :: a1ratio
real                        rat
!--- parameter --
integer,parameter        :: miss_int = -9999

!-- initialize ---
a2in_tmp  = a2in
a1connect = miss_int
a1num     = 0
a2id      = miss_int
id        = 0
!----------------
do iy = 1,ny
  do ix = 1,nx
    if (a2in(ix,iy).ne.miss)then
      !-- 1st search ----
      call mk_8gridsxy(nx, ny, ix, iy, a1x, a1y)
      !*****************************
      ! if (ix,iy) has not been given an ID
      !------------
      if (a2id(ix,iy).eq.miss_int)then
        id = id + 1
        a2id(ix,iy) = id 
      end if
      !***************************
      !----------------------------
      icount_surr = 0
      do ik = 1,8
        iix = a1x(ik)
        iiy = a1y(ik)
        if (a2in_tmp(iix,iiy).ne.miss)then
          icount_surr   = icount_surr + 1
          !--- when adjacent id is -9999
          if (a2id(iix,iiy).eq.miss_int)then
            a2id(iix,iiy) = a2id(ix,iy)
          !--- when adjacent id is not miss
          else if (a2id(iix,iiy).ne.miss_int)then
            if (a2id(ix,iy).eq.a2id(iix,iiy))then
              continue
            else
              idMIN = min(a2id(ix,iy), a2id(iix,iiy))
              idMAX = max(a2id(ix,iy), a2id(iix,iiy))
              MINmin = a1connect(idMIN)
              MAXmin = a1connect(idMAX)
              if ((MINmin.eq.miss_int).and.(MAXmin.eq.miss_int))then
                a1connect(idMAX) = idMIN
              else if ((MINmin.eq.miss_int).and.(MAXmin.ne.miss_int))then
                if (MAXmin.lt.idMIN)then
                  a1connect(idMIN) = MAXmin
                else if (idMIN .le. MAXmin) then
                  a1connect(idMAX) = idMIN
                end if
              else if ((MINmin.ne.miss_int).and.(MAXmin.eq.miss_int))then
                a1connect(idMAX)=MINmin
              else if (MAXmin.eq.MINmin)then
                continue
              else if (MAXmin.gt.MINmin)then
                a1connect(idMAX) = MINmin
              else if (MINmin.gt.MAXmin)then
                a1connect(idMIN) = MAXmin
              end if
            end if
          end if
        end if        
      end do
      if (icount_surr.eq. 0)then
        a2in_tmp(ix,iy) = miss
        a2id(ix,iy)     = miss_int
        id = id -1
      end if
    end if
  end do
end do
!-----
do ik = 1,nx*ny
  if (a1connect(ik).ne.miss_int)then
    stopflag   = 0
    ik_tmp     = ik
    id_min     = a1connect(ik_tmp)
    do while (stopflag .ne.1)
      if (a1connect(ik_tmp).eq.miss_int)then
        stopflag = 1
      else
        if (a1connect(ik_tmp).lt.id_min)then
          id_min = a1connect(ik_tmp)
        else
          continue
        end if
        ik_tmp = a1connect(ik_tmp)
      end if
    end do 
    a1connect(ik) = id_min
  end if
end do
!-----
do iy = 1, ny
  do ix = 1,nx
    if (a2id(ix,iy).ne.miss_int)then
      ik = a2id(ix,iy)
      if (a1connect(ik).ne.miss_int)then
        a2id(ix,iy) = a1connect(ik)
      else
        continue
      end if
    end if
  end do
end do
!********************************
! count id
!-----------------------
CALL mk_a1arearat_sphere(a1lat, ny, a1ratio)
do iy = 1,ny
  rat   = a1ratio(iy)
  do ix = 1,nx
    if (a2id(ix,iy).ne.miss_int)then
      id = a2id(ix,iy)
      a1num(id) = a1num(id) + 1 * rat
    end if
  end do
end do
!********************************
! remove small segments
!-----------------------
a2out = a2in_tmp
do iy = 1,ny
  do ix = 1,nx
    id = a2id(ix,iy)
    if (id.eq.miss_int)then
      cycle
    end if
    if (a1num(id).le.real(thnum))then
      a2out(ix,iy) = miss
    end if
  end do
end do
!--------------------------------
!a2out = real(a2id)
return
END SUBROUTINE del_front_lesseq_ngrids_wgt

!*********************************************************
SUBROUTINE del_front_lesseq_ngrids(a2in, miss, thnum, nx, ny, a2out)
implicit none
!--- in ---------
integer                     nx, ny
integer                     thnum
!f2py intent(in)            thnum
real,dimension(nx,ny)    :: a2in
!f2py intent(in)            a2in
real                        miss
!f2py intent(in)            miss
!--- out --------
real,dimension(nx,ny)    :: a2out
!f2py intent(out)           a2out
!--- calc -------
integer                     ix,iy,ik
integer                     icount_surr
integer                     iix,iiy
integer                     id, id_min, ik_tmp
integer                     stopflag
integer                     idMAX, idMIN, MAXmin, MINmin
integer,dimension(8)     :: a1x, a1y
integer,dimension(nx,ny) :: a2id
integer,dimension(nx*ny) :: a1connect
real,dimension(nx*ny)    :: a1num
real,dimension(nx,ny)    :: a2in_tmp
!--- parameter --
integer,parameter        :: miss_int = -9999

!-- initialize ---
a2in_tmp  = a2in
a1connect = miss_int
a1num     = 0
a2id      = miss_int
id        = 0
!----------------
do iy = 1,ny
  do ix = 1,nx
    if (a2in(ix,iy).ne.miss)then
      !-- 1st search ----
      call mk_8gridsxy(nx, ny, ix, iy, a1x, a1y)
      !*****************************
      ! if (ix,iy) has not been given an ID
      !------------
      if (a2id(ix,iy).eq.miss_int)then
        id = id + 1
        a2id(ix,iy) = id 
      end if
      !***************************
      !----------------------------
      icount_surr = 0
      do ik = 1,8
        iix = a1x(ik)
        iiy = a1y(ik)
        if (a2in_tmp(iix,iiy).ne.miss)then
          icount_surr   = icount_surr + 1
          !--- when adjacent id is -9999
          if (a2id(iix,iiy).eq.miss_int)then
            a2id(iix,iiy) = a2id(ix,iy)
          !--- when adjacent id is not miss
          else if (a2id(iix,iiy).ne.miss_int)then
            if (a2id(ix,iy).eq.a2id(iix,iiy))then
              continue
            else
              idMIN = min(a2id(ix,iy), a2id(iix,iiy))
              idMAX = max(a2id(ix,iy), a2id(iix,iiy))
              MINmin = a1connect(idMIN)
              MAXmin = a1connect(idMAX)
              if ((MINmin.eq.miss_int).and.(MAXmin.eq.miss_int))then
                a1connect(idMAX) = idMIN
              else if ((MINmin.eq.miss_int).and.(MAXmin.ne.miss_int))then
                if (MAXmin.lt.idMIN)then
                  a1connect(idMIN) = MAXmin
                else if (idMIN .le. MAXmin) then
                  a1connect(idMAX) = idMIN
                end if
              else if ((MINmin.ne.miss_int).and.(MAXmin.eq.miss_int))then
                a1connect(idMAX)=MINmin
              else if (MAXmin.eq.MINmin)then
                continue
              else if (MAXmin.gt.MINmin)then
                a1connect(idMAX) = MINmin
              else if (MINmin.gt.MAXmin)then
                a1connect(idMIN) = MAXmin
              end if
            end if
          end if
        end if        
      end do
      if (icount_surr.eq. 0)then
        a2in_tmp(ix,iy) = miss
        a2id(ix,iy)     = miss_int
        id = id -1
      end if
    end if
  end do
end do
!-----
do ik = 1,nx*ny
  if (a1connect(ik).ne.miss_int)then
    stopflag   = 0
    ik_tmp     = ik
    id_min     = a1connect(ik_tmp)
    do while (stopflag .ne.1)
      if (a1connect(ik_tmp).eq.miss_int)then
        stopflag = 1
      else
        if (a1connect(ik_tmp).lt.id_min)then
          id_min = a1connect(ik_tmp)
        else
          continue
        end if
        ik_tmp = a1connect(ik_tmp)
      end if
    end do 
    a1connect(ik) = id_min
  end if
end do
!-----
do iy = 1, ny
  do ix = 1,nx
    if (a2id(ix,iy).ne.miss_int)then
      ik = a2id(ix,iy)
      if (a1connect(ik).ne.miss_int)then
        a2id(ix,iy) = a1connect(ik)
      else
        continue
      end if
    end if
  end do
end do
!!********************************
!! count id
!!-----------------------
!do iy = 1,ny
!  do ix = 1,nx
!    if (a2id(ix,iy).ne.miss)then
!      id = a2id(ix,iy)
!      a1num(id) = a1num(id) + 1
!    end if
!  end do
!end do
!

!********************************
! count id
!-----------------------
do iy = 1,ny
  do ix = 1,nx
    if (a2id(ix,iy).ne.miss_int)then
      id = a2id(ix,iy)
      a1num(id) = a1num(id) + 1
    end if
  end do
end do
!********************************
! remove small segments
!-----------------------
a2out = a2in_tmp
do iy = 1,ny
  do ix = 1,nx
    id = a2id(ix,iy)
    if (id.eq.miss_int)then
      cycle
    end if
    if (a1num(id).le.real(thnum))then
      a2out(ix,iy) = miss
    end if
  end do
end do
!--------------------------------
!a2out = real(a2id)
return
END SUBROUTINE del_front_lesseq_ngrids
!
!*********************************************************



!!*********************************************************
SUBROUTINE mk_a2meanunitaxis_h98_apdx2_2(a2inx, a2iny, miss, nx, ny, a2meanunitaxis_x, a2meanunitaxis_y)
implicit none
!----------------
! return unit vector of angle of mean unit axis
! cos(beta), sin(beta) 
! use the graphical method in Appendix 2.2, Hewson 1998
! For front, a2inx = grad2_x,  a2iny = grad2_y
!
! Theoretically the same as the method of Appendix 2.1
!--- in ---------
integer                  nx, ny
real,dimension(nx,ny) :: a2inx, a2iny
!f2py intent(in)         a2inx, a2iny

real,optional         :: miss
!f2py intent(in)         miss
!!--- out --------
real,dimension(nx,ny) :: a2meanunitaxis_x, a2meanunitaxis_y
!f2py intent(out)        a2meanunitaxis_x, a2meanunitaxis_y
!--- calc -------
real                       vxo, vxn, vxs, vxw, vxe
real                       vyo, vyn, vys, vyw, vye
!
real                       vrxo, vrxn, vrxs, vrxw, vrxe
real                       vryo, vryn, vrys, vryw, vrye
!
real                       vrx, vry, radr
!
real                       coso, cosn, coss, cosw, cose
real                       sino, sinn, sins, sinw, sine
!
real                       cosr, sinr, cosr_hlf, sinr_hlf
!
integer                    ix,  iy
integer                    ixn, ixs, ixw, ixe
integer                    iyn, iys, iyw, iye
!----------------
do iy = 1,ny
  do ix = 1,nx
    !---
    call ixy2iixy(nx, ny, ix, iy+1, ixn, iyn)
    call ixy2iixy(nx, ny, ix, iy-1, ixs, iys)
    call ixy2iixy(nx, ny, ix-1, iy, ixw, iyw)
    call ixy2iixy(nx, ny, ix+1, iy, ixe, iye)
    !---
    vxo = a2inx(ix,  iy)
    vyo = a2iny(ix,  iy)
    vxn = a2inx(ixn, iyn)
    vxs = a2inx(ixs, iys)
    vxw = a2inx(ixw, iyw)
    vxe = a2inx(ixe, iye)
    vyn = a2iny(ixn, iyn)
    vys = a2iny(ixs, iys)
    vyw = a2iny(ixw, iyw)
    vye = a2iny(ixe, iye)
    !---
    if ((vxo .eq. 0.0) .and. (vyo .eq. 0.0))then
      coso = 0.0
      sino = 0.0
    else
      coso = vxo / (vxo**2.0 + vyo**2.0)**0.5
      sino = vyo / (vxo**2.0 + vyo**2.0)**0.5
    endif
    if ((vxn .eq. 0.0) .and. (vyn .eq. 0.0))then
      cosn = 0.0
      sinn = 0.0
    else
      cosn = vxn / (vxn**2.0 + vyn**2.0)**0.5
      sinn = vyn / (vxn**2.0 + vyn**2.0)**0.5
    endif
    if ((vxs .eq. 0.0) .and. (vys .eq. 0.0))then
      coss = 0.0
      sins = 0.0
    else
      coss = vxs / (vxs**2.0 + vys**2.0)**0.5
      sins = vys / (vxs**2.0 + vys**2.0)**0.5
    endif
    if ((vxw .eq. 0.0) .and. (vyw .eq. 0.0))then
      cosw = 0.0
      sinw = 0.0
    else
      cosw = vxw / (vxw**2.0 + vyw**2.0)**0.5
      sinw = vyw / (vxw**2.0 + vyw**2.0)**0.5
    endif
    if ((vxe .eq. 0.0) .and. (vye .eq. 0.0))then
      cose = 0.0
      sine = 0.0
    else
      cose = vxe / (vxe**2.0 + vye**2.0)**0.5
      sine = vye / (vxe**2.0 + vye**2.0)**0.5
    endif

    !---
    vrxo = coso*vxo - sino*vyo
    vrxn = cosn*vxn - sinn*vyn
    vrxs = coss*vxs - sins*vys
    vrxw = cosw*vxw - sinw*vyw
    vrxe = cose*vxe - sine*vye
    !---
    vryo = sino*vxo + coso*vyo
    vryn = sinn*vxn + cosn*vyn
    vrys = sins*vxs + coss*vys
    vryw = sinw*vxw + cosw*vyw
    vrye = sine*vxe + cose*vye
    !---
    vrx  = vrxo + vrxn + vrxs + vrxw + vrxe 
    vry  = vryo + vryn + vrys + vryw + vrye
    !---
    radr = (vrx**2.0 + vry**2.0)**0.5
    !---
    if (radr .gt. 0.0)then
      cosr = vrx / radr 
      sinr = vry / radr
      !---
      cosr_hlf = ((1+cosr)*0.5)**0.5 * sign(1.0, sinr)
      sinr_hlf = ((1-cosr)*0.5)**0.5
      !---
      !a2meanunitaxis_x(ix,iy) = 0.2* radr * cosr_hlf
      !a2meanunitaxis_y(ix,iy) = 0.2* radr * sinr_hlf
      a2meanunitaxis_x(ix,iy) = cosr_hlf
      a2meanunitaxis_y(ix,iy) = sinr_hlf
    else
      a2meanunitaxis_x(ix,iy) = 0.0
      a2meanunitaxis_y(ix,iy) = 0.0
    endif
    !******************
    ! check miss
    !------------------
    if ( present(miss) )then
      if (& 
        & (vxo.eq.miss).or.&
        & (vyo.eq.miss).or.&
        & (vxn.eq.miss).or.&
        & (vxs.eq.miss).or.&
        & (vxw.eq.miss).or.&
        & (vxe.eq.miss).or.&
        & (vyn.eq.miss).or.&
        & (vys.eq.miss).or.&
        & (vyw.eq.miss).or.&
        & (vye.eq.miss)) then
        a2meanunitaxis_x(ix,iy) = miss
        a2meanunitaxis_y(ix,iy) = miss
      end if 
    end if
    !******************
  end do
end do

return
END SUBROUTINE mk_a2meanunitaxis_h98_apdx2_2

!*********************************************************
SUBROUTINE fill_front_gap_trace(a2in, a2trace, miss, nx, ny, a2out)
implicit none
integer                  nx, ny
real,dimension(nx,ny) :: a2in, a2trace
!f2py intent(in)         a2in, a2trace
real                     miss
!f2py intent(in)         miss
!---- out -------
real,dimension(nx,ny) :: a2out
!f2py intent(out)        a2out
!---- calc ------
integer                  ix,iy, i
integer                  icount, icountn, icounts, icountw, icounte
integer                  iix,iiy
integer,dimension(8)  :: a1x, a1y
!---- init ------
a2out     = a2in
!----------------
do iy = 1,ny
  do ix = 2,nx-1
    if (a2trace(ix,iy).ne.miss)then
      if (a2in(ix,iy).eq.miss)then  ! when a2in==miss
        call mk_8gridsxy(nx, ny, ix, iy, a1x, a1y)
        icount = 0
        icountn= 0
        icounts= 0
        icounte= 0
        icountw= 0
        do i = 1,8
          iix = a1x(i)
          iiy = a1y(i)
          if (a2in(iix,iiy).ne.miss)then
            icount = icount + 1
            !--------------
            if (iiy.gt.iy) icountn = icountn+1
            if (iiy.lt.iy) icounts = icounts+1
            if (iix.gt.ix) icounte = icounte+1
            if (iix.lt.ix) icountw = icountw+1
            !--------------
          end if
        end do
        if ((icount.ge.2).and.(icount.le.4))then
          if (max(icountn, icounts, icounte, icountw).lt.3)then
            a2out(ix,iy) = a2trace(ix,iy)
          end if
        end if
      end if 
    end if
  end do
  !-----
  do ix = 1, nx, nx-1
    if (a2trace(ix,iy).ne.miss)then
      if (a2in(ix,iy).eq.miss)then  ! when a2in==miss
        call mk_8gridsxy(nx, ny, ix, iy, a1x, a1y)
        icount = 0
        icountn= 0
        icounts= 0
        icounte= 0
        icountw= 0
        do i = 1,8
          iix = a1x(i)
          iiy = a1y(i)
          if (a2in(iix,iiy).ne.miss)then
            icount = icount + 1
            !--------------
            if (iiy.gt.iy  ) icountn = icountn+1
            if (iiy.lt.iy  ) icounts = icounts+1
            !-- cases for ix=1  --
            if (ix.eq.1) then
              if (iix.eq.2   ) icounte = icounte+1
              if (iix.eq.nx  ) icountw = icountw+1
            end if
            !-- cases for ix=nx --
            if (ix.eq.nx) then
              if (iix.eq.1   ) icounte = icounte+1
              if (iix.eq.nx-1) icountw = icountw+1
            end if
            !--------------
          end if
        end do
        if ((icount.ge.2).and.(icount.le.4))then
          if (max(icountn, icounts, icounte, icountw).lt.3)then
            a2out(ix,iy) = a2trace(ix,iy)
          end if
        end if
      end if 
    end if
  end do
end do
!------------
return
END SUBROUTINE fill_front_gap_trace

!*********************************************************
SUBROUTINE fill_front_gap(a2in, miss, nx, ny, a2out)
implicit none
integer                  nx, ny
real,dimension(nx,ny) :: a2in
!f2py intent(in)         a2in
real                     miss
!f2py intent(in)         miss
!---- out -------
real,dimension(nx,ny) :: a2out
!f2py intent(out)        a2out
!---- calc ------
integer                  ix,iy
integer                  iix,iiy, imx,imy
integer                  dx,dy
!---- init ------
a2out     = a2in
!----------------
do iy = 1,ny
  do ix = 1,nx
    !print '(E10.2)',a2in(ix,iy)
    if (a2in(ix,iy).ne.miss)then
      do dy = -2,2,4
        call ixy2iixy(nx, ny, ix,iy+dy, iix, iiy)
        if (a2in(iix,iiy).ne.miss)then
          imx = int((ix + iix)*0.5)
          imy = int((iy + iiy)*0.5)
          if (a2in(imx,imy).eq.miss)then
            a2out(imx,imy) = (a2in(ix,iy) + a2in(iix,iiy))*0.5
          end if
        end if
      end do
      do dx = -2,2,4
        call ixy2iixy(nx, ny, ix+dx,iy, iix, iiy)
        if (a2in(iix,iiy).ne.miss)then
          !print *,iix,iiy,a2in(iix,iiy)
          imx = int((ix + iix)*0.5)
          imy = int((iy + iiy)*0.5)
          if (a2in(imx,imy).eq.miss)then
            a2out(imx,imy) = (a2in(ix,iy) + a2in(iix,iiy))*0.5
            !print *, imx,imy,a2in(ix,iy), a2out(imx,imy), a2in(iix,iiy)
          end if
        end if
      end do
    end if
  end do
end do
!------------
return
END SUBROUTINE fill_front_gap


!*********************************************************
SUBROUTINE mk_a2contour(a2in, v, vtrue_out, miss, nx, ny, a2contour)
implicit none
!--- in -------
integer                  nx, ny
real,dimension(nx,ny) :: a2in
!f2py intent(in)         a2in
real                     v, vtrue_out, miss
!f2py intent(in)         v, vtrue_out, miss
!--- out ------
real,dimension(nx,ny) :: a2contour
!f2py intent(out)        a2contour
!--- calc -----
integer                  k
integer                  ix, iy
integer                  ixn, ixs, ixw, ixe
integer                  iyn, iys, iyw, iye
real                     vn, vs, vw, ve
real                     vnw, vne, vsw, vse
real                     vo, vmid
real,dimension(8)     :: a1v
!--------------
a2contour = miss
a1v  = miss

!-----------------------
! iy = 1 & ny
!-----------------------
do iy = 1,ny, ny-1
  do ix = 1,nx
    !---
    call ixy2iixy(nx, ny, ix, iy+1, ixn, iyn)
    call ixy2iixy(nx, ny, ix, iy-1, ixs, iys)
    call ixy2iixy(nx, ny, ix-1, iy, ixw, iyw)
    call ixy2iixy(nx, ny, ix+1, iy, ixe, iye)
    !---
    vo  = a2in(ix, iy)
    vn  = a2in(ixn, iyn)
    vs  = a2in(ixs, iys)
    vw  = a2in(ixw, iyw)
    ve  = a2in(ixe, iye)
    vnw = a2in(ixw, iyn)
    vne = a2in(ixe, iyn)
    vsw = a2in(ixw, iys)
    vse = a2in(ixe, iys)
    a1v(1) = vnw
    a1v(2) = vn
    a1v(3) = vne
    a1v(4) = vw
    a1v(5) = ve
    a1v(6) = vsw
    a1v(7) = vs
    a1v(8) = vse
    !---
    if (vo .le. v)then
      do k = 1,8
        if (a1v(k) .ge. v)then
          vmid = (a1v(k) + v)*0.5
          if ( v .lt. vmid) then
            a2contour(ix,iy) = vtrue_out
            exit
          end if
        end if
      end do
    else if (vo .gt. v)then
      do k = 1,8
        if (a1v(k) .ge. v)then
          vmid = (a1v(k) + v)*0.5
          if ( v .ge. vmid) then
            a2contour(ix,iy) = vtrue_out
            exit
          end if
        end if
      end do
    end if
    !---
  end do
end do

!-----------------------
! iy = 2 & ny-1
! ix = 1 & nx
!-----------------------
do iy = 2,ny-1
  do ix = 1,nx, nx-1
    !---
    call ixy2iixy(nx, ny, ix, iy+1, ixn, iyn)
    call ixy2iixy(nx, ny, ix, iy-1, ixs, iys)
    call ixy2iixy(nx, ny, ix-1, iy, ixw, iyw)
    call ixy2iixy(nx, ny, ix+1, iy, ixe, iye)
    !---
    vo  = a2in(ix, iy)
    vn  = a2in(ixn, iyn)
    vs  = a2in(ixs, iys)
    vw  = a2in(ixw, iyw)
    ve  = a2in(ixe, iye)
    vnw = a2in(ixw, iyn)
    vne = a2in(ixe, iyn)
    vsw = a2in(ixw, iys)
    vse = a2in(ixe, iys)
    a1v(1) = vnw
    a1v(2) = vn
    a1v(3) = vne
    a1v(4) = vw
    a1v(5) = ve
    a1v(6) = vsw
    a1v(7) = vs
    a1v(8) = vse
    !---
    if (vo .le. v)then
      do k = 1,8
        if (a1v(k) .ge. v)then
          vmid = (a1v(k) + v)*0.5
          if ( v .lt. vmid) then
            a2contour(ix,iy) = vtrue_out
            exit
          end if
        end if
      end do
    else if (vo .gt. v)then
      do k = 1,8
        if (a1v(k) .ge. v)then
          vmid = (a1v(k) + v)*0.5
          if ( v .ge. vmid) then
            a2contour(ix,iy) = vtrue_out
            exit
          end if
        end if
      end do
    end if
    !---
  end do
end do

!-----------------------
! iy = 2,3, ... , ny-1
! ix = 2,3, ... , nx-1
!-----------------------
do iy = 2,ny-1
  do ix = 2,nx-1
    !---
    vo  = a2in(ix, iy)
    vn  = a2in(ix, iy+1)
    vs  = a2in(ix, iy-1)
    vw  = a2in(ix-1, iy)
    ve  = a2in(ix+1, iy)
    vnw = a2in(ix-1, iy+1)
    vne = a2in(ix+1, iy+1)
    vsw = a2in(ix-1, iy-1)
    vse = a2in(ix+1, iy-1)
    a1v(1) = vnw
    a1v(2) = vn
    a1v(3) = vne
    a1v(4) = vw
    a1v(5) = ve
    a1v(6) = vsw
    a1v(7) = vs
    a1v(8) = vse
    !---
    if (vo .le. v)then
      do k = 1,8
        if (a1v(k) .ge. v)then
          vmid = (a1v(k) + v)*0.5
          if ( v .lt. vmid) then
            a2contour(ix,iy) = vtrue_out
            exit
          end if
        end if
      end do
    else if (vo .gt. v)then
      do k = 1,8
        if (a1v(k) .ge. v)then
          vmid = (a1v(k) + v)*0.5
          if ( v .ge. vmid) then
            a2contour(ix,iy) = vtrue_out
            exit
          end if
        end if
      end do
    end if
    !---
  end do
end do

return
END SUBROUTINE  mk_a2contour

!*********************************************************
SUBROUTINE old_mk_a2contour(a2in, v, vtrue_out, miss, nx, ny, a2contour)
implicit none
!--- in -------
integer                  nx, ny
real,dimension(nx,ny) :: a2in
!f2py intent(in)         a2in
real                     v, vtrue_out, miss
!f2py intent(in)         v, vtrue_out, miss
!--- out ------
real,dimension(nx,ny) :: a2contour
!f2py intent(out)        a2contour
!--- calc -----
integer                  k
integer                  ix, iy
integer                  ixn, ixs, ixw, ixe
integer                  iyn, iys, iyw, iye
real                     vn, vs, vw, ve
real                     vnw, vne, vsw, vse
real                     vo, vmid
real,dimension(8)     :: a1v
!--------------
a2contour = miss
a1v  = miss
!-----------------------
! iy = 1 & ny
! ix = 1 & nx
!-----------------------
do iy = 1,ny, ny-1
  do ix = 1,nx, nx-1
    !---
    call ixy2iixy(nx, ny, ix, iy+1, ixn, iyn)
    call ixy2iixy(nx, ny, ix, iy-1, ixs, iys)
    call ixy2iixy(nx, ny, ix-1, iy, ixw, iyw)
    call ixy2iixy(nx, ny, ix+1, iy, ixe, iye)
    !---
    vo  = a2in(ix, iy)
    vn  = a2in(ixn, iyn)
    vs  = a2in(ixs, iys)
    vw  = a2in(ixw, iyw)
    ve  = a2in(ixe, iye)
    vnw = a2in(ixw, iyn)
    vne = a2in(ixe, iyn)
    vsw = a2in(ixw, iys)
    vse = a2in(ixe, iys)
    a1v(1) = vnw
    a1v(2) = vn
    a1v(3) = vne
    a1v(4) = vw
    a1v(5) = ve
    a1v(6) = vsw
    a1v(7) = vs
    a1v(8) = vse
    !---
    if (vo .le. v)then
      do k = 1,8
        if (a1v(k) .ge. v)then
          vmid = (a1v(k) + v)*0.5
          if ( v .lt. vmid) then
            a2contour(ix,iy) = vtrue_out
            exit
          end if
        end if
      end do
    else if (vo .gt. v)then
      do k = 1,8
        if (a1v(k) .ge. v)then
          vmid = (a1v(k) + v)*0.5
          if ( v .ge. vmid) then
            a2contour(ix,iy) = vtrue_out
            exit
          end if
        end if
      end do
    end if
    !---
  end do
end do

!-----------------------
! iy = 2,3, ... , ny-1
! ix = 2,3, ... , nx-1
!-----------------------
do iy = 2,ny-1
  do ix = 2,nx-1
    !---
    vo  = a2in(ix, iy)
    vn  = a2in(ix, iy+1)
    vs  = a2in(ix, iy-1)
    vw  = a2in(ix-1, iy)
    ve  = a2in(ix+1, iy)
    vnw = a2in(ix-1, iy+1)
    vne = a2in(ix+1, iy+1)
    vsw = a2in(ix-1, iy-1)
    vse = a2in(ix+1, iy-1)
    a1v(1) = vnw
    a1v(2) = vn
    a1v(3) = vne
    a1v(4) = vw
    a1v(5) = ve
    a1v(6) = vsw
    a1v(7) = vs
    a1v(8) = vse
    !---
    if (vo .le. v)then
      do k = 1,8
        if (a1v(k) .ge. v)then
          vmid = (a1v(k) + v)*0.5
          if ( v .lt. vmid) then
            a2contour(ix,iy) = vtrue_out
            exit
          end if
        end if
      end do
    else if (vo .gt. v)then
      do k = 1,8
        if (a1v(k) .ge. v)then
          vmid = (a1v(k) + v)*0.5
          if ( v .ge. vmid) then
            a2contour(ix,iy) = vtrue_out
            exit
          end if
        end if
      end do
    end if
    !---
  end do
end do

return
END SUBROUTINE  old_mk_a2contour

!*********************************************************


!*********************************************************
SUBROUTINE mk_a2meanunitaxis_h98_apdx2_1(a2inx, a2iny, miss, nx, ny, a2meanunitaxis_x, a2meanunitaxis_y)
implicit none
!----------------
! return unit vector of angle of mean unit axis
! cos(beta), sin(beta) 
! use Equation A5 in Appendix 2.1, Hewson 1998
!
! For front, a2inx = grad2_x,  a2iny = grad2_y
!
! Theoretically the same as the method of Appendix 2.2
!--- in ---------
integer                  nx, ny
real,dimension(nx,ny) :: a2inx, a2iny
!f2py intent(in)         a2inx, a2iny

real,optional         :: miss
!f2py intent(in)         miss
!!--- out --------
real,dimension(nx,ny) :: a2meanunitaxis_x, a2meanunitaxis_y
!f2py intent(out)        a2meanunitaxis_x, a2meanunitaxis_y
!--- calc -------
real                       vxo, vxn, vxs, vxw, vxe
real                       vyo, vyn, vys, vyw, vye
!
real                       coso, cosn, coss, cosw, cose, cosm
real                       sino, sinn, sins, sinw, sine, sinm

real                       cos2o, cos2n, cos2s, cos2w, cos2e, cos2m
real                       sin2o, sin2n, sin2s, sin2w, sin2e, sin2m

!
real                       DDo, DDn, DDs, DDw, DDe
real                       P, Q
!
integer                    ix,  iy
integer                    ixn, ixs, ixw, ixe
integer                    iyn, iys, iyw, iye
!----------------
do iy = 1,ny
  do ix = 1,nx
    !---
    call ixy2iixy(nx, ny, ix, iy+1, ixn, iyn)
    call ixy2iixy(nx, ny, ix, iy-1, ixs, iys)
    call ixy2iixy(nx, ny, ix-1, iy, ixw, iyw)
    call ixy2iixy(nx, ny, ix+1, iy, ixe, iye)
    !---
    vxo = a2inx(ix,  iy)
    vyo = a2iny(ix,  iy)
    vxn = a2inx(ixn, iyn)
    vxs = a2inx(ixs, iys)
    vxw = a2inx(ixw, iyw)
    vxe = a2inx(ixe, iye)
    vyn = a2iny(ixn, iyn)
    vys = a2iny(ixs, iys)
    vyw = a2iny(ixw, iyw)
    vye = a2iny(ixe, iye)
    !---
    if ((vxo .eq. 0.0) .and. (vyo .eq. 0.0))then
      DDo  = 0.0
      coso = 0.0
      sino = 0.0
    else
      DDo   = sqrt(vxo**2.0 + vyo**2.0)
      coso = vxo / DDo
      sino = vyo / DDo
    endif
    if ((vxn .eq. 0.0) .and. (vyn .eq. 0.0))then
      DDn  = 0.0
      cosn = 0.0
      sinn = 0.0
    else
      DDn   = sqrt(vxn**2.0 + vyn**2.0)
      cosn = vxn / DDn 
      sinn = vyn / DDn 
    endif
    if ((vxs .eq. 0.0) .and. (vys .eq. 0.0))then
      DDs  = 0.0
      coss = 0.0
      sins = 0.0
    else
      DDs  = sqrt(vxs**2.0 + vys**2.0)
      coss = vxs / DDs
      sins = vys / DDs
    endif
    if ((vxw .eq. 0.0) .and. (vyw .eq. 0.0))then
      DDw  = 0.0
      cosw = 0.0
      sinw = 0.0
    else
      DDw  = sqrt(vxw**2.0 + vyw**2.0)
      cosw = vxw / DDw
      sinw = vyw / DDw
    endif
    if ((vxe .eq. 0.0) .and. (vye .eq. 0.0))then
      DDe  = 0.0
      cose = 0.0
      sine = 0.0
    else
      DDe  = sqrt(vxe**2.0 + vye**2.0)
      cose = vxe / DDe
      sine = vye / DDe
    endif
    !---
    sin2o  = 2.0*sino*coso
    sin2n  = 2.0*sinn*cosn
    sin2s  = 2.0*sins*coss
    sin2w  = 2.0*sinw*cosw
    sin2e  = 2.0*sine*cose

    cos2o  = 2.0*coso**2.0 - 1.0
    cos2n  = 2.0*cosn**2.0 - 1.0
    cos2s  = 2.0*coss**2.0 - 1.0
    cos2w  = 2.0*cosw**2.0 - 1.0
    cos2e  = 2.0*cose**2.0 - 1.0

    P    =   DDo*cos2o&
          & +DDn*cos2n&
          & +DDs*cos2s&
          & +DDw*cos2w&
          & +DDe*cos2e

    Q    =   DDo*sin2o&
          & +DDn*sin2n&
          & +DDs*sin2s&
          & +DDw*sin2w&
          & +DDe*sin2e

    !-----------------------
    ! cos(2*beta) = P/n / Dmean
    ! sin(2*beta) = Q/n / Dmean
    ! where
    !   Dmean = sqrt(P**2.0 + Q**2.0)/n
    ! then
    !   cos(2*beta) = P / sqrt( P**2.0 + Q**2.0 )
    !               = 1 / sqrt( 1 + (Q/P)**2.0 ) * sign(P)
    !
    !   sin(2*beta) = Q / sqrt( P**2.0 + Q**2.0 )
    !               = 1 / sqrt( 1 + (P/Q)**2.0 ) * sign(Q) 
    !-----------------------
    if (Q.eq.0.0)then
      cos2m = sign(1.0, P)
      sin2m = 0.0
    else if (P.eq.0.0)then
      cos2m = 0.0
      sin2m = sign(1.0, Q)
    else
      cos2m = 1.0/sqrt( 1.0+ (Q/P)**2.0 )*sign(1.0, P)
      sin2m = 1.0/sqrt( 1.0+ (P/Q)**2.0 )*sign(1.0, Q)
    end if  

    cosm = sqrt((1.0 + cos2m)*0.5)*sign(1.0, sin2m) 
    sinm = sqrt((1.0 - cos2m)*0.5)

    a2meanunitaxis_x(ix,iy) = cosm
    a2meanunitaxis_y(ix,iy) = sinm

    !******************
    ! check miss
    !------------------
    if ( present(miss) )then
      if (& 
        & (vxo.eq.miss).or.&
        & (vyo.eq.miss).or.&
        & (vxn.eq.miss).or.&
        & (vxs.eq.miss).or.&
        & (vxw.eq.miss).or.&
        & (vxe.eq.miss).or.&
        & (vyn.eq.miss).or.&
        & (vys.eq.miss).or.&
        & (vyw.eq.miss).or.&
        & (vye.eq.miss)) then
        a2meanunitaxis_x(ix,iy) = miss
        a2meanunitaxis_y(ix,iy) = miss
      end if 
    end if
    !******************
  end do
end do

return
END SUBROUTINE mk_a2meanunitaxis_h98_apdx2_1

!*********************************************************
SUBROUTINE mk_a2meanaxisgrad3_h98_eq6(a2tau, a1lon, a1lat, miss, nx, ny, a2meanaxisgrad3)
implicit none
!--------------
! Return Numerator of Hewson's (1998) Equation (6)
!--- in -------
integer                  nx, ny
real,dimension(nx,ny) :: a2tau
!f2py intent(in)         a2tau

real,dimension(nx)    :: a1lon
!f2py intent(in)         a1lon

real,dimension(ny)    :: a1lat
!f2py intent(in)         a1lat

real,optional         :: miss
!--- out ------
real,dimension(nx,ny) :: a2meanaxisgrad3
!f2py intent(out)        a2meanaxisgrad3
!!--- calc -----
integer                  ix,  iy
integer                  ixn, ixs, ixw, ixe
integer                  iyn, iys, iyw, iye
real                     xn, xs, xw, xe   ! resolved to mean axis
real                     yn, ys, yw, ye   ! resolved to mean axis
real                     rn, rs, rw, re
real                     dn, ds, dew, dns
real                     cosb, sinb

real,dimension(nx,ny) :: a2grad_x, a2grad_y
real,dimension(nx,ny) :: a2grad2_x, a2grad2_y
real,dimension(nx,ny) :: a2grad_abs
real,dimension(nx,ny) :: a2meanunitaxis_x, a2meanunitaxis_y
!!---------------
CALL mk_a2grad(a2in=a2tau,      a1lon=a1lon, a1lat=a1lat, miss=miss, nx=nx, ny=ny, a2gradx=a2grad_x, a2grady=a2grad_y)

CALL mk_a2grad_abs(a2in=a2tau,      a1lon=a1lon, a1lat=a1lat, miss=miss, nx=nx, ny=ny, a2gradabs=a2grad_abs)

CALL mk_a2grad(a2in=a2grad_abs, a1lon=a1lon, a1lat=a1lat, miss=miss, nx=nx, ny=ny, a2gradx=a2grad2_x, a2grady=a2grad2_y)
!
CALL mk_a2meanunitaxis_h98_apdx2_2(a2grad_x, a2grad_y, miss, nx, ny, a2meanunitaxis_x, a2meanunitaxis_y)
!!!--------------
do iy = 1,ny
  if ((iy.eq.1).or.(iy.eq.ny))then
    dn  = hubeny_real(a1lat(1), 0.0, a1lat(2), 0.0)
    dns = dn * 2.0
  else
    dn  = hubeny_real(a1lat(iy), 0.0, a1lat(iyn), 0.0)
    ds  = hubeny_real(a1lat(iy), 0.0, a1lat(iys), 0.0)
    dns = (dn + ds)
  end if
  dew = hubeny_real(a1lat(iy), a1lon(1), a1lat(iy), a1lon(2))

  do ix = 1, nx
    !---
    call ixy2iixy(nx, ny, ix, iy+1, ixn, iyn)
    call ixy2iixy(nx, ny, ix, iy-1, ixs, iys)
    call ixy2iixy(nx, ny, ix-1, iy, ixw, iyw)
    call ixy2iixy(nx, ny, ix+1, iy, ixe, iye)
    !---
    cosb = a2meanunitaxis_x(ix,iy)
    sinb = a2meanunitaxis_y(ix,iy)

    rn  = a2grad2_x(ixn,iyn)*cosb + a2grad2_y(ixn,iyn)*sinb
    xn  = rn*cosb 
    yn  = rn*sinb

    rs  = a2grad2_x(ixs,iys)*cosb + a2grad2_y(ixs,iys)*sinb
    xs  = rs*cosb 
    ys  = rs*sinb

    rw  = a2grad2_x(ixw,iyw)*cosb + a2grad2_y(ixw,iyw)*sinb
    xw  = rw*cosb 
    yw  = rw*sinb

    re  = a2grad2_x(ixe,iye)*cosb + a2grad2_y(ixe,iye)*sinb
    xe  = re*cosb 
    ye  = re*sinb

    a2meanaxisgrad3(ix,iy) = (xe - xw)/dew + (yn - ys)/dns
    if (     (a2grad2_x(ixw,iyw).eq.miss)&
        &.or.(a2grad2_x(ixe,iye).eq.miss)&
        &.or.(a2grad2_y(ixn,iyn).eq.miss)&
        &.or.(a2grad2_y(iys,iys).eq.miss)&
        &.or.(a2grad2_x(ix, iy ).eq.miss))then
    
      a2meanaxisgrad3(ix,iy) = miss
    end if
    !---
  end do
end do
!!--------------
!
return
END SUBROUTINE mk_a2meanaxisgrad3_h98_eq6



!*********************************************************
SUBROUTINE mk_a2axisgrad3_h98_eq7(a2tau, a1lon, a1lat, miss, nx, ny, a2axisgrad3)
implicit none
!--------------
! Return Numerator of Hewson's (1998) Equation (7)
!--- in -------
integer                  nx, ny
real,dimension(nx,ny) :: a2tau
!f2py intent(in)         a2tau

real,dimension(nx)    :: a1lon
!f2py intent(in)         a1lon

real,dimension(ny)    :: a1lat
!f2py intent(in)         a1lat

real,optional         :: miss
!--- out ------
real,dimension(nx,ny) :: a2axisgrad3
!f2py intent(out)        a2axisgrad3
!!--- calc -----
integer                  ix,  iy
real,dimension(nx,ny) :: a2grad_abs, a2grad2_abs
real,dimension(nx,ny) :: a2grad2_x, a2grad2_y
real,dimension(nx,ny) :: a2grad3_x, a2grad3_y
!!---------------

CALL mk_a2grad_abs(a2in=a2tau,      a1lon=a1lon, a1lat=a1lat, miss=miss, nx=nx, ny=ny, a2gradabs=a2grad_abs)

CALL mk_a2grad_abs(a2in=a2grad_abs, a1lon=a1lon, a1lat=a1lat, miss=miss, nx=nx, ny=ny, a2gradabs=a2grad2_abs)

CALL mk_a2grad(a2in=a2grad_abs, a1lon=a1lon, a1lat=a1lat, miss=miss, nx=nx, ny=ny, a2gradx=a2grad2_x, a2grady=a2grad2_y)

CALL mk_a2grad(a2in=a2grad2_abs, a1lon=a1lon, a1lat=a1lat, miss=miss, nx=nx, ny=ny, a2gradx=a2grad3_x, a2grady=a2grad3_y)

!!!--------------

a2axisgrad3 = a2grad3_x*a2grad2_x + a2grad3_y*a2grad2_y

do iy = 1,ny
  do ix = 1, nx
    if ((a2grad3_x(ix,iy).eq.miss).or.(a2grad3_y(ix,iy).eq.miss))then
      a2axisgrad3(ix,iy) = miss
    end if
  end do
end do
!!--------------
!
return
END SUBROUTINE mk_a2axisgrad3_h98_eq7

!*********************************************************
SUBROUTINE mk_a2frontmask2(a2tau, a1lon, a1lat, miss, nx, ny, a2frontmask2)
implicit none
!--------------
integer                  nx, ny
!--- in -------
real,dimension(nx,ny) :: a2tau
!f2py intent(in)         a2tau

real,dimension(nx)    :: a1lon
!f2py intent(in)         a1lon

real,dimension(ny)    :: a1lat
!f2py intent(in)         a1lat

real,optional         :: miss
!f2py intent(in)         miss
!--- out ------
real,dimension(nx,ny) :: a2frontmask2
!f2py intent(out)        a2frontmask2
!--- calc -----
integer                  ix,iy
real,dimension(nx,ny) :: a2gradtauabs, a2grad2tauabs
real                     dns, dew, meangridlen
real                     coef
real                     lat, dlon
!-----------------
if (present(miss))then
  CALL mk_a2grad_abs(a2in=a2tau,        a1lon=a1lon, a1lat=a1lat, miss=miss, nx=nx, ny=ny, a2gradabs=a2gradtauabs)
  CALL mk_a2grad_abs(a2in=a2gradtauabs, a1lon=a1lon, a1lat=a1lat, miss=miss, nx=nx, ny=ny, a2gradabs=a2grad2tauabs)
else
  CALL mk_a2grad_abs(a2in=a2tau,        a1lon=a1lon, a1lat=a1lat, nx=nx, ny=ny, a2gradabs=a2gradtauabs)
  CALL mk_a2grad_abs(a2in=a2gradtauabs, a1lon=a1lon, a1lat=a1lat, nx=nx, ny=ny, a2gradabs=a2grad2tauabs)
end if


dlon  = a1lon(2) - a1lon(1)
!-----------------
! iy = 1 & ny
!-----------------
do iy = 1, ny, ny-1
  lat = a1lat(iy)
  dns = hubeny_real(a1lat(1), 0.0, a1lat(2), 0.0)
  dew = hubeny_real(lat, 0.0, lat, dlon)
  meangridlen  = (dns + dew)*0.5
  coef         = meangridlen / (2**0.5)
  a2frontmask2(:,iy) = a2gradtauabs(:,iy) + coef*a2grad2tauabs(:,iy)
 
end do
!-----------------
do iy = 2,ny-1
  lat = a1lat(iy)
  dns = hubeny_real(a1lat(iy-1), 0.0, a1lat(iy+1), 0.0) *0.5
  dew = hubeny_real(lat, 0.0, lat, dlon)
  meangridlen  = (dns + dew)*0.5
  coef         = meangridlen / (2**0.5)
  a2frontmask2(:,iy) = a2gradtauabs(:,iy) + coef*a2grad2tauabs(:,iy)
end do
!-------------
! check miss
!-------------
if (present(miss))then
  do iy=1,ny
    do ix=1,nx
      if ((a2gradtauabs(ix,iy).eq.miss).or.(a2grad2tauabs(ix,iy).eq.miss))then
        a2frontmask2(ix,iy) = miss
      end if
    end do
  end do
end if
!-------------
return
!--------------
END SUBROUTINE mk_a2frontmask2

!*********************************************************
SUBROUTINE mk_a2frontmask1(a2tau, a1lon, a1lat, miss, nx, ny, a2frontmask1)
!--------------
! return left hand side of equation (9),  Hewson 1998
!--------------
implicit none
!--------------
integer                  nx, ny
!--- in -------
real,dimension(nx,ny) :: a2tau
!f2py intent(in)         a2tau

real,dimension(nx)    :: a1lon
!f2py intent(in)         a1lon

real,dimension(ny)    :: a1lat
!f2py intent(in)         a1lat

real,optional         :: miss
!f2py intent(in)         miss
!--- out ------
real,dimension(nx,ny) :: a2frontmask1
!f2py intent(out)        a2frontmask1
!--- calc -----
integer                  ix, iy
integer                  ixn, ixs, ixw, ixe
integer                  iyn, iys, iyw, iye
real                     vn, vs, vw, ve
real,dimension(nx,ny) :: a2gradtaux, a2gradtauy, a2gradtauabs
real,dimension(nx,ny) :: a2grad2taux, a2grad2tauy
real,dimension(nx,ny) :: a2frontmasksingle


!--------------
! For with miss
!--------------
if ( present(miss) )then
  CALL mk_a2grad(a2in=a2tau, a1lon=a1lon, a1lat=a1lat, miss=miss, nx=nx, ny=ny, a2gradx=a2gradtaux, a2grady=a2gradtauy)

!  CALL mk_a2grad_abs(a2in=a2tau, a1lon=a1lon, a1lat=a1lat, miss=miss, nx=nx, ny=ny, a2gradabs=a2gradtauabs)
  a2gradtauabs       = (a2gradtaux**2.0 + a2gradtauy**2.0)**0.5

  CALL mk_a2grad(a2in=a2gradtauabs, a1lon=a1lon, a1lat=a1lat,  miss=miss, nx=nx, ny=ny, a2gradx=a2grad2taux, a2grady=a2grad2tauy)
  a2frontmasksingle  = (-a2grad2taux*a2gradtaux - a2grad2tauy*a2gradtauy)/a2gradtauabs
  do iy=1,ny
    do ix=1,nx
      if ((a2gradtaux(ix,iy).eq.miss).or.(a2gradtauy(ix,iy).eq.miss)&
         &.or.(a2grad2taux(ix,iy).eq.miss).or.(a2grad2tauy(ix,iy).eq.miss))then
        a2frontmasksingle(ix,iy) = miss
      end if
    end do
  end do
!---------------
! For without miss
!---------------
else 
  CALL mk_a2grad(a2in=a2tau, a1lon=a1lon, a1lat=a1lat, nx=nx, ny=ny,  a2gradx=a2gradtaux, a2grady=a2gradtauy)
  a2gradtauabs       = (a2gradtaux**2.0 + a2gradtauy**2.0)**0.5
  CALL mk_a2grad(a2in=a2gradtauabs, a1lon=a1lon, a1lat=a1lat,  nx=nx, ny=ny, a2gradx=a2grad2taux, a2grady=a2grad2tauy)
  a2frontmasksingle  = (-a2grad2taux*a2gradtaux - a2grad2tauy*a2gradtauy)/a2gradtauabs
end if
!---------------
! iy = 1, & ny
!---------------
do iy = 1, ny, ny-1
  do ix = 1, nx
    if (a2frontmasksingle(ix,iy).eq.miss)then
      a2frontmask1(ix,iy) = miss
    else
      !---
      call ixy2iixy(nx, ny, ix, iy+1, ixn, iyn)
      call ixy2iixy(nx, ny, ix, iy-1, ixs, iys)
      call ixy2iixy(nx, ny, ix-1, iy, ixw, iyw)
      call ixy2iixy(nx, ny, ix+1, iy, ixe, iye)
      !---
      vn = a2frontmasksingle(ixn, iyn)
      vs = a2frontmasksingle(ixs, iys)
      vw = a2frontmasksingle(ixw, iyw)
      ve = a2frontmasksingle(ixe, iye)
      !---
      a2frontmask1(ix,iy) = (vn+vs+vw+ve)/4.0
      if ( isnan(a2frontmask1(ix,iy)) )then
        a2frontmask1(ix,iy) = 0.0
      end if
      !--- check miss ------
      if ( present(miss) )then
        if ( (vn.eq.miss).or.(vs.eq.miss).or.(vw.eq.miss).or.(ve.eq.miss))then
          a2frontmask1(ix,iy) = miss
        end if
      end if
      !---------------------
    end if
  end do
end do
!---------------
! iy = 2,3, ... ,ny-1
!---------------
do iy = 2,ny-1
  !-------------
  ! ix = 1 & nx
  !-------------  
  do ix = 1,nx,nx-1
    if (a2frontmasksingle(ix,iy).eq.miss)then
      a2frontmask1(ix,iy) = miss
    else
      !---
      call ixy2iixy(nx, ny, ix-1, iy, ixw, iyw)
      call ixy2iixy(nx, ny, ix+1, iy, ixe, iye)
      !---
      vn = a2frontmasksingle(ix, iy+1)
      vs = a2frontmasksingle(ix, iy-1)
      vw = a2frontmasksingle(ixw, iyw)
      ve = a2frontmasksingle(ixe, iye)
      !---
      a2frontmask1(ix,iy) = (vn+vs+vw+ve)/4.0
      if ( isnan(a2frontmask1(ix,iy)) )then
        a2frontmask1(ix,iy) = 0.0
      end if
      !--- check miss ------
      if ( present(miss) )then
        if ( (vn.eq.miss).or.(vs.eq.miss).or.(vw.eq.miss).or.(ve.eq.miss))then
          a2frontmask1(ix,iy) = miss
        end if
      end if
      !---------------------
    end if
  end do

  !-------------
  ! ix = 2,3, ... , nx-1
  !-------------  
  do ix = 2,nx-1
    if (a2frontmasksingle(ix,iy).eq.miss)then
      a2frontmask1(ix,iy) = miss
    else
      !---
      vn = a2frontmasksingle(ix, iy+1)
      vs = a2frontmasksingle(ix, iy-1)
      vw = a2frontmasksingle(ix-1, iy)
      ve = a2frontmasksingle(ix+1, iy)
      !---
      a2frontmask1(ix,iy) = (vn+vs+vw+ve)/4.0
      if ( isnan(a2frontmask1(ix,iy)) )then
        a2frontmask1(ix,iy) = 0.0
      end if
      !--- check miss ------
      if ( present(miss) )then
        if ( (vn.eq.miss).or.(vs.eq.miss).or.(vw.eq.miss).or.(ve.eq.miss))then
          a2frontmask1(ix,iy) = miss
        end if
      end if
      !---------------------
    end if
  end do
end do

return
END SUBROUTINE mk_a2frontmask1

!*********************************************************
SUBROUTINE mk_a2grad(a2in, a1lon, a1lat, miss, nx, ny, a2gradx, a2grady)
!---------------------------------
! data order should be South -> North, West -> East
! returns two vector map (map of da/dx, map of da/dy)
!---------------------------------
implicit none
!-----------------
integer                    nx
integer                    ny
!--- in ----------
real,dimension(nx,ny)   :: a2in
!f2py intent(in)           a2in

real,dimension(nx)      :: a1lon
!f2py intent(in)           a1lon

real,dimension(ny)      :: a1lat
!f2py intent(in)           a1lon

real,optional           :: miss
!f2py intent(in)           miss
!--- out ---------
real,dimension(nx,ny)   :: a2gradx, a2grady
!f2py intent(out)          a2gradx, a2grady
!--- calc --------
real                       lat_first
real                       dn, ds, dew
real                       vn, vs, vw, ve
real                       lat
integer                    ix,  iy
integer                    ixn, ixs, ixw, ixe
integer                    iyn, iys, iyw, iye
!-----------------
lat_first  = a1lat(1)

!-------------------
! ny = 1 & ny
!-------------------
do iy = 1, ny, ny-1  ! iy = 1 & ny
  call ixy2iixy(nx, ny, ix, iy+1, ixn, iyn)
  call ixy2iixy(nx, ny, ix, iy-1, ixs, iys)
  dn  = hubeny_real(a1lat(iy), 0.0, a1lat(iyn), 0.0)
  ds  = hubeny_real(a1lat(iy), 0.0, a1lat(iys), 0.0)
  dew = hubeny_real(a1lat(iy), a1lon(1), a1lat(iy), a1lon(2))

  do ix = 1, nx
    !---
    if (a2in(ix,iy).eq.miss)then
      a2gradx(ix,iy) = miss 
      a2grady(ix,iy) = miss 
    else
      !---
      call ixy2iixy(nx, ny, ix, iy+1, ixn, iyn)
      call ixy2iixy(nx, ny, ix, iy-1, ixs, iys)
      call ixy2iixy(nx, ny, ix-1, iy, ixw, iyw)
      call ixy2iixy(nx, ny, ix+1, iy, ixe, iye)
      !---
      vn = a2in(ixn, iyn)
      vs = a2in(ixs, iys)
      vw = a2in(ixw, iyw)
      ve = a2in(ixe, iye)
      !---
      a2gradx(ix,iy) = (ve - vw) / (2.0*dew)
      a2grady(ix,iy) = (vn - vs) / (dn + ds)
      !--- check miss ------
      if ((vn.eq.miss).or.(vs.eq.miss).or.(vw.eq.miss).or.(ve.eq.miss))then
        a2gradx(ix,iy) = miss
        a2grady(ix,iy) = miss
      end if
    end if
  end do
end do

!-----------------
! iy = 2,3,..ny-1
!-----------------
do iy = 2, ny-1
  lat = a1lat(iy)
  dn  = hubeny_real(lat, 0.0, a1lat(iy+1), 0.0)
  ds  = hubeny_real(lat, 0.0, a1lat(iy-1), 0.0)
  dew = hubeny_real(lat, a1lon(1), lat, a1lon(2))
  !-------------------
  ! ix = 1 & nx
  !-------------------
  do ix = 1, nx, nx-1  ! ix = 1 & nx
    if (a2in(ix,iy).eq.miss)then
      a2gradx(ix,iy) = miss 
      a2grady(ix,iy) = miss 
    else
      call ixy2iixy(nx, ny, ix-1, iy, ixw, iyw)
      call ixy2iixy(nx, ny, ix+1, iy, ixe, iye)
      !---
      vn = a2in(ix, iy+1)
      vs = a2in(ix, iy-1)
      vw = a2in(ixw, iy)
      ve = a2in(ixe, iy)
      !---
      a2gradx(ix, iy) = (ve - vw) / (2.0*dew)
      a2grady(ix, iy) = (vn - vs) / (dn + ds)
      !--- check miss ------
      if ((vn.eq.miss).or.(vs.eq.miss).or.(vw.eq.miss).or.(ve.eq.miss))then
        a2gradx(ix,iy) = miss
        a2grady(ix,iy) = miss
      end if
    end if
  end do
  !-------------------
  ! ix = 2,3,4,....nx-1
  !-------------------
  do ix = 2, nx-1
    if (a2in(ix,iy).eq.miss)then
      a2gradx(ix,iy) = miss 
      a2grady(ix,iy) = miss 
    else
      !---
      vn = a2in(ix, iy+1)
      vs = a2in(ix, iy-1)
      vw = a2in(ix-1, iy)
      ve = a2in(ix+1, iy)
      !---
      a2gradx(ix,iy) = (ve - vw) / (2.0*dew)
      a2grady(ix,iy) = (vn - vs) / (dn + ds)
      !--- check miss ------
      if ((vn.eq.miss).or.(vs.eq.miss).or.(vw.eq.miss).or.(ve.eq.miss))then
        a2gradx(ix,iy) = miss
        a2grady(ix,iy) = miss
      end if
    end if
  end do
end do
return
END SUBROUTINE mk_a2grad

!*********************************************************
SUBROUTINE mk_a2grad_abs(a2in, a1lon, a1lat, miss, nx, ny, a2gradabs)
!---------------------------------
! data order should be South -> North, West -> East
!---------------------------------
implicit none
!-----------------
integer                    nx
integer                    ny
!--- in ----------
real,dimension(nx,ny)   :: a2in
!f2py intent(in)           a2in

real,dimension(nx)      :: a1lon
!f2py intent(in)           a1lon

real,dimension(ny)      :: a1lat
!f2py intent(in)           a1lon

real,optional           :: miss
!f2py intent(in)           miss
!--- out ---------
real,dimension(nx,ny)   :: a2gradabs
!f2py intent(out)          a2gradabs
!--- calc --------
real                       lat_first
real                       dn, ds, dew
real                       vn, vs, vw, ve
real                       lat
real                       gradx, grady
integer                    ix,  iy
integer                    ixn, ixs, ixw, ixe
integer                    iyn, iys, iyw, iye
!-----------------
lat_first  = a1lat(1)

!-------------------
! ny = 1 & ny
!-------------------
do iy = 1, ny, ny-1  ! iy = 1 & ny
  call ixy2iixy(nx, ny, ix, iy+1, ixn, iyn)
  call ixy2iixy(nx, ny, ix, iy-1, ixs, iys)
  dn  = hubeny_real(a1lat(iy), 0.0, a1lat(iyn), 0.0)
  ds  = hubeny_real(a1lat(iy), 0.0, a1lat(iys), 0.0)
  dew = hubeny_real(a1lat(iy), a1lon(1), a1lat(iy), a1lon(2))

  do ix = 1, nx
    if (a2in(ix,iy).eq.miss)then
      a2gradabs(ix,iy) = miss
    else
      !---
      call ixy2iixy(nx, ny, ix, iy+1, ixn, iyn)
      call ixy2iixy(nx, ny, ix, iy-1, ixs, iys)
      call ixy2iixy(nx, ny, ix-1, iy, ixw, iyw)
      call ixy2iixy(nx, ny, ix+1, iy, ixe, iye)
      !---
      vn = a2in(ixn, iyn)
      vs = a2in(ixs, iys)
      vw = a2in(ixw, iyw)
      ve = a2in(ixe, iye)
      !---
      gradx = (ve - vw) / (2.0*dew)
      grady = (vn - vs) / (dn + ds)
      a2gradabs(ix, iy) = ( gradx**2.0 + grady**2.0)**0.5 
      !--- check miss ------
      if ((vn.eq.miss).or.(vs.eq.miss).or.(vw.eq.miss).or.(ve.eq.miss))then
        a2gradabs(ix,iy) = miss
      end if
    end if
  end do
end do

!-----------------
! iy = 2,3,..ny-1
!-----------------
do iy = 2, ny-1
  lat = a1lat(iy)
  dn  = hubeny_real(lat, 0.0, a1lat(iy+1), 0.0)
  ds  = hubeny_real(lat, 0.0, a1lat(iy-1), 0.0)
  dew = hubeny_real(lat, a1lon(1), lat, a1lon(2))
  !-------------------
  ! ix = 1 & nx
  !-------------------
  do ix = 1, nx, nx-1  ! ix = 1 & nx
    if (a2in(ix,iy).eq.miss)then
      a2gradabs(ix,iy) = miss
    else
      call ixy2iixy(nx, ny, ix-1, iy, ixw, iyw)
      call ixy2iixy(nx, ny, ix+1, iy, ixe, iye)
      !---
      vn = a2in(ix, iy+1)
      vs = a2in(ix, iy-1)
      vw = a2in(ixw, iy)
      ve = a2in(ixe, iy)
      !---
      gradx = (ve - vw) / (2.0*dew)
      grady = (vn - vs) / (dn + ds)
      a2gradabs(ix, iy) = ( gradx**2.0 + grady**2.0)**0.5 
      !--- check miss ------
      if ((vn.eq.miss).or.(vs.eq.miss).or.(vw.eq.miss).or.(ve.eq.miss))then
        a2gradabs(ix,iy) = miss
      end if
    end if
  end do
  !-------------------
  ! ix = 2,3,4,....nx-1
  !-------------------
  do ix = 2, nx-1
    if (a2in(ix,iy).eq.miss)then
      a2gradabs(ix,iy) = miss
    else
      !---
      vn = a2in(ix, iy+1)
      vs = a2in(ix, iy-1)
      vw = a2in(ix-1, iy)
      ve = a2in(ix+1, iy)
      !---
      gradx = (ve - vw) / (2.0*dew)
      grady = (vn - vs) / (dn + ds)
      a2gradabs(ix, iy) = ( gradx**2.0 + grady**2.0)**0.5 
      !--- check miss ------
      if ((vn.eq.miss).or.(vs.eq.miss).or.(vw.eq.miss).or.(ve.eq.miss))then
        a2gradabs(ix,iy) = miss
      end if
    end if
  end do
end do
return
END SUBROUTINE mk_a2grad_abs

!*********************************************************
SUBROUTINE ixy2iixy(nx, ny, ix, iy, iix, iiy)
!--------------------------
! data array order should be "South->North" & "West->East"
!--------------------------
implicit none
!--- input -----------
integer             nx, ny, ix, iy
!f2py intent(in)    nx, ny, ix, iy
!--- output ----------
integer             iix, iiy
!f2py intent(out)   iix, iiy
!--- calc  -----------
!---------------------
if (iy .le. 0)then
  iiy = 1 - iy
  iix = mod(ix + int(nx*0.5), nx)
else if (iy .gt. ny) then
  iiy = 2*ny - iy +1
  iix = mod(ix + int(nx*0.5), nx)
else
  iiy = iy
  iix = ix
end if
!
if (iix .gt. nx) then
  iix = mod(iix, nx)
else if (iix .le. 0) then
  iix = nx - mod(abs(iix), nx)
end if
!
return
END SUBROUTINE ixy2iixy

!**************************************************************
FUNCTION hubeny_real(lat1, lon1, lat2, lon2)
  implicit none
  !-- for input -----------
  real                                  lat1, lon1, lat2, lon2
!f2py intent(in)                        lat1, lon1, lat2, lon2
  !-- for output-----------
  real                                  hubeny_real
!f2py intent(out)                       hubeny_real
  !-- for calc ------------
  real,parameter                     :: pi = atan(1.0)*4.0
  real,parameter                     :: a  = 6378137
  real,parameter                     :: b  = 6356752.314140
  real,parameter                     :: e2 = 0.00669438002301188
  real,parameter                     :: a_1_e2 = 6335439.32708317
  real                                  M, N, W
  real                                  latrad1, latrad2, lonrad1, lonrad2
  real                                  latave, dlat, dlon
  real                                  dlondeg
  !------------------------
  latrad1   = lat1 * pi / 180.0
  latrad2   = lat2 * pi / 180.0
  lonrad1   = lon1 * pi / 180.0
  lonrad2   = lon2 * pi / 180.0
  !
  latave    = (latrad1 + latrad2)/2.0
  dlat      = latrad2 - latrad1
  dlon      = lonrad2 - lonrad1
  !
  dlondeg   = lon2 - lon1
  if ( abs(dlondeg) .gt. 180.0) then
    dlondeg = 180.0 - mod(abs(dlondeg), 180.0)
    dlon    = dlondeg * pi / 180.0
  end if
  !-------
  W  = sqrt(1.0 - e2 * sin(latave)**2.0 )
  M  =  a_1_e2 / (W**3.0)
  N  =  a / W
  hubeny_real  = sqrt( (dlat * M)**2.0 + (dlon * N * cos(latave))**2.0 )
RETURN
END FUNCTION hubeny_real
!!**************************************************************


!!*********************************************************
!SUBROUTINE count_frontlen_grids(a2in, a2mask, miss, nx, ny, a1grids)
!implicit none
!!--- in ---------
!integer                     nx, ny
!real,dimension(nx,ny)    :: a2in, a2mask  ! a2mask: masked grids have miss
!!f2py intent(in)            a2in, a2mask
!real                        miss
!!f2py intent(in)            miss
!!--- out --------
!real,dimension(nx*ny)    :: a1grids
!!f2py intent(out)           a1grids
!!--- calc -------
!integer                     ix,iy,ik
!integer                     icount_surr
!integer                     iix,iiy
!integer                     id, id_min, ik_tmp
!integer                     stopflag
!integer                     idMAX, idMIN, MAXmin, MINmin
!integer,dimension(8)     :: a1x, a1y
!integer,dimension(nx,ny) :: a2id
!integer,dimension(nx*ny) :: a1connect
!integer,dimension(nx*ny) :: a1flag
!real,dimension(nx*ny)    :: a1num
!real,dimension(nx,ny)    :: a2in_tmp
!!--- parameter --
!integer,parameter        :: miss_int = -9999
!!------
!!-- initialize ---
!a2in_tmp  = a2in
!a1connect = miss_int
!a1num     = 0
!a1grids   = 0.0
!a2id      = miss_int
!id        = 0
!a1flag    = 0
!!----------------
!do iy = 1,ny
!  do ix = 1,nx
!    if (a2in(ix,iy).ne.miss)then
!      !-- 1st search ----
!      call mk_8gridsxy(ix,iy, a1x, a1y)
!      !*****************************
!      ! if (ix,iy) has not been given an ID
!      !------------
!      if (a2id(ix,iy).eq.miss_int)then
!        id = id + 1
!        a2id(ix,iy) = id 
!      end if
!      !***************************
!      !----------------------------
!      icount_surr = 0
!      do ik = 1,8
!        iix = a1x(ik)
!        iiy = a1y(ik)
!        if (a2in_tmp(iix,iiy).ne.miss)then
!          icount_surr   = icount_surr + 1
!          !--- when adjacent id is -9999
!          if (a2id(iix,iiy).eq.miss_int)then
!            a2id(iix,iiy) = a2id(ix,iy)
!          !--- when adjacent id is not miss
!          else if (a2id(iix,iiy).ne.miss_int)then
!            if (a2id(ix,iy).eq.a2id(iix,iiy))then
!              continue
!            else
!              idMIN = min(a2id(ix,iy), a2id(iix,iiy))
!              idMAX = max(a2id(ix,iy), a2id(iix,iiy))
!              MINmin = a1connect(idMIN)
!              MAXmin = a1connect(idMAX)
!              if ((MINmin.eq.miss_int).and.(MAXmin.eq.miss_int))then
!                a1connect(idMAX) = idMIN
!              else if ((MINmin.eq.miss_int).and.(MAXmin.ne.miss_int))then
!                if (MAXmin.lt.idMIN)then
!                  a1connect(idMIN) = MAXmin
!                else if (idMIN .le. MAXmin) then
!                  a1connect(idMAX) = idMIN
!                end if
!              else if ((MINmin.ne.miss_int).and.(MAXmin.eq.miss_int))then
!                a1connect(idMAX)=MINmin
!              else if (MAXmin.eq.MINmin)then
!                continue
!              else if (MAXmin.gt.MINmin)then
!                a1connect(idMAX) = MINmin
!              else if (MINmin.gt.MAXmin)then
!                a1connect(idMIN) = MAXmin
!              end if
!            end if
!          end if
!        end if        
!      end do
!    end if
!  end do
!end do
!!-----
!
!do ik = 1,nx*ny
!  if (a1connect(ik).ne.miss_int)then
!    stopflag   = 0
!    ik_tmp     = ik
!    id_min     = a1connect(ik_tmp)
!    do while (stopflag .ne.1)
!      if (a1connect(ik_tmp).eq.miss_int)then
!        stopflag = 1
!      else
!        if (a1connect(ik_tmp).lt.id_min)then
!          id_min = a1connect(ik_tmp)
!        else
!          continue
!        end if
!        ik_tmp = a1connect(ik_tmp)
!      end if
!    end do 
!    a1connect(ik) = id_min
!  end if
!end do
!
!!-----
!do iy = 1, ny
!  do ix = 1,nx
!    if (a2id(ix,iy).ne.miss_int)then
!      ik = a2id(ix,iy)
!      if (a1connect(ik).ne.miss_int)then
!        a2id(ix,iy) = a1connect(ik)
!      else
!        continue
!      end if
!    end if
!  end do
!end do
!
!do iy=1,ny
!  do ix=1,nx
!    id = a2id(ix,iy)
!    a1grids(id) = a1grids(id) + 1.0
!    call mk_8gridsxy(ix,iy, a1x, a1y)
!    !--- check domain edge ----
!    do ik = 1,8
!      iix = a1x(ik)
!      iiy = a1y(ik)
!      if (a2mask(iix,iiy).eq.miss)then
!        a1flag(id) = 1
!      end if
!    end do
!
!
!    !do dy = -2,2
!    !  do dx = -2,2
!    !    iix =ix + dx 
!    !    iiy =iy + dy
!    !    if (a2mask(iix,iiy).eq.miss)then
!    !      a1flag(id) = 1
!    !    end if
!    !  end do
!    !end do
!    !--------------------------
!  end do 
!end do
!
!!-- remove IDs with flag=1 ---
!do id=1,nx*ny
!  if (a1flag(id).eq.1)then
!    a1grids(id) = 0
!  end if
!end do
!!-----------------------------
!return
!END SUBROUTINE count_frontlen_grids
!!
!!*********************************************************
!
!
!!*********************************************************
!SUBROUTINE count_fronts(a2in, miss, nx, ny, nfronts)
!implicit none
!!--- in ---------
!integer                     nx, ny
!real,dimension(nx,ny)    :: a2in
!!f2py intent(in)            a2in
!real                        miss
!!f2py intent(in)            miss
!!--- out --------
!integer                     nfronts
!!f2py intent(out)           nfronts
!!--- calc -------
!integer                     ix,iy,ik
!integer                     icount_surr
!integer                     iix,iiy
!integer                     id, id_min, ik_tmp
!integer                     stopflag
!integer                     idMAX, idMIN, MAXmin, MINmin
!integer,dimension(8)     :: a1x, a1y
!integer,dimension(nx,ny) :: a2id
!integer,dimension(nx*ny) :: a1connect
!integer,dimension(nx*ny) :: a1idflag
!real,dimension(nx*ny)    :: a1num
!real,dimension(nx,ny)    :: a2in_tmp
!!--- parameter --
!integer,parameter        :: miss_int = -9999
!!------
!!-- initialize ---
!a2in_tmp  = a2in
!a1connect = miss_int
!a1num     = 0
!a2id      = miss_int
!id        = 0
!a1idflag  = 0
!!----------------
!do iy = 1,ny
!  do ix = 1,nx
!    if (a2in(ix,iy).ne.miss)then
!      !-- 1st search ----
!      call mk_8gridsxy(ix,iy, a1x, a1y)
!      !*****************************
!      ! if (ix,iy) has not been given an ID
!      !------------
!      if (a2id(ix,iy).eq.miss_int)then
!        id = id + 1
!        a2id(ix,iy) = id 
!      end if
!      !***************************
!      !----------------------------
!      icount_surr = 0
!      do ik = 1,8
!        iix = a1x(ik)
!        iiy = a1y(ik)
!        if (a2in_tmp(iix,iiy).ne.miss)then
!          icount_surr   = icount_surr + 1
!          !--- when adjacent id is -9999
!          if (a2id(iix,iiy).eq.miss_int)then
!            a2id(iix,iiy) = a2id(ix,iy)
!          !--- when adjacent id is not miss
!          else if (a2id(iix,iiy).ne.miss_int)then
!            if (a2id(ix,iy).eq.a2id(iix,iiy))then
!              continue
!            else
!              idMIN = min(a2id(ix,iy), a2id(iix,iiy))
!              idMAX = max(a2id(ix,iy), a2id(iix,iiy))
!              MINmin = a1connect(idMIN)
!              MAXmin = a1connect(idMAX)
!              if ((MINmin.eq.miss_int).and.(MAXmin.eq.miss_int))then
!                a1connect(idMAX) = idMIN
!              else if ((MINmin.eq.miss_int).and.(MAXmin.ne.miss_int))then
!                if (MAXmin.lt.idMIN)then
!                  a1connect(idMIN) = MAXmin
!                else if (idMIN .le. MAXmin) then
!                  a1connect(idMAX) = idMIN
!                end if
!              else if ((MINmin.ne.miss_int).and.(MAXmin.eq.miss_int))then
!                a1connect(idMAX)=MINmin
!              else if (MAXmin.eq.MINmin)then
!                continue
!              else if (MAXmin.gt.MINmin)then
!                a1connect(idMAX) = MINmin
!              else if (MINmin.gt.MAXmin)then
!                a1connect(idMIN) = MAXmin
!              end if
!            end if
!          end if
!        end if        
!      end do
!      !if (icount_surr.eq. 0)then
!      !  a2in_tmp(ix,iy) = miss
!      !  a2id(ix,iy)     = miss_int
!      !  id = id -1
!      !end if
!    end if
!  end do
!end do
!!-----
!
!do ik = 1,nx*ny
!  if (a1connect(ik).ne.miss_int)then
!    stopflag   = 0
!    ik_tmp     = ik
!    id_min     = a1connect(ik_tmp)
!    do while (stopflag .ne.1)
!      if (a1connect(ik_tmp).eq.miss_int)then
!        stopflag = 1
!      else
!        if (a1connect(ik_tmp).lt.id_min)then
!          id_min = a1connect(ik_tmp)
!        else
!          continue
!        end if
!        ik_tmp = a1connect(ik_tmp)
!      end if
!    end do 
!    a1connect(ik) = id_min
!  end if
!end do
!
!!-----
!do iy = 1, ny
!  do ix = 1,nx
!    if (a2id(ix,iy).ne.miss_int)then
!      ik = a2id(ix,iy)
!      if (a1connect(ik).ne.miss_int)then
!        a2id(ix,iy) = a1connect(ik)
!      else
!        continue
!      end if
!    end if
!  end do
!end do
!
!do iy=1,ny
!  do ix=1,nx
!    id = a2id(ix,iy)
!    a1idflag(id) = 1
!  end do 
!end do
!
!nfronts = sum(a1idflag)
!return
!END SUBROUTINE count_fronts 
!!
!!*********************************************************
!!*********************************************************
!SUBROUTINE fill_front_gap(a2in, miss, nx, ny, a2out)
!implicit none
!integer                  nx, ny
!real,dimension(nx,ny) :: a2in
!!f2py intent(in)         a2in
!real                     miss
!!f2py intent(in)         miss
!!---- out -------
!real,dimension(nx,ny) :: a2out
!!f2py intent(out)        a2out
!!---- calc ------
!integer                  ix,iy
!integer                  iix,iiy, imx,imy
!integer                  dx,dy
!!---- init ------
!a2out     = a2in
!!----------------
!do iy = 1,ny
!  do ix = 1,nx
!    !print '(E10.2)',a2in(ix,iy)
!    if (a2in(ix,iy).ne.miss)then
!      do dy = -2,2,4
!        call ixy2iixy_saone(ix,iy+dy, iix, iiy)
!        if (a2in(iix,iiy).ne.miss)then
!          imx = int((ix + iix)*0.5)
!          imy = int((iy + iiy)*0.5)
!          if (a2in(imx,imy).eq.miss)then
!            a2out(imx,imy) = (a2in(ix,iy) + a2in(iix,iiy))*0.5
!          end if
!        end if
!      end do
!      do dx = -2,2,4
!        call ixy2iixy_saone(ix+dx,iy, iix, iiy)
!        if (a2in(iix,iiy).ne.miss)then
!          !print *,iix,iiy,a2in(iix,iiy)
!          imx = int((ix + iix)*0.5)
!          imy = int((iy + iiy)*0.5)
!          if (a2in(imx,imy).eq.miss)then
!            a2out(imx,imy) = (a2in(ix,iy) + a2in(iix,iiy))*0.5
!            !print *, imx,imy,a2in(ix,iy), a2out(imx,imy), a2in(iix,iiy)
!          end if
!        end if
!      end do
!    end if
!  end do
!end do
!!------------
!return
!END SUBROUTINE fill_front_gap
!!*********************************************************
!SUBROUTINE check_exist_3deg_saone(a2in, miss, nx, ny, a2out)
!implicit none
!!------------------------------------------
!! input: 1.0 deg array
!! output: 1.0 deg array, but informations are reduced on 3.0 deg resol
!! check whether any non-miss 1.0deg grid exists in 3.0 deg grid box.
!! if any non-miss value exists in 3.0 degree grid box,
!! fill the 3.0 grid box with 1.0
!!-------------------------------------------
!integer                  nx, ny
!real,dimension(nx,ny) :: a2in
!!f2py intent(in)         a2in
!real                     miss
!!f2py intent(in)         miss
!!--- out --------
!real,dimension(nx,ny) :: a2out
!!f2py intent(out)        a2out
!!--- calc -------
!integer                  ix,iy,ik
!integer                  iix,iiy
!integer                  ixn,ixs,ixe,ixw
!integer                  iyn,iys,iye,iyw
!integer,dimension(8)  :: a1x, a1y
!integer                  icount
!real                     fillv
!!----------------
!a2out = miss
!do iy = 2, ny-1, 3
!  do ix = 2,nx-1, 3
!    a1x    = -9999
!    a1y    = -9999
!    !------
!    call ixy2iixy_saone(ix, iy+1, ixn, iyn)
!    call ixy2iixy_saone(ix, iy-1, ixs, iys)
!    call ixy2iixy_saone(ix-1, iy, ixw, iyw)
!    call ixy2iixy_saone(ix+1, iy, ixe, iye)
!    !------
!    a1x(1) = ixn
!    a1x(2) = ixe
!    a1x(3) = ixe
!    a1x(4) = ixe
!    a1x(5) = ixw
!    a1x(6) = ixw
!    a1x(7) = ixw
!    a1x(8) = ixs     
!    !------
!    a1y(1) = iyn
!    a1y(2) = iyn
!    a1y(3) = iy
!    a1y(4) = iys
!    a1y(5) = iyn
!    a1y(6) = iy
!    a1y(7) = iys
!    a1y(8) = iys     
!    !------
!    icount = 0
!    if (a2in(ix,iy).ne.miss)then
!      icount = 1
!    else
!      do ik = 1,8
!        iix = a1x(ik)
!        iiy = a1y(ik)
!        if (a2in(iix,iiy).ne.miss)then
!          icount = icount + 1
!        end if        
!      end do
!    end if
!    !----
!    if (icount .ge. 1)then
!      fillv = 1.0
!    else
!      fillv = 0.0
!    end if
!    !----
!    a2out(ix,iy) = fillv
!    do ik = 1,8
!      iix = a1x(ik)
!      iiy = a1y(ik)
!      a2out(iix,iiy) = fillv
!    end do
!  end do
!end do
!return
!END SUBROUTINE check_exist_3deg_saone
!!*********************************************************
!SUBROUTINE sum_9grids_saone(a2in, miss, nx, ny, a2out)
!implicit none
!!--- in --------
!integer                  nx, ny
!real,dimension(nx,ny) :: a2in
!!f2py intent(in)         a2in
!real                     miss
!!f2py intent(in)         miss
!!--- out --------
!real,dimension(nx,ny) :: a2out
!!f2py intent(out)        a2out
!!--- calc -------
!integer                  ix,iy,ik
!integer                  iix,iiy
!integer                  ixn,ixs,ixe,ixw
!integer                  iyn,iys,iye,iyw
!integer,dimension(8)  :: a1x, a1y
!integer                  icount
!!----------------
!a2out = miss
!do iy = 1,ny
!  do ix = 1,nx
!    if (a2in(ix,iy).ne.miss)then
!      a2out(ix,iy) = a2in(ix,iy)
!      a1x    = -9999
!      a1y    = -9999
!      !------
!      call ixy2iixy_saone(ix, iy+1, ixn, iyn)
!      call ixy2iixy_saone(ix, iy-1, ixs, iys)
!      call ixy2iixy_saone(ix-1, iy, ixw, iyw)
!      call ixy2iixy_saone(ix+1, iy, ixe, iye)
!      !------
!      a1x(1) = ixn
!      a1x(2) = ixe
!      a1x(3) = ixe
!      a1x(4) = ixe
!      a1x(5) = ixw
!      a1x(6) = ixw
!      a1x(7) = ixw
!      a1x(8) = ixs     
!      !------
!      a1y(1) = iyn
!      a1y(2) = iyn
!      a1y(3) = iy
!      a1y(4) = iys
!      a1y(5) = iyn
!      a1y(6) = iy
!      a1y(7) = iys
!      a1y(8) = iys     
!      !------
!      icount = 1
!      do ik = 1,8
!        iix = a1x(ik)
!        iiy = a1y(ik)
!        if (a2in(iix,iiy).ne.miss)then
!          icount = icount + 1
!          a2out(ix,iy) = a2out(ix,iy) + a2in(iix,iiy)
!        end if        
!      end do
!      a2out(ix,iy) = a2out(ix,iy)
!    end if
!  end do
!end do
!END SUBROUTINE sum_9grids_saone
!!*********************************************************
!SUBROUTINE mean_9grids_saone(a2in, miss, nx, ny, a2out)
!implicit none
!!--- in --------
!integer                  nx, ny
!real,dimension(nx,ny) :: a2in
!!f2py intent(in)         a2in
!real                     miss
!!f2py intent(in)         miss
!!--- out --------
!real,dimension(nx,ny) :: a2out
!!f2py intent(out)        a2out
!!--- calc -------
!integer                  ix,iy,ik
!integer                  iix,iiy
!integer                  ixn,ixs,ixe,ixw
!integer                  iyn,iys,iye,iyw
!integer,dimension(8)  :: a1x, a1y
!integer                  icount
!!----------------
!a2out = miss
!do iy = 1,ny
!  do ix = 1,nx
!    if (a2in(ix,iy).ne.miss)then
!      a2out(ix,iy) = a2in(ix,iy)
!      a1x    = -9999
!      a1y    = -9999
!      !------
!      call ixy2iixy_saone(ix, iy+1, ixn, iyn)
!      call ixy2iixy_saone(ix, iy-1, ixs, iys)
!      call ixy2iixy_saone(ix-1, iy, ixw, iyw)
!      call ixy2iixy_saone(ix+1, iy, ixe, iye)
!      !------
!      a1x(1) = ixn
!      a1x(2) = ixe
!      a1x(3) = ixe
!      a1x(4) = ixe
!      a1x(5) = ixw
!      a1x(6) = ixw
!      a1x(7) = ixw
!      a1x(8) = ixs     
!      !------
!      a1y(1) = iyn
!      a1y(2) = iyn
!      a1y(3) = iy
!      a1y(4) = iys
!      a1y(5) = iyn
!      a1y(6) = iy
!      a1y(7) = iys
!      a1y(8) = iys     
!      !------
!      icount = 1
!      do ik = 1,8
!        iix = a1x(ik)
!        iiy = a1y(ik)
!        if (a2in(iix,iiy).ne.miss)then
!          icount = icount + 1
!          a2out(ix,iy) = a2out(ix,iy) + a2in(iix,iiy)
!        end if        
!      end do
!      a2out(ix,iy) = a2out(ix,iy) / icount
!    end if
!  end do
!end do
!
!
!END SUBROUTINE mean_9grids_saone

!!*********************************************************
!SUBROUTINE del_front_lesseq_ngrids(a2in, wgtflag, miss, thnum, nx, ny, a2out)
!implicit none
!!--- in ---------
!integer                     nx, ny
!integer                     wgtflag
!!f2py intent(in)            wgtflag
!integer                     thnum
!!f2py intent(in)            thnum
!real,dimension(nx,ny)    :: a2in
!!f2py intent(in)            a2in
!real                        miss
!!f2py intent(in)            miss
!!--- out --------
!real,dimension(nx,ny)    :: a2out
!!f2py intent(out)           a2out
!!--- calc -------
!integer                     ix,iy,ik
!integer                     icount_surr
!integer                     iix,iiy
!integer                     id, id_min, ik_tmp
!integer                     stopflag
!integer                     idMAX, idMIN, MAXmin, MINmin
!integer,dimension(8)     :: a1x, a1y
!integer,dimension(nx,ny) :: a2id
!integer,dimension(nx*ny) :: a1connect
!real,dimension(nx*ny)    :: a1num
!real,dimension(nx,ny)    :: a2in_tmp, a2arearat
!real                        arearat
!!--- parameter --
!integer,parameter        :: miss_int = -9999
!!------
!!-- initialize ---
!a2in_tmp  = a2in
!a1connect = miss_int
!a1num     = 0
!a2id      = miss_int
!id        = 0
!!----------------
!do iy = 1,ny
!  do ix = 1,nx
!    if (a2in(ix,iy).ne.miss)then
!      !-- 1st search ----
!      call mk_8gridsxy(ix,iy, a1x, a1y)
!      !*****************************
!      ! if (ix,iy) has not been given an ID
!      !------------
!      if (a2id(ix,iy).eq.miss_int)then
!        id = id + 1
!        a2id(ix,iy) = id 
!      end if
!      !***************************
!      !----------------------------
!      icount_surr = 0
!      do ik = 1,8
!        iix = a1x(ik)
!        iiy = a1y(ik)
!        if (a2in_tmp(iix,iiy).ne.miss)then
!          icount_surr   = icount_surr + 1
!          !--- when adjacent id is -9999
!          if (a2id(iix,iiy).eq.miss_int)then
!            a2id(iix,iiy) = a2id(ix,iy)
!          !--- when adjacent id is not miss
!          else if (a2id(iix,iiy).ne.miss_int)then
!            if (a2id(ix,iy).eq.a2id(iix,iiy))then
!              continue
!            else
!              idMIN = min(a2id(ix,iy), a2id(iix,iiy))
!              idMAX = max(a2id(ix,iy), a2id(iix,iiy))
!              MINmin = a1connect(idMIN)
!              MAXmin = a1connect(idMAX)
!              if ((MINmin.eq.miss_int).and.(MAXmin.eq.miss_int))then
!                a1connect(idMAX) = idMIN
!              else if ((MINmin.eq.miss_int).and.(MAXmin.ne.miss_int))then
!                if (MAXmin.lt.idMIN)then
!                  a1connect(idMIN) = MAXmin
!                else if (idMIN .le. MAXmin) then
!                  a1connect(idMAX) = idMIN
!                end if
!              else if ((MINmin.ne.miss_int).and.(MAXmin.eq.miss_int))then
!                a1connect(idMAX)=MINmin
!              else if (MAXmin.eq.MINmin)then
!                continue
!              else if (MAXmin.gt.MINmin)then
!                a1connect(idMAX) = MINmin
!              else if (MINmin.gt.MAXmin)then
!                a1connect(idMIN) = MAXmin
!              end if
!            end if
!          end if
!        end if        
!      end do
!      if (icount_surr.eq. 0)then
!        a2in_tmp(ix,iy) = miss
!        a2id(ix,iy)     = miss_int
!        id = id -1
!      end if
!    end if
!  end do
!end do
!!-----
!do ik = 1,nx*ny
!  if (a1connect(ik).ne.miss_int)then
!    stopflag   = 0
!    ik_tmp     = ik
!    id_min     = a1connect(ik_tmp)
!    do while (stopflag .ne.1)
!      if (a1connect(ik_tmp).eq.miss_int)then
!        stopflag = 1
!      else
!        if (a1connect(ik_tmp).lt.id_min)then
!          id_min = a1connect(ik_tmp)
!        else
!          continue
!        end if
!        ik_tmp = a1connect(ik_tmp)
!      end if
!    end do 
!    a1connect(ik) = id_min
!  end if
!end do
!!-----
!do iy = 1, ny
!  do ix = 1,nx
!    if (a2id(ix,iy).ne.miss_int)then
!      ik = a2id(ix,iy)
!      if (a1connect(ik).ne.miss_int)then
!        a2id(ix,iy) = a1connect(ik)
!      else
!        continue
!      end if
!    end if
!  end do
!end do
!!!********************************
!!! count id
!!!-----------------------
!!do iy = 1,ny
!!  do ix = 1,nx
!!    if (a2id(ix,iy).ne.miss)then
!!      id = a2id(ix,iy)
!!      a1num(id) = a1num(id) + 1
!!    end if
!!  end do
!!end do
!!
!
!!********************************
!! count id
!!-----------------------
!if (wgtflag .eq. 0)then
!  do iy = 1,ny
!    do ix = 1,nx
!      if (a2id(ix,iy).ne.miss_int)then
!        id = a2id(ix,iy)
!        a1num(id) = a1num(id) + 1
!      end if
!    end do
!  end do
!else if (wgtflag .eq. 1)then
!  do iy = 1,ny
!    CALL mk_a2arearat_sphere_saone(a2arearat)
!    arearat   = a2arearat(1,iy)
!    do ix = 1,nx
!      if (a2id(ix,iy).ne.miss_int)then
!        id = a2id(ix,iy)
!        a1num(id) = a1num(id) + 1*arearat
!      end if
!    end do
!  end do
!end if
!!********************************
!! remove small segments
!!-----------------------
!a2out = a2in_tmp
!do iy = 1,ny
!  do ix = 1,nx
!    id = a2id(ix,iy)
!    if (id.eq.miss_int)then
!      cycle
!    end if
!    if (a1num(id).le.real(thnum))then
!      a2out(ix,iy) = miss
!    end if
!  end do
!end do
!!--------------------------------
!!a2out = real(a2id)
!return
!END SUBROUTINE del_front_lesseq_ngrids
!!
!!*********************************************************
!
!
!!*********************************************************
!SUBROUTINE mk_a2arearat_sphere_saone(a2arearat)
!implicit none
!!--------------------------
!! make ratio of area relative to 1.0x1.0 grid located on the equator
!!-- out -----------------
!real,dimension(360,180)              :: a2arearat
!!f2py intent(out)                       a2arearat
!!-- calc ----------------
!integer                                 iy
!real                                    lats, latn, lonw, lone
!real                                    area_eq
!!------------------------
!area_eq  = cal_area_sphere(-0.5, 0.5, 0.0, 1.0)
!do iy =1,180
!  latn = -90.0 + real(iy)
!  lats = latn - 1.0
!  lonw = 0.0
!  lone = 1.0
!  a2arearat(:,iy) = cal_area_sphere(lats, latn, lonw, lone) / area_eq
!end do
!!------------------------
!return
!END SUBROUTINE mk_a2arearat_sphere_saone
!!*********************************************************
!SUBROUTINE mk_a2area_sphere_saone(a2area)
!implicit none
!!-- out -----------------
!real,dimension(360,180)              :: a2area 
!!f2py intent(out)                       a2area
!!-- calc ----------------
!integer                                 iy
!real                                    lats, latn, lonw, lone
!!------------------------
!do iy =1,180
!  latn = -90.0 + real(iy)
!  lats = latn - 1.0
!  lonw = 0.0
!  lone = 1.0
!  a2area(:,iy) = cal_area_sphere(lats, latn, lonw, lone)
!end do
!!------------------------
!return
!END SUBROUTINE mk_a2area_sphere_saone
!!*********************************************************
!FUNCTION cal_area_sphere(lats, latn, lonw, lone)
!!----------------------
!! estimate area (km2) assuming that the eath is a sphere
!! S = r2 * dlon * pi / 180 * (sin(lat2) - sin(lat1))
!! lat, lon are in degree
!!----------------------
!implicit none
!!--- in -----------
!real                                    lats, latn, lonw, lone  ! (deg.)
!!f2py intent(in)                        lats, latn, lonw, lone
!
!!--- out ----------
!real                                    cal_area_sphere
!!f2py intent(out)                       cal_area_sphere
!!--- para ---------
!real,parameter                       :: r = 6371.012  ! (km)
!real,parameter                       :: pi = 3.1416
!!--- calc ---------
!real                                    dlon
!!----------------------
!dlon = min( abs(lone-lonw), abs(360.0-(lone-lonw)))
!
!cal_area_sphere = r**2.0 *dlon*pi/180.0 *(sin(latn/180.0*pi)-sin(lats/180.0*pi))
!!----------------------
!return
!END FUNCTION cal_area_sphere
!
!
!!*********************************************************
!SUBROUTINE del_front_3grids(a2in, miss, nx, ny, a2out)
!implicit none
!!--- in ---------
!integer                  nx, ny
!real,dimension(nx,ny) :: a2in
!!f2py intent(in)         a2in
!real                     miss
!!f2py intent(in)         miss
!!--- out --------
!real,dimension(nx,ny) :: a2out
!!f2py intent(out)        a2out
!!--- calc -------
!integer                  ix,iy,ik
!integer                  ixc, iyc
!integer                  iix,iiy,iik
!integer,dimension(8)  :: a1x, a1y
!integer,dimension(9)  :: a1xx1, a1yy1
!integer,dimension(9)  :: a1xx2, a1yy2
!
!integer                  fcount1, fcount2, fcount3
!integer                  flag_type2
!!------
!!5!1!2!
!!6! !3!
!!7!8!4l
!!------
!!-- initialize ---
!a2out = a2in
!!----------------
!do iy = 1,ny
!  do ix = 1,nx
!    if (a2in(ix,iy).ne.miss)then
!      !-- init ----------
!      a1xx1 = -9999
!      a1yy1 = -9999
!      a1xx2 = -9999
!      a1yy2 = -9999
!      !-- 1st search ----
!      ixc   = ix
!      iyc   = iy
!      a1xx1(9) = ixc
!      a1yy1(9) = iyc
!      call mk_8gridsxy(ixc,iyc, a1x, a1y)
!      !---
!      fcount1 = 0
!      do ik = 1,8
!        iix = a1x(ik)
!        iiy = a1y(ik)
!        if (a2in(iix,iiy).ne.miss)then
!          fcount1 = fcount1 + 1
!          a1xx1(fcount1) = iix
!          a1yy1(fcount1) = iiy
!        end if        
!      end do
!      if (fcount1.eq.0)then
!        a2out(ix,iy) = miss
!        a2in(ix,iy)  = miss
!      else if (fcount1 .ge. 3)then
!        cycle
!      else if (fcount1.eq.1)then
!        !******************************
!        ! type 1  
!        !        "*"**
!        !
!        !******************************
!        !----- 2nd search -------
!        ixc          = a1xx1(1)
!        iyc          = a1yy1(1)
!        a1xx2(9)     = ixc
!        a1yy2(9)     = iyc
!
!        call mk_8gridsxy(ixc,iyc, a1x, a1y)
!        fcount2 = 0
!        do ik = 1,8
!          iix = a1x(ik)
!          iiy = a1y(ik)
!          if (a2in(iix,iiy).ne.miss)then
!            fcount2 = fcount2 + 1
!            a1xx2(fcount2) = iix
!            a1yy2(fcount2) = iiy
!          end if
!        end do
!        if (fcount2 .eq. 1)then
!          a2out(ix,iy) = miss
!          a2in(ix,iy)  = miss
!          a2out(ixc,iyc) = miss
!          a2in(ixc,iyc)  = miss
!        else if (fcount2 .ge. 3)then
!          cycle
!        else if (fcount2 .eq. 2)then
!          do iik = 1,2
!            if ((a1xx2(iik).ne.a1xx1(9)).or.(a1yy2(iik).ne.a1yy1(9)))then
!              !--- 3rd search for type1 ---
!              ixc          = a1xx2(iik)
!              iyc          = a1yy2(iik)
!              call mk_8gridsxy(ixc,iyc, a1x, a1y)
!              fcount3 = 0
!              do ik = 1,8
!                iix = a1x(ik)
!                iiy = a1y(ik)
!                if (a2in(iix,iiy).ne.miss)then
!                  fcount3 = fcount3 + 1
!                end if
!              end do
!              if (fcount3 .le.1)then
!                a2out(ix,iy)   = miss
!                a2in(ix,iy)    = miss
!                a2out(ixc,iyc) = miss
!                a2in(ixc,iyc)  = miss
!              end if
!            end if
!          end do
!        end if
!      else if (fcount1 .eq. 2)then
!        !******************************
!        ! type 2  
!        !        *"*"*   "*"*     *"*"*
!        !                 *
!        !
!        !******************************
!        flag_type2 = 0
!        do iik = 1,2
!          !-- 2nd and 3rd search ----
!          ixc = a1xx1(iik)
!          iyc = a1yy1(iik)
!          call mk_8gridsxy(ixc,iyc, a1x, a1y)
!          fcount2 = 0
!          do ik = 1,8
!            iix = a1x(ik)
!            iiy = a1y(ik)
!            if (a2in(iix,iiy).ne.miss)then
!              fcount2 = fcount2 + 1
!              a1xx2(fcount2)  = iix
!              a1yy2(fcount2)  = iiy
!            end if
!          end do
!          if (fcount2.ge.3)then
!            flag_type2 = flag_type2+1
!          else if (fcount2.eq.2)then
!            if ( ((a1xx2(1).eq.a1xx2(2)).and.(a1xx2(1).eq.ixc)) &
!               .or. (a1yy2(1).eq.a1yy2(2).and.(a1yy2(1).eq.iyc)) )then
!              flag_type2 = flag_type2 +1
!            end if
!          end if
!          !--------------------------
!        end do
!        if (flag_type2.eq.0)then
!          a2out(ix,iy)   = miss
!          a2in(ix,iy)    = miss
!          do iik =1,2
!            ixc = a1xx1(iik)
!            iyc = a1yy1(iik)
!            a2out(ixc,iyc) = miss
!            a2in(ixc,iyc)  = miss
!          end do
!        end if
!      end if
!    end if
!  end do
!end do
!
!END SUBROUTINE del_front_3grids
!!
!!*********************************************************
!SUBROUTINE del_front_2grids(a2in, miss, nx, ny, a2out)
!implicit none
!!--- in ---------
!integer                  nx, ny
!real,dimension(nx,ny) :: a2in
!!f2py intent(in)         a2in
!real                     miss
!!f2py intent(in)         miss
!!--- out --------
!real,dimension(nx,ny) :: a2out
!!f2py intent(out)        a2out
!!--- calc -------
!integer                  ix,iy,ik
!integer                  iix,iiy,iik
!integer                  iiix,iiiy
!integer                  ixn,ixs,ixe,ixw
!integer                  iyn,iys,iye,iyw
!integer                  iixn,iixs,iixe,iixw
!integer                  iiyn,iiys,iiye,iiyw
!integer,dimension(8)  :: a1x, a1y
!integer,dimension(8)  :: a1xx, a1yy
!integer                  fcount
!!------
!!5!1!2!
!!6! !3!
!!7!8!4l
!!------
!!-- initialize ---
!a2out = a2in
!!----------------
!do iy = 1,ny
!  do ix = 1,nx
!    if (a2in(ix,iy).ne.miss)then
!      a1x    = -9999
!      a1y    = -9999
!      !------
!      call ixy2iixy_saone(ix, iy+1, ixn, iyn)
!      call ixy2iixy_saone(ix, iy-1, ixs, iys)
!      call ixy2iixy_saone(ix-1, iy, ixw, iyw)
!      call ixy2iixy_saone(ix+1, iy, ixe, iye)
!      !------
!      a1x(1) = ixn
!      a1x(2) = ixe
!      a1x(3) = ixe
!      a1x(4) = ixe
!      a1x(5) = ixw
!      a1x(6) = ixw
!      a1x(7) = ixw
!      a1x(8) = ixs     
!      !------
!      a1y(1) = iyn
!      a1y(2) = iyn
!      a1y(3) = iy
!      a1y(4) = iys
!      a1y(5) = iyn
!      a1y(6) = iy
!      a1y(7) = iys
!      a1y(8) = iys     
!      !------
!      fcount = 1
!
!      do ik = 1,8
!        iix = a1x(ik)
!        iiy = a1y(ik)
!        if (a2in(iix,iiy).ne.miss)then
!          fcount = fcount + 1
!        end if        
!      end do
!      if (fcount.eq.1)then
!        a2out(ix,iy) = miss
!      else if (fcount.ge.3)then
!        continue
!      else
!        !qqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqq
!        do ik = 1,8
!          iix = a1x(ik)
!          iiy = a1y(ik)
!          if (a2in(iix,iiy).ne.miss)then
!            call ixy2iixy_saone(iix, iiy+1, iixn, iiyn)
!            call ixy2iixy_saone(iix, iiy-1, iixs, iiys)
!            call ixy2iixy_saone(iix-1, iiy, iixw, iiyw)
!            call ixy2iixy_saone(iix+1, iiy, iixe, iiye)
!            !------
!            a1xx(1) = iixn
!            a1xx(2) = iixe
!            a1xx(3) = iixe
!            a1xx(4) = iixe
!            a1xx(5) = iixw
!            a1xx(6) = iixw
!            a1xx(7) = iixw
!            a1xx(8) = iixs     
!            !------
!            a1yy(1) = iiyn
!            a1yy(2) = iiyn
!            a1yy(3) = iiy
!            a1yy(4) = iiys
!            a1yy(5) = iiyn
!            a1yy(6) = iiy
!            a1yy(7) = iiys
!            a1yy(8) = iiys     
!            !------
!            do iik = 1,8
!              iiix = a1xx(iik)
!              iiiy = a1yy(iik)
!              if (a2in(iiix,iiiy).ne.miss)then
!                if ((iiix.eq.ix).and.(iiiy.eq.iy))then
!                  continue
!                else
!                  fcount = fcount + 1
!                end if
!              end if
!            end do
!            if (fcount.le.2)then
!              a2out(ix,iy)=miss
!            endif
!          end if
!        end do 
!        !qqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqq
!      end if
!    end if
!  end do
!end do
!
!END SUBROUTINE del_front_2grids
!!*********************************************************
!SUBROUTINE mk_a2rh(a2t, a2q, plev, nx, ny, a2rh)
!implicit none
!!--- in ------
!integer                  nx, ny
!real,dimension(nx,ny) :: a2t, a2q
!!f2py intent(in)         a2t, a2q
!real                     plev
!!f2py intent(in)         plev
!!--- out -----
!real,dimension(nx,ny) :: a2rh
!!f2py intent(out)        a2rh
!!--- calc ----
!integer                    ix,  iy
!!-----------------
!do iy = 1,ny
!  do ix = 1, nx
!    a2rh(ix,iy)  = cal_RH(plev, a2t(ix,iy), a2q(ix,iy))
!  end do
!end do
!
!return
!END SUBROUTINE mk_a2rh
!!*********************************************************
!FUNCTION cal_es(rT)
!implicit none
!!-------------------------
!real  rT
!!f2py intent(in)  rT
!
!real  cal_es
!!f2py intent(out) cal_es
!real,parameter            ::  rT0 = 273.16
!real,parameter            ::  res0= 611.73 ![Pa]
!real,parameter            ::  Lv  = 2.5e6  ![J kg-1]
!real,parameter            ::  Rv  = 461.7 ![J K-1 kg -1]
!!
!cal_es = res0 * exp( Lv/Rv *(1.0/rT0 - 1.0/rT))
!RETURN
!END FUNCTION cal_es
!!*********************************************************
!FUNCTION cal_qs(rT, rP)
!implicit none
!real                 rT, rP
!!f2py intent(in)     rT, rP
!real                 res
!real                 cal_qs
!real,parameter    :: repsi = 0.62185
!!
!res = cal_es(rT)
!cal_qs = repsi * res / (rP - res)
!RETURN
!END FUNCTION cal_qs
!!*********************************************************
!FUNCTION cal_RH(rP, rT, rq)
!implicit none
!real                 rT, rP, rq
!!f2py intent(in)     rT, rP, rq
!!
!real                 res, re
!real,parameter    :: repsi = 0.62185
!real                 cal_RH
!!
!res = cal_es(rT)
!re  = rP *  rq /(rq + repsi)
!cal_RH = re / res *100.0
!RETURN
!END FUNCTION cal_RH
!!*********************************************************
!SUBROUTINE mk_a2thermoadv(a2thermo, a2uwind, a2vwind, nx, ny, a2thermoadv)
!implicit none
!!--- in ------
!integer                  nx, ny
!real,dimension(nx,ny) :: a2thermo, a2uwind, a2vwind
!!f2py intent(in)      :: a2thermo, a2uwind, a2vwind
!!--- out -----
!real,dimension(nx,ny) :: a2thermoadv
!!f2py intent(out)     :: a2thermoadv
!!--- cal -----
!real,dimension(nx,ny) :: a2gradx, a2grady
!!-------------
!CALL mk_a2grad_saone(a2in=a2thermo, a2gradx=a2gradx, a2grady=a2grady)
!a2thermoadv  = a2gradx* a2uwind + a2grady*a2vwind
!!-------------
!return
!END SUBROUTINE mk_a2thermoadv
!!*********************************************************
!SUBROUTINE mk_a2frontspeed(a2thermo, a2frontloc, a2uwind, a2vwind, miss, nx, ny, a2frontspeed)
!implicit none
!!--- in ------
!integer                  nx, ny
!real,dimension(nx,ny) :: a2thermo, a2frontloc, a2uwind, a2vwind
!!f2py intent(in)         a2thermo, a2frontloc, a2uwind, a2vwind
!real                     miss
!!f2py intent(in)         miss
!!--- out -----
!real,dimension(nx,ny) :: a2frontspeed
!!f2py intent(out)        a2frontspeed
!!--- calc ----
!real,dimension(nx,ny) :: a2gradabs, a2gradx, a2grady
!real,dimension(nx,ny) :: a2grad2x, a2grad2y, a2grad2abs
!real,dimension(nx,ny) :: a2gradabs_n, a2gradabs_s, a2gradabs_e, a2gradabs_w
!real                       dn, ds, dew
!real                       lat
!real                       speedsign
!integer                    ix,  iy
!integer                    ixn, ixs, ixw, ixe
!integer                    iyn, iys, iyw, iye
!!--- para --------
!real                    :: lat_first = -89.5
!!-----------------
!CALL mk_a2grad_abs_saone(a2in=a2thermo, a2gradabs=a2gradabs)
!CALL mk_a2grad_saone(a2in=a2thermo, a2gradx=a2gradx, a2grady=a2grady)
!do iy = 1,ny
!  lat = lat_first + (iy -1)*1.0
!  dn  = hubeny_real(lat, 0.0, lat+1.0, 0.0)
!  ds  = hubeny_real(lat, 0.0, lat-1.0, 0.0)
!  dew = hubeny_real(lat, 0.0, lat, 1.0)
!  do ix = 1, nx
!    if (a2frontloc(ix,iy) .eq. miss) then
!      a2frontspeed(ix,iy) = miss
!      cycle
!    end if 
!    !---
!    call ixy2iixy_saone(ix, iy+1, ixn, iyn)
!
!    call ixy2iixy_saone(ix, iy-1, ixs, iys)
!    call ixy2iixy_saone(ix-1, iy, ixw, iyw)
!    call ixy2iixy_saone(ix+1, iy, ixe, iye)
!    !---
!    a2gradabs_n = a2gradabs(ixn, iyn)
!    a2gradabs_s = a2gradabs(ixs, iys)
!    a2gradabs_w = a2gradabs(ixw, iyw)
!    a2gradabs_e = a2gradabs(ixe, iye)
!    !---
!    a2grad2x(ix, iy) = (a2gradabs_e(ix,iy) - a2gradabs_w(ix,iy)) / (2.0*dew)
!    a2grad2y(ix, iy) = (a2gradabs_n(ix,iy) - a2gradabs_s(ix,iy)) / (dn + ds)
!    a2grad2abs(ix,iy)= ( (a2grad2x(ix,iy))**2.0 + (a2grad2y(ix,iy))**2.0 )**0.5
!    !---
!    speedsign           = -sign(1.0, a2uwind(ix,iy)*a2gradx(ix,iy)   &
!                            +a2vwind(ix,iy)*a2grady(ix,iy))
!
!    a2frontspeed(ix,iy) = abs( a2uwind(ix,iy)*a2grad2x(ix,iy)   &
!                           +a2vwind(ix,iy)*a2grad2y(ix,iy) ) &
!                           / a2grad2abs(ix,iy)               &
!                          * speedsign
!    !---       
!  end do
!end do
!
!return
!END SUBROUTINE mk_a2frontspeed
!!*********************************************************
!SUBROUTINE mk_a2frontogen_def(a2thermo, a2uwind, a2vwind, nx, ny, a2frontogen_def)
!implicit none
!!--- in -------
!integer                  nx, ny
!real,dimension(nx,ny) :: a2thermo
!!f2py intent(in)         a2thermo
!real,dimension(nx,ny) :: a2uwind, a2vwind
!!f2py intent(in)         a2uwind, a2vwind
!!--- out ------
!real,dimension(nx,ny) :: a2frontogen_def
!!f2py intent(out)        a2frontogen_def
!!--- calc -----
!real                     lat, dns, dew
!integer                  ix, iy
!integer                  ixn, ixs, ixw, ixe
!integer                  iyn, iys, iyw, iye
!real                     un, us, uw, ue, vn, vs, vw, ve
!real                     thermon, thermos, thermow, thermoe
!real                     abs_nabla_thermo
!real                     D1, D2
!!----parameter -
!real                  :: lat_first = -89.5
!!--------------
!do iy = 1,ny
!  lat = lat_first + (iy -1)*1.0
!  dns = hubeny_real(lat-1.0, 0.0, lat+1.0, 0.0)
!  dew = hubeny_real(lat, 0.0, lat, 1.0)
!  do ix = 1,nx
!    !---
!    call ixy2iixy_saone(ix, iy+1, ixn, iyn)
!    call ixy2iixy_saone(ix, iy-1, ixs, iys)
!    call ixy2iixy_saone(ix-1, iy, ixw, iyw)
!    call ixy2iixy_saone(ix+1, iy, ixe, iye)
!    !---
!    thermon  = a2thermo(ixn, iyn)
!    thermos  = a2thermo(ixs, iys)
!    thermow  = a2thermo(ixw, iyw)
!    thermoe  = a2thermo(ixe, iye)
!    !
!    un       = a2uwind(ixn,iyn)
!    us       = a2uwind(ixs,iys)
!    uw       = a2uwind(ixw,iyw)
!    ue       = a2uwind(ixe,iye)
!    vn       = a2vwind(ixn,iyn)
!    vs       = a2vwind(ixs,iys)
!    vw       = a2vwind(ixw,iyw)
!    ve       = a2vwind(ixe,iye)
!    !
!    abs_nabla_thermo  =( ((thermow-thermoe)/dew)**2.0 + ((thermon-thermos)/dns)**2.0 )**0.5
!    !
!    D1       = (ue-uw)/dew - (vn-vs)/dns
!    D2       = (ve-vw)/dew + (un-us)/dns
!    !---
!    a2frontogen_def(ix,iy) =-0.5/abs_nabla_thermo    &
!           *(   D1*( ( (thermoe-thermow)/dew )**2.0  &
!                    -( (thermon-thermos)/dns )**2.0 )&
!             +2*D2*(thermoe-thermow)/dew*(thermon-thermos)/dns )
!    !---
!  end do
!end do
!return
!END SUBROUTINE mk_a2frontogen_def
!
!!*********************************************************
!SUBROUTINE mk_a2frontogen_div(a2thermo, a2uwind, a2vwind, nx, ny, a2frontogen_div)
!implicit none
!!--- in -------
!integer                  nx, ny
!real,dimension(nx,ny) :: a2thermo
!!f2py intent(in)         a2thermo
!real,dimension(nx,ny) :: a2uwind, a2vwind
!!f2py intent(in)         a2uwind, a2vwind
!!--- out ------
!real,dimension(nx,ny) :: a2frontogen_div
!!f2py intent(out)        a2frontogen_div
!!--- calc -----
!real                     lat, dns, dew
!integer                  ix, iy
!integer                  ixn, ixs, ixw, ixe
!integer                  iyn, iys, iyw, iye
!real                     uw, ue, vn, vs
!real                     thermon, thermos, thermow, thermoe
!real                     abs_nabla_thermo, div
!!----parameter -
!real                  :: lat_first = -89.5
!!--------------
!do iy = 1,ny
!  lat = lat_first + (iy -1)*1.0
!  dns = hubeny_real(lat-1.0, 0.0, lat+1.0, 0.0)
!  dew = hubeny_real(lat, 0.0, lat, 1.0)
!  do ix = 1,nx
!    !---
!    call ixy2iixy_saone(ix, iy+1, ixn, iyn)
!    call ixy2iixy_saone(ix, iy-1, ixs, iys)
!    call ixy2iixy_saone(ix-1, iy, ixw, iyw)
!    call ixy2iixy_saone(ix+1, iy, ixe, iye)
!    !---
!    thermon  = a2thermo(ixn, iyn)
!    thermos  = a2thermo(ixs, iys)
!    thermow  = a2thermo(ixw, iyw)
!    thermoe  = a2thermo(ixe, iye)
!    !
!    uw       = a2uwind(ixw,iyw)
!    ue       = a2uwind(ixe,iye)
!    vn       = a2vwind(ixn,iyn)
!    vs       = a2vwind(ixs,iys)
!    !
!    abs_nabla_thermo  =( ((thermow-thermoe)/dew)**2.0 + ((thermon-thermos)/dns)**2.0 )**0.5
!    !
!    div      = (ue-uw)/dew + (vn-vs)/dns
!    !---
!    a2frontogen_div(ix,iy)  = -0.5*abs_nabla_thermo*div
!!    a2frontogen_div(ix,iy)  = a2uwind(ix,iy)
!!    !---
!  end do
!end do
!return
!END SUBROUTINE mk_a2frontogen_div
!!*********************************************************
!SUBROUTINE mk_a2relative_vorticity(a2uwind, a2vwind, nx, ny, a2relative_vort)
!implicit none
!!--- in -------
!integer                  nx, ny
!real,dimension(nx,ny) :: a2uwind, a2vwind
!!f2py intent(in)         a2uwind, a2vwind
!!--- out ------
!real,dimension(nx,ny) :: a2relative_vort
!!f2py intent(out)        a2relative_vort
!!--- calc -----
!real                     lat, dns, dew
!integer                  ix, iy
!integer                  ixn, ixs, ixw, ixe
!integer                  iyn, iys, iyw, iye
!real                     vw, ve, un, us
!!----parameter -
!real                  :: lat_first = -89.5
!!--------------
!do iy = 1,ny
!  lat = lat_first + (iy -1)*1.0
!  dns = hubeny_real(lat-1.0, 0.0, lat+1.0, 0.0)
!  dew = hubeny_real(lat, 0.0, lat, 1.0)
!  do ix = 1,nx
!    !---
!    call ixy2iixy_saone(ix, iy+1, ixn, iyn)
!    call ixy2iixy_saone(ix, iy-1, ixs, iys)
!    call ixy2iixy_saone(ix-1, iy, ixw, iyw)
!    call ixy2iixy_saone(ix+1, iy, ixe, iye)
!    !---
!    un = a2uwind(ixn, iyn)
!    us = a2uwind(ixs, iys)
!    vw = a2vwind(ixw, iyw)
!    ve = a2vwind(ixe, iye)
!    !---
!    a2relative_vort(ix,iy) = (ve-vw)/dew - (un-us)/dns
!    !---
!  end do
!end do
!
!return
!
!END SUBROUTINE mk_a2relative_vorticity
!!*********************************************************
!SUBROUTINE mk_a2contour_regional(a2in, v, vtrue_out, miss, nx, ny, a2contour)
!implicit none
!!--- in -------
!integer                  nx, ny
!real,dimension(nx,ny) :: a2in
!!f2py intent(in)         a2in
!real                     v, vtrue_out, miss
!!f2py intent(in)         v, vtrue_out, miss
!!--- out ------
!real,dimension(nx,ny) :: a2contour
!!f2py intent(out)        a2contour
!!--- calc -----
!integer                  k
!integer                  ix, iy
!real                     vn, vs, vw, ve
!real                     vnw, vne, vsw, vse
!real                     vo, vmid
!real,dimension(8)     :: a1v
!!--------------
!a2contour = miss
!a1v  = miss
!vn   = miss
!vs   = miss
!vw   = miss
!ve   = miss
!vnw  = miss
!vne  = miss
!vsw  = miss
!vse  = miss
!do iy = 1,ny
!  do ix = 1,nx
!    !---
!    !---
!    vo  = a2in(ix, iy)
!    if (iy .eq. 1) then
!      vs  = miss
!      vsw = miss
!      vse = miss
!    else if (iy .eq. ny) then
!      vn  = miss
!      vnw = miss
!      vne = miss
!    else if (ix .eq. 1) then
!      vw  = miss
!      vnw = miss
!      vsw = miss
!    else if (ix .eq. nx) then
!      ve  = miss
!      vne = miss
!      vse = miss
!    else
!      vn  = a2in(ix, iy+1)
!      vs  = a2in(ix, iy-1)
!      vw  = a2in(ix-1, iy)
!      ve  = a2in(ix+1, iy)
!      vnw = a2in(ix-1, iy+1)
!      vne = a2in(ix+1, iy+1)
!      vsw = a2in(ix-1, iy-1)
!      vse = a2in(ix+1, iy-1)
!    end if
!    a1v(1) = vnw
!    a1v(2) = vn
!    a1v(3) = vne
!    a1v(4) = vw
!    a1v(5) = ve
!    a1v(6) = vsw
!    a1v(7) = vs
!    a1v(8) = vse
!    !---
!    if (vo .le. v)then
!      do k = 1,8
!        !print *,"vo, a1v(k)",vo, a1v(k)
!        if (a1v(k) .eq. miss) then
!          cycle
!        else if (a1v(k) .ge. v)then
!          vmid = (a1v(k) + v)*0.5
!          if ( v .lt. vmid) then
!            a2contour(ix,iy) = vtrue_out
!            exit
!          end if
!        end if
!      end do
!    else if (vo .gt. v)then
!      do k = 1,8
!        !print *,"vo, a1v(k)",vo, a1v(k)
!        if (a1v(k) .eq. miss) then
!          cycle
!        else if (a1v(k) .ge. v)then
!          vmid = (a1v(k) + v)*0.5
!          if ( v .ge. vmid) then
!            a2contour(ix,iy) = vtrue_out
!            exit
!          end if
!        end if
!      end do
!    end if
!    !---
!  end do
!end do
!
!return
!END SUBROUTINE  mk_a2contour_regional
!!*********************************************************
!
!
!!*********************************************************
!SUBROUTINE mk_a2divergence(a2u, a2v, nx, ny, a2div)
!implicit none
!!--- in -------
!integer                  nx, ny
!real,dimension(nx,ny) :: a2u, a2v
!!f2py intent(in)         a2u, a2v
!!--- out ------
!real,dimension(nx,ny) :: a2div
!!f2py intent(out)        a2div
!!--- calc -----
!real                     lat, dns, dew
!integer                  ix, iy
!integer                  ixn, ixs, ixw, ixe
!integer                  iyn, iys, iyw, iye
!real                     vn, vs, vw, ve
!!----parameter -
!real                  :: lat_first = -89.5
!!--------------
!do iy = 1,ny
!  lat = lat_first + (iy -1)*1.0
!  dns = hubeny_real(lat-1.0, 0.0, lat+1.0, 0.0)
!  dew = hubeny_real(lat, 0.0, lat, 1.0)
!  do ix = 1,nx
!    !---
!    call ixy2iixy_saone(ix, iy+1, ixn, iyn)
!    call ixy2iixy_saone(ix, iy-1, ixs, iys)
!    call ixy2iixy_saone(ix-1, iy, ixw, iyw)
!    call ixy2iixy_saone(ix+1, iy, ixe, iye)
!    !---
!    vn = a2v(ixn, iyn)
!    vs = a2v(ixs, iys)
!    vw = a2u(ixw, iyw)
!    ve = a2u(ixe, iye)
!    !---
!    a2div(ix,iy) = (vn-vs)/dns + (ve-vw)/dew
!    !a2div(ix,iy) = -(vn-vs)/dns - (vw-ve)/dew
!    !---
!  end do
!end do
!
!return
!END SUBROUTINE mk_a2divergence
!!*********************************************************
!SUBROUTINE mk_a2adj_multigrids(a2in, grids, nx, ny, a2adj)
!implicit none
!!--- in -------
!integer                  nx, ny
!real,dimension(nx,ny) :: a2in
!!f2py intent(in)         a2in
!integer                  grids
!!f2py intent(in)         grids
!!--- out ------
!real,dimension(nx,ny) :: a2adj
!!f2py intent(out)        a2adj
!!--- calc -----
!integer                  iy
!real,dimension(nx,ny) :: a2gradinabs
!real                     dns, dew, meangridlen
!real                     coef
!real                     lat
!!--- para --------
!real                    :: lat_first = -89.5
!!-----------------
!CALL mk_a2grad_abs_saone(a2in=a2in, a2gradabs=a2gradinabs)
!!-----------------
!do iy = 1,ny
!  lat = lat_first + (iy -1)*1.0
!  dns = hubeny_real(lat-1.0, 0.0, lat+1.0, 0.0) *0.5
!  dew = hubeny_real(lat, 0.0, lat, 1.0)
!  meangridlen  = (dns + dew)*0.5
!  !coef         = meangridlen / (2**0.5)
!  coef         = meangridlen *grids
!  a2adj(:,iy) = a2in(:,iy) &
!                + coef * a2gradinabs(:,iy)
!
!end do
!
!return
!!--------------
!END SUBROUTINE mk_a2adj_multigrids
!!*********************************************************
!SUBROUTINE mk_a2adj(a2in, nx, ny, a2adj)
!implicit none
!!--- in -------
!integer                  nx, ny
!real,dimension(nx,ny) :: a2in
!!f2py intent(in)         a2in
!!--- out ------
!real,dimension(nx,ny) :: a2adj
!!f2py intent(out)        a2adj
!!--- calc -----
!integer                  iy
!real,dimension(nx,ny) :: a2gradinabs
!real                     dns, dew, meangridlen
!real                     coef
!real                     lat
!!--- para --------
!real                    :: lat_first = -89.5
!!-----------------
!CALL mk_a2grad_abs_saone(a2in=a2in, a2gradabs=a2gradinabs)
!!-----------------
!do iy = 1,ny
!  lat = lat_first + (iy -1)*1.0
!  dns = hubeny_real(lat-1.0, 0.0, lat+1.0, 0.0) *0.5
!  dew = hubeny_real(lat, 0.0, lat, 1.0)
!  meangridlen  = (dns + dew)*0.5
!  !coef         = meangridlen / (2**0.5)
!  coef         = meangridlen 
!  a2adj(:,iy) = a2in(:,iy) &
!                + coef * a2gradinabs(:,iy)
!
!end do
!
!return
!!--------------
!END SUBROUTINE mk_a2adj
!!*********************************************************
!SUBROUTINE mk_a2axisgrad(a2inx, a2iny, miss, nx, ny, a2axisgrad)
!implicit none
!!--- in -------
!integer                  nx, ny
!real,dimension(nx,ny) :: a2inx, a2iny
!!f2py intent(in)         a2inx, a2iny
!real,optional         :: miss
!!--- out ------
!real,dimension(nx,ny) :: a2axisgrad
!!f2py intent(out)        a2axisgrad
!!!--- calc -----
!integer                  ix,  iy
!integer                  ixn, ixs, ixw, ixe
!integer                  iyn, iys, iyw, iye
!real,dimension(nx,ny) :: a2meanunitaxis_x, a2meanunitaxis_y
!real                     lat, dns, dew
!real                     resolvabsn, resolvabss, resolvabsw, resolvabse
!real                     resolvxw, resolvxe
!real                     resolvyn, resolvys
!real                     vxn, vxs, vxw, vxe
!real                     vyn, vys, vyw, vye
!
!!!--- parameter -
!real                  :: lat_first = -89.5
!!!---------------
!!
!if ( present(miss) )then
!  CALL mk_a2meanunitaxis(a2inx=a2inx, a2iny=a2iny, miss=miss, nx=nx, ny=ny&
!       &, a2meanunitaxis_x=a2meanunitaxis_x, a2meanunitaxis_y=a2meanunitaxis_y)
!else
!  CALL mk_a2meanunitaxis(a2inx=a2inx, a2iny=a2iny, nx=nx, ny=ny&
!       &, a2meanunitaxis_x=a2meanunitaxis_x, a2meanunitaxis_y=a2meanunitaxis_y)
!
!end if
!!!!--------------
!do iy = 1,ny
!  lat = lat_first + (iy -1)*1.0
!  dns = hubeny_real(lat-1.0, 0.0, lat+1.0, 0.0)
!  dew = hubeny_real(lat, 0.0, lat, 1.0)
!  do ix = 1, nx
!    !---
!    call ixy2iixy_saone(ix, iy+1, ixn, iyn)
!    call ixy2iixy_saone(ix, iy-1, ixs, iys)
!    call ixy2iixy_saone(ix-1, iy, ixw, iyw)
!    call ixy2iixy_saone(ix+1, iy, ixe, iye)
!    !---
!    vxn = a2inx(ixn, iyn)
!    vxs = a2inx(ixs, iys)
!    vxw = a2inx(ixw, iyw)
!    vxe = a2inx(ixe, iye)
!    vyn = a2iny(ixn, iyn)
!    vys = a2iny(ixs, iys)
!    vyw = a2iny(ixw, iyw)
!    vye = a2iny(ixe, iye)
!    !---
!    resolvabsn = vxn * a2meanunitaxis_x(ix,iy) + vyn * a2meanunitaxis_y(ix,iy)
!    resolvabss = vxs * a2meanunitaxis_x(ix,iy) + vys * a2meanunitaxis_y(ix,iy)
!    resolvabsw = vxw * a2meanunitaxis_x(ix,iy) + vyw * a2meanunitaxis_y(ix,iy)
!    resolvabse = vxe * a2meanunitaxis_x(ix,iy) + vye * a2meanunitaxis_y(ix,iy)
!    !---
!    resolvxw   = resolvabsw * a2meanunitaxis_x(ix,iy)
!    resolvxe   = resolvabse * a2meanunitaxis_x(ix,iy)
!    resolvyn   = resolvabsn * a2meanunitaxis_y(ix,iy)
!    resolvys   = resolvabsw * a2meanunitaxis_y(ix,iy)
!    !---
!    a2axisgrad(ix,iy) = (resolvyn - resolvys)/dns + (resolvxe - resolvxw)/dew
!    !*********************
!    ! check miss
!    !---------------------
!    if ( present(miss) )then
!      if ((a2meanunitaxis_x(ix,iy).eq.miss).or.(a2meanunitaxis_y(ix,iy).eq.miss))then
!        !!!! checking a2meanunitaxis is enough !!!!
!        !!!! see their subroutine
!        a2axisgrad(ix,iy) = miss
!      end if 
!    end if
!    !*********************
!  end do
!end do
!!!--------------
!!
!return
!END SUBROUTINE mk_a2axisgrad
!
!!!*********************************************************
!SUBROUTINE mk_a2meanunitaxis(a2inx, a2iny, miss, nx, ny, a2meanunitaxis_x, a2meanunitaxis_y)
!implicit none
!!--- in ---------
!integer                  nx, ny
!real,dimension(nx,ny) :: a2inx, a2iny
!!f2py intent(in)         a2inx, a2iny
!real,optional         :: miss
!!f2py intent(in)         miss
!!!--- out --------
!real,dimension(nx,ny) :: a2meanunitaxis_x, a2meanunitaxis_y
!!f2py intent(out)        a2meanunitaxis_x, a2meanunitaxis_y
!!--- calc -------
!real                       vxo, vxn, vxs, vxw, vxe
!real                       vyo, vyn, vys, vyw, vye
!!
!real                       vrxo, vrxn, vrxs, vrxw, vrxe
!real                       vryo, vryn, vrys, vryw, vrye
!!
!real                       vrx, vry, radr
!!
!real                       coso, cosn, coss, cosw, cose
!real                       sino, sinn, sins, sinw, sine
!!
!real                       cosr, sinr, cosr_hlf, sinr_hlf
!!
!integer                    ix,  iy
!integer                    ixn, ixs, ixw, ixe
!integer                    iyn, iys, iyw, iye
!!----------------
!do iy = 1,ny
!  do ix = 1,nx
!    !---
!    call ixy2iixy_saone(ix, iy+1, ixn, iyn)
!    call ixy2iixy_saone(ix, iy-1, ixs, iys)
!    call ixy2iixy_saone(ix-1, iy, ixw, iyw)
!    call ixy2iixy_saone(ix+1, iy, ixe, iye)
!    !---
!    vxo = a2inx(ix,  iy)
!    vyo = a2iny(ix,  iy)
!    vxn = a2inx(ixn, iyn)
!    vxs = a2inx(ixs, iys)
!    vxw = a2inx(ixw, iyw)
!    vxe = a2inx(ixe, iye)
!    vyn = a2iny(ixn, iyn)
!    vys = a2iny(ixs, iys)
!    vyw = a2iny(ixw, iyw)
!    vye = a2iny(ixe, iye)
!    !---
!    if ((vxo .eq. 0.0) .and. (vyo .eq. 0.0))then
!      coso = 0.0
!      sino = 0.0
!    else
!      coso = vxo / (vxo**2.0 + vyo**2.0)**0.5
!      sino = vyo / (vxo**2.0 + vyo**2.0)**0.5
!    endif
!    if ((vxn .eq. 0.0) .and. (vyn .eq. 0.0))then
!      cosn = 0.0
!      sinn = 0.0
!    else
!      cosn = vxn / (vxn**2.0 + vyn**2.0)**0.5
!      sinn = vyn / (vxn**2.0 + vyn**2.0)**0.5
!    endif
!    if ((vxs .eq. 0.0) .and. (vys .eq. 0.0))then
!      coss = 0.0
!      sins = 0.0
!    else
!      coss = vxs / (vxs**2.0 + vys**2.0)**0.5
!      sins = vys / (vxs**2.0 + vys**2.0)**0.5
!    endif
!    if ((vxw .eq. 0.0) .and. (vyw .eq. 0.0))then
!      cosw = 0.0
!      sinw = 0.0
!    else
!      cosw = vxw / (vxw**2.0 + vyw**2.0)**0.5
!      sinw = vyw / (vxw**2.0 + vyw**2.0)**0.5
!    endif
!    if ((vxe .eq. 0.0) .and. (vye .eq. 0.0))then
!      cose = 0.0
!      sine = 0.0
!    else
!      cose = vxe / (vxe**2.0 + vye**2.0)**0.5
!      sine = vye / (vxe**2.0 + vye**2.0)**0.5
!    endif
!
!    !---
!    vrxo = coso*vxo - sino*vyo
!    vrxn = cosn*vxn - sinn*vyn
!    vrxs = coss*vxs - sins*vys
!    vrxw = cosw*vxw - sinw*vyw
!    vrxe = cose*vxe - sine*vye
!    !---
!    vryo = sino*vxo + coso*vyo
!    vryn = sinn*vxn + cosn*vyn
!    vrys = sins*vxs + coss*vys
!    vryw = sinw*vxw + cosw*vyw
!    vrye = sine*vxe + cose*vye
!    !---
!    vrx  = vrxo + vrxn + vrxs + vrxw + vrxe 
!    vry  = vryo + vryn + vrys + vryw + vrye
!    !---
!    radr = (vrx**2.0 + vry**2.0)**0.5
!    !---
!    if (radr .gt. 0.0)then
!      cosr = vrx / radr 
!      sinr = vry / radr
!      !---
!      cosr_hlf = ((1+cosr)*0.5)**0.5 * sign(1.0, sinr)
!      sinr_hlf = ((1-cosr)*0.5)**0.5
!      !---
!      !a2meanunitaxis_x(ix,iy) = 0.2* radr * cosr_hlf
!      !a2meanunitaxis_y(ix,iy) = 0.2* radr * sinr_hlf
!      a2meanunitaxis_x(ix,iy) = cosr_hlf
!      a2meanunitaxis_y(ix,iy) = sinr_hlf
!    else
!      a2meanunitaxis_x(ix,iy) = 0.0
!      a2meanunitaxis_y(ix,iy) = 0.0
!    endif
!    !******************
!    ! check miss
!    !------------------
!    if ( present(miss) )then
!      if (& 
!        & (vxo.eq.miss).or.&
!        & (vyo.eq.miss).or.&
!        & (vxn.eq.miss).or.&
!        & (vxs.eq.miss).or.&
!        & (vxw.eq.miss).or.&
!        & (vxe.eq.miss).or.&
!        & (vyn.eq.miss).or.&
!        & (vys.eq.miss).or.&
!        & (vyw.eq.miss).or.&
!        & (vye.eq.miss)) then
!        a2meanunitaxis_x(ix,iy) = miss
!        a2meanunitaxis_y(ix,iy) = miss
!      end if 
!    end if
!    !******************
!  end do
!end do
!
!return
!END SUBROUTINE mk_a2meanunitaxis
!
!!*********************************************************
!SUBROUTINE mk_a2theta(plev, a2T, nx, ny, a2theta)
!implicit none
!!--- in ---------
!integer                  nx, ny
!real,dimension(nx,ny) :: a2T
!!f2py intent(in)         a2T
!real                     plev       ! (Pa)
!!f2py intent(in)         plev
!!--- out --------
!real,dimension(nx,ny) :: a2theta
!!f2py intent(out)        a2theta
!!--- para -------------
!real              :: P1000 = 1000.0*100.0
!real              :: chi   = 0.286  ! Rd/cp
!!--- calc -------
!integer                  ix, iy
!real                     t
!!----------------
!do iy = 1, ny
!  do ix = 1, nx
!    t  = a2T(ix,iy)
!    a2theta(ix,iy) = t*(P1000/plev)**(chi)
!  end do
!end do
!!
!return
!END SUBROUTINE mk_a2theta
!!*********************************************************
!SUBROUTINE mk_a2theta_e(plev, a2T, a2q, nx, ny, a2theta_e)
!implicit none
!!--- in ---------
!integer                  nx, ny
!real,dimension(nx,ny) :: a2T, a2q
!!f2py intent(in)         a2T, a2q
!real                     plev       ! (Pa)
!!f2py intent(in)         plev
!real,dimension(nx,ny) :: a2theta_e
!!f2py intent(out)        a2theta_e
!!--- calc -------
!integer                  ix, iy
!real                     t,  q
!!----------------
!do iy = 1, ny
!  do ix = 1, nx
!    t  = a2T(ix,iy)
!    q  = a2q(ix,iy) 
!    a2theta_e(ix,iy) = cal_theta_e(Plev, t, q)
!  end do
!end do
!!
!return
!END SUBROUTINE mk_a2theta_e
!!*********************************************************
!SUBROUTINE mk_a2wetbulbtheta(plev, a2T, a2q, nx, ny, a2wetbulbtheta)
!implicit none
!!--- in ---------
!integer                  nx, ny
!real,dimension(nx,ny) :: a2T, a2q
!!f2py intent(in)         a2T, a2q
!real                     plev       ! (Pa)
!!f2py intent(in)         plev
!!--- out --------
!real,dimension(nx,ny) :: a2wetbulbtheta
!!f2py intent(out)        a2wetbulbtheta
!!--- calc -------
!integer                  ix, iy
!real                     t, q
!!---------------------------------
!do iy = 1, ny
!  do ix = 1, nx
!    t  = a2T(ix,iy)
!    q  = a2q(ix,iy) 
!    a2wetbulbtheta(ix,iy) = cal_wetbulbtheta(Plev, t, q)
!  end do
!end do
!return
!END SUBROUTINE mk_a2wetbulbtheta
!
!!*********************************************************
!SUBROUTINE mk_a2slope_abs_saone(a2in, a2slopeabs)
!!---------------------------------
!! data order should be South -> North, West -> East
!! returns slope, not gradient
!!---------------------------------
!implicit none
!!--- in ----------
!integer                 :: ny = 180
!integer                 :: nx = 360
!real,dimension(360,180) :: a2in
!!f2py intent(in)           a2in
!!--- out ---------
!real,dimension(360,180) :: a2slopeabs
!!f2py intent(out)          a2slopeabs
!!--- para --------
!real                    :: lat_first = -89.5
!!--- calc --------
!real                       dn, ds, dew
!real                       vn, vs, vw, ve
!real                       lat
!real                       slopen, slopes, slopee, slopew
!integer                    ix,  iy
!integer                    ixn, ixs, ixw, ixe
!integer                    iyn, iys, iyw, iye
!!-----------------
!do iy = 1, ny
!  lat = lat_first + (iy -1)*1.0
!  dn  = hubeny_real(lat, 0.0, lat+1.0, 0.0)
!  ds  = hubeny_real(lat, 0.0, lat-1.0, 0.0)
!  dew = hubeny_real(lat, 0.0, lat, 1.0)
!  do ix = 1, nx
!    !---
!    call ixy2iixy_saone(ix, iy+1, ixn, iyn)
!    call ixy2iixy_saone(ix, iy-1, ixs, iys)
!    call ixy2iixy_saone(ix-1, iy, ixw, iyw)
!    call ixy2iixy_saone(ix+1, iy, ixe, iye)
!    !---
!    vn = a2in(ixn, iyn)
!    vs = a2in(ixs, iys)
!    vw = a2in(ixw, iyw)
!    ve = a2in(ixe, iye)
!    !---
!    slopen= abs((vn - a2in(ix,iy))/dn)
!    slopes= abs((vs - a2in(ix,iy))/ds)
!    slopee= abs((ve - a2in(ix,iy))/dew)
!    slopew= abs((vw - a2in(ix,iy))/dew)
!    a2slopeabs(ix, iy) = (slopen + slopes + slopew + slopee)/4.0
!    !---
!  end do
!end do
!!-----------------
!return
!END SUBROUTINE mk_a2slope_abs_saone

!!*********************************************************
!FUNCTION cal_theta_e(P, T, q)
!implicit none
!!----------------------------
!! estimate equivalent potential temperature
!!----------------------------
!!--- in ---------------
!real                 P, T, q  ! (Pa), (K), (kg/kg: mixing ratio)
!!f2py intent(in)     P, T, q
!!--- out --------------
!real                 cal_theta_e
!!f2py intent(out)    cal_theta_e
!!--- para -------------
!real              :: P1000 = 1000.0*100.0
!real              :: chi   = 0.286  ! Rd/cp
!real              :: Lv    = 2.5e6  ![J kg-1]
!real              :: cp    = 1004.0 ![J kg-1 K-1]
!!--- calc -------------
!real                 Tlcl, theta
!!----------------------
!theta   = T * (P1000/P)**chi
!Tlcl    = cal_tlcl(P, T, q)
!cal_theta_e = theta * exp(Lv*q/(cp*Tlcl))
!!
!return
!END FUNCTION cal_theta_e
!!*********************************************************
!FUNCTION cal_wetbulbtheta(P, T, q)
!implicit none
!!----------------------------
!! estimate wet-bulb potential temperature
!!----------------------------
!!--- in -------------
!real                  P, T, q  ! (Pa), (K), (kg/kg: mixing ratio)
!!f2py intent(in)      P, T, q
!!--- out ------------
!real                  cal_wetbulbtheta
!!f2py intent(out)     cal_wetbulbtheta
!!--- para -----------
!real               :: P1000 = 1000.0*100.0
!!--- calc -----------
!real                  Plcl, Tlcl
!!--------------------
!Plcl  = cal_plcl(P, T, q)
!Tlcl  = cal_tlcl(P, T, q)
!!
!cal_wetbulbtheta = t1_to_t2_moistadia(Plcl, P1000, Tlcl)
!!
!return
!END FUNCTION cal_wetbulbtheta
!!*********************************************************
!SUBROUTINE a2t1_to_a2t2(P1,P2, a2t1, a2q1, miss, nx, ny, a2t2)
!implicit none
!integer                   nx, ny
!!--- in -------------
!real                      P1, P2
!!f2py intent(in)          P1, P2
!real,dimension(nx,ny)  :: a2t1, a2q1
!!f2py intent(in)          a2t1, a2q1
!real                      miss
!!f2py intent(in)          miss
!!--- out ------------
!real,dimension(nx,ny)  :: a2t2
!!f2py intent(out)         a2t2
!!--- para -----------
!
!!--- calc -----------
!integer                   ix,iy
!real                      t1,q1
!!--------------------
!a2t2  = miss
!do iy = 1,ny
!  do ix = 1,nx
!    t1 = a2t1(ix,iy)
!    q1 = a2q1(ix,iy)
!    if ((t1.eq.miss).or.(q1.eq.miss))then
!      cycle
!    else
!      a2t2(ix,iy) = t1_to_t2(P1,P2, t1, q1)
!    end if
!  end do
!end do
!!--------------------
!return
!END SUBROUTINE a2t1_to_a2t2
!!*********************************************************
!FUNCTION t1_to_t2(P1,P2, T1, q1)
!implicit none
!!--- in -------------
!real                 P1, P2, T1, q1  ! (Pa), (K), (kg/kg)
!!f2py intent(in)     P1, P2, T1, q1
!!--- out ------------
!real                 t1_to_t2
!!f2py intent(out)    t1_to_t2
!!--- para -----------
!real               :: Cpd   = 1004.0 !(J kg^-1 K^-1)
!real               :: Rd    = 287.04  !(J kg^-1 K^-1)
!!--- calc -----------
!real                 Plcl, Tlcl 
!!--------------------
!t1_to_t2 = 0.0
!
!Tlcl  = cal_tlcl(p1,T1,q1)
!Plcl  = P1*(Tlcl/T1)**(Cpd/Rd)
!if (Plcl.gt.P2)then
!  t1_to_t2  = t1_to_t2_moistadia(Plcl, P2, Tlcl)
!else if (Plcl.le.P2)then
!  t1_to_t2  = T1*(P2/P1)**(Rd/Cpd)
!end if
!!--------------------
!return
!END FUNCTION t1_to_t2
!!!!*********************************************************
!FUNCTION t1_to_t2_moistadia(P1,P2, T1)
!implicit none
!!--- in -------------
!real                  P1, P2, T1  ! pressure:(Pa), temperature:(K)
!!f2py intent(in)      P1, P2, T1
!!--- out ------------
!real                  t1_to_t2_moistadia
!!f2py intent(out)     t1_to_t2_moistadia
!!--- para -----------
!real               :: thres = 0.01
!integer            :: imax  = 50
!!--- calc -----------
!real                  Tnow, Tnext, d
!integer               i
!!--------------------
!Tnow = T1
!Tnext = 0.0
!do i = 1, imax
!  call t1_to_t2_moistadia_sub(P1, P2, T1, Tnow, Tnext, d)
!  if ( abs(d) .le. thres )then
!    exit
!  end if
!  Tnow = Tnext
!end do
!t1_to_t2_moistadia  = Tnext
!!
!return
!END FUNCTION t1_to_t2_moistadia
!!!*********************************************************
!SUBROUTINE t1_to_t2_moistadia_sub(P1,P2, T1, Tnow, Tnext, d)
!implicit none
!!--- in ----------
!real                  P1, P2, T1, Tnow  ! pressure:(Pa), temerature:(K)
!!f2py intent(in)      P1, P2, T1, Tnow
!!--- out ---------
!real                  Tnext, d
!!f2py intent(out)     Tnext, d
!!--- para --------
!real               :: Cpd   = 1004.0 !(J kg^-1 K^-1)
!real               :: Rd    = 287.04  !(J kg^-1 K^-1)
!real               :: epsi  = 0.62185 ! = Rd/Rv = Mv/Md
!real               :: Lv    = 2.5e6   ! for vaporization
!real               :: P1000 = 1000.0*100.0
!!--- calc --------
!real                  f, df, dqs
!real                  rA
!real                  es, qs
!!-----------------
!es           = cal_es_water(Tnow)
!qs           = cal_qs_water(Tnow, P2)
!f            = cal_sateqtheta(P2, Tnow) - cal_sateqtheta(P1, T1)
!rA           = (P1000/(P2-es))**(Rd/Cpd)
!dqs          = epsi**2.0 * Lv*P2*es / (Rd *(P2-es)**2.0 *Tnow**2.0)
!df           = rA*exp(Lv*qs/(Cpd*Tnow)) &
!               + rA*Tnow* (Lv*(dqs*Tnow - qs)/(Cpd*Tnow**2.0)) &
!                 *exp(Lv*qs/(Cpd*Tnow))
!!
!Tnext        = Tnow - f/df
!d            = cal_sateqtheta(P2, Tnext) - cal_sateqtheta(P1, T1)
!!-----------------
!return
!END SUBROUTINE t1_to_t2_moistadia_sub
!!!*********************************************************
!FUNCTION cal_sateqtheta(P, T)
!implicit none
!!--- in --------
!real                  P, T  ! (Pa), (K)
!!f2py intent(in)      P, T
!!--- out -------
!real                  cal_sateqtheta
!!f2py intent(out)     cal_sateqtheta
!!--- para ------
!real               :: Cpd   = 1004.0 !(J kg^-1 K^-1)
!real               :: Rd    = 287.04  !(J kg^-1 K^-1)             
!real               :: Lv    = 2.5e6   ! for vaporization 
!real               :: P1000 = 1000.0*100.0  !(Pa)
!!--- calc ------
!real                  drytheta
!real                  es, qs
!!---------------
!es         = cal_es_water(T)
!qs         = cal_qs_water(T, P)
!drytheta   = T*( P1000/(P-es))**(Rd/Cpd)
!cal_sateqtheta  = drytheta * exp((Lv*qs)/(Cpd*T))
!
!
!END FUNCTION cal_sateqtheta
!!!*********************************************************
!
!FUNCTION cal_plcl(P0, T0, q0)
!!!------------------------------------
!!! Following the procedure from Yoshizaki and Kato (2007)
!!! "Go-u Go-setsu no Kisho-gaku", A-3
!!!------------------------------------
!implicit none
!!---- in ----------------
!real                  P0, T0, q0  ! P:(Pa)  T:(K)  q:(kg/kg)
!!f2py intent(in)      P0, T0, q0
!!---- out ---------------
!real                  cal_plcl    ! (Pa)
!!---- calc --------------
!real                  Plcl, Tlcl
!!---- para --------------
!real               :: Cpd   = 1004.0 !(J kg^-1 K^-1)
!real               :: Rd    = 287.04  !(J kg^-1 K^-1)
!!------------------------
!Tlcl      = cal_tlcl(P0, T0, q0)
!Plcl      = P0*(Tlcl/T0)**(Cpd/Rd)
!cal_plcl  = Plcl
!!------------------------
!END FUNCTION cal_plcl
!
!!**************************************************************
!FUNCTION cal_tlcl(P0, T0, q0)
!!!------------------------------------
!!! Following the procedure from Yoshizaki and Kato (2007)
!!! "Go-u Go-setsu no Kisho-gaku", A-3
!!!------------------------------------
!implicit none
!!---- in  ---------------
!real                  P0, T0, q0   ! P:(Pa)  T:(K)  q:(kg/kg)
!!f2py intent(in)      P0, T0, q0
!!---- out ---------------
!real                  cal_tlcl     ! (K)
!!f2py intent(out)     cal_tlcl
!!---- calc --------------
!real                  qs0
!real                  Plcl, Tlcl
!real                  Ttemp, Ptemp
!!---- para --------------
!real               :: Cpd   = 1004.0 !(J kg^-1 K^-1)
!real               :: Rd    = 287.04  !(J kg^-1 K^-1)
!!------------------------
!!-- check saturation level ---
!qs0   = cal_qs_water(T0, P0)
!if (qs0 .le. q0)then
!  cal_tlcl = T0
!  return
!endif
!!-- FIRST guess ------
!Tlcl  = cal_tlcl_sub(P0, T0, q0)
!Plcl  = P0*(Tlcl/T0)**(Cpd/Rd)
!
!!-- SECOND guess ------
!Ttemp = Tlcl
!Ptemp = Plcl
!
!Tlcl  = cal_tlcl_sub(Ptemp, Ttemp, q0)
!
!cal_tlcl = Tlcl
!return
!END FUNCTION cal_tlcl
!
!!**************************************************************
!FUNCTION cal_tlcl_sub(P0, T0, q0)
!implicit none
!!---- input  -------------
!real                  P0, T0, q0   ! (Pa), (K), (kg/kg)
!!---- output -------------
!real                  cal_tlcl_sub
!!---- calc   -------------
!real                  qs0
!real                  rA, rB, rC
!real                  bunshi, bunbo
!!---- para   -------------
!real               :: Cpd   = 1004.0 !(J kg^-1 K^-1
!real               :: Rd    = 287.04  !(J kg^-1 K^-1)
!real               :: Lv    = 2.5e6   ! for vaporization
!real               :: epsi  = 0.62185 ! = Rd/Rv = Mv/Md
!!-------------------------
!qs0        = cal_qs_water(T0, P0)
!
!rA         = 2.0 * (qs0 - q0)/ (qs0 + q0)
!rB         = Cpd * (epsi + (qs0 + q0)*0.5 )/(epsi*Rd)
!rC         = Lv  * (epsi + (qs0 + q0)*0.5 )/ Rd
!
!bunshi     = rB*T0 + rC - ( (rB*T0 - rC)**2.0 + 2.0*rA*rC*T0)**0.5
!bunbo      = 2.0*rB - rA
!
!cal_tlcl_sub= 2.0*bunshi/bunbo - T0
!
!return
!END FUNCTION cal_tlcl_sub  
!
!!*********************************************************************
!FUNCTION cal_qs_water(rT, rP)
!  implicit none
!  real                 rT, rP
!!f2py intent(in)       rT, rP
!  real                 res
!  real                 cal_qs_water
!!f2py intent(out)      cal_qs_water
!  real,parameter    :: repsi = 0.62185
!!
!res = cal_es_water(rT)
!cal_qs_water = repsi * res / (rP - res)
!RETURN
!END FUNCTION cal_qs_water
!!*********************************************************************
!FUNCTION cal_es_water(rT)
!  real rT
!  real cal_es_water
!!
!  real,parameter            ::  rT0 = 273.16
!  real,parameter            ::  res0= 611.73 ![Pa]
!  real,parameter            ::  Lv  = 2.5e6  ![J kg-1]
!  real,parameter            ::  Rv  = 461.7 ![J K-1 kg -1]
!!
!cal_es_water = res0 * exp( Lv/Rv *(1.0/rT0 - 1.0/rT))
!RETURN
!END FUNCTION cal_es_water
!!*********************************************************************
!
!!!**************************************************************
!!FUNCTION lcl_old(rPsfc, rTsfc, rqsfc)
!!!###########################################################
!!! original code was obtained from
!!! http://www1.doshisha.ac.jp/~jmizushi/program/Fortran90/4.2.txt
!!! modified by: N.Utsumi
!!! f(x)=x**3+6*x**2+21*x+32
!!!###########################################################
!!implicit none
!!
!!real                  rPsfc, rTsfc, rqsfc   ! rPsfc:[Pa]
!!!f2py intent(in)      rPsfc, rTsfc, rqsfc
!!real                  lcl_old
!!!f2py intent(out)     lcl_old
!!
!!double precision      dPsfc_hPa, dTsfc, dq
!!double precision      x, xk, fx
!!double precision      delta
!!integer               k
!!INTEGER,PARAMETER :: KMAX=200
!!!-------------
!!!Psfc = 1000   !(hPa)
!!!Tsfc = 293.15 !(K)
!!!q    = 0.0087268029 !(kg/kg)
!!!-------------
!!dPsfc_hPa = dble(rPsfc)*0.01d0  ! Pa -> hPa
!!dTsfc = dble(rTsfc)
!!dq    = dble(rqsfc)
!!!-------------
!!x=1000.d0
!!delta=1.d-10
!!!-------------
!!fx=func(x, dPsfc_hPa, dTsfc, dq)
!!k=0
!!!WRITE(*,"('x(',i2,')=',1PE15.8,', f(',i2,')=',1PE15.8)") k,x,k,fx
!!!WRITE(*,*)
!!
!!DO k=1,KMAX
!!
!!xk=fnewton(x, dPsfc_hPa, dTsfc, dq)
!!fx=func(xk, dPsfc_hPa, dTsfc, dq)
!!
!!!WRITE(*,"('x(',i2,')=',1PE15.8,', f(',i2,')=',1PE15.8)") k,xk,k,fx
!!
!!    IF(abs(fx)<delta)GOTO 100
!!
!!x=xk    ! LCL [hPa]
!!
!!END DO
!!
!!WRITE(*,*) 'could not solve.'
!!print *, "Psfc=",dPsfc_hPa
!!print *, "Tsfc=",dTsfc
!!print *, "q=",dq
!!print *, "fx=",fx
!!if (.not.isnan(x)) then
!!  STOP
!!endif
!!
!!100 CONTINUE
!!!
!!if (isnan(x) ) then
!!  lcl_old = x    ! lcl = nan
!!else
!!  lcl_old = real(x) *100.0  ! [hPa] -> [Pa]
!!endif
!!!-----------------
!!! for the case: lcl is lower than the surface (RH > 100%)
!!!-----------------
!!if (-lcl_old .lt. -rPsfc) then
!!  lcl_old = rPsfc
!!endif
!!!-----------------
!!return
!!END FUNCTION lcl_old
!!**************************************************************
!FUNCTION func(P, Psfc, Tsfc, q)
!  implicit none
!  double precision      P, Psfc, Tsfc, q
!  double precision      f1, f2, func
!  double precision      L
!!
!  double precision :: T0    = 273.16  !(K)
!  double precision :: e0    = 6.1173  !(hPa)
!  double precision :: Rv    = 461.7   !(J kg^-1 K^-1)
!  !double precision :: Lv    = 2.500d6 !(J kg^-1)
!  double precision :: epsi  = 0.62185 !(-)
!  double precision :: Rd    = 287.04  !(J kg^-1 K^-1)
!  double precision :: Cpd   = 1004.0 !(J kg^-1 K^-1)
!!
!L = dble(cal_latentheat( real(Tsfc) ))
!f1 = (1/T0 - Rv/L *log( q * P /( e0*(epsi + q) ) ) )**(-1)
!f2 = Tsfc * ( P / Psfc )**(Rd/Cpd)
!func = f1 - f2
!RETURN
!END FUNCTION func
!!**************************************************************
!FUNCTION fnewton(P, Psfc, Tsfc, q)
!  implicit none
!  double precision       P, Psfc, Tsfc, q
!  double precision       f1, f2, func
!  double precision       df1_P, df2_P, df_P
!  double precision       fnewton
!
!!
!  double precision    L
!  double precision :: T0    = 273.16  !(K)
!  double precision :: e0    = 6.1173  !(hPa)
!  double precision :: Rv    = 461.7   !(J kg^-1 K^-1)
!  !double precision :: Lv    = 2.500d6 !(J kg^-1)
!  double precision :: epsi  = 0.62185 !(-)
!  double precision :: Rd    = 287.04  !(J kg^-1 K^-1)
!  double precision :: Cpd   = 1004.0 !(J kg^-1 K^-1)
!!
!L = dble(cal_latentheat( real(Tsfc) ))
!f1 = (1/T0 - Rv/L *log( q * P /( e0*(epsi + q) ) ) )**(-1)
!f2 = Tsfc * ( P / Psfc )**(Rd/Cpd)
!func = f1 - f2
!!
!df1_P = 1/P * Rv/L *(1/T0 - Rv/L*log( q*P /(e0*(epsi + q)) ) )**(-2)
!df2_P = Tsfc* (1/Psfc)**(Rd/Cpd) * Rd/Cpd * (P **(Rd/Cpd -1))
!df_P  = df1_P - df2_P
!!
!fnewton = P - func / df_P
!RETURN
!END FUNCTION fnewton
!!**************************************************************
!FUNCTION cal_latentheat(rT)
!  implicit none
!  real                  rT
!  real,parameter     :: Lv = 2.5e6  ! for vaporization
!  real,parameter     :: Ld = 2.834e6 ! for sublimation
!  real,parameter     :: rTliq = 273.15  !   0 deg.C
!  real,parameter     :: rTice = 250.15   ! -23 deg.C
!  real               cal_latentheat
!!
!if ( rT .ge. rTliq) then
!  cal_latentheat = Lv
!else if ( rT .le. rTice ) then
!  cal_latentheat = Ld
!else
!  cal_latentheat = ((rT - rTice)*Lv + (rTliq - rT)*Ld)/(rTliq - rTice)
!end if
!RETURN
!END FUNCTION cal_latentheat
!
!!**************************************************************
!!*********************************************************
END MODULE front_fsub
