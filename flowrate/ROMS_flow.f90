PROGRAM ROMS_flow
    USE netcdf
    implicit none
    integer :: ncid, varid, dimid, nx, ny, nz,nc
    integer :: i,j,k,ui,ujstr(2),ujend(2),vid(5),n,cnt
    integer :: xdimid,ydimid,zdimid,tdimid,dimids2(2)
    integer :: date(3),datestp(3),nrec,pts(2,2) 
    real(kind=8) :: hc, dy, myfx, a, b, S, time, tmptime(4), dif, flow(4)
    real(kind=8), allocatable :: h(:,:),cs_w(:),zeta(:,:),depth(:,:,:) 
    real(kind=8), allocatable :: mask(:,:),s_w(:),depth_psi(:,:,:),mask_u(:,:)
    real(kind=8), allocatable :: u(:,:,:), pn(:,:), myflux(:,:)
    real(kind=8), allocatable :: tmpzeta(:,:,:),tmpu(:,:,:,:)
    character(len=200) :: filename,path
    character(len=50) :: tunits,tcal

    100 format(a,'/ocean_avg_',i4,'-',i2.2,'-',i2.2,'.nc')
    path = '/users/work/mmuzyka/CSDIR/metro_560x600_era2004v12/run/baltic'
!    path = '/users/magazyn/mmuzyka/ROMS_CASES/metro_560x600_era2004'

    ui = 159
    ujstr(1) = 76
    ujend(1) = 120
    ujstr(2) = 20
    ujend(2) = 73

    nrec = 4

    date = (/2003,11,1/)
    datestp = (/2003,11,25/)

    write(filename,100) trim(path),date(1),date(2),date(3)
!READ GRID
    CALL check(nf90_open(trim(filename),NF90_NOWRITE,ncid),310)

    CALL check(nf90_inq_dimid(ncid, "xi_rho", dimid),311)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nx),312)

    CALL check(nf90_inq_dimid(ncid, "eta_rho", dimid),313)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=ny),314)

    CALL check(nf90_inq_dimid(ncid, "s_rho", dimid),315)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nz),316)

    allocate( h(nx,ny), cs_w(nz+1), zeta(nx,ny), depth(nx,ny,nz+1) )
    allocate( mask(nx,ny), s_w(nz+1), depth_psi(nx-1,ny-1,nz+1) )
    allocate( u(nx-1,ny,nz), pn(nx,ny), myflux(nx-2,ny-1), mask_u(nx-1,ny) )
    allocate( tmpzeta(nx,ny,nrec), tmpu(nx-1,ny,nz,nrec) )

    CALL check(nf90_inq_varid(ncid,"hc",varid),317)
    CALL check(nf90_get_var(ncid,varid,hc),318)

    CALL check(nf90_inq_varid(ncid,"Cs_w",varid),319)
    CALL check(nf90_get_var(ncid,varid,cs_w),320)

    CALL check(nf90_inq_varid(ncid,"s_w",varid),319)
    CALL check(nf90_get_var(ncid,varid,s_w),320)

    CALL check(nf90_inq_varid(ncid,"h",varid),319)
    CALL check(nf90_get_var(ncid,varid,h),320)

    CALL check(nf90_inq_varid(ncid,"pn",varid),319)
    CALL check(nf90_get_var(ncid,varid,pn),320)

    CALL check(nf90_inq_varid(ncid,"zeta",varid),319)
    CALL check(nf90_get_var(ncid,varid,zeta,start=(/1,1,2/)),320)

    CALL check(nf90_inq_varid(ncid,"u",varid),319)
    CALL check(nf90_get_var(ncid,varid,u,start=(/1,1,1,2/)),320)

    CALL check(nf90_inq_varid(ncid,"mask_rho",varid),319)
    CALL check(nf90_get_var(ncid,varid,mask),320)

    CALL check(nf90_inq_varid(ncid,"mask_u",varid),319)
    CALL check(nf90_get_var(ncid,varid,mask_u),320)

    CALL check(nf90_inq_varid(ncid,"ocean_time",varid),319)
    CALL check(nf90_get_att(ncid,varid,'units',tunits),320)
    CALL check(nf90_get_att(ncid,varid,'calendar',tcal),320)

    CALL check(nf90_close(ncid),360)
!END READ GRID

    CALL check( nf90_create( 'ROMS_flow_era2004v12.nc',NF90_NETCDF4,nc ), 500 )

    CALL check( nf90_def_dim( nc, 'time', NF90_UNLIMITED, tdimid ), 501 )

    CALL check(nf90_def_var( nc, 'time', NF90_DOUBLE,&
        (/tdimid/), vid(1) ), 508)
    CALL check( nf90_put_att( nc, vid(1), 'units',&
        'seconds since 1968-05-23 00:00:00' ), 305 )
    CALL check( nf90_put_att( nc, vid(1), 'standard_name',&
        'time' ), 306 )

    CALL check(nf90_def_var( nc, 'north', NF90_DOUBLE,&
        (/tdimid/), vid(2) ), 508)
    CALL check( nf90_put_att( nc, vid(2), 'longname',&
        'net flow on north section' ), 305 )
    CALL check( nf90_put_att( nc, vid(2), 'units',&
        'km**3'), 305 )

    CALL check(nf90_def_var( nc, 'south', NF90_DOUBLE,&
        (/tdimid/), vid(3) ), 508)
    CALL check( nf90_put_att( nc, vid(3), 'longname',&
        'net flow on south section' ), 305 )
    CALL check( nf90_put_att( nc, vid(3), 'units',&
        'km**3'), 305 )

    CALL check(nf90_def_var( nc, 'north_abs', NF90_DOUBLE,&
        (/tdimid/), vid(4) ), 508)
    CALL check( nf90_put_att( nc, vid(4), 'longname',&
        'total flow on north section' ), 305 )
    CALL check( nf90_put_att( nc, vid(4), 'units',&
        'km**3'), 305 )

    CALL check(nf90_def_var( nc, 'south_abs', NF90_DOUBLE,&
        (/tdimid/), vid(5) ), 508)
    CALL check( nf90_put_att( nc, vid(5), 'longname',&
        'total flow on south section' ), 305 )
    CALL check( nf90_put_att( nc, vid(5), 'units',&
        'km**3'), 305 )

    cnt = 1
    do
        write(filename,100) trim(path),date(1),date(2),date(3)
        write(*,*) trim(filename)

        CALL check(nf90_open(trim(filename),NF90_NOWRITE,ncid),310)

        CALL check(nf90_inq_varid(ncid,"ocean_time",varid),319)
        CALL check(nf90_get_var(ncid,varid,tmptime),320)

        time = 0.5*(tmptime(2)+tmptime(3))

        CALL check(nf90_inq_varid(ncid,"u",varid),319)
        CALL check(nf90_get_var(ncid,varid,tmpu),320)

        u = sum(tmpu,4)/real(nrec)

        CALL check(nf90_inq_varid(ncid,"zeta",varid),319)
        CALL check(nf90_get_var(ncid,varid,tmpzeta),320)

        zeta = sum(tmpzeta,3)/real(nrec)

        CALL check(nf90_close(ncid),360)

    do k = 1,nz
      where(mask_u.eq.0) u(:,:,k)=0.0
    end do
    
    where(mask.lt.0.5) zeta = 0.0
    do i = 1, nx   
      do j = 1, ny   
        do k = 1, nz+1   
          depth(i,j,k) = zeta(i,j)+(zeta(i,j)+h(i,j))*(hc*s_w(k)+h(i,j)*cs_w(k))/(hc+h(i,j))
        end do
      end do
    end do

    do i = 1, nx-1
      do j = 1, ny-1
        depth_psi(i,j,:) = 0.25*(depth(i,j,:)+depth(i+1,j,:)+&
            depth(i,j+1,:)+depth(i+1,j+1,:))
      end do
    end do

    flow = 0.0
    do j = ujstr(1), ujend(1)
      dy = 0.5*(1/pn(ui,j)+1/pn(ui+1,j))
      do k = 1, nz
        a = depth_psi(ui,j-1,k+1)-depth_psi(ui,j-1,k)
        b = depth_psi(ui,j,k+1)-depth_psi(ui,j,k)
        S = 0.5*(a+b)*dy
        flow(1) = flow(1)+S*u(ui,j,k)
        flow(3) = flow(3)+S*abs(u(ui,j,k))
      end do 
    end do
    do j = ujstr(2), ujend(2)
      dy = 0.5*(1/pn(ui,j)+1/pn(ui+1,j))
      do k = 1, nz
        a = depth_psi(ui,j-1,k+1)-depth_psi(ui,j-1,k)
        b = depth_psi(ui,j,k+1)-depth_psi(ui,j,k)
        S = 0.5*(a+b)*dy
        flow(2) = flow(2)+S*u(ui,j,k)
        flow(4) = flow(4)+S*abs(u(ui,j,k))
      end do
    end do 

!    write(*,*) v(v_istr-1,v_j,1),v(v_istr,v_j,1),v(v_iend,v_j,1),v(v_iend+1,v_j,1)


!    myflux = 0.0
!    do i = 1, nx-2
!      do j = 1, ny-1
!        dx = 0.5*(1/pm(i+1,j)+1/pm(i+1,j+1))
!        do k = 1, nz
!          a = depth_psi(i,j,k+1)-depth_psi(i,j,k)
!          b = depth_psi(i+1,j,k+1)-depth_psi(i+1,j,k)
!          S = 0.5*(a+b)*dx
!          myflux(i,j) = myflux(i,j) + S*v(i+1,j,k)
!        end do
!      end do
!    end do

!    CALL check( nf90_create( 'flux.nc',NF90_NETCDF4,ncid ), 500 )
!    CALL check( nf90_def_dim( ncid, 'nx', nx-2, xdimid ), 501 )
!    CALL check( nf90_def_dim( ncid, 'ny', ny-1, ydimid ), 502 )
!    dimids2 = (/ xdimid, ydimid /)
!    CALL check(nf90_def_var( ncid, 'nflux', NF90_DOUBLE,&
!        dimids2, varid ), 508)
!    CALL check( nf90_put_att( ncid, varid, 'long_name',&
!        'flux towards the north' ), 305 )
!    CALL check( nf90_put_att( ncid, varid, 'units',&
!        'm**3/s' ), 306 )
!    CALL check( nf90_enddef(ncid), 516 )
!    CALL check( nf90_put_var( ncid, varid, myflux ), 517 )
!    CALL check( nf90_close(ncid), 600 )


      write(*,*) time
        
      CALL check( nf90_put_var( nc, vid(1), time, start=(/cnt/)), 517 )
      CALL check( nf90_put_var( nc, vid(2), flow(1)*1e-9, start=(/cnt/)), 517 )
      CALL check( nf90_put_var( nc, vid(3), flow(2)*1e-9, start=(/cnt/)), 517 )
      CALL check( nf90_put_var( nc, vid(4), flow(3)*1e-9, start=(/cnt/)), 517 )
      CALL check( nf90_put_var( nc, vid(5), flow(4)*1e-9, start=(/cnt/)), 517 )       
 
        cnt = cnt+1


        if (date(1).ge.datestp(1).and.&
            date(2).ge.datestp(2).and.&
            date(3).ge.datestp(3)) EXIT
        CALL add_day(date,.true.)
    end do
    CALL check( nf90_close(nc), 600 ) 


    deallocate( h, cs_w, zeta, depth, mask, s_w, depth_psi, u, pn, myflux, mask_u )
    deallocate( tmpu, tmpzeta )

END PROGRAM ROMS_flow

SUBROUTINE check(status,label)
    USE NETCDF
    implicit none
    INTEGER, INTENT(in) :: status , label

    IF (status .ne. nf90_noerr) then
        print*, trim(nf90_strerror(status)),' label= ',label
        stop "Stopped"
    END IF
END SUBROUTINE check

SUBROUTINE add_day(date,use_leap)
    implicit none
    integer, intent(inout) :: date(3)
    logical, intent(in) :: use_leap

    logical :: leap
    integer :: days_in_month(12)

    days_in_month = (/31,28,31,30,31,30,31,31,30,31,30,31/)

    if ( date(2).eq.2 .and. date(3).eq.28 .and. use_leap) then
        leap = .false.
        if ( mod(date(1),4).eq.0 ) leap = .true.
        if ( mod(date(1),100).eq.0 ) leap = .false.
        if ( mod(date(1),400).eq.0 ) leap = .true.

        if ( leap ) days_in_month(2) = 29
    end if

    date(3) = date(3)+1
    if ( date(3).gt.days_in_month(date(2)) ) then
        date(3) = 1
        date(2) = date(2)+1
        if ( date(2).gt.12 ) then
            date(2) = 1
            date(1) = date(1)+1
        end if
    end if

END SUBROUTINE add_day

SUBROUTINE dens(temp,salt,depth,den)
    implicit none
    real(kind=8), intent(in) :: temp,salt,depth
    real(kind=8), intent(out) :: den

    integer :: i,j,k
    real(kind=8) :: Tp, Tpr10, Ts, Tt, sqrtTs, cff
    real(kind=8), dimension(0:9) :: C
    real(kind=8) :: den1,bulk,bulk0,bulk1,bulk2

    real(kind=8), parameter :: A00 = +1.909256e+04
    real(kind=8), parameter :: A01 = +2.098925e+02
    real(kind=8), parameter :: A02 = -3.041638e+00
    real(kind=8), parameter :: A03 = -1.852732e-03
    real(kind=8), parameter :: A04 = -1.361629e-05
    real(kind=8), parameter :: B00 = +1.044077e+02
    real(kind=8), parameter :: B01 = -6.500517e+00
    real(kind=8), parameter :: B02 = +1.553190e-01
    real(kind=8), parameter :: B03 = +2.326469e-04
    real(kind=8), parameter :: D00 = -5.587545e+00
    real(kind=8), parameter :: D01 = +7.390729e-01
    real(kind=8), parameter :: D02 = -1.909078e-02
    real(kind=8), parameter :: E00 = +4.721788e-01
    real(kind=8), parameter :: E01 = +1.028859e-02
    real(kind=8), parameter :: E02 = -2.512549e-04
    real(kind=8), parameter :: E03 = -5.939910e-07
    real(kind=8), parameter :: F00 = -1.571896e-02
    real(kind=8), parameter :: F01 = -2.598241e-04
    real(kind=8), parameter :: F02 = +7.267926e-06
    real(kind=8), parameter :: G00 = +2.042967e-03
    real(kind=8), parameter :: G01 = +1.045941e-05
    real(kind=8), parameter :: G02 = -5.782165e-10
    real(kind=8), parameter :: G03 = +1.296821e-07
    real(kind=8), parameter :: H00 = -2.595994e-07
    real(kind=8), parameter :: H01 = -1.248266e-09
    real(kind=8), parameter :: H02 = -3.508914e-09
    real(kind=8), parameter :: Q00 = +9.99842594e+02
    real(kind=8), parameter :: Q01 = +6.793952e-02
    real(kind=8), parameter :: Q02 = -9.095290e-03
    real(kind=8), parameter :: Q03 = +1.001685e-04
    real(kind=8), parameter :: Q04 = -1.120083e-06
    real(kind=8), parameter :: Q05 = +6.536332e-09
    real(kind=8), parameter :: U00 = +8.24493e-01
    real(kind=8), parameter :: U01 = -4.08990e-03
    real(kind=8), parameter :: U02 = +7.64380e-05
    real(kind=8), parameter :: U03 = -8.24670e-07
    real(kind=8), parameter :: U04 = +5.38750e-09
    real(kind=8), parameter :: V00 = -5.72466e-03
    real(kind=8), parameter :: V01 = +1.02270e-04
    real(kind=8), parameter :: V02 = -1.65460e-06
    real(kind=8), parameter :: W00 = +4.8314e-04


          Tt=MAX(-2.0,temp)
          Ts=MAX(0.0,salt)
          sqrtTs=SQRT(Ts)
          Tp=depth
          Tpr10=0.1*Tp


!
!-----------------------------------------------------------------------
!  Compute density (kg/m3) at standard one atmosphere pressure.
!-----------------------------------------------------------------------
!
          C(0)=Q00+Tt*(Q01+Tt*(Q02+Tt*(Q03+Tt*(Q04+Tt*Q05))))
          C(1)=U00+Tt*(U01+Tt*(U02+Tt*(U03+Tt*U04)))
          C(2)=V00+Tt*(V01+Tt*V02)
          den1=C(0)+Ts*(C(1)+sqrtTs*C(2)+Ts*W00)

!
!-----------------------------------------------------------------------
!  Compute secant bulk modulus.
!-----------------------------------------------------------------------
!
          C(3)=A00+Tt*(A01+Tt*(A02+Tt*(A03+Tt*A04)))
          C(4)=B00+Tt*(B01+Tt*(B02+Tt*B03))
          C(5)=D00+Tt*(D01+Tt*D02)
          C(6)=E00+Tt*(E01+Tt*(E02+Tt*E03))
          C(7)=F00+Tt*(F01+Tt*F02)
          C(8)=G01+Tt*(G02+Tt*G03)
          C(9)=H00+Tt*(H01+Tt*H02)

          bulk0=C(3)+Ts*(C(4)+sqrtTs*C(5))
          bulk1=C(6)+Ts*(C(7)+sqrtTs*G00)
          bulk2=C(8)+Ts*C(9)
          bulk =bulk0-Tp*(bulk1-Tp*bulk2)
!
!-----------------------------------------------------------------------
!  Compute local "in situ" density anomaly (kg/m3 - 1000).
!-----------------------------------------------------------------------
!
          cff=1.0/(bulk+Tpr10)
          den=den1*bulk*cff

END SUBROUTINE dens

