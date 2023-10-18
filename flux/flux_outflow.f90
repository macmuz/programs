PROGRAM flux_outflow
    USE netcdf
    implicit none
    integer :: ncid, varid, dimid, nx, ny, nz, varid2
    integer :: i,j,k,v_j,v_istr,v_iend,day
    integer :: u_i,u_jstr,u_jend,varid3
    integer :: tdimid,d,n,cnt,ntimes,win,nr
    integer :: nt,date(3),enddate(3),maxrec,start(3)
    real(kind=8) :: hc, dx, myfx, a, b, S, time
    real(kind=8), allocatable :: h(:,:),cs_w(:),zeta(:,:,:),depth(:,:,:) 
    real(kind=8), allocatable :: mask(:,:),s_w(:),depth_psi(:,:,:),pn(:,:)
    real(kind=8), allocatable :: v(:,:,:,:),pm(:,:),myflux(:,:),flux(:)
    real(kind=8), allocatable :: flux_sum(:),rflux(:),rflux_sum(:),rtmp(:)
    real(kind=8), allocatable :: u(:,:,:,:),flux2(:),flux2_sum(:)
    character(len=100) :: path,rpath
    character(len=200) :: filename
    logical :: leapy

    101 format(a,'ocean_avg_',i4.4,'-',i2.2,'-',i2.2,'.nc')
    102 format(a,'rivers_560x600_',i4.4,'_copy.nc')

    path = '/users/magazyn/mmuzyka/ROMS_CASES/new_metro_560x600_02/run/baltic/'
    rpath = '/users/work/mmuzyka/CSDIR/input_560x600/'

    nt = 2

    v_j = 80
    v_istr = 1
    v_iend = 120

    u_i = 120
    u_jstr = 50
    u_jend = 90

    time = 86400.0/nt

    date = (/2016,7,1/)
    enddate = (/2019,6,29/)
!    enddate = (/2016,7,30/)
    win = 365
    cnt = 0
    do
      cnt = cnt+1
      if(all(date.eq.enddate)) EXIT
      CALL add_day(date,.true.) 
    end do
    write(*,*) cnt
    ntimes = cnt+1-win
    if (ntimes.lt.1) ntimes = 1

    allocate(flux(cnt),flux_sum(ntimes))
    allocate(flux2(cnt),flux2_sum(ntimes))
    allocate(rflux(cnt),rflux_sum(ntimes))

    date = (/2016,7,1/)
    write(filename,101) trim(path),date(1),date(2),date(3)
    write(*,*) trim(filename)

!READ GRID
    CALL check(nf90_open(trim(filename),NF90_NOWRITE,ncid),310)

    CALL check(nf90_inq_dimid(ncid, "xi_rho", dimid),311)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nx),312)

    CALL check(nf90_inq_dimid(ncid, "eta_rho", dimid),313)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=ny),314)

    CALL check(nf90_inq_dimid(ncid, "s_rho", dimid),315)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nz),316)

    allocate( h(nx,ny), cs_w(nz+1), zeta(nx,ny,nt), depth(nx,ny,nz+1) )
    allocate( mask(nx,ny), s_w(nz+1), depth_psi(nx-1,ny-1,nz+1) )
    allocate( v(nx,ny-1,nz,nt), pm(nx,ny), pn(nx,ny), myflux(nx-2,ny-1) )
    allocate( u(nx-1,ny,nz,nt) )

    CALL check(nf90_inq_varid(ncid,"hc",varid),317)
    CALL check(nf90_get_var(ncid,varid,hc),318)

    CALL check(nf90_inq_varid(ncid,"Cs_w",varid),319)
    CALL check(nf90_get_var(ncid,varid,cs_w),321)

    CALL check(nf90_inq_varid(ncid,"s_w",varid),319)
    CALL check(nf90_get_var(ncid,varid,s_w),322)

    CALL check(nf90_inq_varid(ncid,"h",varid),319)
    CALL check(nf90_get_var(ncid,varid,h),323)

    CALL check(nf90_inq_varid(ncid,"pm",varid),319)
    CALL check(nf90_get_var(ncid,varid,pm),324)

    CALL check(nf90_inq_varid(ncid,"pn",varid),319)
    CALL check(nf90_get_var(ncid,varid,pn),325)

    CALL check(nf90_inq_varid(ncid,"zeta",varid),319)
    CALL check(nf90_get_var(ncid,varid,zeta),326)

    CALL check(nf90_inq_varid(ncid,"u",varid),319)
    CALL check(nf90_get_var(ncid,varid,u),327)

    CALL check(nf90_inq_varid(ncid,"v",varid),319)
    CALL check(nf90_get_var(ncid,varid,v),328)

    CALL check(nf90_inq_varid(ncid,"mask_rho",varid),319)
    CALL check(nf90_get_var(ncid,varid,mask),329)

    CALL check(nf90_close(ncid),360)
!END READ GRID


    do d = 1,cnt

    CALL check(nf90_open(trim(filename),NF90_NOWRITE,ncid),310)
    CALL check(nf90_inq_varid(ncid,"zeta",varid),319)
    CALL check(nf90_get_var(ncid,varid,zeta),320)
    CALL check(nf90_inq_varid(ncid,"v",varid),319)
    CALL check(nf90_get_var(ncid,varid,v),320)
    CALL check(nf90_inq_varid(ncid,"u",varid),319)
    CALL check(nf90_get_var(ncid,varid,u),320)
    CALL check(nf90_close(ncid),360)

    myfx = 0.0
    flux(d) = 0.0
    flux2(d) = 0.0
    do n = 1,nt

    where(mask.lt.0.5) zeta(:,:,n) = 0.0
    CALL psi_lvl(nx,ny,nz,zeta(:,:,n),h,hc,s_w,cs_w,depth_psi)

    where(v(:,:,:,n).gt.100) v(:,:,:,n)=0.0
    CALL vflux(nx,ny,nz,v_istr,v_iend,v_j,pm,depth_psi,v(:,:,:,n),0,myfx)
    flux(d) = flux(d) + myfx*time
    
    where(u(:,:,:,n).gt.100) u(:,:,:,n)=0.0
    CALL uflux(nx,ny,nz,u_jstr,u_jend,u_i,pn,depth_psi,u(:,:,:,n),0,myfx)
    flux2(d) = flux2(d) + myfx*time

    end do

    flux(d) = flux(d)*1e-9
    flux2(d) = flux2(d)*1e-9

    CALL add_day(date,.true.)
    write(filename,101) trim(path),date(1),date(2),date(3)
    write(*,*) trim(filename)

    end do

    date = (/2016,7,1/)
    start = (/date(1),1,1/)
    day = 0
    do
      day = day+1
      if(all(date.eq.start)) EXIT
      CALL add_day(start,.true.)
    end do

    write(filename,102) trim(rpath),date(1)
    date(1) = date(1)+1
    write(*,*) trim(filename)
    CALL check(nf90_open(trim(filename),NF90_NOWRITE,ncid),410)
    
    CALL check(nf90_inq_dimid(ncid, "river", dimid),311)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nr),312)

    CALL check(nf90_inq_dimid(ncid, "river_time", dimid),313)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=maxrec),314)

    allocate(rtmp(nr))

    do i = 1,cnt 
      if (day.gt.maxrec) then
        day = 1
        CALL check(nf90_close(ncid),460)
        write(filename,102) trim(rpath),date(1)
        date(1) = date(1)+1
        write(*,*) trim(filename)
        CALL check(nf90_open(trim(filename),NF90_NOWRITE,ncid),410)
        CALL check(nf90_inq_dimid(ncid, "river_time", dimid),313)
        CALL check(nf90_inquire_dimension(ncid, dimid, len=maxrec),314)
      end if 
      
      CALL check(nf90_inq_varid(ncid,"river_transport",varid),319)
      CALL check(nf90_get_var(ncid,varid,rtmp,&
        start=(/1,day/),count=(/nr,1/)),320)
      rflux(i) = sum(rtmp)*86400*1e-9 

      write(*,*) day 
      day = day+1
    end do    

    CALL check(nf90_close(ncid),460) 


    if (ntimes.eq.1) then
      flux_sum(1) = sum(flux)
      flux2_sum(1) = sum(flux2)
      rflux_sum(1) = sum(rflux)
    else
      do i = 1,ntimes
        flux_sum(i) = sum(flux(i:i+win-1)) 
        flux2_sum(i) = sum(flux2(i:i+win-1)) 
        rflux_sum(i) = sum(rflux(i:i+win-1)) 
      end do
    end if
    

    CALL check( nf90_create( 'outflow.nc',NF90_NETCDF4,ncid ), 500 )
    CALL check( nf90_def_dim( ncid, 'nt', ntimes, tdimid ), 501 )
    CALL check(nf90_def_var( ncid, 'outflow', NF90_DOUBLE,&
        (/tdimid/), varid ), 508)
    CALL check( nf90_put_att( ncid, varid, 'long_name',&
        'outflow' ), 305 )
    CALL check( nf90_put_att( ncid, varid, 'units',&
        'km**3' ), 306 )
    CALL check(nf90_def_var( ncid, 'rivers', NF90_DOUBLE,&
        (/tdimid/), varid2 ), 508)
    CALL check( nf90_put_att( ncid, varid2, 'long_name',&
        'river run-off' ), 305 )
    CALL check( nf90_put_att( ncid, varid2, 'units',&
        'km**3' ), 306 )
    CALL check(nf90_def_var( ncid, 'outflow2', NF90_DOUBLE,&
        (/tdimid/), varid3 ), 508)
    CALL check( nf90_put_att( ncid, varid3, 'long_name',&
        'outflow2' ), 305 )
    CALL check( nf90_put_att( ncid, varid3, 'units',&
        'km**3' ), 306 )
    CALL check( nf90_enddef(ncid), 516 )
    CALL check( nf90_put_var( ncid, varid, flux_sum ), 517 )
    CALL check( nf90_put_var( ncid, varid2, rflux_sum ), 517 )
    CALL check( nf90_put_var( ncid, varid3, flux2_sum ), 517 )
    CALL check( nf90_close(ncid), 600 ) 


    deallocate( h, cs_w, zeta, depth, mask, s_w, depth_psi, u, v, myflux )
    deallocate( flux, flux_sum, pm, pn, rflux, rflux_sum, rtmp )
    deallocate( flux2, flux2_sum )
END PROGRAM flux_outflow

SUBROUTINE check(status,label)
    USE NETCDF
    implicit none
    INTEGER, INTENT(in) :: status , label

    IF (status .ne. nf90_noerr) then
        print*, trim(nf90_strerror(status)),' label= ',label
        stop "Stopped"
    END IF
END SUBROUTINE check

SUBROUTINE leap(year,answer)
    implicit none
    integer, intent(in) :: year
    logical, intent(out) :: answer

    answer = .false.
    if( mod(year,4).eq.0 ) answer = .true.
    if( mod(year,4).eq.100 ) answer = .false.
    if( mod(year,4).eq.400 ) answer = .true.
END SUBROUTINE leap

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

SUBROUTINE psi_lvl(nx,ny,nz,zeta,h,hc,s_w,cs_w,depth_psi)
    implicit none
    integer, intent(in) :: nx,ny,nz
    real(kind=8), intent(in) :: zeta(nx,ny),h(nx,ny)
    real(kind=8), intent(in) :: hc,s_w(nz+1),cs_w(nz+1)
    real(kind=8), intent(out) :: depth_psi(nx-1,ny-1,nz+1)

    integer :: i,j,k
    real(kind=8) :: depth(nx,ny,nz+1)
        
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
END SUBROUTINE psi_lvl

SUBROUTINE uflux(nx,ny,nz,jstr,jend,i,pn,depth_psi,u,dir,flx)
    implicit none
    integer, intent(in) :: nx,ny,nz,jstr,jend,i,dir
    real(kind=8), intent(in) :: pn(nx,ny),depth_psi(nx-1,ny-1,nz+1)
    real(kind=8), intent(in) :: u(nx-1,ny,nz)
    real(kind=8), intent(out) :: flx

    integer :: j,k
    real(kind=8) :: dy,a,b,S

    flx = 0.0

    do j = jstr, jend
      dy = 0.5*(1/pn(i,j)+1/pn(i+1,j))
      do k = 1, nz
        a = depth_psi(i,j-1,k+1)-depth_psi(i,j-1,k)
        b = depth_psi(i,j,k+1)-depth_psi(i,j,k)
        S = 0.5*(a+b)*dy
    
        if (dir.eq.0) then
          flx = flx + S*u(i,j,k)
        elseif (dir.gt.0) then
          if ( u(i,j,k).gt.0.0 ) flx = flx + S*u(i,j,k)
        elseif (dir.lt.0) then
          if ( u(i,j,k).lt.0.0 ) flx = flx + S*u(i,j,k)
        end if
      end do
    end do
END SUBROUTINE uflux

SUBROUTINE vflux(nx,ny,nz,istr,iend,j,pm,depth_psi,v,dir,flx)
    implicit none
    integer, intent(in) :: nx,ny,nz,istr,iend,j,dir
    real(kind=8), intent(in) :: pm(nx,ny),depth_psi(nx-1,ny-1,nz+1)
    real(kind=8), intent(in) :: v(nx,ny-1,nz)
    real(kind=8), intent(out) :: flx 

    integer :: i,k
    real(kind=8) :: dx,a,b,S

    flx = 0.0

    do i = istr, iend
      dx = 0.5*(1/pm(i,j)+1/pm(i,j+1))
      do k = 1, nz
        a = depth_psi(i-1,j,k+1)-depth_psi(i-1,j,k)
        b = depth_psi(i,j,k+1)-depth_psi(i,j,k)
        S = 0.5*(a+b)*dx
    
        if (dir.eq.0) then
          flx = flx + S*v(i,j,k)
        elseif (dir.gt.0) then
          if ( v(i,j,k).gt.0.0 ) flx = flx + S*v(i,j,k)
        elseif (dir.lt.0) then
          if ( v(i,j,k).lt.0.0 ) flx = flx + S*v(i,j,k)
        end if
      end do
    end do
END SUBROUTINE vflux
