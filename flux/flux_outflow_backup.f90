PROGRAM flux_outflow
    USE netcdf
    implicit none
    integer :: ncid, varid, dimid, nx, ny, nz
    integer :: i,j,k,v_j,v_istr,v_iend
    integer :: tdimid,d,n
    integer :: year,nt,days,date(3) 
    real(kind=8) :: hc, dx, myfx, a, b, S
    real(kind=8), allocatable :: h(:,:),cs_w(:),zeta(:,:,:),depth(:,:,:) 
    real(kind=8), allocatable :: mask(:,:),s_w(:),depth_psi(:,:,:)
    real(kind=8), allocatable :: v(:,:,:,:),pm(:,:),myflux(:,:),flux(:)
    character(len=100) :: path
    character(len=200) :: filename
    logical :: leapy

    101 format(a,'ocean_avg_',i4.4,'-',i2.2,'-',i2.2,'.nc')

    path = '/users/magazyn/mmuzyka/ROMS_CASES/new_metro_560x600_02/run/baltic/'

    year = 2018
    nt = 2

!    v_j = 201
!    v_istr = 57
!    v_iend = 87
    v_j = 80
    v_istr = 1
    v_iend = 120

    CALL leap(year,leapy) 
    days = 365
    if (leapy) days = 366

    allocate(flux(days))

    date = (/year,1,1/)
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
    allocate( v(nx,ny-1,nz,nt), pm(nx,ny), myflux(nx-2,ny-1) )

    CALL check(nf90_inq_varid(ncid,"hc",varid),317)
    CALL check(nf90_get_var(ncid,varid,hc),318)

    CALL check(nf90_inq_varid(ncid,"Cs_w",varid),319)
    CALL check(nf90_get_var(ncid,varid,cs_w),320)

    CALL check(nf90_inq_varid(ncid,"s_w",varid),319)
    CALL check(nf90_get_var(ncid,varid,s_w),320)

    CALL check(nf90_inq_varid(ncid,"h",varid),319)
    CALL check(nf90_get_var(ncid,varid,h),320)

    CALL check(nf90_inq_varid(ncid,"pm",varid),319)
    CALL check(nf90_get_var(ncid,varid,pm),320)

    CALL check(nf90_inq_varid(ncid,"zeta",varid),319)
    CALL check(nf90_get_var(ncid,varid,zeta),320)

    CALL check(nf90_inq_varid(ncid,"v",varid),319)
    CALL check(nf90_get_var(ncid,varid,v),320)

    CALL check(nf90_inq_varid(ncid,"mask_rho",varid),319)
    CALL check(nf90_get_var(ncid,varid,mask),320)

    CALL check(nf90_close(ncid),360)
!END READ GRID


    do d = 1,days

    CALL check(nf90_open(trim(filename),NF90_NOWRITE,ncid),310)
    CALL check(nf90_inq_varid(ncid,"zeta",varid),319)
    CALL check(nf90_get_var(ncid,varid,zeta),320)
    CALL check(nf90_inq_varid(ncid,"v",varid),319)
    CALL check(nf90_get_var(ncid,varid,v),320)
    CALL check(nf90_close(ncid),360)

    myfx = 0.0
    do n = 1,nt

    where(v(:,:,:,n).gt.100) v(:,:,:,n)=0.0
    where(mask.lt.0.5) zeta(:,:,n) = 0.0

    do i = 1, nx   
      do j = 1, ny   
        do k = 1, nz+1   
          depth(i,j,k) = zeta(i,j,n)+(zeta(i,j,n)+h(i,j))*(hc*s_w(k)+h(i,j)*cs_w(k))/(hc+h(i,j))
        end do
      end do
    end do

    do i = 1, nx-1
      do j = 1, ny-1
        depth_psi(i,j,:) = 0.25*(depth(i,j,:)+depth(i+1,j,:)+&
            depth(i,j+1,:)+depth(i+1,j+1,:))
      end do
    end do

    do i = v_istr, v_iend
      dx = 0.5*(1/pm(i,v_j)+1/pm(i,v_j+1))
      do k = 1, nz
        a = depth_psi(i-1,v_j,k+1)-depth_psi(i-1,v_j,k) 
        b = depth_psi(i,v_j,k+1)-depth_psi(i,v_j,k) 
        S = 0.5*(a+b)*dx
       
        if (v(i,v_j,k,n).gt.0.0) myfx = myfx + S*v(i,v_j,k,n)
      end do
    end do
    
    end do

    flux(d) = myfx/nt

    CALL add_day(date,.true.)
    write(filename,101) trim(path),date(1),date(2),date(3)
    write(*,*) trim(filename)

    end do
    

    CALL check( nf90_create( 'outflow.nc',NF90_NETCDF4,ncid ), 500 )
    CALL check( nf90_def_dim( ncid, 'nt', days, tdimid ), 501 )
    CALL check(nf90_def_var( ncid, 'outflow', NF90_DOUBLE,&
        (/tdimid/), varid ), 508)
    CALL check( nf90_put_att( ncid, varid, 'long_name',&
        'outflow' ), 305 )
    CALL check( nf90_put_att( ncid, varid, 'units',&
        'm**3/s' ), 306 )
    CALL check( nf90_enddef(ncid), 516 )
    CALL check( nf90_put_var( ncid, varid, flux ), 517 )
    CALL check( nf90_close(ncid), 600 ) 


    deallocate( h, cs_w, zeta, depth, mask, s_w, depth_psi, v, pm, myflux )
    deallocate( flux )
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
