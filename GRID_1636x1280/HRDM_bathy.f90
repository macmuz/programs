PROGRAM HRDM_bathy
    USE NETCDF
    implicit none
    integer, parameter :: imt=5001
    integer :: nx,ny,ncid,dimid,vid,dx,dy,i,j,m
    integer :: xdim,ydim,dimids2(2),varid(3)
    integer :: istr(1),iend(1),jstr(1),jend(1)
    real(kind=8), allocatable :: y1d(:),y21d(:),selev(:,:),slon(:),slat(:)
    real(kind=8), allocatable :: elev(:,:),hi_lon(:),hi_lat(:),y2(:,:)
    real(kind=8) :: lon(imt), lat(imt), h(imt,imt)
    real(kind=8) :: splint
    character(len=20) :: f_elev,f_lonlat

    f_elev = 'D5_2020.nc'
    f_lonlat = 'lonlat.cdf'

    CALL check(nf90_open(trim(f_elev),NF90_NOWRITE,ncid),10)
    CALL check(nf90_inq_dimid(ncid, "lon", dimid),11)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nx),12)
    CALL check(nf90_inq_dimid(ncid, "lat", dimid),13)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=ny),14)
    allocate( elev(nx,ny), hi_lon(nx), hi_lat(ny) )
    CALL check(nf90_inq_varid(ncid,"elevation",vid),17)
    CALL check(nf90_get_var(ncid,vid,elev),18)
    CALL check(nf90_inq_varid(ncid,"lon",vid),19)
    CALL check(nf90_get_var(ncid,vid,hi_lon),20)
    CALL check(nf90_inq_varid(ncid,"lat",vid),21)
    CALL check(nf90_get_var(ncid,vid,hi_lat),22)
    CALL check(nf90_close(ncid),35)

    CALL check(nf90_open(trim(f_lonlat),NF90_NOWRITE,ncid),10)
    CALL check(nf90_inq_varid(ncid,"LON",vid),19)
    CALL check(nf90_get_var(ncid,vid,lon),20)
    CALL check(nf90_inq_varid(ncid,"LAT",vid),21)
    CALL check(nf90_get_var(ncid,vid,lat),22)
    CALL check(nf90_close(ncid),35)

    write(*,*) minval(lon),maxval(lon)
    write(*,*) minval(lat),maxval(lat)

    istr = minloc(abs(hi_lon-lon(1)))-5
    iend = minloc(abs(hi_lon-lon(imt)))+5
    jstr = minloc(abs(hi_lat-lat(1)))-5
    jend = minloc(abs(hi_lat-lat(imt)))+5

    dx = 1+iend(1)-istr(1)
    dy = 1+jend(1)-jstr(1)

    write(*,*) istr,iend,dx
    write(*,*) jstr,jend,dy
   
    write(*,*) hi_lon(istr(1)),hi_lon(iend(1)) 
    write(*,*) hi_lat(jstr(1)),hi_lat(jend(1))

    allocate( y2(dx,dy), y1d(dx), y21d(dx) )
    allocate( selev(dx,dy), slon(dx), slat(dy) )

    selev = elev(istr(1):iend(1),jstr(1):jend(1))
    slon = hi_lon(istr(1):iend(1))
    slat = hi_lat(jstr(1):jend(1))

    write(*,*) minval(slon),maxval(slon)
    write(*,*) minval(slat),maxval(slat)
    write(*,*) minval(selev),maxval(selev)

    CALL splie2( slat, dx, dy, selev, y2)

    do j = 1,imt
      do m = 1,dx
        y1d(m) = splint( slat, selev(m,:), y2(m,:), dy, lat(j) )
        CALL spline( slon, y1d, dx, y21d )
        do i = 1,imt
          h(i,j) = splint( slon, y1d, y21d, dx, lon(i) )
        end do 
      end do
      write(*,*) j 
    end do 

    CALL check( nf90_create( 'HRDM_elev.nc',NF90_NETCDF4,ncid ), 300 )
    CALL check( nf90_def_dim( ncid, 'lon', imt, xdim ), 301 )
    CALL check( nf90_def_dim( ncid, 'lat', imt, ydim ), 301 )
    dimids2 = (/ xdim, ydim /)
    CALL check(nf90_def_var( ncid, 'h', NF90_DOUBLE,&
        dimids2, varid(1) ), 304)
    CALL check( nf90_put_att( ncid, varid(1), 'units',&
        'meter' ), 314 )
    CALL check(nf90_def_var( ncid, 'lon', NF90_DOUBLE,&
        (/xdim/), varid(2) ), 304)
    CALL check( nf90_put_att( ncid, varid(2), 'long_name',&
        'longitude' ), 314 )
    CALL check( nf90_put_att( ncid, varid(2), 'units',&
        'degree_east' ), 314 )
    CALL check(nf90_def_var( ncid, 'lat', NF90_DOUBLE,&
        (/ydim/), varid(3) ), 304)
    CALL check( nf90_put_att( ncid, varid(3), 'long_name',&
        'latitude' ), 314 )
    CALL check( nf90_put_att( ncid, varid(3), 'units',&
        'degree_north' ), 314 )
    CALL check( nf90_enddef(ncid), 341 )
    CALL check( nf90_put_var( ncid, varid(1), h ), 352 )
    CALL check( nf90_put_var( ncid, varid(2), lon ), 352 )
    CALL check( nf90_put_var( ncid, varid(3), lat ), 352 )
    CALL check( nf90_close(ncid), 353 )

    deallocate( elev, hi_lon, hi_lat, y2, y1d, y21d )
    deallocate( selev, slon, slat )

END PROGRAM HRDM_bathy

SUBROUTINE check(status,label)
    USE NETCDF
    implicit none
    INTEGER, INTENT(in) :: status , label

    IF (status .ne. nf90_noerr) then
        print*, trim(nf90_strerror(status)),' label=',label
        stop "Stopped"
    END IF
END SUBROUTINE check

SUBROUTINE spline(x,y,n,y2)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    REAL(kind=8), DIMENSION(n), INTENT(IN) :: x,y
    REAL(kind=8), DIMENSION(n), INTENT(OUT) :: y2

    REAL(kind=8), DIMENSION(n) :: a,b,c,r

    c(1:n-1)=x(2:n)-x(1:n-1)
    r(1:n-1)=6.0*((y(2:n)-y(1:n-1))/c(1:n-1))
    r(2:n-1)=r(2:n-1)-r(1:n-2)
    a(2:n-1)=c(1:n-2)
    b(2:n-1)=2.0*(c(2:n-1)+a(2:n-1))
    b(1)=1.0
    b(n)=1.0
    r(1)=0.0
    c(1)=0.0
    r(n)=0.0
    a(n)=0.0
    call tridag(n,a(2:n),b(1:n),c(1:n-1),r(1:n),y2(1:n))
END SUBROUTINE spline

FUNCTION splint(xa,ya,y2a,n,x)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    REAL(kind=8), DIMENSION(n),  INTENT(IN) :: xa,ya,y2a
    REAL(kind=8), INTENT(IN) ::  x

    INTEGER :: khi,klo,locate,jl,jm,ju
    REAL(kind=8) :: a,b,h,splint
    LOGICAL :: ascnd

    ascnd=(xa(n)>=xa(1))
    jl=0
    ju=n+1
    do
        if (ju-jl<=1) exit
        jm=(ju+jl)/2
        if (ascnd .eqv. (x>=xa(jm))) then
            jl=jm
        else
            ju=jm
        end if
    end do

    if (x==xa(1)) then
        locate=1
    else if (x==xa(n)) then
        locate=n-1
    else
        locate=jl
    end if

    klo=max(min(locate,n-1),1)
    khi=klo+1
    h=xa(khi)-xa(klo)
    if (h == 0.0) then
        write(*,*) 'bad xa input in splint'
        STOP
    end if
    a=(xa(khi)-x)/h
    b=(x-xa(klo))/h
    splint=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0
END FUNCTION splint

SUBROUTINE splie2(x2a,m,n,ya,y2a)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: m,n
    REAL(kind=8), DIMENSION(n), INTENT(IN) :: x2a
    REAL(kind=8), DIMENSION(m,n), INTENT(IN) :: ya
    REAL(kind=8), DIMENSION(m,n), INTENT(OUT) :: y2a

    INTEGER :: j

    DO j=1,m
        call spline(x2a,ya(j,:),n,y2a(j,:))
    END DO
END SUBROUTINE splie2

SUBROUTINE tridag(n,a,b,c,r,u)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    REAL(kind=8), DIMENSION(n), INTENT(IN) :: b,r
    REAL(kind=8), DIMENSION(n-1), INTENT(IN) :: a,c
    REAL(kind=8), DIMENSION(n), INTENT(OUT) :: u

    REAL(kind=8) :: gam(n),bet
    INTEGER :: j

    bet=b(1)

    if (bet==0.0) then
        write(*,*) 'tridag: Error at code stage 1'
        STOP
    end if

    u(1)=r(1)/bet
    do j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j-1)*gam(j)
        if (bet==0.0) then
            write(*,*) 'tridag: Error at code stage 2'
            STOP
        end if
        u(j)=(r(j)-a(j-1)*u(j-1))/bet
    end do
    do j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
    end do
END SUBROUTINE tridag
