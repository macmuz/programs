PROGRAM getWAVE
    USE netcdf
    implicit none
    integer :: idx(2),x,y
    integer :: ncid,dimid,ncid2,tdimid
    integer :: ni,nj,nt,vid,varid(2),p,t
    character(len=200) :: wavefile
    character(len=30) :: coord
    real(kind=8) :: W(2,2),pt(2)
    real(kind=8), allocatable :: lon(:),lat(:)
    real(kind=8), allocatable :: data3d(:,:,:),time(:),output(:)

    101 format(f7.4,'N, ',f8.4,'E')
    wavefile = 'WAVE.nc'
    
    
    pt = (/16.5,55.05/)
    write(coord,101) pt(2),pt(1) 


    CALL check(nf90_open(trim(wavefile),NF90_NOWRITE,ncid),310)

    CALL check(nf90_inq_dimid(ncid, "lon", dimid),311)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=ni),312)

    CALL check(nf90_inq_dimid(ncid, "lat", dimid),313)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nj),314)

    CALL check(nf90_inq_dimid(ncid, "time", dimid),315)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nt),316)

    allocate( lon(ni), lat(nj) )

    CALL check(nf90_inq_varid(ncid,'lon',vid),119)
    CALL check(nf90_get_var(ncid,vid,lon),18)
    CALL check(nf90_inq_varid(ncid,'lat',vid),121)
    CALL check(nf90_get_var(ncid,vid,lat),18)

    do p = 1,ni
     if (pt(1).lt.lon(p)) then
        idx(1) = p-1
        exit
      end if
    end do
    do p = 1,nj
      if (pt(2).lt.lat(p)) then
        idx(2) = p-1
        exit
      end if
    end do
    x = idx(1)
    y = idx(2)
    write(*,*) pt
    write(*,*) x,y
    call calc_w( lon(x:x+1), lat(y:y+1), pt, W )
    write(*,*) lon(x),lat(y),W(1,1)
    write(*,*) lon(x+1),lat(y),W(2,1)
    write(*,*) lon(x),lat(y+1),W(1,2)
    write(*,*) lon(x+1),lat(y+1),W(2,2)
    write(*,*) sum(W)

    CALL check(nf90_create('vhm0_pt.nc',NF90_NETCDF4,ncid2),310)
    CALL check(nf90_inq_varid(ncid,'time',vid),122)

    CALL check( nf90_def_dim( ncid2, 'time', NF90_UNLIMITED, tdimid ),503)
    CALL check( nf90_def_var( ncid2, 'time', NF90_DOUBLE, (/tdimid/), varid(1) ), 504)
    CALL check( nf90_copy_att(ncid, vid, 'units', ncid2, varid(1)), 504)
    CALL check( nf90_copy_att(ncid, vid, 'calendar', ncid2, varid(1)), 504)

    CALL check( nf90_def_var( ncid2, 'VHM0', NF90_FLOAT, (/tdimid/),varid(2) ), 504)
    CALL check( nf90_put_att( ncid2, varid(2), 'long_name', &
        'Spectral significant wave height (Hm0)'), 504)
    CALL check( nf90_put_att( ncid2, varid(2), 'units', &
        'm'), 504)
    CALL check( nf90_put_att( ncid2, varid(2), 'coordinates', &
        trim(coord) ), 504)

    CALL check( nf90_enddef(ncid2), 516 )

    allocate( time(nt), data3d(ni,nj,nt), output(nt) )

    CALL check(nf90_inq_varid(ncid,'time',vid),123)
    CALL check(nf90_get_var(ncid,vid,time),18)

    CALL check(nf90_inq_varid(ncid,'VHM0',vid),124)
    CALL check(nf90_get_var(ncid,vid,data3d),18)
    CALL check(nf90_close(ncid),360)

    do t = 1,nt
      x = idx(1)
      y = idx(2)
      output(t) = data3d(x,y,t)*W(1,1)+&
                  data3d(x+1,y,t)*W(2,1)+&
                  data3d(x,y+1,t)*W(1,2)+&
                  data3d(x+1,y+1,t)*W(2,2)
    end do


    CALL check( nf90_put_var( ncid2, varid(1), time ), 517 )
    CALL check( nf90_put_var( ncid2, varid(2), output ), 517 )
    CALL check(nf90_close(ncid2),360)
   
    deallocate(time,data3d,output)
    deallocate(lon,lat) 
    
END PROGRAM getWAVE

SUBROUTINE check(status,label)
    USE NETCDF
    implicit none
    INTEGER, INTENT(in) :: status , label

    IF (status .ne. nf90_noerr) then
        print*, trim(nf90_strerror(status)),' label= ',label
        stop "Stopped"
    END IF
END SUBROUTINE check

SUBROUTINE calc_w(lon,lat,point_xy,w)
    implicit none
    real(kind=8), intent(in) :: lon(2),lat(2)
    real(kind=8), intent(in) :: point_xy(2)
    real(kind=8), intent(out) :: w(2,2)

    real(kind=8) :: x,x1,x2,y,y1,y2

    x = point_xy(1)
    y = point_xy(2)
    x1 = lon(1)
    x2 = lon(2)
    y1 = lat(1)
    y2 = lat(2)

    w(1,1) = (y2-y)*(x2-x)/(y2-y1)/(x2-x1)
    w(2,1) = (y2-y)*(x-x1)/(y2-y1)/(x2-x1)
    w(1,2) = (y-y1)*(x2-x)/(y2-y1)/(x2-x1)
    w(2,2) = (y-y1)*(x-x1)/(y2-y1)/(x2-x1)
END SUBROUTINE calc_w
