PROGRAM windpts
    USE netcdf
    implicit none
    integer :: cnt,idx(2,2),x,y
    integer :: ncid,dimid,ncid2,tdimid
    integer :: ni,nj,nt,vid,varid(5),i,n,p,t
    character(len=200) :: filename,path
    real(kind=8) :: W(2,2,2),pt(2,2),ao,sf
    integer, allocatable :: time(:)
    integer(kind=2), allocatable :: tmp(:,:,:)
    real(kind=8), allocatable :: lon(:),lat(:),tmp1d(:)
    real(kind=8), allocatable :: data3d(:,:,:),uvel(:,:),vvel(:,:)
    logical :: first,score,fexit

    100 format(a,'/ERA5_',i4.4,'.nc')
    path = '/users/work/mmuzyka/programs/ERA5'
    first = .true.
    cnt = 0
    
    pt(1,:) = (/18.81,54.6/)
    pt(2,:) = (/15.1,55.05/)

    CALL check(nf90_create('windpts_ERA.nc',NF90_NETCDF4,ncid2),310)

    do i = 1990,2020!2019
        write(filename,100) trim(path),i

        write(*,*) trim(filename)

        CALL check(nf90_open(trim(filename),NF90_NOWRITE,ncid),310)

        CALL check(nf90_inq_dimid(ncid, "time", dimid),315)
        CALL check(nf90_inquire_dimension(ncid, dimid, len=nt),316)

        if (first) then
            CALL check(nf90_inq_dimid(ncid, "longitude", dimid),311)
            CALL check(nf90_inquire_dimension(ncid, dimid, len=ni),312)

            CALL check(nf90_inq_dimid(ncid, "latitude", dimid),313)
            CALL check(nf90_inquire_dimension(ncid, dimid, len=nj),314)

            allocate( lon(ni), lat(nj), tmp1d(nj) )

            CALL check(nf90_inq_varid(ncid,'longitude',vid),119)
            CALL check(nf90_get_var(ncid,vid,lon),18)
            CALL check(nf90_inq_varid(ncid,'latitude',vid),121)
            CALL check(nf90_get_var(ncid,vid,tmp1d),18)

            do p = 1,nj
              lat(p) = tmp1d(nj+1-p)
            end do
    
            do n = 1,2
              do p = 1,ni
                if (pt(n,1).lt.lon(p)) then
                  idx(n,1) = p-1
                  exit
                end if
              end do
              do p = 1,nj
                if (pt(n,2).lt.lat(p)) then
                  idx(n,2) = p-1
                  exit
                end if
              end do
              x = idx(n,1)
              y = idx(n,2)
              write(*,*) pt(n,:)
              write(*,*) x,y
              write(*,*) lon(x:x+1) 
              write(*,*) lat(y:y+1)
              call calc_w( lon(x:x+1), lat(y:y+1), pt(n,:), W(n,:,:) )
              write(*,*) sum(W(n,:,:))
            end do

            CALL check(nf90_inq_varid(ncid,'time',vid),122)

            CALL check( nf90_def_dim( ncid2, 'time', NF90_UNLIMITED, tdimid ),503)
            CALL check( nf90_def_var( ncid2, 'time', NF90_INT, (/tdimid/), varid(1) ), 504)
            CALL check( nf90_copy_att(ncid, vid, 'units', ncid2, varid(1)), 504)
            CALL check( nf90_copy_att(ncid, vid, 'calendar', ncid2, varid(1)), 504)

            CALL check( nf90_def_var( ncid2, 'upt1', NF90_FLOAT, (/tdimid/),varid(2) ), 504)
            CALL check( nf90_put_att( ncid2, varid(2), 'long_name', &
                '10-meter u-wind component for pt1(54.6 N, 18.81 E)'), 504)
            CALL check( nf90_put_att( ncid2, varid(2), 'units', &
                'meter second-1'), 504)
            CALL check( nf90_def_var( ncid2, 'upt2', NF90_FLOAT, (/tdimid/),varid(3) ), 505)
            CALL check( nf90_put_att( ncid2, varid(3), 'long_name', &
                '10-meter u-wind component for pt2(55.05 N, 15.1 E)'), 504)
            CALL check( nf90_put_att( ncid2, varid(3), 'units', &
                'meter second-1'), 504)
            
            CALL check( nf90_def_var( ncid2, 'vpt1', NF90_FLOAT, (/tdimid/),varid(4) ), 504)
            CALL check( nf90_put_att( ncid2, varid(4), 'long_name', &
                '10-meter v-wind component for pt1(54.6 N, 18.81 E)'), 504)
            CALL check( nf90_put_att( ncid2, varid(4), 'units', &
                'meter second-1'), 504)
            CALL check( nf90_def_var( ncid2, 'vpt2', NF90_FLOAT, (/tdimid/),varid(5) ), 505)
            CALL check( nf90_put_att( ncid2, varid(5), 'long_name', &
                '10-meter v-wind component for pt2(55.05 N, 15.1 E)'), 504)
            CALL check( nf90_put_att( ncid2, varid(5), 'units', &
                'meter second-1'), 504)

            CALL check( nf90_enddef(ncid2), 516 )
            first = .false.
        end if

        allocate( time(nt), tmp(ni,nj,nt), data3d(ni,nj,nt), uvel(2,nt), vvel(2,nt) )

        CALL check(nf90_inq_varid(ncid,'time',vid),123)
        CALL check(nf90_get_var(ncid,vid,time),18)

        CALL check(nf90_inq_varid(ncid,'u10',vid),124)
        CALL check(nf90_get_var(ncid,vid,tmp(:,:,:)),18)
        CALL check(nf90_get_att(ncid,vid,'scale_factor',sf),18)
        CALL check(nf90_get_att(ncid,vid,'add_offset',ao),18)

        data3d = real(tmp,8)*sf+ao
    
        do n = 1,2
          do t = 1,nt
            x = idx(n,1)
            y = idx(n,2)
            uvel(n,t) = data3d(x,y,t)*W(n,1,1)+&
                        data3d(x+1,y,t)*W(n,2,1)+&
                        data3d(x,y+1,t)*W(n,1,2)+&
                        data3d(x+1,y+1,t)*W(n,2,2)
          end do
        end do

        CALL check(nf90_inq_varid(ncid,'v10',vid),124)
        CALL check(nf90_get_var(ncid,vid,tmp(:,:,:)),18)
        CALL check(nf90_get_att(ncid,vid,'scale_factor',sf),18)
        CALL check(nf90_get_att(ncid,vid,'add_offset',ao),18)

        data3d = real(tmp,8)*sf+ao

        do n = 1,2
          do t = 1,nt
            x = idx(n,1)
            y = idx(n,2)
            vvel(n,t) = data3d(x,y,t)*W(n,1,1)+&
                        data3d(x+1,y,t)*W(n,2,1)+&
                        data3d(x,y+1,t)*W(n,1,2)+&
                        data3d(x+1,y+1,t)*W(n,2,2)
          end do
        end do

        CALL check(nf90_close(ncid),360)

        CALL check( nf90_put_var( ncid2, varid(1), time, start=(/cnt+1/) ), 517 )
        CALL check( nf90_put_var( ncid2, varid(2), uvel(1,:),&
            start=(/cnt+1/) ), 517 )
        CALL check( nf90_put_var( ncid2, varid(3), uvel(2,:),&
            start=(/cnt+1/) ), 517 )
        CALL check( nf90_put_var( ncid2, varid(4), vvel(1,:),&
            start=(/cnt+1/) ), 517 )
        CALL check( nf90_put_var( ncid2, varid(5), vvel(2,:),&
            start=(/cnt+1/) ), 517 )
        
        cnt = cnt+nt
    
        deallocate(time,tmp,data3d,uvel,vvel)

    end do

    CALL check(nf90_close(ncid2),360)
   
    deallocate(lon,lat,tmp1d) 
    
END PROGRAM windpts

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

SUBROUTINE in_convex_polygon(corners,point,score)
    implicit none
    real(kind=8), intent(in) :: corners(4,2)
    real(kind=8), intent(in) :: point(2)
    logical, intent(out) :: score
    logical :: state,side
    real(kind=8) :: xi,yi,xj,yj,d
    integer*1 :: i

    score=.TRUE.
    state=.TRUE.

    do i = 1,4
      xi = corners(i,1)
      yi = corners(i,2)
      xj = corners(modulo(i,4)+1,1)
      yj = corners(modulo(i,4)+1,2)
      d = (point(1)-xi)*(yj-yi)-(point(2)-yi)*(xj-xi)
      if ( d==0.0 ) then
        go to 10
      else
        if (state) then
          state=.FALSE.
          side=(d>0.0)
        else if ((d>0.0)/=side) then
          score=.FALSE.
          exit
        end if
      end if
      10 continue
    end do
END SUBROUTINE in_convex_polygon

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

SUBROUTINE wind(nt,udir,vspe)
    implicit none
    integer, intent(in) :: nt
    real(kind=8), intent(inout) :: udir(nt),vspe(nt)

    real(kind=8) :: wdir(nt),wspeed(nt),pi
    integer :: i

    pi = 4.0*ATAN(1.0)

    wdir = (270.0-udir(:))*pi/180.0
    wspeed = vspe(:)

    udir(:) = wspeed(:)*cos(wdir(:))
    vspe(:) = wspeed(:)*sin(wdir(:))
END SUBROUTINE wind
