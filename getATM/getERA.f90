PROGRAM getERA
    USE netcdf
    implicit none
    integer :: cnt,idx(2),x,y
    integer :: ncid,dimid,ncid2,tdimid
    integer :: ni,nj,nt,vid,varid(7),i,n,p,t
    character(len=200) :: filename,path
    character(len=30) :: coord
    real(kind=8) :: W(2,2),pt(2),ao,sf
    integer, allocatable :: time(:)
    integer(kind=2), allocatable :: tmp(:,:,:),reverse(:,:,:)
    real(kind=8), allocatable :: lon(:),lat(:),tmp1d(:),RH(:)
    real(kind=8), allocatable :: data3d(:,:,:),uvel(:),vvel(:)
    real(kind=8), allocatable :: dew(:),tair(:),pair(:),qair(:)
    logical :: first,score,fexit

    100 format(a,'/ERA5_',i4.4,'.nc')
    101 format(f7.4,'N, ',f8.4,'E')
    path = '/users/work/mmuzyka/programs/ERA5'
    first = .true.
    cnt = 0
    
    pt = (/16.5,55.05/)
    write(coord,101) pt(2),pt(1) 

    CALL check(nf90_create('atm_ERA.nc',NF90_NETCDF4,ncid2),310)

    do i = 1990,2020
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

            CALL check(nf90_inq_varid(ncid,'time',vid),122)

            CALL check( nf90_def_dim( ncid2, 'time', NF90_UNLIMITED, tdimid ),503)
            CALL check( nf90_def_var( ncid2, 'time', NF90_INT, (/tdimid/), varid(1) ), 504)
            CALL check( nf90_copy_att(ncid, vid, 'units', ncid2, varid(1)), 504)
            CALL check( nf90_copy_att(ncid, vid, 'calendar', ncid2, varid(1)), 504)

            CALL check( nf90_def_var( ncid2, 'u10', NF90_FLOAT, (/tdimid/),varid(2) ), 504)
            CALL check( nf90_put_att( ncid2, varid(2), 'long_name', &
                '10 metre U wind component'), 504)
            CALL check( nf90_put_att( ncid2, varid(2), 'units', &
                'meter second-1'), 504)
            CALL check( nf90_put_att( ncid2, varid(2), 'coordinates', &
                trim(coord) ), 504)
            CALL check( nf90_def_var( ncid2, 'v10', NF90_FLOAT, (/tdimid/),varid(3) ), 505)
            CALL check( nf90_put_att( ncid2, varid(3), 'long_name', &
                '10 metre V wind component'), 504)
            CALL check( nf90_put_att( ncid2, varid(3), 'units', &
                'meter second-1'), 504)
            CALL check( nf90_put_att( ncid2, varid(3), 'coordinates', &
                trim(coord) ), 504)
            CALL check( nf90_def_var( ncid2, 'msl', NF90_FLOAT, (/tdimid/),varid(4) ), 505)
            CALL check( nf90_put_att( ncid2, varid(4), 'long_name', &
                'Mean sea level pressure'), 504)
            CALL check( nf90_put_att( ncid2, varid(4), 'units', &
                'hPa'), 504)
            CALL check( nf90_put_att( ncid2, varid(4), 'coordinates', &
                trim(coord) ), 504)
            CALL check( nf90_def_var( ncid2, 't2m', NF90_FLOAT, (/tdimid/),varid(5) ), 505)
            CALL check( nf90_put_att( ncid2, varid(5), 'long_name', &
                '2 metre temperature'), 504)
            CALL check( nf90_put_att( ncid2, varid(5), 'units', &
                'K'), 504)
            CALL check( nf90_put_att( ncid2, varid(5), 'coordinates', &
                trim(coord) ), 504)
            CALL check( nf90_def_var( ncid2, 'qair', NF90_FLOAT, (/tdimid/),varid(6) ), 505)
            CALL check( nf90_put_att( ncid2, varid(6), 'long_name', &
                'surface air specific humidity'), 504)
            CALL check( nf90_put_att( ncid2, varid(6), 'units', &
                'g/kg'), 504)
            CALL check( nf90_put_att( ncid2, varid(6), 'coordinates', &
                trim(coord) ), 504)
            CALL check( nf90_def_var( ncid2, 'RH', NF90_FLOAT, (/tdimid/),varid(7) ), 505)
            CALL check( nf90_put_att( ncid2, varid(7), 'long_name', &
                'surface air relative humidity'), 504)
            CALL check( nf90_put_att( ncid2, varid(7), 'units', &
                '%'), 504)
            CALL check( nf90_put_att( ncid2, varid(7), 'coordinates', &
                trim(coord) ), 504)

            CALL check( nf90_enddef(ncid2), 516 )
            first = .false.
        end if

        allocate( time(nt), data3d(ni,nj,nt), uvel(nt), vvel(nt) )
        allocate( tmp(ni,nj,nt), reverse(ni,nj,nt), RH(nt) )
        allocate( dew(nt), tair(nt), pair(nt), qair(nt) ) 

        CALL check(nf90_inq_varid(ncid,'time',vid),123)
        CALL check(nf90_get_var(ncid,vid,time),18)

        CALL check(nf90_inq_varid(ncid,'u10',vid),124)
        CALL check(nf90_get_var(ncid,vid,tmp),18)
        CALL check(nf90_get_att(ncid,vid,'scale_factor',sf),18)
        CALL check(nf90_get_att(ncid,vid,'add_offset',ao),18)

        do p = 1,nj
          reverse(:,p,:) = tmp(:,nj+1-p,:)
        end do
        data3d = real(reverse,8)*sf+ao
    
        do t = 1,nt
          x = idx(1)
          y = idx(2)
          uvel(t) = data3d(x,y,t)*W(1,1)+&
                      data3d(x+1,y,t)*W(2,1)+&
                      data3d(x,y+1,t)*W(1,2)+&
                      data3d(x+1,y+1,t)*W(2,2)
        end do

        CALL check(nf90_inq_varid(ncid,'v10',vid),124)
        CALL check(nf90_get_var(ncid,vid,tmp),18)
        CALL check(nf90_get_att(ncid,vid,'scale_factor',sf),18)
        CALL check(nf90_get_att(ncid,vid,'add_offset',ao),18)

        do p = 1,nj
          reverse(:,p,:) = tmp(:,nj+1-p,:)
        end do
        data3d = real(reverse,8)*sf+ao

        do t = 1,nt
          x = idx(1)
          y = idx(2)
          vvel(t) = data3d(x,y,t)*W(1,1)+&
                      data3d(x+1,y,t)*W(2,1)+&
                      data3d(x,y+1,t)*W(1,2)+&
                      data3d(x+1,y+1,t)*W(2,2)
        end do

        CALL check(nf90_inq_varid(ncid,'t2m',vid),124)
        CALL check(nf90_get_var(ncid,vid,tmp),18)
        CALL check(nf90_get_att(ncid,vid,'scale_factor',sf),18)
        CALL check(nf90_get_att(ncid,vid,'add_offset',ao),18)

        do p = 1,nj
          reverse(:,p,:) = tmp(:,nj+1-p,:)
        end do
        data3d = real(reverse,8)*sf+ao

        do t = 1,nt
          x = idx(1)
          y = idx(2)
          tair(t) = data3d(x,y,t)*W(1,1)+&
                      data3d(x+1,y,t)*W(2,1)+&
                      data3d(x,y+1,t)*W(1,2)+&
                      data3d(x+1,y+1,t)*W(2,2)
        end do

        CALL check(nf90_inq_varid(ncid,'msl',vid),124)
        CALL check(nf90_get_var(ncid,vid,tmp),18)
        CALL check(nf90_get_att(ncid,vid,'scale_factor',sf),18)
        CALL check(nf90_get_att(ncid,vid,'add_offset',ao),18)

        do p = 1,nj
          reverse(:,p,:) = tmp(:,nj+1-p,:)
        end do
        data3d = real(reverse,8)*sf+ao

        do t = 1,nt
          x = idx(1)
          y = idx(2)
          pair(t) = data3d(x,y,t)*W(1,1)+&
                      data3d(x+1,y,t)*W(2,1)+&
                      data3d(x,y+1,t)*W(1,2)+&
                      data3d(x+1,y+1,t)*W(2,2)
        end do

        CALL check(nf90_inq_varid(ncid,'d2m',vid),124)
        CALL check(nf90_get_var(ncid,vid,tmp),18)
        CALL check(nf90_get_att(ncid,vid,'scale_factor',sf),18)
        CALL check(nf90_get_att(ncid,vid,'add_offset',ao),18)

        do p = 1,nj
          reverse(:,p,:) = tmp(:,nj+1-p,:)
        end do
        data3d = real(reverse,8)*sf+ao

        do t = 1,nt
          x = idx(1)
          y = idx(2)
          dew(t) = data3d(x,y,t)*W(1,1)+&
                      data3d(x+1,y,t)*W(2,1)+&
                      data3d(x,y+1,t)*W(1,2)+&
                      data3d(x+1,y+1,t)*W(2,2)
        end do
        CALL spec_humid(nt,dew,pair,qair)
        CALL rel_humid(nt,dew,tair,RH)

        CALL check(nf90_close(ncid),360)

        CALL check( nf90_put_var( ncid2, varid(1), time, start=(/cnt+1/) ), 517 )
        CALL check( nf90_put_var( ncid2, varid(2), uvel,&
            start=(/cnt+1/) ), 517 )
        CALL check( nf90_put_var( ncid2, varid(3), vvel,&
            start=(/cnt+1/) ), 517 )
        CALL check( nf90_put_var( ncid2, varid(4), pair/100,&
            start=(/cnt+1/) ), 517 )
        CALL check( nf90_put_var( ncid2, varid(5), tair,&
            start=(/cnt+1/) ), 517 )
        CALL check( nf90_put_var( ncid2, varid(6), qair*1000,&
            start=(/cnt+1/) ), 517 )
        CALL check( nf90_put_var( ncid2, varid(7), RH,&
            start=(/cnt+1/) ), 517 )
        
        cnt = cnt+nt
    
        deallocate(time,tmp,reverse,data3d,uvel,vvel)
        deallocate(dew,tair,pair,qair,RH)

    end do

    CALL check(nf90_close(ncid2),360)
   
    deallocate(lon,lat,tmp1d) 
    
END PROGRAM getERA

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

SUBROUTINE spec_humid(nt,dew,pair,output)
    implicit none
    integer, intent(in) :: nt
    real(kind=8), intent(in) :: dew(nt),pair(nt)
    real(kind=8), intent(out) :: output(nt)

    integer :: t
    real(kind=8) :: Rdry,Rvap,a1,a3,a4,T0
    real(kind=8) :: esat(nt)

    Rdry = 287.0597
    Rvap = 461.5250
    a1 = 611.21
    a3 = 17.502
    a4 = 32.19
    T0 = 273.16

    do t = 1, nt
      esat(t) = a1*exp(a3*((dew(t)-T0)/(dew(t)-a4)))
    end do

    do t = 1, nt
      output(t) = (Rdry*esat(t)/Rvap)/(pair(t)-(1-Rdry/Rvap)*esat(t))
    end do
END SUBROUTINE spec_humid

SUBROUTINE rel_humid(nt,dew,tair,output)
    implicit none
    integer, intent(in) :: nt
    real(kind=8), intent(in) :: dew(nt),tair(nt)
    real(kind=8), intent(out) :: output(nt)

    real(kind=8) :: dewC(nt),tairC(nt),E(nt),Es(nt)

    dewC = dew-273.15
    tairC = tair-273.15

    E = exp( (17.625*dewC)/(243.04+dewC) ) 
    Es = exp( (17.625*tairC)/(243.04+tairC) )

    output = 100*(E/Es) 

END SUBROUTINE rel_humid
