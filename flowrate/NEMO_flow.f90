PROGRAM NEMO_flow
    USE NETCDF
    implicit none
    integer, parameter :: xin=763,yin=774,zin=56
    integer :: ncid, varid, i, j, k, nc, tdimid, vid(5), cnt, tmp(1)
    integer :: date(3), datestp(3), jstr(2), jend(2), istr(2), iend(2)
    real(kind=8) :: lon_tr(xin), lat_tr(yin), depth_tr(zin), hz_tr(zin)
    real(kind=8) :: hlat(2,2),hlon(2),fv,d,dx(xin),flux(4),time,dy(2,yin) 
    real(kind=8), allocatable :: uvelin(:,:,:), vvelin(:,:,:)
    character(len=150) :: trfile
    character(len=200) :: path
    logical :: first
  
    path = '/users/work/mmuzyka/programs/ROMS_bc/new_files/' 

    hlon(1) = 14.75 
    hlon(2) = 14.75 

    hlat(1,:) = (/55.2,56.2/)
    hlat(2,:) = (/54.0,55.1/)


    allocate( uvelin(xin,yin,zin), vvelin(xin,yin,zin) )


    CALL check( nf90_create( 'NEMO_flow.nc',NF90_NETCDF4,nc ), 500 )

    CALL check( nf90_def_dim( nc, 'time', NF90_UNLIMITED, tdimid ), 501 )

    CALL check(nf90_def_var( nc, 'time', NF90_DOUBLE,&
        (/tdimid/), vid(1) ), 508)
    CALL check( nf90_put_att( nc, vid(1), 'units',&
        'days since 1900-01-01 00:00:00' ), 305 )
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

    CALL check( nf90_enddef(nc), 516 )


    date = (/2003,11,1/)
    datestp = (/2004,12,30/)
!READ TRACERS FROM FILES
    100 format(a,i4.4,'/',i2.2,'/BAL-MYP-NEMO_PHY-DailyMeans-',i4,i2.2,i2.2,'.nc')
    first = .true.
    cnt = 1
    do
      write(trfile,100) trim(path),date(1),date(2),date(1),date(2),date(3)
      write(*,*) trim(trfile) 
      CALL check(nf90_open(trim(trfile),NF90_NOWRITE,ncid),1016)
      
      if (first) then  
        CALL check(nf90_inq_varid(ncid,"depth",varid),17)
        CALL check(nf90_get_var(ncid,varid,depth_tr),18)

        hz_tr = 0.0
        do k = 1,zin    
          hz_tr(k) = (depth_tr(k)-sum(hz_tr))*2
        end do

        CALL check(nf90_inq_varid(ncid,"lon",varid),19)
        CALL check(nf90_get_var(ncid,varid,lon_tr),18)


        CALL check(nf90_inq_varid(ncid,"lat",varid),20)
        CALL check(nf90_get_var(ncid,varid,lat_tr),18)

        tmp = minloc(abs(lon_tr-hlon(1)))
        istr(1) = tmp(1)
        tmp = minloc(abs(lon_tr-hlon(2)))
        istr(2) = tmp(1)
        tmp = minloc(abs(lat_tr-hlat(1,1)))
        jstr(1) = tmp(1)
        tmp = minloc(abs(lat_tr-hlat(2,1)))
        jstr(2) = tmp(1)
        tmp = minloc(abs(lat_tr-hlat(1,2)))
        jend(1) = tmp(1)
        tmp = minloc(abs(lat_tr-hlat(2,2)))
        jend(2) = tmp(1)


      end if

      CALL check(nf90_inq_varid(ncid,"time",varid),17)
      CALL check(nf90_get_var(ncid,varid,time),18)

!      CALL check(nf90_inq_varid(ncid,"thetao",vid),17)
!      CALL check(nf90_get_att(ncid,vid,"_FillValue",fv),471)
!      CALL check(nf90_get_var(ncid,vid,tempin),18)

!      CALL check(nf90_inq_varid(ncid,"so",vid),17)
!      CALL check(nf90_get_var(ncid,vid,saltin),18)

!      CALL check(nf90_inq_varid(ncid,"uo",vid),17)
!      CALL check(nf90_get_var(ncid,vid,uvelin),18)

!      CALL check(nf90_inq_varid(ncid,"vo",varid),17)
!      CALL check(nf90_get_att(ncid,varid,"_FillValue",fv),471)
!      CALL check(nf90_get_var(ncid,varid,vvelin),18)

      CALL check(nf90_inq_varid(ncid,"uo",varid),17)
      CALL check(nf90_get_att(ncid,varid,"_FillValue",fv),471)
      CALL check(nf90_get_var(ncid,varid,uvelin),18)

      where(uvelin.eq.fv) uvelin=0.0

      if (first) then
        do j = jstr(1),jend(1)
          CALL haversine((/lon_tr(istr(1)),lon_tr(istr(1))/),&
              (/lat_tr(j-1),lat_tr(j)/),d)
          dy(1,j) = d*0.5
          CALL haversine((/lon_tr(istr(1)),lon_tr(istr(1))/),&
              (/lat_tr(j),lat_tr(j+1)/),d)
          dy(1,j) = dy(1,j)+d*0.5
          write(*,*) 'j=',j,'dy=',dy(1,j),'m'
        end do
        do j = jstr(2),jend(2)
          CALL haversine((/lon_tr(istr(2)),lon_tr(istr(2))/),&
              (/lat_tr(j-1),lat_tr(j)/),d)
          dy(2,j) = d*0.5
          CALL haversine((/lon_tr(istr(2)),lon_tr(istr(2))/),&
              (/lat_tr(j),lat_tr(j+1)/),d)
          dy(2,j) = dy(2,j)+d*0.5
          write(*,*) 'j=',j,'dy=',dy(2,j),'m'
        end do

        first = .false.
      end if

      flux = 0.0
      do j = jstr(1),jend(1)
      do k = 1, zin
        flux(1) = flux(1)+uvelin(istr(1),j,k)*dy(1,j)*hz_tr(k)
        flux(3) = flux(3)+abs(uvelin(istr(1),j,k))*dy(1,j)*hz_tr(k)
      end do
      end do
      do j = jstr(2),jend(2)
      do k = 1, zin
        flux(2) = flux(2)+uvelin(istr(2),j,k)*dy(2,j)*hz_tr(k)
        flux(4) = flux(4)+abs(uvelin(istr(2),j,k))*dy(2,j)*hz_tr(k)
      end do
      end do

      CALL check( nf90_put_var( nc, vid(1), time, start=(/cnt/)), 517 )
      CALL check( nf90_put_var( nc, vid(2), flux(1)*1e-9, start=(/cnt/)), 517 )
      CALL check( nf90_put_var( nc, vid(3), flux(2)*1e-9, start=(/cnt/)), 517 )
      CALL check( nf90_put_var( nc, vid(4), flux(3)*1e-9, start=(/cnt/)), 517 )
      CALL check( nf90_put_var( nc, vid(5), flux(4)*1e-9, start=(/cnt/)), 517 )
      cnt = cnt+1

      if (date(1).ge.datestp(1).and.&
            date(2).ge.datestp(2).and.&
            date(3).ge.datestp(3)) EXIT
      CALL add_day(date,.true.)
    end do

    CALL check( nf90_close(nc), 600 )

    deallocate( uvelin, vvelin ) 
!END READ TRACERS

END PROGRAM NEMO_flow

SUBROUTINE check(status,label)
    USE NETCDF
    implicit none
    INTEGER, INTENT(in) :: status , label

    IF (status .ne. nf90_noerr) then
        print*, trim(nf90_strerror(status)),' label= ',label
        stop "Stopped"
    END IF
END SUBROUTINE check

SUBROUTINE haversine(lon,lat,d)
    implicit none
    real(kind=8), intent(in) :: lon(2),lat(2)
    real(kind=8), intent(out) :: d

    real(kind=8) :: R,a,c,pi,fi(2),lambda(2)

    pi = 4.D0*DATAN(1.D0)

    lambda = lon*pi/180.0
    fi = lat*pi/180.0

    R = 6.371e6

    a = sin(abs(fi(2)-fi(1))*0.5)**2+cos(fi(1))*cos(fi(2))*&
        sin(abs(lambda(2)-lambda(1))*0.5)**2

    c = 2*atan2(sqrt(a),sqrt(1-a))

    d = R*c

END SUBROUTINE haversine

SUBROUTINE leap(year,answer)
    implicit none
    integer, intent(in) :: year
    logical, intent(out) :: answer

    answer = .false.
    if( mod(year,4).eq.0 ) answer = .true.
    if( mod(year,4).eq.100 ) answer = .false.
    if( mod(year,4).eq.400 ) answer = .true.
END SUBROUTINE leap

SUBROUTINE extrap(a,mask,lon,lat,maxscn,met)
    implicit none
    integer, intent(in) :: lon,lat,maxscn,met
    real(kind=8), intent(inout) :: a(lon,lat)
    logical, intent(in) :: mask(lon,lat)

    integer :: i,j,n,cnt,overall
    real(kind=8) :: relc,ave
    real(kind=8), dimension(lon,lat) :: sor,res
    logical :: mask_tmp(lon,lat),mask_tmp2(lon,lat)

    relc=1.0
    sor = 0.0
    where(.not.mask) sor=relc

    select case(met)
    case(0)
      where(.not.mask) a=0.0
    case(1)
      cnt = 0
      ave = 0.0
      do i=1,lon
      do j=1,lat
        if (mask(i,j)) then
            ave=ave+a(i,j)
            cnt=cnt+1
        end if
      end do
      end do
      if ( cnt.GT.0 ) ave = ave/real(cnt,8)
      where(.not.mask) a=ave
    case(2)
      mask_tmp2 = mask
      mask_tmp = mask_tmp2
      do

      overall = 0
      do i = 1, lon
      do j = 1, lat

        if (.not.mask_tmp(i,j)) then

          cnt = 0
          ave = 0.0

          if ( i.gt.1 ) then
            if ( mask_tmp(i-1,j) ) then
              ave = ave+a(i-1,j)
              cnt = cnt+1
            end if
          end if

          if ( j.gt.1 ) then
            if ( mask_tmp(i,j-1) ) then
              ave = ave+a(i,j-1)
              cnt = cnt+1
            end if
          end if

          if ( i.lt.lon ) then
            if ( mask_tmp(i+1,j) ) then
              ave = ave+a(i+1,j)
              cnt = cnt+1
            end if
          end if

          if ( j.lt.lat ) then
            if ( mask_tmp(i,j+1) ) then
              ave = ave+a(i,j+1)
              cnt = cnt+1
            end if
          end if

          if ( cnt.gt.0 ) then
            a(i,j) = ave/(real(cnt,8))
            overall = overall+cnt
            mask_tmp2(i,j) = .true.
          end if

        end if

      end do
      end do

      mask_tmp = mask_tmp2
      if ( overall.eq.0 ) EXIT

      end do
    end select

    do n=1,maxscn

!!$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j),SCHEDULE(DYNAMIC)
        do i=2,lon-1
          do j=2,lat-1
            res(i,j)=0.25*(a(i-1,j)+a(i+1,j)+a(i,j-1)+a(i,j+1))-a(i,j)
          end do
        end do
!!$OMP END PARALLEL DO

        do i=2,lon-1
          res(i,1)=0.3333*(a(i-1,1)+a(i+1,1)+a(i,2))-a(i,1)
          res(i,lat)=0.3333*(a(i-1,lat)+a(i+1,lat)+a(i,lat-1))-a(i,lat)
        end do

        do j=2,lat-1
          res(1,j)=0.3333*(a(1,j-1)+a(1,j+1)+a(2,j))-a(1,j)
          res(lon,j)=0.3333*(a(lon,j-1)+a(lon,j+1)+a(lon-1,j))-a(lon,j)
        end do

        res(1,1)=0.5*(a(1,2)+a(2,1))-a(1,1)
        res(lon,1)=0.5*(a(lon,2)+a(lon-1,1))-a(lon,1)
        res(1,lat)=0.5*(a(1,lat-1)+a(2,lat))-a(1,lat)
        res(lon,lat)=0.5*(a(lon,lat-1)+a(lon-1,lat))-a(lon,lat)

        res=res*sor
        a=a+res
    end do
END SUBROUTINE extrap

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

SUBROUTINE calc_w(lonin,latin,point_xy,w)
    implicit none
    real(kind=8), intent(in) :: lonin(2),latin(2)
    real(kind=8), intent(in) :: point_xy(2)
    real(kind=8), intent(out) :: w(4)

    real(kind=8) :: x, y, lon(2,2), lat(2,2)
    real(kind=8) :: A, B, C, s, t
    real(kind=8) :: y1,y2,x1,x2

    x = point_xy(1)
    y = point_xy(2)
    lon(1,:) = lonin
    lon(2,:) = lonin
    lat(:,1) = latin
    lat(:,2) = latin

    A = ( lon(1,1)-lon(1,2) )*( lat(2,1)-lat(2,2) )-&
        ( lat(1,1)-lat(1,2) )*( lon(2,1)-lon(2,2) )

    B = y*( ( lon(2,1)-lon(2,2) )-( lon(1,1)-lon(1,2) ) )-&
        x*( ( lat(2,1)-lat(2,2) )-( lat(1,1)-lat(1,2) ) )+&
        ( lon(1,1)-lon(1,2) )*lat(2,2)-&
        ( lat(1,1)-lat(1,2) )*lon(2,2)+&
        ( lat(2,1)-lat(2,2) )*lon(1,2)-&
        ( lon(2,1)-lon(2,2) )*lat(1,2)

    C = y*( lon(2,2)-lon(1,2) )-&
        x*( lat(2,2)-lat(1,2) )+&
        lon(1,2)*lat(2,2)-&
        lat(1,2)*lon(2,2)

    if ( A.ne.0 ) then 
      t = ( -B+sqrt(B**2-4*A*C) )/(2*A)
      s = (  y-lat(1,2)-( lat(1,1)-lat(1,2) )*t  )/&
        (  lat(2,2)+( lat(2,1)-lat(2,2) )*t-&
        lat(1,2)-( lat(1,1)-lat(1,2) )*t  )

      if ( (t<0).OR.(t>1) ) then
        t = ( -B-sqrt(B**2-4*A*C) )/(2*A)
        s = (  y-lat(1,2)-( lat(1,1)-lat(1,2) )*t  )/&
          (  lat(2,2)+( lat(2,1)-lat(2,2) )*t-&
          lat(1,2)-( lat(1,1)-lat(1,2) )*t  )
      end if

      w(1) = (1-s)*t
      w(2) = (1-s)*(1-t)
      w(3) = s*(1-t)
      w(4) = s*t
    else
      x1 = lonin(1)
      x2 = lonin(2)
      y1 = latin(1)
      y2 = latin(2)
      w(1) = ((x2-x)/(x2-x1))*((y2-y)/(y2-y1)) 
      w(2) = ((x2-x)/(x2-x1))*((y-y1)/(y2-y1)) 
      w(3) = ((x-x1)/(x2-x1))*((y-y1)/(y2-y1)) 
      w(4) = ((x-x1)/(x2-x1))*((y2-y)/(y2-y1)) 
    end if

END SUBROUTINE calc_w

SUBROUTINE calcme(nxin,nyin,nxout,inarray,outarray,W,idx)
    implicit none
    integer, intent(in) :: nxin,nyin,nxout
    real(kind=8), intent(in) :: inarray(nxin,nxin),W(nxout,4)
    real(kind=8), intent(out) :: outarray(nxout)
    integer, intent(in) :: idx(nxout,4,2)

    real(kind=8) :: tmp
    integer :: i,k,x,y

    do i = 1,nxout
      if (idx(i,1,1).ne.0) then
        tmp = 0
        do k = 1,4
          x = idx(i,k,1)
          y = idx(i,k,2)
          if (x.gt.0 .and. y.gt.0) then
            tmp = tmp+inarray(x,y)*W(i,k)
          end if
        end do
        outarray(i) = tmp
      end if
    end do
END SUBROUTINE calcme

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

SUBROUTINE elapsed(date,date_0,days)
    implicit none
    integer, intent(in) :: date(3),date_0(3)
    integer, intent(out) :: days

    integer :: date_tmp(3)

    date_tmp(1) = date_0(1)
    date_tmp(2) = date_0(2)
    date_tmp(3) = date_0(3)

    days = 0

    do
        if( date_tmp(1).eq.date(1) .and.&
            date_tmp(2).eq.date(2) .and.&
            date_tmp(3).eq.date(3) ) EXIT
        days = days+1
        CALL add_day(date_tmp,.true.)
    end do
END SUBROUTINE elapsed

SUBROUTINE extrap3d(a,mask,lon,lat,n,maxscn,met)
    USE omp_lib
    implicit none
    integer, intent(in) :: lon,lat,n,maxscn,met
    real(kind=8), intent(inout) :: a(lon,lat,n)
    logical, intent(in) :: mask(lon,lat,n)

    integer :: i,j,k,l,cnt,overall
    real(kind=8) :: relc,ave
    real(kind=8), dimension(lon,lat,n) :: sor,res
    logical :: mask_tmp(lon,lat,n),mask_tmp2(lon,lat,n)

    relc = 0.6
    sor = 0.0
    where(.not.mask) sor=relc

    !FILLING LAND WITH LAYER AVERAGE
    select case(met)
    case(0)
      where(.not.mask) a=0.0
    case(1)
      do k = 1,n

      cnt = 0
      ave = 0.0

      do i=1,lon
      do j=1,lat
        if (mask(i,j,k)) then
            ave=ave+a(i,j,k)
            cnt=cnt+1
        end if
      end do
      end do
      if ( cnt.GT.0 ) ave = ave/real(cnt,8)
      where(.not.mask(:,:,k)) a(:,:,k)=ave

      end do
    case(2)
      mask_tmp2 = mask
      mask_tmp = mask_tmp2
      do

      overall = 0
      do i = 1, lon
      do j = 1, lat
      do k = 1, n

        if (.not.mask_tmp(i,j,k)) then

          cnt = 0
          ave = 0.0

          if ( i.gt.1 ) then
            if ( mask_tmp(i-1,j,k) ) then
              ave = ave+a(i-1,j,k)
              cnt = cnt+1
            end if
          end if

          if ( j.gt.1 ) then
            if ( mask_tmp(i,j-1,k) ) then
              ave = ave+a(i,j-1,k)
              cnt = cnt+1
            end if
          end if

          if ( k.gt.1 ) then
            if ( mask_tmp(i,j,k-1) ) then
              ave = ave+a(i,j,k-1)
              cnt = cnt+1
            end if
          end if

          if ( i.lt.lon ) then
            if ( mask_tmp(i+1,j,k) ) then
              ave = ave+a(i+1,j,k)
              cnt = cnt+1
            end if
          end if

          if ( j.lt.lat ) then
            if ( mask_tmp(i,j+1,k) ) then
              ave = ave+a(i,j+1,k)
              cnt = cnt+1
            end if
          end if

          if ( k.lt.n ) then
            if ( mask_tmp(i,j,k+1) ) then
              ave = ave+a(i,j,k+1)
              cnt = cnt+1
            end if
          end if

          if ( cnt.gt.0 ) then
            a(i,j,k) = ave/(real(cnt,8))
            overall = overall+cnt
            mask_tmp2(i,j,k) = .true.
          end if

        end if

      end do
      end do
      end do
      mask_tmp = mask_tmp2
      if ( overall.eq.0 ) EXIT

      end do
    end select
    !END FILLING

    do l = 1,maxscn

!$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k),SCHEDULE(DYNAMIC)
      do i=2,lon-1
        do j=2,lat-1
          res(i,j,n)=0.2*(a(i-1,j,n)+a(i+1,j,n)+a(i,j-1,n)+a(i,j+1,n)+&
            a(i,j,n-1))-a(i,j,n)
          res(i,j,1)=0.2*(a(i-1,j,1)+a(i+1,j,1)+a(i,j-1,1)+a(i,j+1,1)+&
            a(i,j,2))-a(i,j,1)
        do k = 2,n-1
          res(i,j,k)=0.1667*(a(i-1,j,k)+a(i+1,j,k)+a(i,j-1,k)+a(i,j+1,k)+&
            a(i,j,k+1)+a(i,j,k-1))-a(i,j,k)
        end do
        end do
      end do
!$OMP END PARALLEL DO

      do k = 2,n-1
      do i=2,lon-1
        res(i,1,k)=0.2*(a(i-1,1,k)+a(i+1,1,k)+a(i,1,k+1)+a(i,1,k-1)+a(i,2,k))-a(i,1,k)
        res(i,lat,k)=0.2*(a(i-1,lat,k)+a(i+1,lat,k)+a(i,lat,k+1)+a(i,lat,k-1)+a(i,lat-1,k))-a(i,lat,k)
      end do
      do j=2,lat-1
        res(1,j,k)=0.2*(a(1,j-1,k)+a(1,j+1,k)+a(1,j,k+1)+a(1,j,k-1)+a(2,j,k))-a(1,j,k)
        res(lon,j,k)=0.2*(a(lon,j-1,k)+a(lon,j+1,k)+a(lon,j,k+1)+a(lon,j,k-1)+a(lon-1,j,k))-a(lon,j,k)
      end do
      end do

      do i=2,lon-1
        res(i,1,1)=0.25*(a(i-1,1,1)+a(i+1,1,1)+a(i,1,2)+a(i,2,1))-a(i,1,1)
        res(i,lat,1)=0.25*(a(i-1,lat,1)+a(i+1,lat,1)+a(i,lat,2)+a(i,lat-1,1))-a(i,lat,1)
        res(i,1,n)=0.25*(a(i-1,1,n)+a(i+1,1,n)+a(i,1,n-1)+a(i,2,n))-a(i,1,n)
        res(i,lat,n)=0.25*(a(i-1,lat,n)+a(i+1,lat,n)+a(i,lat-1,n)+a(i,lat,n-1))-a(i,lat,n)
      end do

      do j=2,lat-1
        res(1,j,1)=0.25*(a(1,j-1,1)+a(1,j+1,1)+a(1,j,2)+a(2,j,1))-a(1,j,1)
        res(lon,j,1)=0.25*(a(lon,j-1,1)+a(lon,j+1,1)+a(lon,j,2)+a(lon-1,j,1))-a(lon,j,1)
        res(1,j,n)=0.25*(a(1,j-1,n)+a(1,j+1,n)+a(1,j,n-1)+a(2,j,n))-a(1,j,n)
        res(lon,j,n)=0.25*(a(lon,j-1,n)+a(lon,j+1,n)+a(lon,j,n-1)+a(lon-1,j,n))-a(lon,j,n)
      end do

      do k = 2,n-1
       res(1,1,k)=0.25*(a(1,1,k-1)+a(1,1,k+1)+a(2,1,k)+a(1,2,k))-a(1,1,k)
       res(lon,1,k)=0.25*(a(lon,1,k-1)+a(lon,1,k+1)+a(lon-1,1,k)+a(lon,2,k))-a(lon,1,k)
       res(1,lat,k)=0.25*(a(1,lat,k-1)+a(1,lat,k+1)+a(2,lat,k)+a(1,lat-1,k))-a(1,lat,k)
       res(lon,lat,k)=0.25*(a(lon,lat,k-1)+a(lon,lat,k+1)+a(lon-1,lat,k)+a(lon,lat-1,k))-a(lon,lat,k)
      end do

      res(1,1,1)=0.3333*(a(2,1,1)+a(1,2,1)+a(1,1,2))-a(1,1,1)
      res(lon,1,1)=0.3333*(a(lon-1,1,1)+a(lon,2,1)+a(lon,1,2))-a(lon,1,1)
      res(1,lat,1)=0.3333*(a(2,lat,1)+a(1,lat-1,1)+a(1,lat,2))-a(1,lat,1)
      res(lon,lat,1)=0.3333*(a(lon-1,lat,1)+a(lon,lat-1,1)+a(lon,lat,2))-a(lon,lat,1)
      res(1,1,n)=0.3333*(a(2,1,n)+a(1,2,n)+a(1,1,n-1))-a(1,1,n)
      res(lon,1,n)=0.3333*(a(lon-1,1,n)+a(lon,2,n)+a(lon,1,n-1))-a(lon,1,n)
      res(1,lat,n)=0.3333*(a(2,lat,n)+a(1,lat-1,n)+a(1,lat,n-1))-a(1,lat,n)
      res(lon,lat,n)=0.3333*(a(lon-1,lat,n)+a(lon,lat-1,n)+a(lon,lat,n-1))-a(lon,lat,n)


      res=res*sor
      a=a+res
    end do

END SUBROUTINE extrap3d

SUBROUTINE vertical_interp(nx,Nlvl,zin,depth,z_rho,input,output)
    implicit none
    integer, intent(in) :: nx,Nlvl,zin
    real(kind=8), intent(in) :: depth(zin),z_rho(nx,Nlvl)
    real(kind=8), intent(in) :: input(nx,zin)
    real(kind=8), intent(out) :: output(nx,Nlvl)

    integer :: i,k,m
    real(kind=8) :: y0,y1,x0,x1,x

    do i = 1, nx
    do k = Nlvl,1,-1
      do m = zin,2,-1
        if (z_rho(i,k).le.depth(m) .and. z_rho(i,k).gt.depth(m-1)) then
          x = z_rho(i,k)
          x0 = depth(m-1)
          x1 = depth(m)
          y0 = input(i,m-1)
          y1 = input(i,m)
          output(i,k) = y0+(y1-y0)*(x-x0)/(x1-x0)
          EXIT
        else
           continue
        end if
      end do
    end do
    end do
END SUBROUTINE vertical_interp

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


!spline3d(xin,yin,zin+1,lon_tr,lat_tr,depth_tr,&
!   span,nz,lon(istr:iend,jbc),lat(istr:iend,jbc),
!   depth(istr:iend,:),tempin(:,:,:,l),temp_day(:,:,l))
SUBROUTINE spline3d(nxin,nyin,nzin,inlon,inlat,indp,&
            nxout,nzout,outlon,outlat,outdp,&
            input, output)
    implicit none
    integer, intent(in) :: nxin, nyin, nzin, nxout, nzout
    real(kind=8), intent(in) :: inlon(nxin), inlat(nyin),indp(nzin)
    real(kind=8), intent(in) :: outlon(nxout), outlat(nxout),outdp(nxout,nzout)
    real(kind=8), intent(in) :: input(nxin,nyin,nzin)
    real(kind=8), intent(out) :: output(nxout,nzout)

    integer :: i,j,k,x,z
    real(kind=8) :: splint
    real(kind=8) :: y2a(nxin,nyin,nzin)
    real(kind=8) :: slice(nxin,nzin), y2_slice(nxin,nzin)
    real(kind=8) :: line(nzin), y2_line(nzin)

    do i = 1, nxin
      do k = 1, nzin
        call spline(inlat,input(i,:,k),nyin,y2a(i,:,k)) 
      end do
    end do 

    do x = 1, nxout

      do i = 1, nxin
        do k = 1, nzin
          slice(i,k) = splint(inlat,input(i,:,k),y2a(i,:,k),nyin,outlat(x)) 
        end do
      end do
      
      do k = 1, nzin
        call spline(inlon,slice(:,k),nxin,y2_slice(:,k))
      end do

      do k = 1, nzin
        line(k) = splint(inlon,slice(:,k),y2_slice(:,k),nxin,outlon(x))
      end do

      call spline(indp,line,nzin,y2_line)
      do z = 1, nzout
        output(x,z) = splint(indp,line,y2_line,nzin,outdp(x,z))
      end do

    end do

END SUBROUTINE spline3d
