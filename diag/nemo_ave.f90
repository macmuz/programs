PROGRAM nemo_temp
    USE netcdf
    implicit none
    integer :: date(3),datestp(3),cnt,i,j,t,k,l,zmax
    integer :: ncid,dimid,ncid2,tdimid,tmp(2)
    integer :: ni,nj,nk,nt,vid,varid(3)
    character(len=200) :: filename,path
    character(len=30) :: vname(2),varn
    character(len=60) :: nfile
    integer, allocatable :: ones(:,:,:)
    real(kind=8) :: pt(2),W(2,2),r,dplim(2),npt,fv
    real(kind=8), allocatable :: time(:),lon(:),lat(:)
    real(kind=8), allocatable :: dataout(:),dp(:),depth(:)
    real(kind=8), allocatable :: datatmp(:,:,:,:),dist(:,:)
    logical :: first,score
    logical, allocatable :: mask3d(:,:,:)

    100 format(a,'/',i4.4,'/',i2.2,'/',a)
    101 format('BAL-MYP-NEMO_PHY-DailyMeans-',i4.4,i2.2,i2.2,'.nc')
!    path = '/users/work/mmuzyka/CSDIR/metro_05NM_era5/run/baltic'
!    path = '/users/magazyn/mmuzyka/ROMS_CASES/metro_560x600_era5v2'
!    path = '/users/magazyn/mmuzyka/ROMS_CASES/metro_560x600_era5v2'
!    path = '/users/work/mmuzyka/CSDIR/metro_560x600_era5test4/run/baltic'
!    path = '/users/magazyn/mmuzyka/ROMS_CASES/metro_560x600_era5test4'
    path = '/users/work/mmuzyka/programs/ROMS_bc/new_files'
!    path = '/users/work/mmuzyka/CSDIR/metro_560x600_era5test5/run/baltic'
    
    date = (/1993,1,2/)
    datestp = (/2004,12,31/)

    pt = (/15.8,55.32/)
    dplim(:) = (/60,120/)
    r = 0.15

    vname(1) = 'thetao'
    vname(2) = 'so'

    first = .true.
    cnt = 0
    CALL check(nf90_create('nemo_pt_15_55_60_120.nc',NF90_NETCDF4,ncid2),310)

    do
        write(nfile,101) date(1),date(2),date(3)
        write(filename,100) trim(path),date(1),date(2),trim(nfile)
        write(*,*) trim(filename)        

        CALL check(nf90_open(trim(filename),NF90_NOWRITE,ncid),310)
        if (first) then
            CALL check(nf90_inq_dimid(ncid, "lon", dimid),311)
            CALL check(nf90_inquire_dimension(ncid, dimid, len=ni),312)

            CALL check(nf90_inq_dimid(ncid, "lat", dimid),313)
            CALL check(nf90_inquire_dimension(ncid, dimid, len=nj),314)
        
            CALL check(nf90_inq_dimid(ncid, "depth", dimid),313)
            CALL check(nf90_inquire_dimension(ncid, dimid, len=nk),314)

            CALL check(nf90_inq_dimid(ncid, "time", dimid),315)
            CALL check(nf90_inquire_dimension(ncid, dimid, len=nt),316)

            allocate( time(nt), datatmp(ni,nj,nk,nt), dataout(nt) )
            allocate( dp(nk), depth(nk), lon(ni), lat(nj) )
            allocate( mask3d(ni,nj,nk), dist(ni,nj), ones(ni,nj,nk) )


            ones = 1
            CALL check(nf90_inq_varid(ncid,"lon",vid),29)
            CALL check(nf90_get_var(ncid,vid,lon),30)
            CALL check(nf90_inq_varid(ncid,"lat",vid),31)
            CALL check(nf90_get_var(ncid,vid,lat),32)
            CALL check(nf90_inq_varid(ncid,"depth",vid),31)
            CALL check(nf90_get_var(ncid,vid,depth),32)

!            write(*,*) lon(1),lon(ni),lat(1),lat(nj)
!            write(*,*) dp(1),dp(nk)
!            CALL calc_dp(nk,dp)

!            do k = ni,1,-1
!              if (pt(1).ge.lon(k)) EXIT
!            end do
!            do l = nj,1,-1
!              if (pt(2).ge.lat(l)) EXIT
!            end do
!            write(*,*) pt
!            write(*,*) k,l
!            write(*,*) lon(k),lon(k+1)
!            write(*,*) lat(l),lat(l+1)

!            CALL coef(lon(k:k+1),lat(l:l+1),pt,W)
!            write(*,*) W(1,1),W(1,2),W(2,2),W(2,1)

            write(*,*) lon(1),lat(1),lon(ni),lat(nj),pt
            mask3d = .false.
            do i = 1,ni
            do j = 1,nj
              dist(i,j) = sqrt( (lon(i)-pt(1))**2+(lat(j)-pt(2))**2 )
!              if (dist(i,j).lt.r) write(*,*) 'idx',i,j,dist(i,j)
            end do
            end do
!            write(*,*) dist(244,140)
            do k = 1,nk
              where(dist.lt.r) mask3d(:,:,k) = .true.
              if (depth(k).lt.dplim(1) .or. depth(k).gt.dplim(2)) mask3d(:,:,k) = .false.
            end do
            npt = real(sum(ones,mask=mask3d),8)
            write(*,*) 'ONES=',npt 

            CALL check(nf90_inq_varid(ncid,'time',vid),120)

            CALL check( nf90_def_dim( ncid2, 'time', NF90_UNLIMITED, tdimid ),503)
            CALL check( nf90_def_var( ncid2, 'time', NF90_DOUBLE, (/tdimid/), varid(1) ), 504)
            CALL check( nf90_copy_att(ncid, vid, 'long_name', ncid2, varid(1)), 504)
            CALL check( nf90_copy_att(ncid, vid, 'units', ncid2, varid(1)), 504)
            CALL check( nf90_copy_att(ncid, vid, 'axis', ncid2, varid(1)), 504)

            do i = 1, 2
            CALL check( nf90_def_var( ncid2, trim(vname(i)), NF90_DOUBLE,&
              (/tdimid/),varid(1+i) ), 504)
            end do

            CALL check(nf90_inq_varid(ncid,vname(1),vid),120)
            CALL check(nf90_get_var(ncid,vid,datatmp),18)
            CALL check(nf90_get_att(ncid,vid,'_FillValue',fv),120)
       
            where( abs(datatmp(:,:,:,1)-fv).lt.1 ) mask3d = .false. 
            npt = real(sum(ones,mask=mask3d),8)
            write(*,*) 'ONES=',npt 

            CALL check( nf90_enddef(ncid2), 516 )
            first = .false.
        end if

        CALL check(nf90_inq_varid(ncid,'time',vid),120)
        CALL check(nf90_get_var(ncid,vid,time),18)
        CALL check( nf90_put_var( ncid2, varid(1), time, start=(/cnt+1/) ), 517 )


!        write(*,*) 'after dp'
        do i = 1,2
          CALL check(nf90_inq_varid(ncid,vname(i),vid),120)
          CALL check(nf90_get_var(ncid,vid,datatmp),18)
!          CALL ave(zmax,nt,datatmp(:,:,1:zmax,:),W,dp(1:zmax),dataout)
!          write(*,*) 'after ave'
          do t = 1,nt
            dataout(t) = sum(datatmp(:,:,:,t),mask=mask3d)/npt
          end do
          CALL check( nf90_put_var( ncid2, varid(1+i),&
            dataout, start=(/cnt+1/) ), 518 )
        end do
        
        CALL check(nf90_close(ncid),360)

        cnt = cnt+nt

        if (date(1).ge.datestp(1).and.&
            date(2).ge.datestp(2).and.&
            date(3).ge.datestp(3)) EXIT
        CALL add_day(date,.true.)
    end do

    CALL check(nf90_close(ncid2),360)
   
    deallocate(time,datatmp,dp,depth,lon,lat,dataout,mask3d,dist,ones)
    
END PROGRAM nemo_temp

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

SUBROUTINE calc_dp(nz,dp)
    implicit none
    integer, intent(in) :: nz
    real(kind=8), intent(inout) :: dp(nz)

    integer :: k
    real(kind=8) :: tmp(nz),prev

    tmp = dp
    prev = 0.0
    do k = 1,nz
      dp(k) = (tmp(k)-prev)*2
      prev = prev+dp(k)
!      write(*,*) dp(k) 
    end do
END SUBROUTINE calc_dp

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

SUBROUTINE coef(lon,lat,point,W)
    implicit none
    real(kind=8), intent(in) :: lon(2),lat(2)
    real(kind=8), intent(in) :: point(2)
    real(kind=8), intent(out) :: W(2,2)

    W(1,1) = (lon(2)-point(1))*(lat(2)-point(2))/(lon(2)-lon(1))/&
            (lat(2)-lat(1)) 
    W(1,2) = (lon(2)-point(1))*(point(2)-lat(1))/(lon(2)-lon(1))/&
            (lat(2)-lat(1)) 
    W(2,1) = (point(1)-lon(1))*(lat(2)-point(2))/(lon(2)-lon(1))/&
            (lat(2)-lat(1)) 
    W(2,2) = (point(1)-lon(1))*(point(2)-lat(1))/(lon(2)-lon(1))/&
            (lat(2)-lat(1)) 
END SUBROUTINE coef

SUBROUTINE ave(nk,nt,tmp,W,dp,output)
    implicit none
    integer, intent(in) :: nk,nt
    real(kind=8), intent(in) :: tmp(2,2,nk,nt), W(2,2), dp(nk)
    real(kind=8), intent(out) :: output(nt)

    integer :: z,t
    
    output = 0
    do t = 1,nt
    do z = 1,nk
      output(t) = output(t)+dp(z)*(tmp(1,1,z,t)*W(1,1)+tmp(1,2,z,t)*W(1,2)+&
            &   tmp(2,2,z,t)*W(2,2)+tmp(2,1,z,t)*W(2,1))
    end do
    output(t) = output(t)/real(sum(dp))
    write(*,*) output(t)
    end do
END SUBROUTINE ave
