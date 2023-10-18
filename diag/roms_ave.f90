PROGRAM roms_ave
    USE netcdf
    implicit none
    integer, parameter :: npt4d=1,nvar4d=2
    integer :: date(3),datestp(3),cnt,i,j,t
    integer :: ncid,dimid,ncid2,tdimid
    integer :: ni,nj,nk,nt,vid,varid(1+npt4d*nvar4d)
    integer, allocatable :: mask(:,:),ones(:,:,:)
    character(len=200) :: filename,path
    character(len=30) :: vname4d(nvar4d),varn
    real(kind=8) :: dplim(2),r,hc,pt4d(npt4d,2),npt
    real(kind=8), allocatable :: time(:),Cs_r(:),s_r(:),lon(:,:),lat(:,:)
    real(kind=8), allocatable :: h(:,:),zeta(:,:,:),dataout(:),dist(:,:)
    real(kind=4), allocatable :: datatmp(:,:,:,:)
    logical :: first
    logical, allocatable :: mask3d(:,:,:)

    100 format(a,'/ocean_avg_',i4,'-',i2.2,'-',i2.2,'.nc')
    101 format(a,'_pt',i2.2)
!    path = '/users/work/mmuzyka/CSDIR/metro_05NM_era5/run/baltic'
!    path = '/users/magazyn/mmuzyka/ROMS_CASES/metro_560x600_era5v2'
!    path = '/users/magazyn/mmuzyka/ROMS_CASES/metro_560x600_era5v2'
!    path = '/users/work/mmuzyka/CSDIR/metro_560x600_era5test4/run/baltic'
!    path = '/users/magazyn/mmuzyka/ROMS_CASES/metro_560x600_era5test4'
    path = '/users/magazyn/mmuzyka/ROMS_CASES/metro_560x600_uerraMY'
!    path = '/users/work/mmuzyka/CSDIR/metro_560x600_era2004v14/run/baltic'
    
    date = (/2007,1,1/)
    datestp = (/2007,12,31/)

    pt4d(1,:) = (/15.8,55.32/)
    dplim(:) = (/60,120/)
    r = 0.15
!    pt4d(1,:) = (/190,83,1/)

!    vname4d(1) = 'temp'
!    vname4d(2) = 'salt'
    vname4d(1) = 'salt'
    vname4d(2) = 'temp'

    first = .true.
    cnt = 0
    CALL check(nf90_create('uerraMY_salt_temp_08.nc',NF90_NETCDF4,ncid2),310)

    do
        write(filename,100) trim(path),date(1),date(2),date(3)
        write(*,*) trim(filename)        

        CALL check(nf90_open(trim(filename),NF90_NOWRITE,ncid),310)
        if (first) then
            CALL check(nf90_inq_dimid(ncid, "xi_rho", dimid),311)
            CALL check(nf90_inquire_dimension(ncid, dimid, len=ni),312)

            CALL check(nf90_inq_dimid(ncid, "eta_rho", dimid),313)
            CALL check(nf90_inquire_dimension(ncid, dimid, len=nj),314)
        
            CALL check(nf90_inq_dimid(ncid, "s_rho", dimid),313)
            CALL check(nf90_inquire_dimension(ncid, dimid, len=nk),314)

            CALL check(nf90_inq_dimid(ncid, "ocean_time", dimid),315)
            CALL check(nf90_inquire_dimension(ncid, dimid, len=nt),316)

            allocate( time(nt), datatmp(ni,nj,nk,nt), dataout(nt) )
            allocate( Cs_r(nk), s_r(nk), lon(ni,nj), lat(ni,nj) )
            allocate( h(ni,nj), zeta(ni,nj,nt), mask(ni,nj), dist(ni,nj) )
            allocate( mask3d(ni,nj,nk), ones(ni,nj,nk) )


            CALL check(nf90_inq_varid(ncid,"mask_rho",vid),27)
            CALL check(nf90_get_var(ncid,vid,ones(:,:,1)),28)
            CALL check(nf90_inq_varid(ncid,"lon_rho",vid),29)
            CALL check(nf90_get_var(ncid,vid,lon),30)
            CALL check(nf90_inq_varid(ncid,"lat_rho",vid),31)
            CALL check(nf90_get_var(ncid,vid,lat),32)
            CALL check(nf90_inq_varid(ncid,"h",vid),33)
            CALL check(nf90_get_var(ncid,vid,h),34)
            CALL check(nf90_inq_varid(ncid,"hc",vid),35)
            CALL check(nf90_get_var(ncid,vid,hc),36)
            CALL check(nf90_inq_varid(ncid,"Cs_r",vid),41)
            CALL check(nf90_get_var(ncid,vid,Cs_r),42)
            CALL check(nf90_inq_varid(ncid,"s_rho",vid),43)
            CALL check(nf90_get_var(ncid,vid,s_r),44)

            dist = sqrt((lon-pt4d(1,1))**2+(lat-pt4d(1,2))**2)
            mask = 0
            where(ones(:,:,1).eq.1 .and. dist.lt.r) mask=1
            write(*,*) pt4d(1,1),pt4d(1,2)
            write(*,*) sum(mask)
!            do i = 1,ni
!            do j = 1,nj
!              if (mask(i,j).eq.1) then
!                write(*,*) i,j
!              end if
!            end do
!            end do
            ones = 1

            CALL check(nf90_inq_varid(ncid,'ocean_time',vid),120)
            CALL check(nf90_get_var(ncid,vid,time),18)

            CALL check(nf90_inq_varid(ncid,'ocean_time',vid),120)

            CALL check( nf90_def_dim( ncid2, 'time', NF90_UNLIMITED, tdimid ),503)
            CALL check( nf90_def_var( ncid2, 'time', NF90_DOUBLE, (/tdimid/), varid(1) ), 504)
            CALL check( nf90_copy_att(ncid, vid, 'long_name', ncid2, varid(1)), 504)
            CALL check( nf90_copy_att(ncid, vid, 'units', ncid2, varid(1)), 504)
            CALL check( nf90_copy_att(ncid, vid, 'calendar', ncid2, varid(1)), 504)
            CALL check( nf90_copy_att(ncid, vid, 'field', ncid2, varid(1)), 504)

            do i = 1, nvar4d
            do j = 1, npt4d
            write(varn,101) trim(vname4d(i)),j
            CALL check( nf90_def_var( ncid2, trim(varn), NF90_FLOAT,&
              (/tdimid/),varid(1+j+(i-1)*npt4d) ), 504)
            end do
            end do

            CALL check( nf90_enddef(ncid2), 516 )
            first = .false.
        end if

        CALL check(nf90_inq_varid(ncid,'ocean_time',vid),120)
        CALL check(nf90_get_var(ncid,vid,time),18)
        CALL check( nf90_put_var( ncid2, varid(1), sum(time)/real(nt,8),&
             start=(/cnt+1/) ), 517 )

        CALL check(nf90_inq_varid(ncid,'zeta',vid),120)
        CALL check(nf90_get_var(ncid,vid,zeta),18)
        do i = 1, nvar4d
        CALL check(nf90_inq_varid(ncid,vname4d(i),vid),120)
        CALL check(nf90_get_var(ncid,vid,datatmp),18)
        do j = 1, npt4d
        do t = 1,nt
          CALL calc_dp(ni,nj,nk,zeta(:,:,t),h,hc,s_r,Cs_r,dplim,mask,mask3d)
!          write(*,*) 'mask sum=',sum(ones,mask=mask3d)
          npt = real(sum(ones,mask=mask3d),8) 
          dataout(t) = sum(datatmp(:,:,:,t),mask=mask3d)/npt
        end do
        CALL check( nf90_put_var( ncid2, varid(1+j+(i-1)*npt4d),&
            sum(dataout)/real(nt,8), start=(/cnt+1/) ), 517 )
        end do
        end do
        
        CALL check(nf90_close(ncid),360)

        cnt = cnt+1

        if (all(date.eq.datestp)) EXIT
        CALL add_day(date,.true.)
    end do

    CALL check(nf90_close(ncid2),360)
   
    deallocate(time,datatmp,Cs_r,s_r,lon,lat,dataout)
    deallocate(h,zeta,mask,dist,mask3d) 
    
END PROGRAM roms_ave

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

SUBROUTINE calc_dp(nx,ny,nz,zeta,h,hc,s_r,Cs_r,lim,mask,mask3d)
    implicit none
    integer, intent(in) :: nx,ny,nz
    real(kind=8), intent(in) :: zeta(nx,ny),h(nx,ny),hc,s_r(nz),Cs_r(nz)
    real(kind=8), intent(in) :: lim(2)
    integer, intent(in) :: mask(nx,ny)
    logical, intent(out) :: mask3d(nx,ny,nz)

    integer :: i,j,k
    real(kind=8) :: dp(nx,ny,nz)

    dp = 0.0

    do i = 1,nx
    do j = 1,ny
    if (mask(i,j).eq.1) then
      do k = 1,nz
        dp(i,j,k) = -1*(zeta(i,j)+(zeta(i,j)+h(i,j))*(hc*s_r(k)+h(i,j)*Cs_r(k))/&
            (hc+h(i,j)))
      end do
    end if
    end do
    end do

    mask3d = .false.
    where(dp.ge.lim(1) .and. dp.le.lim(2)) mask3d = .true.

!    write(*,*) 'minval=',minval(dp)
!    write(*,*) 'maxval=',maxval(dp)

END SUBROUTINE calc_dp
