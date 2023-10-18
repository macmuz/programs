PROGRAM read_files
    USE NETCDF
    implicit none 
    real :: r
    integer :: ncid2,tdimid,cnt,idx(2),minidx(2)
    integer :: fstart(3),fstop(3),cnt2
    integer :: date0(3),ncid,dimid,nt,vid,date(3),date1(3)
    integer :: i,reason,NstationFiles,iStation,days,t
    integer :: tmask(560,600),dstart,dstop
    integer, allocatable :: varid(:),coord(:,:)
    real(kind=8) :: tlon(560,600),tlat(560,600),dis_array(560,600)
    real(kind=8), allocatable :: time(:),myzeta(:,:),mytime(:)
    real(kind=4) :: mylon,mylat,zeta(560,600,24),timetmp(24) 
    real(kind=4), allocatable :: lon(:),lat(:)
    character(LEN=100), dimension(:), allocatable :: stationFileNames,shortNames
    character(len=100) :: fname,path

    101 format(a,'/ocean_his_',i4.4,'-',i2.2,'-',i2.2,'.nc')

    date0 = (/1950,1,1/)
    fname = '/users/work/mmuzyka/CSDIR/input_560x600/ROMS_grid_125NM_may.nc'
    path = '/users/magazyn/mmuzyka/ROMS_CASES/metro_560x600_uerraMY'

    fstart = (/1993,1,3/)
    fstop = (/2007,12,31/)

    cnt = 0

    CALL elapsed(fstart,date0,dstart)
    CALL elapsed(fstop,date0,dstop)

    write(*,*) 'dstart=',dstart
    write(*,*) 'dstop=',dstop

    CALL check(nf90_open(trim(fname),NF90_NOWRITE,ncid),310)
    CALL check(nf90_inq_varid(ncid,"lon_rho",vid),27)
    CALL check(nf90_get_var(ncid,vid,tlon),28)
    CALL check(nf90_inq_varid(ncid,"lat_rho",vid),27)
    CALL check(nf90_get_var(ncid,vid,tlat),28)
    CALL check(nf90_inq_varid(ncid,"mask_rho",vid),27)
    CALL check(nf90_get_var(ncid,vid,tmask),28)
    CALL check(nf90_close(ncid),360)

!    write(*,*) tlon(1,1),tlat(1,1),tmask(1,1) 

    ! get the files
    call system('ls ./*TG*.nc > fileContents.txt')
    open(31,FILE='fileContents.txt',action="read")
    !how many
    i = 0
    do
     read(31,FMT='(a)',iostat=reason) r
     if (reason/=0) EXIT
     i = i+1
    end do
    NstationFiles = i
    write(*,'(a,I0)') "Number of station files: " , NstationFiles
    allocate(stationFileNames(NstationFiles),shortNames(NstationFiles))
    allocate(varid(1+NstationFiles),coord(NstationFiles,2))
    allocate(mytime((1+dstop-dstart)*24),myzeta(NstationFiles,(1+dstop-dstart)*24))
    rewind(31)

    
    CALL check(nf90_create('zeta_uerraMY.nc',NF90_NETCDF4,ncid2),310)
    CALL check( nf90_def_dim( ncid2, 'time', NF90_UNLIMITED, tdimid ),503)
    CALL check( nf90_def_var( ncid2, 'time', NF90_DOUBLE, (/tdimid/), varid(1) ), 504)
    CALL check( nf90_put_att(ncid2, varid(1), 'long_name', 'Time'), 504)
    CALL check( nf90_put_att(ncid2, varid(1), 'units', &
    'days since 1950-01-01T00:00:00Z'), 504)
    CALL check( nf90_put_att(ncid2, varid(1), 'calendar', 'standard'), 504)
    CALL check( nf90_put_att(ncid2, varid(1), 'axis', 'T'), 504)


    do i = 1,NstationFiles
     read(31,'(a)') stationFileNames(i)
     CALL check(nf90_open(trim(stationFileNames(i)),NF90_NOWRITE,ncid),310)
     CALL check(nf90_inq_dimid(ncid, "TIME", dimid),311)
     CALL check(nf90_inquire_dimension(ncid, dimid, len=nt),312)
     allocate(time(nt),lon(nt),lat(nt))
     CALL check(nf90_inq_varid(ncid,"TIME",vid),27)
     CALL check(nf90_get_var(ncid,vid,time),28)
     CALL check(nf90_inq_varid(ncid,"LONGITUDE",vid),27)
     CALL check(nf90_get_var(ncid,vid,lon),28)
     CALL check(nf90_inq_varid(ncid,"LATITUDE",vid),27)
     CALL check(nf90_get_var(ncid,vid,lat),28)

     date = date0
     days = 0
     if (int(time(1)).lt.0) then
       do
         CALL sub_day(date,.true.)
         days = days-1 
         if (days.eq.int(time(1))) then
           date1 = date
           exit
         end if
       enddo
     else
       do
         CALL add_day(date,.true.)
         days = days+1
         if (days.eq.int(time(1))) then
           date1 = date
           exit
         end if
       enddo
     end if

     do 
        CALL add_day(date,.true.) 
        days = days+1 
        if (days.eq.int(time(nt))) then
          EXIT
        end if
     end do
        
!     if (date1(1).le.2003 .and. date(1).ge.2000) then
     if ( (dstart.ge.int(time(1)).and.dstart.le.int(time(nt))).or.&
          (dstop.ge.int(time(1)).and.dstop.le.int(time(nt))).or.&
          (dstart.lt.int(time(1)).and.dstart.gt.int(time(nt))) ) then
       cnt = cnt+1 
       write(*,*) trim(stationFileNames(i))
       write(*,*) int(time(1)),int(time(nt))

       write(*,*) date1(1),date1(2),date1(3)
       write(*,*) date(1),date(2),date(3)
       idx(1) = INDEX(stationFileNames(i), '_', back=.true.)
       idx(2) = INDEX(stationFileNames(i), '.nc')
       shortNames(cnt) = stationFileNames(i)(idx(1)+1:idx(2)-1) 
       mylon = sum(lon)/real(nt) 
       mylat = sum(lat)/real(nt) 
       write(*,*) trim(shortNames(cnt)), mylon, mylat 
      
       dis_array = sqrt((tlon-mylon)**2+(tlat-mylat)**2)
       do 
         minidx = minloc(dis_array)
         if(tmask(minidx(1),minidx(2)).eq.1) then
           coord(cnt,:) = minidx
           write(*,*) coord(cnt,1),coord(cnt,2)
           EXIT
         else
           dis_array(minidx(1),minidx(2))=maxval(dis_array)
         end if
       end do  

       CALL check( nf90_def_var( ncid2, trim(shortNames(cnt)), &
         NF90_DOUBLE, (/tdimid/), varid(1+cnt) ), 504)
       CALL check( nf90_put_att(ncid2, varid(1+cnt), 'long_name',&
         'Water surface height'), 504)     
       CALL check( nf90_put_att(ncid2, varid(1+cnt), 'units',&
         'm'), 504)     
       CALL check( nf90_put_att(ncid2, varid(1+cnt), 'lon',&
         coord(cnt,1)), 504)     
       CALL check( nf90_put_att(ncid2, varid(1+cnt), 'lat',&
         coord(cnt,2)), 504)     

     end if

     deallocate(time,lon,lat) 
     CALL check(nf90_close(ncid),360)
    end do 
    close(31)

    CALL check( nf90_enddef(ncid2), 516 )

    cnt2 = 0
    date = fstart
    CALL elapsed(date,date0,days)
    do
      do t = 1,24
        timetmp(t) = real(days,8)+real(t,8)/24.0
      end do
!      CALL check( nf90_put_var(ncid2, varid(1), mytime, start=(/cnt2+1/)), 504)
      mytime(cnt2+1:cnt2+24) = timetmp

      write(fname,101) trim(path),date(1),date(2),date(3)
      write(*,*) trim(fname)

      CALL check(nf90_open(trim(fname),NF90_NOWRITE,ncid),310)
      CALL check(nf90_inq_varid(ncid,"zeta",vid),27)
      CALL check(nf90_get_var(ncid,vid,zeta),28)
      CALL check(nf90_close(ncid),360)

      do t = 1,cnt
        myzeta(t,cnt2+1:cnt2+24) = zeta(coord(t,1),coord(t,2),:)
!        CALL check( nf90_put_var(ncid2, varid(1+t), myzeta, start=(/cnt2+1/)),504) 
      end do

      if (all(date.eq.fstop)) EXIT
      CALL add_day(date,.true.)
      days = days+1
      cnt2 = cnt2+24 
    end do

    CALL check( nf90_put_var(ncid2, varid(1), mytime), 504)
    do t = 1,cnt
      myzeta(t,:) = myzeta(t,:)-sum(myzeta(t,:))/real((1+dstop-dstart)*24)
      CALL check( nf90_put_var(ncid2, varid(1+t), myzeta(t,:)),504)
    end do

    deallocate(stationFileNames,shortNames,varid,coord,mytime,myzeta)

    CALL check(nf90_close(ncid2),360)
END PROGRAM read_files

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
   
SUBROUTINE sub_day(date,use_leap)
    implicit none
    integer, intent(inout) :: date(3)
    logical, intent(in) :: use_leap

    logical :: leap
    integer :: days_in_month(12)

    days_in_month = (/31,28,31,30,31,30,31,31,30,31,30,31/)

    if ( date(2).eq.3 .and. date(3).eq.1 .and. use_leap) then
        leap = .false.
        if ( mod(date(1),4).eq.0 ) leap = .true.
        if ( mod(date(1),100).eq.0 ) leap = .false.
        if ( mod(date(1),400).eq.0 ) leap = .true.

        if ( leap ) days_in_month(2) = 29
    end if

    date(3) = date(3)-1
    if ( date(3).eq.0 ) then
        date(2) = date(2)-1
        if ( date(2).eq.0 ) then
            date(2) = 12
            date(1) = date(1)-1
        end if
        date(3) = days_in_month(date(2))
    end if
END SUBROUTINE sub_day

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
