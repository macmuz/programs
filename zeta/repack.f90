PROGRAM repack
    USE NETCDF
    implicit none 
    real :: r
    integer :: ncid2,tdimid,cnt,idx(2),minidx(2)
    integer :: fstart(3),fstop(3),cnt2
    integer :: ncid3,varid3
    integer :: date0(3),ncid,dimid,nt,vid,date(3),date1(3)
    integer :: i,reason,NstationFiles,iStation,days,t,t2
    integer :: dstart,dstop,fv,myloc(1)
    integer, allocatable :: varid(:),tmp(:)
    real(kind=8), allocatable :: mytime(:),myzeta(:,:)
    real(kind=8), allocatable :: time(:),zeta(:),dist(:)
    real(kind=8), allocatable :: modeltime(:),modelzeta(:)
    real(kind=8) :: s1,s2,ave1,ave2,nom,myfv
    real(kind=4) :: mylon,mylat,af,sf,dt,suma
    real(kind=4), allocatable :: lon(:),lat(:),coord(:,:)
    character(LEN=100), dimension(:), allocatable :: stationFileNames,shortNames
    character(len=100) :: fname,path,modelfile,statfile,output

    date0 = (/1950,1,1/)

    fstart = (/1993,1,3/)
    fstop = (/2007,12,31/)

    cnt = 0
    dt = 0.007
    myfv = -999.0

    modelfile = 'zeta_uerraMY.nc' 
    statfile = 'correlation_uerraMY.txt'
    output = 'zeta_copernicus_uerraMY.nc'

    CALL elapsed(fstart,date0,dstart) 
    CALL elapsed(fstop,date0,dstop)

    write(*,*) 'dstart=',dstart 
    write(*,*) 'dstop=',dstop 
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
    allocate(modeltime((1+dstop-dstart)*24),modelzeta((1+dstop-dstart)*24))
    rewind(31)

    do t = dstart,dstop
    do i = 1,24
       cnt = cnt+1
       mytime(cnt) = t+real(i,8)/24.0
    end do
    end do
    myzeta = myfv
   
    cnt = 0 
    CALL check(nf90_create(trim(output),NF90_NETCDF4,ncid2),310)
    CALL check( nf90_def_dim( ncid2, 'time', NF90_UNLIMITED, tdimid ),503)
    CALL check( nf90_def_var( ncid2, 'time', NF90_DOUBLE, (/tdimid/), varid(1) ), 504)
    CALL check( nf90_put_att(ncid2, varid(1), 'long_name', 'Time'), 504)
    CALL check( nf90_put_att(ncid2, varid(1), 'units', &
    'days since 1950-01-01T00:00:00Z'), 504)
    CALL check( nf90_put_att(ncid2, varid(1), 'calendar', 'standard'), 504)
    CALL check( nf90_put_att(ncid2, varid(1), 'axis', 'T'), 504)


    do i = 1,NstationFiles
     read(31,'(a)') stationFileNames(i)
!     write(*,*) trim(stationFileNames(i))
     CALL check(nf90_open(trim(stationFileNames(i)),NF90_NOWRITE,ncid),310)
     CALL check(nf90_inq_dimid(ncid, "TIME", dimid),311)
     CALL check(nf90_inquire_dimension(ncid, dimid, len=nt),312)
     allocate(time(nt),zeta(nt),lon(nt),lat(nt),tmp(nt),dist(nt))
     CALL check(nf90_inq_varid(ncid,"TIME",vid),27)
     CALL check(nf90_get_var(ncid,vid,time),28)
     CALL check(nf90_inq_varid(ncid,"LONGITUDE",vid),29)
     CALL check(nf90_get_var(ncid,vid,lon),30)
     CALL check(nf90_inq_varid(ncid,"LATITUDE",vid),31)
     CALL check(nf90_get_var(ncid,vid,lat),32)
     reason = nf90_inq_varid(ncid,"SLEV",vid)
     if (reason.eq.NF90_NOERR) then
     CALL check(nf90_get_var(ncid,vid,tmp,start=(/1,1/),&
          count=(/1,nt/)),34)
     CALL check(nf90_get_att(ncid,vid,'_FillValue',fv),35)
     CALL check(nf90_get_att(ncid,vid,'scale_factor',sf),36)
     CALL check(nf90_get_att(ncid,vid,'add_offset',af),37)
     zeta = real(tmp)*sf+af
     where(tmp.eq.fv) zeta = myfv
     CALL check(nf90_inq_varid(ncid,"SLEV_QC",vid),38)
     CALL check(nf90_get_var(ncid,vid,tmp,start=(/1,1/),&
          count=(/1,nt/)),39)
     where(tmp.ne.1) zeta = myfv  
     CALL check(nf90_inq_varid(ncid,"TIME_QC",vid),40)
     CALL check(nf90_get_var(ncid,vid,tmp),41)
     where(tmp.ne.1) zeta = myfv
     where(zeta.lt.-900.0) zeta=myfv
     elseif(reason.eq.NF90_ENOTVAR) then
        write(*,*) 'SLEV variable not found'
        deallocate(time,lon,lat,zeta,tmp,dist)
        CYCLE
     end if  
     CALL check(nf90_close(ncid),360)
!     write(*,*) zeta(1:5)

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
       coord(cnt,1) = mylon 
       coord(cnt,2) = mylat 
      
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
       CALL check( nf90_put_att(ncid2, varid(1+cnt), '_FillValue',&
         myfv), 504)     

       cnt2 = 1
       do t = 1,(1+dstop-dstart)*24
         do t2 = cnt2,nt-1
            if ( (time(t2).gt.mytime(t)) .and. ((time(t2)-mytime(t)).gt.dt) ) EXIT
            if ( (mytime(t).ge.time(t2)) .and. (mytime(t).lt.time(t2+1)) ) then
              if ( mytime(t)-time(t2).lt.dt ) then
                myzeta(cnt,t) = zeta(t2) 
                cnt2 = t2
                EXIT
              elseif ( time(t2+1)-mytime(t).lt.dt ) then
                myzeta(cnt,t) = zeta(t2+1) 
                cnt2 = t2
                EXIT
              endif
              cnt2 = t2
              EXIT
            end if
         end do
!         dist = abs(time-mytime(t))
!         myloc = minloc(dist)
!         myzeta(cnt,t) = zeta(myloc(1))
       end do
     end if
!     if (cnt.eq.5) EXIT

     deallocate(time,lon,lat,zeta,tmp,dist) 
    end do 
    close(31)

    CALL check( nf90_enddef(ncid2), 516 )

    CALL check( nf90_put_var(ncid2, varid(1), mytime), 504)

!    cnt2 = 0
!    date = fstart
!    CALL elapsed(date,date0,days)
!    do
!      do t = 1,24
!        mytime(t) = real(days,8)+real(t,8)/24.0
!      end do
!      CALL check( nf90_put_var(ncid2, varid(1), mytime, start=(/cnt2+1/)), 504)
!
!      write(fname,101) trim(path),date(1),date(2),date(3)
!      write(*,*) trim(fname)

!      CALL check(nf90_open(trim(fname),NF90_NOWRITE,ncid),310)
!      CALL check(nf90_inq_varid(ncid,"zeta",vid),27)
!      CALL check(nf90_get_var(ncid,vid,zeta),28)
!      CALL check(nf90_close(ncid),360)

      do t = 1,cnt
        cnt2 = 0
        suma = 0.0
        do t2 = 1,(1+dstop-dstart)*24
          if (myzeta(t,t2).ne.myfv) then
            cnt2 = cnt2+1
            suma = suma+myzeta(t,t2)
          end if  
        end do
        myzeta(t,:) = myzeta(t,:)-suma/real(cnt2)
      end do


    CALL check(nf90_open(trim(modelfile),NF90_NOWRITE,ncid3),310)
    CALL check(nf90_inq_varid(ncid3,"time",varid3),27)
    CALL check(nf90_get_var(ncid3,varid3,modeltime),28)

    if (all(int(modeltime).eq.int(mytime))) write(*,*) 'TIME is OK'
      ! write to ascii file
    open(unit=1,file=trim(statfile),status='replace',form='formatted')
    11 format(A,F7.3,F7.3,I,F7.3)

    where(myzeta.lt.-900.0) myzeta=myfv

    do t = 1,cnt
      reason = nf90_inq_varid(ncid3,trim(shortNames(t)),varid3)
      if (reason.eq.NF90_NOERR) then
        CALL check(nf90_get_var(ncid3,varid3,modelzeta),28)
        cnt2 = 0
        ave1 = 0.0
        ave2 = 0.0
        s1 = 0.0
        s2 = 0.0
        nom = 0.0
        do i = 1,(1+dstop-dstart)*24
          if (myzeta(t,i).gt.-900.0) then
            cnt2 = cnt2+1
            ave1 = ave1+myzeta(t,i)
            ave2 = ave2+modelzeta(i)
          end if 
        end do
        ave1 = ave1/real(cnt2)
        ave2 = ave2/real(cnt2)
        do i = 1,(1+dstop-dstart)*24
          if (myzeta(t,i).gt.-900.0) then
            nom = nom+(myzeta(t,i)-ave1)*(modelzeta(i)-ave2)
            s1 = s1+(myzeta(t,i)-ave1)**2
            s2 = s2+(modelzeta(i)-ave2)**2
          end if
        end do
        write(1,11)  trim(shortNames(t)),coord(t,2),coord(t,1),cnt2,&
            nom/sqrt(s1*s2)
        write(1,*) 
      end if
!        myzeta = zeta(coord(t,1),coord(t,2),:)
      CALL check( nf90_put_var(ncid2, varid(1+t), myzeta(t,:)),504) 
    end do

!      if (all(date.eq.fstop)) EXIT
!      CALL add_day(date,.true.)
!      days = days+1
!      cnt2 = cnt2+24 
!    end do

    deallocate(stationFileNames,shortNames,varid,coord,mytime,myzeta)
    deallocate(modeltime,modelzeta)

    close(1)
    CALL check(nf90_close(ncid2),360)
    CALL check(nf90_close(ncid3),361)
END PROGRAM repack

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
