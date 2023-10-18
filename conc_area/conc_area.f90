PROGRAM conc_area
    USE netcdf
    implicit none
    integer, parameter :: tmax = 168, nx = 5001, dt = 4
    integer :: date0(2),date(2),edate(2),dlevel(3),i,k
    integer :: cnt,curr_date(2),reclen,greclen
    integer :: xdimid,ydimid,tdimid,dimids3(3)
    integer :: ncid2,varid2(3),tdimid2
    integer :: ncid,varid(5),days,days2
    integer, allocatable :: mask(:,:)
    real(kind=8), allocatable :: lon(:),lat(:)
    real(kind=4), allocatable :: field2d(:,:)
    character(len=200) :: fname,gridname
    character(len=100) :: path,gridout,fnameout
    character(len=10) :: subpath(3),agent(3)
    logical :: cormix,ex

    cormix = .false.
    curr_date = (/2021,7/)

100 format(a,a,'/result_',a,'_',i4.4,i2.2,'01.bin')
101 format('currents_',i4.4,i2.2,'.nc')
102 format(a,a,'/lonlat_',i4.4,i2.2,'01.bin')
103 format(a,'_conc_area_',i4.4,i2.2,'01.nc')

    path = '/users/work/jakacki/SKAGERRAK_JOHN/vers04'
    subpath(1) = '/Tabun'
    subpath(2) = '/mustard'
    subpath(3) = '/Clark'
    agent(1) = 'tabun'
    agent(2) = 'mustard'
    agent(3) = 'clark'
    date0 = (/2018,12/)
    edate = (/2021,8/)
    
    if (cormix) then
      dlevel(1) = 1300.e-6
      dlevel(2) = 25000.e-6
      dlevel(3) = 290.e-6
    else
      dlevel(1) = 1.e-6
      dlevel(2) = 100.e-6
      dlevel(3) = 10.e-6
    end if
    
    allocate( lon(nx), lat(nx), field2d(nx,nx), mask(nx,nx) )

    cnt = 0
    do i = 1,tmax,dt
      cnt = cnt + 1
    end do
    write(*,*) cnt

    CALL unix((/curr_date(1),curr_date(2),1/),days)

    write(gridout,101) curr_date(1),curr_date(2) 
    CALL check( nf90_create( trim(gridout),NF90_NETCDF4,ncid ), 300 )
    CALL check( nf90_def_dim( ncid, 'longitude', nx, xdimid ), 301 )
    CALL check( nf90_def_dim( ncid, 'latitude', nx, ydimid ), 302 )
    CALL check( nf90_def_dim( ncid, 'time', cnt, tdimid ), 303 )
    dimids3 = (/xdimid,ydimid,tdimid/)

    CALL check(nf90_def_var( ncid, 'longitude', NF90_DOUBLE,&
        (/xdimid/), varid(1) ), 304)
    CALL check( nf90_put_att( ncid, varid(1), 'long_name',&
        'longitude' ), 305 )
    CALL check( nf90_put_att( ncid, varid(1), 'units',&
        'degree_east' ), 306 )

    CALL check(nf90_def_var( ncid, 'latitude', NF90_DOUBLE,&
        (/ydimid/), varid(2) ), 307)
    CALL check( nf90_put_att( ncid, varid(2), 'long_name',&
        'latitude' ), 308 )
    CALL check( nf90_put_att( ncid, varid(2), 'units',&
        'degree_north' ), 309 )

    CALL check(nf90_def_var( ncid, 'time', NF90_DOUBLE,&
        (/tdimid/), varid(3) ), 310)
    CALL check( nf90_put_att( ncid, varid(3), 'long_name',&
        'time' ), 311 )
    CALL check( nf90_put_att( ncid, varid(3), 'units',&
        'days since 1970-01-01 00:00:00 GMT' ), 312 )
    CALL check( nf90_put_att( ncid, varid(3), 'calendar',&
        'gregorian' ), 313 )

    CALL check(nf90_def_var( ncid, 'u', NF90_FLOAT,&
        dimids3, varid(4) ), 314)
    CALL check( nf90_put_att( ncid, varid(4), 'long_name',&
        'u component of surface currents' ), 315 )
    CALL check( nf90_put_att( ncid, varid(4), 'units',&
        'm/s' ), 316 )

    CALL check(nf90_def_var( ncid, 'v', NF90_FLOAT,&
        dimids3, varid(5) ), 317)
    CALL check( nf90_put_att( ncid, varid(5), 'long_name',&
        'v component of surface currents' ), 315 )
    CALL check( nf90_put_att( ncid, varid(5), 'units',&
        'm/s' ), 316 )

    CALL check( nf90_enddef( ncid ), 340 )

    inquire (iolength=reclen) field2d
    do i = 1, 3
      date = date0
      do
        write(fname,100) trim(path),trim(subpath(i)),trim(agent(i)),&
            date(1),date(2)
        INQUIRE(FILE=trim(fname), EXIST=ex)

        if (ex) then
       
          CALL unix((/date(1),date(2),1/),days2)        
 
          open(2,file=trim(fname),access='direct',&
            form='unformatted',recl=reclen)

        write(fnameout,103) trim(agent(i)),date(1),date(2)
        CALL check( nf90_create( trim(fnameout),NF90_NETCDF4,ncid2 ), 600 )
        CALL check( nf90_def_dim( ncid2, 'time', tmax, tdimid2 ), 601 )
        
        CALL check(nf90_def_var( ncid2, 'time', NF90_DOUBLE,&
          (/tdimid2/), varid2(1) ), 602)
        CALL check( nf90_put_att( ncid2, varid2(1), 'long_name',&
          'time' ), 603 )
        CALL check( nf90_put_att( ncid2, varid2(1), 'units',&
          'days since 1970-01-01 00:00:00 GMT' ), 604 )
        CALL check( nf90_put_att( ncid2, varid2(1), 'calendar',&
          'gregorian' ), 605 )

        CALL check(nf90_def_var( ncid2, 'conc', NF90_INT,&
          (/tdimid2/), varid2(2) ), 606)
        CALL check( nf90_put_att( ncid2, varid2(2), 'long_name',&
          'area with concentration above treshold' ), 607 )
        CALL check( nf90_put_att( ncid2, varid2(2), 'units',&
          'm**2' ), 608 )

        CALL check(nf90_def_var( ncid2, 'conchf', NF90_INT,&
          (/tdimid2/), varid2(3) ), 609)
        CALL check( nf90_put_att( ncid2, varid2(3), 'long_name',&
          'area with concentration (+hl) above treshold' ), 610 )
        CALL check( nf90_put_att( ncid2, varid2(3), 'units',&
          'm**2' ), 611 )
        
        CALL check( nf90_enddef( ncid2 ), 640 )

        do k = 1, tmax
          CALL check( nf90_put_var( ncid2, varid2(1),&
            real(days2,8)+real((k-1)/24.0,8), start=(/k/) ), 641 )
    
          mask = 0
          read(2,rec=2*tmax+k) field2d
          where(field2d.gt.dlevel(i)) mask = 1
          CALL check( nf90_put_var( ncid2, varid2(2),&
              100*sum(mask), start=(/k/) ), 642 )
    
          mask = 0
          read(2,rec=3*tmax+k) field2d
          where(field2d.gt.dlevel(i)) mask = 1
          CALL check( nf90_put_var( ncid2, varid2(3),&
              100*sum(mask), start=(/k/) ), 643 )
        end do


        CALL check( nf90_close( ncid2 ), 650 )

        if (date(1).eq.curr_date(1).and.date(2).eq.curr_date(2).and.&
            i.eq.1) then
          inquire (iolength=greclen) lon
          write(gridname,102) trim(path),trim(subpath(i)),date(1),date(2)
          write(*,*) trim(gridname)
          open(1,file=trim(gridname),access='direct',&
            form='unformatted',recl=greclen)
          read(1,rec=1) lon
          read(1,rec=2) lat
          close(1)
          CALL check( nf90_put_var( ncid, varid(1), lon ), 341 )
          CALL check( nf90_put_var( ncid, varid(2), lat ), 342 )

          cnt = 1
          do k = 1, tmax,dt
            CALL check( nf90_put_var( ncid, varid(3),&
              real(days,8)+real((k-1)/24.0,8), start=(/cnt/) ), 343 )

            read(2,rec=k) field2d
            CALL check( nf90_put_var( ncid, varid(4),&
              field2d, start=(/1,1,cnt/) ), 343 ) 

            read(2,rec=tmax+k) field2d
            CALL check( nf90_put_var( ncid, varid(5),&
              field2d, start=(/1,1,cnt/) ), 343 )

            cnt = cnt+1
          end do

        end if
          close(2)
        end if

        if (date(1).eq.edate(1).and.date(2).eq.edate(2)) EXIT
        date(2) = date(2)+1
        if (date(2).gt.12) then
          date(1) = date(1)+1
          date(2) = 1
        end if
      end do
    end do

    CALL check( nf90_close( ncid ), 350 )
    

    deallocate( lon, lat, field2d, mask )

END PROGRAM conc_area

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

SUBROUTINE unix(date,days)
    implicit none
    integer, intent(in) :: date(3)
    integer, intent(out) :: days

    integer :: date_0(3),date_tmp(3)

    date_0 = (/1970,1,1/)
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
END SUBROUTINE unix
