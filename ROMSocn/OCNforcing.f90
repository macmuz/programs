PROGRAM OCNforcing
    USE NETCDF
    implicit none
    integer :: date(3),date0(3),enddate(3),date2(3)
    integer :: days,xdimid,ydimid,tdimid,dimids3(3)
    integer :: ncid,oncid,ovarid(7),datetmp(3)
    integer :: ni,nj,dimid,vid,k,cnt,cnt2
    character(len=20) :: oname
    character(len=100) :: ROMSpath,CICEpath
    character(len=200) :: fname
    real(kind=8) :: time
    real :: fv
    real(kind=4), allocatable :: tmp2d(:,:)
    logical :: closed,first
   
    100 format(a,'/ocean_his_',i4.4,'-',i2.2,'-',i2.2,'.nc')    
    101 format(a,'/iceh_01h.',i4.4,'-',i2.2,'-',i2.2,'-',i5.5,'.nc')
    102 format('ocean_',i4.4,'.nc')   
 
    ROMSpath = '/users/magazyn/mmuzyka/ROMS_CASES/metro_560x600_era5MYtest'
    CICEpath = '/users/magazyn/mmuzyka/ROMS_CASES/metro_560x600_era5MYtest/cice'

    date0=(/1968,5,23/)
    date=(/2001,1,1/)
    enddate=(/2001,12,31/)
    
    closed = .true.
    first = .true.

    write(fname,100) trim(ROMSpath),date(1),date(2),date(3)
!    write(*,*) trim(fname)
    CALL check(nf90_open(trim(fname),NF90_NOWRITE,ncid),310)
    CALL check(nf90_inq_dimid(ncid, "xi_rho", dimid),311)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=ni),312)

    CALL check(nf90_inq_dimid(ncid, "eta_rho", dimid),313)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nj),314)
    CALL check(nf90_close(ncid),360)

    allocate( tmp2d(ni,nj) )
   
    cnt = 0
    cnt2 = 0
    do

      if (closed) then
        write(oname,102) date(1)
        CALL check( nf90_create(trim(oname),NF90_NETCDF4,oncid), 310 )
        CALL check( nf90_def_dim(oncid,'ni',ni,xdimid ),503)
        CALL check( nf90_def_dim(oncid,'nj',nj,ydimid ),503)
        CALL check( nf90_def_dim(oncid,'time',NF90_UNLIMITED,tdimid ),503)
        dimids3 = (/xdimid,ydimid,tdimid/)

        CALL check( nf90_def_var(oncid, 'time', NF90_DOUBLE, &
         (/tdimid/), ovarid(1) ), 504)
        CALL check( nf90_put_att(oncid, ovarid(1), 'long_name',&
        'time since initialization' ), 305 )
        CALL check( nf90_put_att(oncid, ovarid(1), 'units',&
        'days since 1968-05-23 00:00:00' ), 305 )
        CALL check( nf90_put_att(oncid, ovarid(1), 'calendar',&
        'proleptic_gregorian' ), 305 )
        CALL check( nf90_put_att(oncid, ovarid(1), 'field',&
        'time, scalar, series' ), 305 )
        
        CALL check( nf90_def_var(oncid, 'zeta', NF90_FLOAT, &
         dimids3, ovarid(2) ), 504)       
        CALL check( nf90_put_att(oncid, ovarid(2), 'long_name',&
        'free-surface' ), 305 )        
        CALL check( nf90_put_att(oncid, ovarid(2), 'units',&
        'meter' ), 305 )
        write(fname,100) trim(ROMSpath),date(1),date(2),date(3)
        CALL check(nf90_open(trim(fname),NF90_NOWRITE,ncid),310)
        CALL check(nf90_inq_varid(ncid,"zeta",vid),17) 
        CALL check(nf90_get_att(ncid,vid,"_FillValue",fv),17)
        CALL check(nf90_close(ncid),360)
        CALL check( nf90_put_att(oncid, ovarid(2), '_FillValue',&
         fv), 305 )

        CALL check( nf90_def_var(oncid, 'sst', NF90_FLOAT, &
         dimids3, ovarid(3) ), 504)
        CALL check( nf90_put_att(oncid, ovarid(3), 'long_name',&
        'sea surface temperature' ), 305 )
        CALL check( nf90_put_att(oncid, ovarid(3), 'units',&
        'C' ), 305 )
        write(fname,101) trim(CICEpath),date(1),date(2),date(3),3600
        CALL check(nf90_open(trim(fname),NF90_NOWRITE,ncid),310)
        CALL check(nf90_inq_varid(ncid,"sst_h",vid),17)
        CALL check(nf90_get_att(ncid,vid,"_FillValue",fv),17)
        CALL check(nf90_close(ncid),360)
        CALL check( nf90_put_att(oncid, ovarid(3), '_FillValue',&
         fv), 305 )

        CALL check( nf90_def_var(oncid, 'sss', NF90_FLOAT, &
         dimids3, ovarid(4) ), 504)
        CALL check( nf90_put_att(oncid, ovarid(4), 'long_name',&
        'sea surface salinity' ), 305 )
        CALL check( nf90_put_att(oncid, ovarid(4), 'units',&
        'ppt' ), 305 )
        CALL check( nf90_put_att(oncid, ovarid(4), '_FillValue',&
         fv), 305 )

        CALL check( nf90_def_var(oncid, 'uocn', NF90_FLOAT, &
         dimids3, ovarid(5) ), 504)
        CALL check( nf90_put_att(oncid, ovarid(5), 'long_name',&
        'ocean current (x)' ), 305 )
        CALL check( nf90_put_att(oncid, ovarid(5), 'units',&
        'm/s' ), 305 )
        CALL check( nf90_put_att(oncid, ovarid(5), '_FillValue',&
         fv), 305 )

        CALL check( nf90_def_var(oncid, 'vocn', NF90_FLOAT, &
         dimids3, ovarid(6) ), 504)
        CALL check( nf90_put_att(oncid, ovarid(6), 'long_name',&
        'ocean current (y)' ), 305 )
        CALL check( nf90_put_att(oncid, ovarid(6), 'units',&
        'm/s' ), 305 )
        CALL check( nf90_put_att(oncid, ovarid(6), '_FillValue',&
         fv), 305 )

        CALL check( nf90_def_var(oncid, 'frzmlt', NF90_FLOAT, &
         dimids3, ovarid(7) ), 504)
        CALL check( nf90_put_att(oncid, ovarid(7), 'long_name',&
        'freeze/melt potential' ), 305 )
        CALL check( nf90_put_att(oncid, ovarid(7), 'units',&
        'W/m^2' ), 305 )
        CALL check( nf90_put_att(oncid, ovarid(7), '_FillValue',&
         fv), 305 )

        CALL check( nf90_enddef( oncid ), 207 )  
        
 
        closed = .false.
      end if

      if (first) then
      if ( date(2).ne.1 .or. date(3).ne.1 ) then
        datetmp = (/date(1),1,1/)
        do 
          write(*,*) datetmp
          CALL elapsed(datetmp,date0,days)
          do k = 1,24
            cnt = cnt+1
            cnt2 = cnt2+1
            tmp2d = 0.0
            time = real(days,8)+real(k,8)/24.0
            CALL check(nf90_put_var(oncid,ovarid(1),time,&
             start=(/cnt/)),18)
            CALL check(nf90_put_var(oncid,ovarid(2),tmp2d,&
             start=(/1,1,cnt/)),19)
            CALL check(nf90_put_var(oncid,ovarid(3),tmp2d,&
             start=(/1,1,cnt/)),19) 
            CALL check(nf90_put_var(oncid,ovarid(4),tmp2d,&
             start=(/1,1,cnt/)),19) 
            CALL check(nf90_put_var(oncid,ovarid(5),tmp2d,&
             start=(/1,1,cnt/)),19) 
            CALL check(nf90_put_var(oncid,ovarid(6),tmp2d,&
             start=(/1,1,cnt/)),19) 
            CALL check(nf90_put_var(oncid,ovarid(7),tmp2d,&
             start=(/1,1,cnt/)),19) 
          enddo
          CALL add_day(datetmp,.true.)
          if (all(datetmp==date)) exit
        enddo
      end if
      first = .false.
      end if


      CALL elapsed(date,date0,days)
      write(fname,100) trim(ROMSpath),date(1),date(2),date(3)
      CALL check(nf90_open(trim(fname),NF90_NOWRITE,ncid),310)
      CALL check(nf90_inq_varid(ncid,"zeta",vid),17)
      do k = 1,24
        cnt = cnt+1
        time = real(days,8)+real(k,8)/24.0
        CALL check(nf90_get_var(ncid,vid,tmp2d,&
         start=(/1,1,k/)),17)
        CALL check(nf90_put_var(oncid,ovarid(1),time,&
         start=(/cnt/)),18) 
        CALL check(nf90_put_var(oncid,ovarid(2),tmp2d,&
         start=(/1,1,cnt/)),19) 
      end do 
      CALL check(nf90_close(ncid),360)

      do k = 1,23
        cnt2 = cnt2+1
        write(fname,101) trim(CICEpath),date(1),date(2),date(3),3600*k
        CALL check(nf90_open(trim(fname),NF90_NOWRITE,ncid),310)

        CALL check(nf90_inq_varid(ncid,"sst_h",vid),17)
        CALL check(nf90_get_var(ncid,vid,tmp2d),17)
        CALL check(nf90_put_var(oncid,ovarid(3),tmp2d,&
         start=(/1,1,cnt2/)),19)

        CALL check(nf90_inq_varid(ncid,"sss_h",vid),17)
        CALL check(nf90_get_var(ncid,vid,tmp2d),17)
        CALL check(nf90_put_var(oncid,ovarid(4),tmp2d,&
         start=(/1,1,cnt2/)),19)

        CALL check(nf90_inq_varid(ncid,"uocn_h",vid),17)
        CALL check(nf90_get_var(ncid,vid,tmp2d),17)
        CALL check(nf90_put_var(oncid,ovarid(5),tmp2d,&
         start=(/1,1,cnt2/)),19)

        CALL check(nf90_inq_varid(ncid,"vocn_h",vid),17)
        CALL check(nf90_get_var(ncid,vid,tmp2d),17)
        CALL check(nf90_put_var(oncid,ovarid(6),tmp2d,&
         start=(/1,1,cnt2/)),19)

        CALL check(nf90_inq_varid(ncid,"frzmlt_h",vid),17)
        CALL check(nf90_get_var(ncid,vid,tmp2d),17)
        CALL check(nf90_put_var(oncid,ovarid(7),tmp2d,&
         start=(/1,1,cnt2/)),19)

        CALL check(nf90_close(ncid),360)
      end do

      date2 = date
      CALL add_day(date2,.true.)
      cnt2 = cnt2+1
      write(fname,101) trim(CICEpath),date2(1),date2(2),date2(3),0
      CALL check(nf90_open(trim(fname),NF90_NOWRITE,ncid),310)

      CALL check(nf90_inq_varid(ncid,"sst_h",vid),17)
      CALL check(nf90_get_var(ncid,vid,tmp2d),17)
      CALL check(nf90_put_var(oncid,ovarid(3),tmp2d,&
       start=(/1,1,cnt2/)),19)

      CALL check(nf90_inq_varid(ncid,"sss_h",vid),17)
      CALL check(nf90_get_var(ncid,vid,tmp2d),17)
      CALL check(nf90_put_var(oncid,ovarid(4),tmp2d,&
       start=(/1,1,cnt2/)),19)

      CALL check(nf90_inq_varid(ncid,"uocn_h",vid),17)
      CALL check(nf90_get_var(ncid,vid,tmp2d),17)
      CALL check(nf90_put_var(oncid,ovarid(5),tmp2d,&
       start=(/1,1,cnt2/)),19)

      CALL check(nf90_inq_varid(ncid,"vocn_h",vid),17)
      CALL check(nf90_get_var(ncid,vid,tmp2d),17)
      CALL check(nf90_put_var(oncid,ovarid(6),tmp2d,&
       start=(/1,1,cnt2/)),19)

      CALL check(nf90_inq_varid(ncid,"frzmlt_h",vid),17)
      CALL check(nf90_get_var(ncid,vid,tmp2d),17)
      CALL check(nf90_put_var(oncid,ovarid(7),tmp2d,&
       start=(/1,1,cnt2/)),19)

      CALL check(nf90_close(ncid),360)


      write(*,*) date(1),date(2),date(3),'elapsed:',days
      if (all(date==enddate)) then
        if ( date(2).ne.12 .or. date(3).ne.31 ) then
          datetmp = enddate
          do
            CALL add_day(datetmp,.true.)
            write(*,*) datetmp

            CALL elapsed(datetmp,date0,days)
            do k = 1,24
              cnt = cnt+1
              cnt2 = cnt2+1
              tmp2d = 0.0
              time = real(days,8)+real(k,8)/24.0
              CALL check(nf90_put_var(oncid,ovarid(1),time,&
                start=(/cnt/)),18)
              CALL check(nf90_put_var(oncid,ovarid(2),tmp2d,&
                start=(/1,1,cnt/)),19)
              CALL check(nf90_put_var(oncid,ovarid(3),tmp2d,&
                start=(/1,1,cnt/)),19)
              CALL check(nf90_put_var(oncid,ovarid(4),tmp2d,&
                start=(/1,1,cnt/)),19)
              CALL check(nf90_put_var(oncid,ovarid(5),tmp2d,&
                start=(/1,1,cnt/)),19)
              CALL check(nf90_put_var(oncid,ovarid(6),tmp2d,&
                start=(/1,1,cnt/)),19)
              CALL check(nf90_put_var(oncid,ovarid(7),tmp2d,&
                start=(/1,1,cnt/)),19)
            enddo
            if ( all(datetmp(2:3)==(/12,31/)) ) exit
          enddo 
        endif
        CALL check(nf90_close(oncid),360)
        EXIT
      end if
      if (all(date(2:3)==(/12,31/))) then
        CALL check(nf90_close(oncid),360)
        closed = .true.
        cnt = 0
        cnt2 = 0
      end if
      CALL add_day(date,.true.)
    end do

    deallocate(tmp2d)
END PROGRAM OCNforcing

SUBROUTINE check(status,label)
    USE NETCDF
    implicit none
    INTEGER, INTENT(in) :: status , label

    IF (status .ne. nf90_noerr) then
        print*, trim(nf90_strerror(status)),' label=',label
        stop "Stopped"
    END IF
END SUBROUTINE check

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
