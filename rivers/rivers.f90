PROGRAM rivers
    USE NETCDF
    implicit none
    integer :: year0,year1,year2,i,t,z
    integer :: oncid,nd,zd,td,ovid(9)
    integer :: ncid,dimid,nin,zin,tin,vid
    integer :: date0(3),date(3),julian,cnt,nt
    real(kind=8), allocatable :: time(:),tmp1d(:),tmp2d(:,:)
    real(kind=8), allocatable :: tmp3d(:,:,:),transport(:,:)
    real(kind=8), allocatable :: temp(:,:,:),trans(:,:),y0(:),y1(:),y(:)
    real(kind=8) :: x0,x1,x,xyear(365),xyearl(366)
    character(len=200) :: pref,path,outname
    logical :: yleap,first

    100 FORMAT(a,'/rivers_125NM_may_vertical_',i4.4,'_Lw.nc')
    pref = '/users/work/mmuzyka/CSDIR/input_560x600'
    
    year0 = 1992
    year1 = 2020

    year2 = 2023

    write(outname,100) '.',year2
    write(*,*) trim(outname)

    CALL leap(year2,yleap)
    nt = 365
    if (yleap) nt = 366

    xyear(1) = 0.0
    xyearl(1) = 0.0
    do t = 2,365
      xyear(t) = real(t-1,8)/real(364,8)
    end do
    do t = 2,366
      xyearl(t) = real(t-1,8)/real(365,8)
    end do
    

    CALL check(nf90_create( trim(outname), NF90_NETCDF4, oncid ), 200)
    CALL check(nf90_def_dim( oncid, 'river_time', NF90_UNLIMITED, td ), 201)

    first = .true.
    cnt = 0
    do i = year0,year1
      write(path,100) trim(pref),i
      write(*,*) trim(path)

      cnt = cnt+1

      CALL leap(i,yleap)

      CALL check(nf90_open(trim(path),NF90_NOWRITE,ncid),10)
      if (first) then
        CALL check(nf90_inq_dimid(ncid, "river", dimid),11)
        CALL check(nf90_inquire_dimension(ncid, dimid, len=nin),12)
        CALL check(nf90_inq_dimid(ncid, "s_rho", dimid),11)
        CALL check(nf90_inquire_dimension(ncid, dimid, len=zin),12)
        CALL check(nf90_inq_dimid(ncid, "river_time", dimid),11)
        CALL check(nf90_inquire_dimension(ncid, dimid, len=tin),12)

        allocate(time(nt),tmp1d(nin),tmp2d(nin,zin))
        allocate(tmp3d(nin,zin,nt),transport(nin,nt))
        allocate(temp(nin,zin,366),trans(nin,366))
        allocate(y0(nin),y1(nin),y(nin))

        CALL check(nf90_def_dim( oncid, 'river', nin, nd ), 201)
        CALL check(nf90_def_dim( oncid, 's_rho', zin, zd ), 201)

        CALL check(nf90_inq_varid(ncid,"river",vid),27)
        CALL check(nf90_def_var( oncid, 'river', NF90_DOUBLE, (/nd/), ovid(1)), 217)
        CALL check(nf90_copy_att(ncid, vid, "long_name", oncid, ovid(1)),28)

        CALL check(nf90_inq_varid(ncid,"river_direction",vid),29)
        CALL check(nf90_def_var( oncid, 'river_direction', NF90_DOUBLE, (/nd/), ovid(2)), 217)
        CALL check(nf90_copy_att(ncid, vid, "long_name", oncid, ovid(2)),30)
        CALL check(nf90_copy_att(ncid, vid, "flag_values", oncid, ovid(2)),31)
        CALL check(nf90_copy_att(ncid, vid, "flag_meanings", oncid, ovid(2)),32)

        CALL check(nf90_inq_varid(ncid,"river_Xposition",vid),33)
        CALL check(nf90_def_var( oncid, 'river_Xposition', NF90_DOUBLE, (/nd/), ovid(3)), 217)
        CALL check(nf90_copy_att(ncid, vid, "long_name", oncid, ovid(3)),34)
        CALL check(nf90_copy_att(ncid, vid, "LuvSrc_meaning", oncid, ovid(3)),35)
        CALL check(nf90_copy_att(ncid, vid, "LwSrc_meaning", oncid, ovid(3)),36)

        CALL check(nf90_inq_varid(ncid,"river_Eposition",vid),37)
        CALL check(nf90_def_var( oncid, 'river_Eposition', NF90_DOUBLE, (/nd/), ovid(4)), 217)
        CALL check(nf90_copy_att(ncid, vid, "long_name", oncid, ovid(4)),38)
        CALL check(nf90_copy_att(ncid, vid, "LuvSrc_True_meaning", oncid, ovid(4)),39)
        CALL check(nf90_copy_att(ncid, vid, "LwSrc_True_meaning", oncid, ovid(4)),40)

        CALL check(nf90_inq_varid(ncid,"river_Vshape",vid),27)
        CALL check(nf90_def_var( oncid, 'river_Vshape', NF90_DOUBLE, (/nd,zd/), ovid(5)), 217)
        CALL check(nf90_copy_att(ncid, vid, "long_name", oncid, ovid(5)),27)
        CALL check(nf90_copy_att(ncid, vid, "requires", oncid, ovid(5)),27)

        CALL check(nf90_inq_varid(ncid,"river_time",vid),27)
        CALL check(nf90_def_var( oncid, 'river_time', NF90_DOUBLE, (/td/), ovid(6)), 217)
        CALL check(nf90_copy_att(ncid, vid, "long_name", oncid, ovid(6)),27)
        CALL check(nf90_copy_att(ncid, vid, "units", oncid, ovid(6)),27)

        CALL check(nf90_inq_varid(ncid,"river_transport",vid),27)
        CALL check(nf90_def_var( oncid, 'river_transport', NF90_DOUBLE, (/nd,td/), ovid(7)), 217)
        CALL check(nf90_copy_att(ncid, vid, "long_name", oncid, ovid(7)),27)
        CALL check(nf90_copy_att(ncid, vid, "units", oncid, ovid(7)),27)
        CALL check(nf90_copy_att(ncid, vid, "positive", oncid, ovid(7)),27)
        CALL check(nf90_copy_att(ncid, vid, "negative", oncid, ovid(7)),27)
        CALL check(nf90_copy_att(ncid, vid, "time", oncid, ovid(7)),27)

        CALL check(nf90_inq_varid(ncid,"river_temp",vid),27)
        CALL check(nf90_def_var( oncid, 'river_temp', NF90_DOUBLE, (/nd,zd,td/), ovid(8)), 217)
        CALL check(nf90_copy_att(ncid, vid, "long_name", oncid, ovid(8)),27)
        CALL check(nf90_copy_att(ncid, vid, "units", oncid, ovid(8)),27)
        CALL check(nf90_copy_att(ncid, vid, "time", oncid, ovid(8)),27)

        CALL check(nf90_inq_varid(ncid,"river_salt",vid),27)
        CALL check(nf90_def_var( oncid, 'river_salt', NF90_DOUBLE, (/nd,zd,td/), ovid(9)), 217)
        CALL check(nf90_copy_att(ncid, vid, "long_name", oncid, ovid(9)),27)
        CALL check(nf90_copy_att(ncid, vid, "time", oncid, ovid(9)),27)

        CALL check( nf90_enddef( oncid ), 208 )

        date0=(/1968,5,23/)
        date=(/year2,1,1/)
        CALL elapsed(date,date0,julian)
        time(1) = real(julian)
        do t = 2,nt
          time(t) = time(t-1)+1.0

        end do

        CALL check( nf90_put_var( oncid, ovid(6), time ), 211 )

        CALL check(nf90_inq_varid(ncid,"river",vid),27)
        CALL check(nf90_get_var(ncid,vid,tmp1d),28)
        CALL check(nf90_put_var( oncid, ovid(1), tmp1d ), 211 )
        
        CALL check(nf90_inq_varid(ncid,"river_direction",vid),27)
        CALL check(nf90_get_var(ncid,vid,tmp1d),28)
        CALL check(nf90_put_var( oncid, ovid(2), tmp1d ), 211 )
        
        CALL check(nf90_inq_varid(ncid,"river_Xposition",vid),27)
        CALL check(nf90_get_var(ncid,vid,tmp1d),28)
        CALL check(nf90_put_var( oncid, ovid(3), tmp1d ), 211 )
        
        CALL check(nf90_inq_varid(ncid,"river_Eposition",vid),27)
        CALL check(nf90_get_var(ncid,vid,tmp1d),28)
        CALL check(nf90_put_var( oncid, ovid(4), tmp1d ), 211 )
        
        CALL check(nf90_inq_varid(ncid,"river_Vshape",vid),27)
        CALL check(nf90_get_var(ncid,vid,tmp2d),28)
        CALL check(nf90_put_var( oncid, ovid(5), tmp2d ), 211 )
       
        transport = 0.0
        tmp3d = 0.0 
        CALL check(nf90_put_var( oncid, ovid(9), tmp3d ), 211 )

!        CALL check(nf90_put_var( oncid, ovid(8), tmp3d ), 211 )
!        transport = 0.0
!        CALL check(nf90_put_var( oncid, ovid(7), transport ), 211 )
        first = .false.
      end if
      if (yleap) then
        CALL check(nf90_inq_varid(ncid,"river_transport",vid),27)
        CALL check(nf90_get_var(ncid,vid,trans),28)
        CALL check(nf90_inq_varid(ncid,"river_temp",vid),27)
        CALL check(nf90_get_var(ncid,vid,temp),28)
        if (nt.eq.366) then
          transport = transport+trans
          tmp3d = tmp3d+temp
        else
          transport(:,1) = transport(:,1)+trans(:,1)
          transport(:,365) = transport(:,365)+trans(:,366)
          tmp3d(:,:,1) = tmp3d(:,:,1)+temp(:,:,1)
          tmp3d(:,:,365) = tmp3d(:,:,365)+temp(:,:,366)
            
          do t = 2,364
            x0 = xyearl(t)
            x1 = xyearl(t+1)
            x = xyear(t)
            y0 = trans(:,t)
            y1 = trans(:,t+1)
            CALL linear(nin,x0,y0,x1,y1,x,y)
            transport(:,t) = transport(:,t)+y
            do z = 1,zin
              y0 = temp(:,z,t)
              y1 = temp(:,z,t+1)
              CALL linear(nin,x0,y0,x1,y1,x,y)
              tmp3d(:,z,t) = tmp3d(:,z,t)+y
            end do
!            write(*,*) x0,x1,x,x.gt.x0.and.x.le.x1 
          end do
        end if
      else
        CALL check(nf90_inq_varid(ncid,"river_transport",vid),27)
        CALL check(nf90_get_var(ncid,vid,trans(:,1:365)),28)
        CALL check(nf90_inq_varid(ncid,"river_temp",vid),27)
        CALL check(nf90_get_var(ncid,vid,temp(:,:,1:365)),28)
        if (nt.eq.365) then
          transport = transport+trans(:,1:365)
          tmp3d = tmp3d+temp(:,:,1:365)
        else
          transport(:,1) = transport(:,1)+trans(:,1)
          transport(:,366) = transport(:,366)+trans(:,365)
          tmp3d(:,:,1) = tmp3d(:,:,1)+temp(:,:,1)
          tmp3d(:,:,366) = tmp3d(:,:,366)+temp(:,:,365)
        
          do t = 2,365
            x0 = xyear(t-1)
            x1 = xyear(t)
            x = xyearl(t)
            y0 = trans(:,t-1)
            y1 = trans(:,t)
            CALL linear(nin,x0,y0,x1,y1,x,y)
            transport(:,t) = transport(:,t)+y
            do z = 1,zin
              y0 = temp(:,z,t-1)
              y1 = temp(:,z,t)
              CALL linear(nin,x0,y0,x1,y1,x,y)
              tmp3d(:,z,t) = tmp3d(:,z,t)+y
            end do
!            write(*,*) x0,x1,x,x.gt.x0.and.x.le.x1 
          end do
        end if
      end if
      CALL check(nf90_close(ncid),35)

    end do

    transport = transport/real(cnt,8)
    tmp3d = tmp3d/real(cnt,8)
    CALL check(nf90_put_var( oncid, ovid(7), transport ), 211 )
    CALL check(nf90_put_var( oncid, ovid(8), tmp3d ), 211 )
    CALL check(nf90_close(oncid),35)

    deallocate(time,tmp1d,tmp2d,tmp3d,transport,temp,trans)
    deallocate(y0,y1,y)
    
END PROGRAM rivers

SUBROUTINE linear(nr,x0,y0,x1,y1,x,y)
    implicit none
    integer, intent(in) :: nr
    real(kind=8), intent(in) :: x0,y0(nr),x1,y1(nr),x
    real(kind=8), intent(out) :: y(nr)
    
    y = y0*(x1-x)/(x1-x0)+y1*(x-x0)/(x1-x0)
END SUBROUTINE linear

SUBROUTINE check(status,label)
    USE NETCDF
    implicit none
    INTEGER, INTENT(in) :: status , label

    IF (status .ne. nf90_noerr) then
        print*, trim(nf90_strerror(status)),' label=',label
        stop "Stopped"
    END IF
END SUBROUTINE check

SUBROUTINE leap(year,answer)
    implicit none
    integer, intent(in) :: year
    logical, intent(out) :: answer

    answer = .false.
    if( mod(year,4).eq.0 ) answer = .true.
    if( mod(year,4).eq.100 ) answer = .false.
    if( mod(year,4).eq.400 ) answer = .true.
END SUBROUTINE leap

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
