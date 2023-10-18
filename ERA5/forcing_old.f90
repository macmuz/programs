PROGRAM forcing
    USE NETCDF
    implicit none
    integer :: i,j,k,ncid,nx,ny,nz,nt,dimid,varid
    integer :: date(3),date0(3),dayinmonth(12)
    integer :: narg,year,julian,cnt
    real, allocatable :: lon(:),lat(:)
    integer(kind=2), allocatable :: input(:,:,:,:)
    character(len=4) :: varname(8)
    character(len=10) :: buffer
    character(len=20) :: filename
    logical :: yleap,ex1

101 FORMAT('ERA5_',i4.4,'_',i2.2,'.nc')

    varname(1) = "u10"
    varname(2) = "v10"
    varname(3) = "d2m"
    varname(4) = "t2m"
    varname(5) = "msl"
    varname(6) = "ssr"
    varname(7) = "strd"
    varname(8) = "tp"

!READ command line arg: YEAR
    narg = command_argument_count()
    if (narg.ne.1) then
      write(*,*) "Program must have YEAR as argument"
      stop
    end if
    call get_command_argument(1,buffer)
    read(buffer, *) year
!END READ

    date0=(/1968,5,23/)
    date=(/1968,5,24/)
    dayinmonth = (/31,28,31,30,31,30,31,31,30,31,30,31/)

    CALL leap(year,yleap)
    nz = 365*24
    if (yleap) then 
      dayinmonth(2)=29
      nz = 366*24
    end if
   
    do i = 1,13
      if (i.eq.13) then
        write(filename,101) year+1,1
        INQUIRE(FILE=trim(filename),EXIST=ex1)
      else
        write(filename,101) year,i
      end if
      write(*,*) trim(filename)

      if (i.le.12 .or. ex1) then
      CALL check(nf90_open(trim(filename),NF90_NOWRITE,ncid),10)
      CALL check(nf90_inq_dimid(ncid, "time", dimid),8)
      CALL check(nf90_inquire_dimension(ncid, dimid, len=nt),9)
      if (i.eq.1) then
        CALL check(nf90_inq_dimid(ncid, "longitude", dimid),11)
        CALL check(nf90_inquire_dimension(ncid, dimid, len=nx),12)
        CALL check(nf90_inq_dimid(ncid, "latitude", dimid),13)
        CALL check(nf90_inquire_dimension(ncid, dimid, len=ny),14)
        allocate( lon(nx), lat(ny), input(nx,ny,nz,8) )
        CALL check(nf90_inq_varid(ncid,"longitude",varid),21)
        CALL check(nf90_get_var(ncid,varid,lon),22)
        CALL check(nf90_inq_varid(ncid,"latitude",varid),23)
        CALL check(nf90_get_var(ncid,varid,lat),24)
        do j = 1,8
          CALL check(nf90_inq_varid(ncid,varname(j),varid),23)
          CALL check(nf90_get_var(ncid,varid,input(:,:,1:nt-1,j), &
            start=(/1,1,2/),count=(/nx,ny,nt-1/) ),24)
        end do
        cnt = nt-1
      elseif (i.eq.13) then
        do j = 1,8
          CALL check(nf90_inq_varid(ncid,varname(j),varid),23)
          CALL check(nf90_get_var(ncid,varid,input(:,:,nz,j), &
            start=(/1,1,1/),count=(/nx,ny,1/) ),24)
        end do
      end if
      CALL check(nf90_close(ncid),35)
      end if
    end do

!    do i = 1, nx
!        write(*,*) lon(i)
!    end do
!    do i = 1, ny
!        write(*,*) lat(i)
!    end do
    

!    CALL elapsed(date,date0,julian)
!    write(*,*) julian

    deallocate( lon, lat, input )

END PROGRAM forcing

SUBROUTINE check(status,label)
    USE NETCDF
    implicit none
    INTEGER, INTENT(in) :: status , label

    IF (status .ne. nf90_noerr) then
        print*, trim(nf90_strerror(status)),' label= ',label
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

SUBROUTINE leap(year,answer)
    implicit none
    integer, intent(in) :: year
    logical, intent(out) :: answer

    answer = .false.
    if( mod(year,4).eq.0 ) answer = .true.
    if( mod(year,4).eq.100 ) answer = .false.
    if( mod(year,4).eq.400 ) answer = .true.
END SUBROUTINE leap
