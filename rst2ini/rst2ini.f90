PROGRAM rst2ini
    USE netcdf
    implicit none
    integer :: narg,id,julian,nx,ny,nz,dimid
    integer :: ncid,ncid2,varid,varid2,ymd(3),date0(3)
    real(kind=8) :: time
    real, allocatable :: a2d(:,:),a3d(:,:,:),u(:,:,:),v(:,:,:)
    real, allocatable :: ubar(:,:),vbar(:,:)
    character(len=100) :: tmpfile,outfile
    character(len=250) :: infile,cmd
    character(len=10) :: date

    tmpfile = 'ROMS_grid_2_3km_560x600_NetCDF4_initial_template.nc'
    date0=(/1968,5,23/)
  
    narg = command_argument_count()
    if (narg.ne.1) then
      write(*,*) "Program must have 1 argument - rst filename"
      stop
    end if
    call get_command_argument(1,infile)

    id = INDEX(infile, '.nc', .true.)
    date = infile(id-10:id-1)
    outfile = 'ROMS_grid_2_3km_560x600_NetCDF4_initial_'//&
        date//'.nc'
    cmd = 'cp '//trim(tmpfile)//' '//trim(outfile)
    CALL EXECUTE_COMMAND_LINE(trim(cmd))

    read(date(1:4),*) ymd(1)
    read(date(6:7),*) ymd(2)
    read(date(9:10),*) ymd(3)
!    ymd(1) = 2014
!    ymd(2) = 6
!    ymd(3) = 26
    write(*,*) ymd(1),ymd(2),ymd(3)
    CALL elapsed(ymd,date0,julian)
    write(*,*) julian
    
    CALL check(nf90_open(trim(infile),NF90_NOWRITE,ncid),310)
    CALL check(nf90_open(trim(outfile),NF90_WRITE,ncid2),311)

    CALL check(nf90_inq_dimid(ncid, "xi_rho", dimid),11)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nx),12)
    CALL check(nf90_inq_dimid(ncid, "eta_rho", dimid),13)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=ny),14)
    CALL check(nf90_inq_dimid(ncid, "s_rho", dimid),15)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nz),16)
    allocate( a2d(nx,ny),a3d(nx,ny,nz),u(nx-1,ny,nz),v(nx,ny-1,nz) )
    allocate( ubar(nx-1,ny),vbar(nx,ny-1) )

    CALL check(nf90_inq_varid(ncid2,"ocean_time",varid2),319)
    CALL check(nf90_put_var(ncid2,varid2,real(julian,8)),420)
 
    CALL check(nf90_inq_varid(ncid,"zeta",varid),319) 
    CALL check(nf90_get_var(ncid,varid,a2d,start=(/1,1,3,1/),&
        count=(/nx,ny,1,1/)),320)
    CALL check(nf90_inq_varid(ncid2,"zeta",varid2),319)
    CALL check(nf90_put_var(ncid2,varid2,a2d),421)

    CALL check(nf90_inq_varid(ncid,"ubar",varid),319)
    CALL check(nf90_get_var(ncid,varid,ubar,start=(/1,1,3,1/),&
        count=(/nx-1,ny,1,1/)),321)
    CALL check(nf90_inq_varid(ncid2,"ubar",varid2),319)
    CALL check(nf90_put_var(ncid2,varid2,ubar),421)

    CALL check(nf90_inq_varid(ncid,"vbar",varid),319)
    CALL check(nf90_get_var(ncid,varid,vbar,start=(/1,1,3,1/),&
        count=(/nx,ny-1,1,1/)),322)
    CALL check(nf90_inq_varid(ncid2,"vbar",varid2),319)
    CALL check(nf90_put_var(ncid2,varid2,vbar),422) 

    CALL check(nf90_inq_varid(ncid,"u",varid),319)
    CALL check(nf90_get_var(ncid,varid,u,start=(/1,1,1,2,1/),&
        count=(/nx-1,ny,nz,1,1/)),323)
    CALL check(nf90_inq_varid(ncid2,"u",varid2),319)
    CALL check(nf90_put_var(ncid2,varid2,u),423)

    CALL check(nf90_inq_varid(ncid,"v",varid),319)
    CALL check(nf90_get_var(ncid,varid,v,start=(/1,1,1,2,1/),&
        count=(/nx,ny-1,nz,1,1/)),324)
    CALL check(nf90_inq_varid(ncid2,"v",varid2),319)
    CALL check(nf90_put_var(ncid2,varid2,v),424)

    CALL check(nf90_inq_varid(ncid,"temp",varid),319)
    CALL check(nf90_get_var(ncid,varid,a3d,start=(/1,1,1,2,1/),&
        count=(/nx,ny,nz,1,1/)),325)
    CALL check(nf90_inq_varid(ncid2,"temp",varid2),319)
    CALL check(nf90_put_var(ncid2,varid2,a3d),425)

    CALL check(nf90_inq_varid(ncid,"salt",varid),319)
    CALL check(nf90_get_var(ncid,varid,a3d,start=(/1,1,1,2,1/),&
        count=(/nx,ny,nz,1,1/)),326)
    CALL check(nf90_inq_varid(ncid2,"salt",varid2),319)
    CALL check(nf90_put_var(ncid2,varid2,a3d),426)

    deallocate( a2d,a3d,u,v,ubar,vbar )
    CALL check(nf90_close(ncid2),360)
    CALL check(nf90_close(ncid),361)

END PROGRAm rst2ini

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
