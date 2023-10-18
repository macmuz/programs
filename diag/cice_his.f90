PROGRAM cice_his
    USE netcdf
    implicit none
    integer, parameter :: n=2,npts=8
    integer :: ni,nj,ntime,pts(npts,2)
    integer :: date(3),datestp(3)
    integer :: ncid, varid, dimid, i, j, t
    real, allocatable :: var3d(:,:,:),varpts(:,:)
    character(len=200) :: filename,path
    character(len=30) :: vname(n),outname

    100 format(a,'/iceh.',i4,'-',i2.2,'-',i2.2,'.nc')
    path = '/users/work/mmuzyka/CSDIR/metro_atm05NM_era5cice/run/baltic/cice/rundir/history'
    
    vname(1) = 'hi'
    vname(2) = 'frzmlt'

    pts(1,:) = (/832,1389/)
    pts(2,:) = (/830,1389/)
    pts(3,:) = (/831,1389/)
    pts(4,:) = (/833,1389/)
    pts(5,:) = (/834,1389/)
    pts(6,:) = (/832,1388/)
    pts(7,:) = (/832,1390/)
    pts(8,:) = (/832,1391/)

    date = (/2016,11,8/)
    datestp = (/2016,12,17/)   

    write(filename,100) trim(path),date(1),date(2),date(3) 
    write(*,*) trim(filename)

    CALL check(nf90_open(trim(filename),NF90_NOWRITE,ncid),310)

    CALL check(nf90_inq_dimid(ncid, "ni", dimid),311)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=ni),312)

    CALL check(nf90_inq_dimid(ncid, "nj", dimid),313)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nj),314)

    CALL check(nf90_inq_dimid(ncid, "time", dimid),315)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=ntime),316)

    CALL check(nf90_close(ncid),360)

    write(*,*) ni,nj,ntime
    allocate( var3d(ni,nj,ntime), varpts(npts,ntime) )

    do i = 1,n
      write(outname,"(a,'.txt')") trim(vname(i)) 
      open (unit = 500+i, file = trim(outname), form="FORMATTED")
    end do

    open (unit = 400, file = "dates.txt", form="FORMATTED")

    do

        write(filename,100) trim(path),date(1),date(2),date(3)
        CALL check(nf90_open(trim(filename),NF90_NOWRITE,ncid),310)

        do i = 1,n
          CALL check(nf90_inq_varid(ncid,trim(vname(i)),varid),319)
          CALL check(nf90_get_var(ncid,varid,var3d),320)
          do j = 1,npts
            do t = 1, ntime
              varpts(j,t) = var3d(pts(j,1),pts(j,2),t)
            end do 
          end do

          do t = 1, ntime
            write(500+i,"(*(f,a))") (varpts(j,t),achar(9),j=1,npts)
          end do
        end do

        CALL check(nf90_close(ncid),360)
   
        write(400,"(3(i,a))") (date(j),achar(9),j=1,3)
 
        if (date(1).ge.datestp(1).and.&
            date(2).ge.datestp(2).and.&
            date(3).ge.datestp(3)) EXIT
        CALL add_day(date,.true.)
    end do

    do i = 1,n
      close (500+i)
    end do
    close(400)

    deallocate( var3d, varpts )
        
END PROGRAM cice_his

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
