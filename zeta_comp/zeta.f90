PROGRAM zeta
    USE NETCDF
    implicit none
    integer, parameter :: istr=124,iend=149,jend=236
    integer :: date(3),date_start(3),date_end(3)
    integer :: dx,dt,i,j,k,ncid,varid,cnt
    character(len=150) :: bcpath,modelpath,filename
    character(len=50) :: outfile
    logical :: isleap
    real(kind=4), allocatable :: bc(:,:),bcave(:),model(:,:,:)

    101 format(a,'/baltic_bry_nemo_', i4.4, '.nc')
    102 format('zeta', i3.3, '.txt')
    103 format(a,'/ocean_avg_', i4.4, '-', i2.2, '-', i2.2, '.nc')
    1200 format(i4.4,3(a,i2.2),*(a,f13.6))
    bcpath = '/users/work/mmuzyka/CSDIR/metro_files'
    modelpath = '/users/work/mmuzyka/CSDIR/new_metro_test01/run/baltic'

    date_start = (/2016,7,1/)
    date_end = (/2017,6,30/)

!    write(*,*) '2016:',isleap(2016)
!    write(*,*) '2017:',isleap(2017)
!    write(*,*) '2000:',isleap(2000)
!    write(*,*) '1900:',isleap(1900)

    dx = iend-istr+1
    dt = 0    
    do i = date_start(1),date_end(1)
        if (isleap(i)) then
            dt = dt + 366*24 
        else
            dt = dt + 365*24 
        end if 
    end do
    write(*,*) dx,dt
    allocate(bc(dx,dt),bcave(dx),model(dx,5,4))

    cnt = 0
    do i = date_start(1),date_end(1)
        write(filename,101) trim(bcpath),i
        if (isleap(i)) then
            dt = 366*24 
        else
            dt = 365*24 
        end if 
        CALL check(nf90_open(trim(filename),NF90_NOWRITE,ncid),100)
        CALL check(nf90_inq_varid(ncid, "zeta_north", varid),101)
        CALL check(nf90_get_var(ncid, varid, bc(:,cnt+1:cnt+dt), &
                    start = (/ istr, 1 /), &
                    count = (/ dx, dt /)), 102)
        CALL check(nf90_close(ncid),105)
        cnt = cnt+dt
    end do

    date(1) = date_start(1)
    date(2) = 1
    date(3) = 1
    cnt = 1
    do
        if( date(1).eq.date_start(1) .and.&
            date(2).eq.date_start(2) .and.&
            date(3).eq.date_start(3) ) EXIT
        cnt = cnt+24
        CALL add_day(date,.true.)
    end do
 
    open( 200, file='zetabc.txt', status='replace' )
    do i = 1,5
        write(outfile,102) jend+1-i
        open( 200+i, file=trim(outfile), status='replace' )
    end do

    date(1) = date_start(1)
    date(2) = date_start(2)
    date(3) = date_start(3)
    do
        write(*,*) date(1), date(2), date(3)
        write(*,*) cnt
        write(filename,103) trim(modelpath),date(1),date(2),date(3)
        write(*,*) trim(filename)
        CALL check(nf90_open(trim(filename),NF90_NOWRITE,ncid),200)
        CALL check(nf90_inq_varid(ncid, "zeta", varid),201)
        CALL check(nf90_get_var(ncid, varid, model, &
                    start = (/ istr, jend-4, 1 /), &
                    count = (/ dx, 5, 4 /)), 202)
        write(*,*) model(1,5,4)
        CALL check(nf90_close(ncid),105)
        where(model>100) model=-999.0
        bcave = 0
        do i = 1, 4
            bcave = sum(bc(:,cnt:cnt+6),2)/7
            cnt = cnt+6
            write(200,1200) date(1), ( achar(9), date(k), k=2,3 ),&
                            achar(9), i*6-3, ( achar(9), bcave(k), k=1,dx )
            do j = 1,5
            write(200+j,1200) date(1), ( achar(9), date(k), k=2,3 ),&
                            achar(9), i*6-3, ( achar(9), model(k,6-j,i), k=1,dx )
            end do
        end do

        if( date(1).eq.date_end(1) .and.&
            date(2).eq.date_end(2) .and.&
            date(3).eq.date_end(3) ) EXIT
        CALL add_day(date,.true.)
    end do  
    deallocate(bc,bcave,model)

    close(200)
    do i = 1,5
        close(200+i)
    end do

END PROGRAM zeta

SUBROUTINE check(status,label)
    USE NETCDF
    implicit none
    INTEGER, INTENT(in) :: status , label

    IF (status .ne. nf90_noerr) then
        print*, trim(nf90_strerror(status)),'label=',label
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

LOGICAL FUNCTION isleap(year)
    implicit none
    integer, intent(in) :: year
    
    isleap = .false.
    if ( mod(year,4).eq.0 ) isleap = .true.
    if ( mod(year,100).eq.0 ) isleap = .false.
    if ( mod(year,400).eq.0 ) isleap = .true.
END FUNCTION
