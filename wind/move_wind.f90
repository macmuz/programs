PROGRAM move_wind
    USE NETCDF
    implicit none
    integer :: ncid,varid,ncid2,varid2,i
    real :: tmp(560,600)
    character(len=250) :: input,output    

    input = "/users/work/mmuzyka/CSDIR/forcing_560x600/baltic_wind_1998.nc"
    output = "/users/work/mmuzyka/CSDIR/forcing_560x600/baltic_wind_1999.nc"

    CALL check(nf90_open(trim(input),NF90_NOWRITE,ncid),310)
    CALL check(nf90_open(trim(output),NF90_WRITE,ncid2),410)

    do i = 1,8760

    CALL check(nf90_inq_varid(ncid,"Uwind",varid),319)
    CALL check(nf90_get_var(ncid,varid,tmp,&
         start=(/1,1,i/) ),320)

    CALL check(nf90_inq_varid(ncid2,"Uwind",varid2),419)
    CALL check(nf90_put_var(ncid2,varid2,tmp,&
        start=(/1,1,i/)),420)
    
    CALL check(nf90_inq_varid(ncid,"Vwind",varid),319)
    CALL check(nf90_get_var(ncid,varid,tmp,&
         start=(/1,1,i/) ),320)

    CALL check(nf90_inq_varid(ncid2,"Vwind",varid2),419)
    CALL check(nf90_put_var(ncid2,varid2,tmp,&
        start=(/1,1,i/)),420)
    
    write(*,*) i
    end do

    CALL check(nf90_close(ncid),360)
    CALL check(nf90_close(ncid2),460)

END PROGRAM move_wind

SUBROUTINE check(status,label)
    USE NETCDF
    implicit none
    INTEGER, INTENT(in) :: status , label

    IF (status .ne. nf90_noerr) then
        print*, trim(nf90_strerror(status)),' label= ',label
        stop "Stopped"
    END IF
END SUBROUTINE check
