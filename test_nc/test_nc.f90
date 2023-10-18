PROGRAM test_nc
    USE netcdf
    implicit none
    integer :: ncid, varid
    real(kind=4) :: var(600*640+2)

    character(len=150) :: filename

    filename = '/users/work/mmuzyka/CICE_forcing/cice.grid.nc'

    CALL check(nf90_open(trim(filename), NF90_NOWRITE, ncid),110)
    CALL check(nf90_inq_varid(ncid, 'ulon', varid ),111)
    CALL check(nf90_get_var(ncid, varid, var, &
        start = (/1,1/), count=(/600,640/) ),112)
    CALL check(nf90_close( ncid ), 115 )


    write(*,*) 'var(1)=',var(1)
    write(*,*) 'var(600)=',var(600)
    write(*,*) 'var(601)=',var(601)
    write(*,*) 'var(640)=',var(640)
END PROGRAM test_nc

SUBROUTINE check(status,label)
      USE NETCDF, only: nf90_noerr, nf90_strerror
      implicit none
      INTEGER, INTENT(in) :: status , label

      IF (status .ne. nf90_noerr) then
          print*, trim(nf90_strerror(status)),'label=',label
          stop "Stopped"
      END IF
END SUBROUTINE check
