PROGRAM down2net
    USE NETCDF
    implicit none
    integer :: i, ncid,dimid,nx,ny,nt,vid,varid(2)
    integer :: xdimid,ydimid,tdimid,dimids3(3)
    real(kind=4) :: albedo
    real(kind=4), allocatable :: tmp(:,:,:)
    real(kind=8), allocatable :: time(:)
    character(50) :: path
    character(100) :: fname,oname

    100 format(a,'/baltic_swrad_down_',i4,'.nc')
    101 format(a,'/baltic_swrad_',i4,'.nc')
    path = "/users/magazyn/mmuzyka/UERRA"

    albedo = 0.07

    do i = 1991,2019
      write(fname,100) trim(path),i
      write(oname,101) trim(path),i
      write(*,*) trim(fname)
      write(*,*) trim(oname)

      CALL check(nf90_open(trim(fname),NF90_NOWRITE,ncid),10)
      CALL check(nf90_inq_dimid(ncid, "xi_rho", dimid),11)
      CALL check(nf90_inquire_dimension(ncid, dimid, len=nx),12)
      CALL check(nf90_inq_dimid(ncid, "eta_rho", dimid),13)
      CALL check(nf90_inquire_dimension(ncid, dimid, len=ny),14)
      CALL check(nf90_inq_dimid(ncid, "srf_time", dimid),15)
      CALL check(nf90_inquire_dimension(ncid, dimid, len=nt),16)

      allocate( time(nt), tmp(nx,ny,nt) )

      CALL check(nf90_inq_varid(ncid,"srf_time",vid),33)
      CALL check(nf90_get_var(ncid,vid,time),34) 
      CALL check(nf90_inq_varid(ncid,"swrad_down",vid),33)
      CALL check(nf90_get_var(ncid,vid,tmp),34) 
 
      CALL check(nf90_close(ncid),35)

      tmp = tmp*(1-albedo)

      CALL check(nf90_create( trim(oname), NF90_NETCDF4, ncid ), 200)
      CALL check(nf90_def_dim( ncid, 'xi_rho', nx, xdimid ), 201)
      CALL check(nf90_def_dim( ncid, 'eta_rho', ny, ydimid ), 202)
      CALL check(nf90_def_dim( ncid, 'srf_time', NF90_UNLIMITED, tdimid ), 203)

      dimids3 = (/ xdimid, ydimid, tdimid /)

      CALL check(nf90_def_var( ncid, 'srf_time', NF90_DOUBLE, (/tdimid/), varid(1)), 204)
      CALL check( nf90_put_att( ncid, varid(1), 'long_name',&
        'time for heat flux' ), 205 )
      CALL check( nf90_put_att( ncid, varid(1), 'units',&
        'days since 1968-05-23 00:00:00 GMT' ), 206 )
      CALL check( nf90_put_att( ncid, varid(1), 'calendar',&
        'gregorian' ), 207 )

      CALL check(nf90_def_var( ncid, 'swrad', NF90_FLOAT, dimids3, varid(2)), 208)
      CALL check( nf90_put_att( ncid, varid(2), 'long_name',&
        'solar shortwave radiation' ), 209 )
      CALL check( nf90_put_att( ncid, varid(2), 'units',&
        'Watts meter-2' ), 210 )
      CALL check( nf90_put_att( ncid, varid(2), 'time',&
        'srf_time' ), 211 )

      CALL check(nf90_enddef( ncid ), 212)
        
      CALL check( nf90_put_var( ncid, varid(1), time ), 213 )
      CALL check( nf90_put_var( ncid, varid(2), tmp ), 214 )

      CALL check(nf90_close(ncid),215)
      

      deallocate( time, tmp )
    end do
    
END PROGRAM down2net

SUBROUTINE check(status,label)
    USE NETCDF
    implicit none
    INTEGER, INTENT(in) :: status , label

    IF (status .ne. nf90_noerr) then
        print*, trim(nf90_strerror(status)),' label= ',label
        stop "Stopped"
    END IF
END SUBROUTINE check
