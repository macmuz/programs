PROGRAM calc    
    USE netcdf
    implicit none
    integer :: my_task, size_Of_Cluster, ierror
    integer, parameter :: master_task = 0
    integer :: ncid,ncid_in,varid,dimids2(2),dimids3(3)
    integer :: dim1, dim2, x_dimid, y_dimid, t_dimid, vid(5)
    integer :: i, j, dim1in, dim2in, timeid, record(2), dt
    integer :: varid2, ios
    integer, parameter :: fh=15
    integer(kind=2), allocatable :: idx(:,:,:,:)
    real(kind=4), allocatable :: W(:,:,:)
    real(kind=8), allocatable :: datain(:,:), dataout(:,:)
    real(kind=8), allocatable :: lon(:,:), lat(:,:)
    real(kind=8) :: time(1)
    character(len=250) :: input, output, outgrid, grid, buffer
    character(len=25) :: var_name,time_name,var2_name
    character(len=50) :: time_att(3),var_att(4),lon_att(3),lat_att(3)
    character(len=50) :: var2_att(4),paramfile
    logical :: var2

    INTERFACE
      SUBROUTINE calcme(inarray,outarray,W,idx)
        real(kind=8), dimension(:,:) :: inarray
        real(kind=8), dimension(:,:) :: outarray
        real(kind=4), dimension(:,:,:) :: W
        integer(kind=2), dimension(:,:,:,:) :: idx
      END SUBROUTINE
    END INTERFACE

    call get_command_argument(1,paramfile)
    write(*,*) trim(paramfile)
    open(fh, file=trim(paramfile))

     read(fh, '(A)', iostat=ios) buffer
     if (ios == 0) then
        read(buffer, *, iostat=ios) input
        write(*,*) trim(input)
     end if
     read(fh, '(A)', iostat=ios) buffer
     if (ios == 0) then
        read(buffer, *, iostat=ios) output
        write(*,*) trim(output)
     end if
     read(fh, '(A)', iostat=ios) buffer
     if (ios == 0) then
        read(buffer, *, iostat=ios) time_name
        write(*,*) trim(time_name)
     end if
     read(fh, '(A)', iostat=ios) buffer
     if (ios == 0) then
        read(buffer, *, iostat=ios) var_name
        write(*,*) trim(var_name)
     end if
     read(fh, '(A)', iostat=ios) buffer
     if (ios == 0) then
        read(buffer, *, iostat=ios) record(1)
        write(*,*) record(1)
     end if
     read(fh, '(A)', iostat=ios) buffer
     if (ios == 0) then
        read(buffer, *, iostat=ios) record(2)
        write(*,*) record(2)
     end if
     read(fh, '(A)', iostat=ios) buffer
     if (ios == 0) then
        read(buffer, *, iostat=ios) var2
        write(*,*) var2
     end if
     read(fh, '(A)', iostat=ios) buffer
     if (ios == 0) then
        read(buffer, *, iostat=ios) var2_name
        write(*,*) trim(var2_name)
     end if

    close(fh)

!    input = '/users/work/mmuzyka/CSDIR/forcing/baltic_wind_2016.nc'
!    output = 'baltic_wind_2016_2700x3200.nc'
!    var_name = 'Uwind'
!    time_name = 'wind_time'
!    record = (/2567,3313/)
!    record = (/2567,2570/)

!    var2 = .true.
!    var2_name = 'Vwind'

    outgrid = '/users/work/mmuzyka/programs/interpolation/out_025NM.nc'
    grid = '/users/work/mmuzyka/CSDIR/input_025NM/ROMS_grid_025NM.nc'

    CALL check(nf90_open(trim(outgrid), NF90_NOWRITE, ncid),110)
    CALL check(nf90_inq_varid(ncid, 'W', varid),111)
    CALL check(nf90_inquire_variable(ncid, varid, dimids=dimids3),112)
    CALL check(nf90_inquire_dimension(ncid, dimids3(1), len=dim1),113)
    CALL check(nf90_inquire_dimension(ncid, dimids3(2), len=dim2),114)
    allocate( idx(dim1,dim2,4,2), W(dim1,dim2,4) )
    allocate( dataout(dim1,dim2) )
    allocate( lon(dim1,dim2), lat(dim1,dim2) )
    CALL check(nf90_get_var(ncid, varid, W),115)
    CALL check(nf90_inq_varid(ncid, 'idx', varid),116)
    CALL check(nf90_get_var(ncid, varid, idx),117)
    CALL check(nf90_close( ncid ), 118 )

    CALL check(nf90_open(trim(grid), NF90_NOWRITE, ncid), 110 )
    CALL check(nf90_inq_varid(ncid, 'lon_rho', varid), 111 )
    CALL check(nf90_get_var(ncid, varid, lon), 115 )
    CALL check(nf90_inq_varid(ncid, 'lat_rho', varid), 111 )
    CALL check(nf90_get_var(ncid, varid, lat), 115 )
    CALL check(nf90_close( ncid ), 118 )


    CALL check(nf90_open(trim(input), NF90_NOWRITE, ncid_in),100)

    CALL check(nf90_inq_varid(ncid_in, trim(time_name), timeid),101)
    CALL check(nf90_get_att(ncid_in, timeid, "long_name", time_att(1)),101)
    CALL check(nf90_get_att(ncid_in, timeid, "units", time_att(2)),101)
    CALL check(nf90_get_att(ncid_in, timeid, "calendar", time_att(3)),101)

    CALL check(nf90_inq_varid(ncid_in, "lon", varid),101)
    CALL check(nf90_get_att(ncid_in, varid, "long_name", lon_att(1)),101)
    CALL check(nf90_get_att(ncid_in, varid, "units", lon_att(2)),101)
    CALL check(nf90_get_att(ncid_in, varid, "standard_name", lon_att(3)),101)

    CALL check(nf90_inq_varid(ncid_in, "lat", varid),101)
    CALL check(nf90_get_att(ncid_in, varid, "long_name", lat_att(1)),101)
    CALL check(nf90_get_att(ncid_in, varid, "units", lat_att(2)),101)
    CALL check(nf90_get_att(ncid_in, varid, "standard_name", lat_att(3)),101)

    if (var2) then
    CALL check(nf90_inq_varid(ncid_in, trim(var2_name), varid2),101)
    CALL check(nf90_get_att(ncid_in, varid2, "long_name", var2_att(1)),101)
    CALL check(nf90_get_att(ncid_in, varid2, "units", var2_att(2)),101)
    CALL check(nf90_get_att(ncid_in, varid2, "time", var2_att(3)),101)
    CALL check(nf90_get_att(ncid_in, varid2, "coordinates", var2_att(4)),101)
    end if

    CALL check(nf90_inq_varid(ncid_in, trim(var_name), varid),101)
    CALL check(nf90_get_att(ncid_in, varid, "long_name", var_att(1)),101)
    CALL check(nf90_get_att(ncid_in, varid, "units", var_att(2)),101)
    CALL check(nf90_get_att(ncid_in, varid, "time", var_att(3)),101)
    CALL check(nf90_get_att(ncid_in, varid, "coordinates", var_att(4)),101)
    CALL check(nf90_inquire_variable(ncid_in, varid, dimids=dimids3),102)
    CALL check(nf90_inquire_dimension(ncid_in, dimids3(1), len=dim1in),103)
    CALL check(nf90_inquire_dimension(ncid_in, dimids3(2), len=dim2in),104)
    allocate( datain(dim1in,dim2in) )


    CALL check(nf90_create( trim(output), NF90_CLOBBER, ncid ), 200)
    CALL check(nf90_def_dim( ncid, 'xi_rho', dim1, x_dimid ), 201)
    CALL check(nf90_def_dim( ncid, 'eta_rho', dim2, y_dimid ), 202)
    CALL check(nf90_def_dim( ncid, trim(time_name), NF90_UNLIMITED, t_dimid ), 202)
    dimids2 = (/ x_dimid, y_dimid /)
    dimids3 = (/ x_dimid, y_dimid, t_dimid /)

    CALL check(nf90_def_var( ncid, trim(time_name), NF90_DOUBLE, (/t_dimid/), vid(1)), 205)
    CALL check(nf90_put_att( ncid, vid(1), "long_name", trim(time_att(1))), 205)
    CALL check(nf90_put_att( ncid, vid(1), "units", trim(time_att(2))), 205)
    CALL check(nf90_put_att( ncid, vid(1), "calendar", trim(time_att(3))), 205)

    CALL check(nf90_def_var( ncid, trim(var_name), NF90_FLOAT, dimids3, vid(2)), 205)
    CALL check(nf90_put_att( ncid, vid(2), "long_name", trim(var_att(1))), 205)
    CALL check(nf90_put_att( ncid, vid(2), "units", trim(var_att(2))), 205)
    CALL check(nf90_put_att( ncid, vid(2), "time", trim(var_att(3))), 205)
    CALL check(nf90_put_att( ncid, vid(2), "coordinates", trim(var_att(4))), 205)

    if (var2) then
    CALL check(nf90_def_var( ncid, trim(var2_name), NF90_FLOAT, dimids3, vid(5)), 205)
    CALL check(nf90_put_att( ncid, vid(5), "long_name", trim(var2_att(1))), 205)
    CALL check(nf90_put_att( ncid, vid(5), "units", trim(var2_att(2))), 205)
    CALL check(nf90_put_att( ncid, vid(5), "time", trim(var2_att(3))), 205)
    CALL check(nf90_put_att( ncid, vid(5), "coordinates", trim(var2_att(4))), 205)
    end if

    CALL check(nf90_def_var( ncid, "lon", NF90_FLOAT, dimids2, vid(3)), 205)
    CALL check(nf90_put_att( ncid, vid(3), "long_name", trim(lon_att(1))), 205)
    CALL check(nf90_put_att( ncid, vid(3), "units", trim(lon_att(2))), 205)
    CALL check(nf90_put_att( ncid, vid(3), "standard_name", trim(lon_att(3))), 205)

    CALL check(nf90_def_var( ncid, "lat", NF90_FLOAT, dimids2, vid(4)), 205)
    CALL check(nf90_put_att( ncid, vid(4), "long_name", trim(lat_att(1))), 205)
    CALL check(nf90_put_att( ncid, vid(4), "units", trim(lat_att(2))), 205)
    CALL check(nf90_put_att( ncid, vid(4), "standard_name", trim(lat_att(3))), 205)

    CALL check( nf90_enddef( ncid ), 207 )
    CALL check( nf90_put_var( ncid, vid(3), lon ), 209 )
    CALL check( nf90_put_var( ncid, vid(4), lat ), 209 )


    dt = 1
    do i = record(1),record(2)
    write(*,*) 'i=',i
    CALL check(nf90_get_var(ncid_in, timeid, time, start = (/i/), count = (/1/)),106)
    CALL check(nf90_get_var(ncid_in, varid, datain,&
                start = (/1,1,i/),&
                count = (/dim1in,dim2in,1/) ),107)

    CALL calcme(datain,dataout,W,idx)

    CALL check( nf90_put_var( ncid, vid(1), time, start = (/dt/) ), 209 )
    CALL check( nf90_put_var( ncid, vid(2), dataout, start = (/1,1,dt/) ), 209 )

    if (var2) then
        CALL check(nf90_get_var(ncid_in, varid2, datain,&
                start = (/1,1,i/),&
                count = (/dim1in,dim2in,1/) ),107)
        CALL calcme(datain,dataout,W,idx)
        CALL check( nf90_put_var( ncid, vid(5), dataout,&
                start = (/1,1,dt/) ), 209 )
    end if

    dt = dt+1
    end do
    

    CALL check(nf90_close( ncid ), 210 )
    CALL check(nf90_close( ncid_in ), 108 )
    deallocate(idx,W)
    deallocate(lon,lat)
    deallocate(datain,dataout)
END PROGRAM calc

SUBROUTINE check(status,label)
    USE NETCDF
    implicit none
    INTEGER, INTENT(in) :: status , label

    IF (status .ne. nf90_noerr) then
        print*, trim(nf90_strerror(status)),'label=',label
        stop "Stopped"
    END IF
END SUBROUTINE check

SUBROUTINE calcme(inarray,outarray,W,idx)
    implicit none
    real(kind=8), dimension(:,:), intent(in) :: inarray
    real(kind=8), dimension(:,:), intent(out) :: outarray
    real(kind=4), dimension(:,:,:), intent(in) :: W
    integer(kind=2), dimension(:,:,:,:), intent(in) :: idx

    real(kind=8) :: tmp
    integer :: i,j,k,x,y
   
    do i = 1,size(outarray,1) 
    do j = 1,size(outarray,2)
      tmp = 0
      do k = 1,4
        x = idx(i,j,k,1)
        y = idx(i,j,k,2)
        if (x.gt.0 .and. y.gt.0) then
          tmp = tmp+inarray(x,y)*W(i,j,k)
        end if
      end do 
      outarray(i,j) = tmp
    end do
    end do
END SUBROUTINE calcme

