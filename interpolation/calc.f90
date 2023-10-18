PROGRAM calc    
    USE netcdf
    implicit none
    include 'mpif.h'
    integer :: my_task, size_Of_Cluster, ierror
    integer, parameter :: master_task = 0
    integer :: ncid,varid,dimids2(2),dimids3(3)
    integer :: dim1, dim2, x_dimid, y_dimid, vid(3)
    integer :: i, j
    integer(kind=2), allocatable :: idx(:,:,:,:)
    real(kind=4), allocatable :: W(:,:,:)
    real(kind=8), allocatable :: inTair(:,:), outTair(:,:)
    real(kind=8), allocatable :: lon(:,:), lat(:,:)
    character(len=250) :: ingrid, outgrid, grid

    INTERFACE
      SUBROUTINE calcme(inarray,outarray,W,idx)
        real(kind=8), dimension(:,:) :: inarray
        real(kind=8), dimension(:,:) :: outarray
        real(kind=4), dimension(:,:,:) :: W
        integer(kind=2), dimension(:,:,:,:) :: idx
      END SUBROUTINE
    END INTERFACE

    call MPI_INIT(ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,size_Of_Cluster,ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD,my_task,ierror)


    ingrid = '/users/work/mmuzyka/CSDIR/forcing_560x600/baltic_Tair_1999.nc'
    outgrid = 'out_05NM.nc'
    grid = 'ROMS_grid_05NM_015.nc'

    CALL check(nf90_open(trim(ingrid), NF90_NOWRITE, ncid),100)
    CALL check(nf90_inq_varid(ncid, 'Tair', varid),101)
    CALL check(nf90_inquire_variable(ncid, varid, dimids=dimids3),102)
    CALL check(nf90_inquire_dimension(ncid, dimids3(1), len=dim1),103)
    CALL check(nf90_inquire_dimension(ncid, dimids3(2), len=dim2),104)
    allocate( inTair(dim1,dim2) )
    CALL check(nf90_get_var(ncid, varid, inTair,&
                start = (/1,1,1/),&
                count = (/dim1,dim2,1/) ),107)
    CALL check(nf90_close( ncid ), 108 )


    CALL check(nf90_open(trim(outgrid), NF90_NOWRITE, ncid),110)
    CALL check(nf90_inq_varid(ncid, 'W', varid),111)
    CALL check(nf90_inquire_variable(ncid, varid, dimids=dimids3),112)
    CALL check(nf90_inquire_dimension(ncid, dimids3(1), len=dim1),113)
    CALL check(nf90_inquire_dimension(ncid, dimids3(2), len=dim2),114)
    allocate( idx(dim1,dim2,4,2), W(dim1,dim2,4) )
    allocate( outTair(dim1,dim2) )
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

    CALL calcme(inTair,outTair,W,idx)

    if (my_task==master_task) then
      CALL check(nf90_create( 'after_interp_05NM.nc', NF90_CLOBBER, ncid ), 200)
      CALL check(nf90_def_dim( ncid, 'nx', dim1, x_dimid ), 201)
      CALL check(nf90_def_dim( ncid, 'ny', dim2, y_dimid ), 202)
      dimids2 = (/ x_dimid, y_dimid /)
      CALL check(nf90_def_var( ncid, 'Tair', NF90_DOUBLE, dimids2, vid(1)), 205)
      CALL check(nf90_def_var( ncid, 'lon', NF90_DOUBLE, dimids2, vid(2)), 205)
      CALL check(nf90_def_var( ncid, 'lat', NF90_DOUBLE, dimids2, vid(3)), 205)
      CALL check( nf90_enddef( ncid ), 207 )
      CALL check( nf90_put_var( ncid, vid(1), outTair ), 209 )
      CALL check( nf90_put_var( ncid, vid(2), lon ), 209 )
      CALL check( nf90_put_var( ncid, vid(3), lat ), 209 )
      CALL check(nf90_close( ncid ), 210 )
    end if
    

    deallocate(idx,W)
    deallocate(lon,lat)
    deallocate(inTair,outTair)
    call MPI_FINALIZE(ierror)
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

