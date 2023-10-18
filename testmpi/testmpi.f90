PROGRAM testmpi
    USE netcdf
    implicit none
    include 'mpif.h'
    integer :: my_task, size_Of_Cluster, ierror
    integer, parameter :: master_task = 0
    integer :: i,j
    integer :: ncid,varid
    integer :: x_dimid,y_dimid,dimids2(2)
    real(kind=8) :: mydata(100,200), array100(100),array2(2)

    call MPI_INIT(ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,size_Of_Cluster,ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD,my_task,ierror)

    array100 = 0

    do j = 1,200
      if (mod(j,size_Of_Cluster).eq.my_task) then
        write(*,*) "Process ", my_task, " is processing j=:",j
          mydata(:,j) = j
      end if
    end do

    call MPI_BARRIER(MPI_COMM_WORLD,ierror)

    if (my_task==master_task) then

      do j = 1,200
        if (mod(j,size_Of_Cluster).ne.master_task) then
          call MPI_RECV(mydata(:,j), 100, MPI_DOUBLE_PRECISION, &
             mod(j,size_Of_Cluster), j, MPI_COMM_WORLD, &
             MPI_STATUS_IGNORE, ierror) 
        end if
      end do

      CALL check(nf90_create( 'out.nc', NF90_CLOBBER, ncid ), 200)
      CALL check(nf90_def_dim( ncid, 'xi_rho', 100, x_dimid ), 201)
      CALL check(nf90_def_dim( ncid, 'eta_rho', 200, y_dimid ), 202)
      dimids2 = (/ x_dimid, y_dimid /)
      CALL check(nf90_def_var( ncid, 'data', NF90_DOUBLE, &
                            dimids2, varid), 203)
      CALL check( nf90_enddef( ncid ), 205 )
      CALL check( nf90_put_var( ncid, varid, mydata ), 206)
      CALL check(nf90_close( ncid ), 207)
    else
      do j = 1,200
        if (mod(j,size_Of_Cluster).eq.my_task) then
          call MPI_SEND(mydata(:,j), 100, MPI_DOUBLE_PRECISION, &
                master_task, j, MPI_COMM_WORLD, ierror)
        end if
      end do
    end if

    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  
    call MPI_Scatter(array100, 2, MPI_DOUBLE_PRECISION, array2, 2, &
        MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror) 
    
    array2 = my_task

    call MPI_Gather(array2, 2, MPI_DOUBLE_PRECISION, &
     &                array100, 2, MPI_DOUBLE_PRECISION, 0, &
     &                MPI_COMM_WORLD, ierror)

    if (my_task==master_task) then
        do i = 1,100
            write(*,*) array100(i)
        end do
    end if

    call MPI_FINALIZE(ierror)
END PROGRAM testmpi

SUBROUTINE check(status,label)
    USE NETCDF
    implicit none
    INTEGER, INTENT(in) :: status , label

    IF (status .ne. nf90_noerr) then
        print*, trim(nf90_strerror(status)),'label=',label
        stop "Stopped"
    END IF
END SUBROUTINE check
