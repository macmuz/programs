PROGRAM closing
    USE netcdf
    implicit none
    integer :: ncid, varid, x_dimid, y_dimid, dimids2(2)
    character(len=250) :: grid,outgrid
    character(len=20) :: buffer
    integer :: mask(2700,3200), narg, norm_no, diag_no
    integer :: i, j, k, l, n, cnt, wet_cell(2),overall
    logical :: diag

    outgrid = 'ROMC_grid_025NM_lowcoast.nc'

    narg = command_argument_count()
    if (narg.lt.2 .or. narg .gt.3) then
        write(*,*) "Program must have 2 or 3 command line arguments"
        stop
    end if

    call get_command_argument(1,buffer)
    read(buffer, *) norm_no
    call get_command_argument(2,buffer)
    read(buffer, *) diag
    if (diag .and. narg.eq.2) then
        write(*,*) "Program needs 3rd argument, because 2nd one was TRUE"
        stop
    end if
    if (diag) then
      call get_command_argument(3,buffer)
      read(buffer, *) diag_no
    end if

    wet_cell = (1400,500)
    grid = '/users/work/mmuzyka/CSDIR/input_025NM/ROMS_grid_025NM.nc'
    
    CALL check(nf90_open(trim(grid), NF90_NOWRITE, ncid), 110 )
    CALL check(nf90_inq_varid(ncid, 'mask_rho', varid), 111 )
    CALL check(nf90_get_var(ncid, varid, mask), 115 )
    CALL check(nf90_close( ncid ), 118 )


    do 

    cnt = 0
    if (diag) then

    do i = 1,2699-diag_no
    do j = 1,3199-diag_no
      if (mask(i,j).eq.0 .and. mask(i+1,j+1).eq.1) then
        do k = 1,diag_no
          if ( mask(i+1+k,j+1+k).eq.0) then
            cnt = cnt+1
            do l = k,1,-1
              mask(i+l,j+l) = 0
            end do
            exit
          end if
        end do
      end if
    end do
    end do

    do i = 2+diag_no,2700
    do j = 1,3199-diag_no
      if (mask(i,j).eq.0 .and. mask(i-1,j+1).eq.1) then
        do k = 1,diag_no
          if ( mask(i-1-k,j+1+k).eq.0) then
            cnt = cnt+1
            do l = k,1,-1
              mask(i-l,j+l) = 0
            end do
            exit
          end if
        end do
      end if
    end do
    end do

    end if

    do i = 1,2700 
    do j = 1,3199-norm_no
      if (mask(i,j).eq.0 .and. mask(i,j+1).eq.1) then
        do k = 1,norm_no
          if ( mask(i,j+1+k).eq.0) then
            cnt = cnt+1
            do l = k,1,-1
              mask(i,j+l) = 0
            end do 
            exit
          end if
        end do    
      end if 
    end do 
    end do

    do j = 1,3200
    do i = 1,2699-norm_no
      if (mask(i,j).eq.0 .and. mask(i+1,j).eq.1) then
        do k = 1,norm_no
          if ( mask(i+1+k,j).eq.0) then
            cnt = cnt+1
            do l = k,1,-1
              mask(i+l,j) = 0
            end do
            exit
          end if
        end do
      end if
    end do
    end do 

    write(*,*) 'cnt=',cnt
    if (cnt.eq.0) exit
    end do

    overall = 0
    mask(wet_cell(1),wet_cell(2)) = 2
    do
      cnt = 0
      do i = 2,2699
      do j = 2,3199
        if (mask(i,j).eq.2 .and. mask(i+1,j).eq.1) then
           mask(i+1,j) = 2
           cnt = cnt+1
        end if
        if (mask(i,j).eq.2 .and. mask(i-1,j).eq.1) then
           mask(i-1,j) = 2
           cnt = cnt+1
        end if
        if (mask(i,j).eq.2 .and. mask(i,j+1).eq.1) then
           mask(i,j+1) = 2
           cnt = cnt+1
        end if
        if (mask(i,j).eq.2 .and. mask(i,j-1).eq.1) then
           mask(i,j-1) = 2
           cnt = cnt+1
        end if
      end do
      end do

      overall = overall+cnt
      write(*,*) 'cnt=',cnt
      write(*,*) 'overall=',overall
      if (cnt.eq.0) exit
    end do

    where(mask.ne.2) mask=0
    where(mask.eq.2) mask=1

    CALL check(nf90_create( "mask_closed.nc", NF90_CLOBBER, ncid ), 200)
    CALL check(nf90_def_dim( ncid, 'xi_rho', 2700, x_dimid ), 201)
    CALL check(nf90_def_dim( ncid, 'eta_rho', 3200, y_dimid ), 202)
    dimids2 = (/ x_dimid, y_dimid /)
    CALL check(nf90_def_var( ncid, 'mask_rho', NF90_INT, dimids2, varid), 205)
    CALL check( nf90_enddef( ncid ), 207 )
    CALL check( nf90_put_var( ncid, varid, mask ), 209 )
    CALL check(nf90_close( ncid ), 210 )

    CALL check(nf90_open( trim(outgrid), nf90_write, ncid ), 200)
    CALL check(nf90_inq_varid(ncid, "mask_rho", varid),201)
    CALL check(nf90_put_var(ncid, varid,real(mask,8)),202)
    CALL check(nf90_close( ncid ), 210 )
    
END PROGRAM closing

SUBROUTINE check(status,label)
    USE NETCDF
    implicit none
    INTEGER, INTENT(in) :: status , label

    IF (status .ne. nf90_noerr) then
        print*, trim(nf90_strerror(status)),'label=',label
        stop "Stopped"
    END IF
END SUBROUTINE check
