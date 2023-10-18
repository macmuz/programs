PROGRAM new_flat
    USE NETCDF
    implicit none
    integer :: nx,ny,ncid,dimid,varid
    integer :: ip(2),jp(2)
    real(kind=8) :: tr,rx0,pt
    real(kind=8), allocatable :: h(:,:)
    integer, allocatable :: mask(:,:)
    character(len=100) :: filein,fileout
    character(len=200) :: command

    filein = 'ROMS_grid_05NM_015.nc'
    fileout = 'ROMS_grid_05NM_015_flat.nc'
    tr = 0.152
    pt = 45.0
    ip = (/100,230/)
    jp = (/400,500/)

    write(command,'("cp ",a," ",a )'),trim(filein),trim(fileout)
    write(*,*) trim(command)
    CALL EXECUTE_COMMAND_LINE(trim(command))

    CALL check(nf90_open(trim(filein), NF90_NOWRITE, ncid), 110 )
    CALL check(nf90_inq_dimid(ncid, "xi_rho", dimid), 200)
    CALL check(nf90_inquire_dimension(ncid, dimid, len = nx), 201)
    CALL check(nf90_inq_dimid(ncid, "eta_rho", dimid), 202)
    CALL check(nf90_inquire_dimension(ncid, dimid, len = ny), 203)
    allocate( mask(nx,ny), h(nx,ny) )
    CALL check(nf90_inq_varid(ncid, 'mask_rho', varid), 111 )
    CALL check(nf90_get_var(ncid, varid, mask), 112 )
    CALL check(nf90_inq_varid(ncid, 'h', varid), 113 )
    CALL check(nf90_get_var(ncid, varid, h), 114 )
    CALL check(nf90_close( ncid ), 115 )

    CALL calc_rx0(nx,ny,mask,h,rx0)
    write(*,*) maxval(h),minval(h)
    write(*,*) 'rx0= ',rx0
!    h(ip(1):ip(2),jp(1):jp(2)) = pt
    where(h(ip(1):ip(2),jp(1):jp(2)).gt.pt) &
        h(ip(1):ip(2),jp(1):jp(2))=pt
    CALL dc_scheme(nx,ny,mask,h,tr) 

    CALL check(nf90_open(trim(fileout), NF90_WRITE, ncid), 301 )
    CALL check(nf90_inq_varid(ncid, 'h', varid), 302 )
    CALL check(nf90_put_var(ncid, varid, h), 303 )
    CALL check(nf90_close( ncid ), 304 )
    
 
    deallocate( mask, h )


END PROGRAM new_flat

SUBROUTINE check(status,label)
    USE NETCDF
    implicit none
    INTEGER, INTENT(in) :: status , label

    IF (status .ne. nf90_noerr) then
        print*, trim(nf90_strerror(status)),'label=',label
        stop "Stopped"
    END IF
END SUBROUTINE check

SUBROUTINE calc_rx0(nx,ny,mask,h,rx0max)
    implicit none
    integer, intent(in) :: nx,ny,mask(nx,ny)
    real(kind=8), intent(in) :: h(nx,ny)
    real(kind=8), intent(out) :: rx0max

    integer :: i,j
    real(kind=8) :: rx0(nx,ny),rx0_tmp

    rx0 = 0

    do i = 2,nx-1
    do j = 2,ny-1
      if (mask(i,j).eq.1) then
        if(mask(i+1,j).eq.1) then
          rx0_tmp = abs(h(i,j)-h(i+1,j))/(h(i,j)+h(i+1,j))
          rx0(i,j) = max(rx0_tmp,rx0(i,j))
        end if
        if(mask(i-1,j).eq.1) then
          rx0_tmp = abs(h(i,j)-h(i-1,j))/(h(i,j)+h(i-1,j))
          rx0(i,j) = max(rx0_tmp,rx0(i,j))
        end if
        if(mask(i,j+1).eq.1) then
          rx0_tmp = abs(h(i,j)-h(i,j+1))/(h(i,j)+h(i,j+1))
          rx0(i,j) = max(rx0_tmp,rx0(i,j))
        end if
        if(mask(i,j-1).eq.1) then
          rx0_tmp = abs(h(i,j)-h(i,j-1))/(h(i,j)+h(i,j-1))
          rx0(i,j) = max(rx0_tmp,rx0(i,j))
        end if
      end if
    end do
    end do

    rx0max = maxval(rx0)
END SUBROUTINE calc_rx0

SUBROUTINE dc_scheme(nx,ny,mask,h,target_rx0)
    implicit none
    integer, intent(in) :: nx,ny,mask(nx,ny)
    real(kind=8), intent(inout) :: h(nx,ny)
    real(kind=8), intent(in) :: target_rx0

    integer :: i,j
    real(kind=8) :: rx0

    do
      CALL calc_rx0(nx,ny,mask,h,rx0)
      write(*,*) 'rx0=',rx0
      if ( rx0.le.(target_rx0+0.001) ) EXIT

      do j = 2,nx-1
      do i = 2,ny-1
      if (mask(i,j).eq.1) then
        if(mask(i+1,j).eq.1) then
          if ( (h(i,j)-h(i+1,j))/(h(i,j)+h(i+1,j)).gt.target_rx0 ) then
            h(i,j) = h(i+1,j)*(1+target_rx0)/(1-target_rx0)
          end if
        end if
        if(mask(i-1,j).eq.1) then
          if ( (h(i,j)-h(i-1,j))/(h(i,j)+h(i-1,j)).gt.target_rx0 ) then
            h(i,j) = h(i-1,j)*(1+target_rx0)/(1-target_rx0)
          end if
        end if
        if(mask(i,j+1).eq.1) then
          if ( (h(i,j)-h(i,j+1))/(h(i,j)+h(i,j+1)).gt.target_rx0 ) then
            h(i,j) = h(i,j+1)*(1+target_rx0)/(1-target_rx0)
          end if
        end if
        if(mask(i,j-1).eq.1) then
          if ( (h(i,j)-h(i,j-1))/(h(i,j)+h(i,j-1)).gt.target_rx0 ) then
            h(i,j) = h(i,j-1)*(1+target_rx0)/(1-target_rx0)
          end if
        end if
      end if
      end do
      end do
    end do
END SUBROUTINE dc_scheme
