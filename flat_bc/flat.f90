PROGRAM flat
    USE NETCDF
    implicit none
    integer, parameter :: istr=220,iend=400,jend=1091,dy=30
    integer :: mask(2700,3200), ncid, varid
    real(kind=8) :: h(2700,3200),trx0=0.09,hin=30.0
    character(len=100) :: filein,fileout

    filein = '/users/work/mmuzyka/CSDIR/input_025NM/ROMS_grid_025NM.nc'
    fileout = 'ROMS_grid_025NM_orig_flat.nc'


    CALL check(nf90_open(trim(filein), NF90_NOWRITE, ncid), 110 )
    CALL check(nf90_inq_varid(ncid, 'mask_rho', varid), 111 )
    CALL check(nf90_get_var(ncid, varid, mask), 112 )
    CALL check(nf90_inq_varid(ncid, 'h', varid), 113 )
    CALL check(nf90_get_var(ncid, varid, h), 114 )
    CALL check(nf90_close( ncid ), 115 )

    where( mask(istr:iend,jend+1-dy:jend).eq.1 ) &
        h(istr:iend,jend+1-dy:jend) = hin
    !CALL dc_scheme(2700,3200,mask,h,trx0)

    CALL check(nf90_open( trim(fileout), nf90_write, ncid ), 200)
    CALL check(nf90_inq_varid(ncid, "h", varid),201)
    CALL check(nf90_put_var(ncid, varid, h),202)
    CALL check(nf90_close( ncid ), 203 )
END PROGRAM flat

SUBROUTINE check(status,label)
    USE NETCDF
    implicit none
    INTEGER, INTENT(in) :: status , label

    IF (status .ne. nf90_noerr) then
        print*, trim(nf90_strerror(status)),'label=',label
        stop "Stopped"
    END IF
END SUBROUTINE check

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
