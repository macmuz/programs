PROGRAM hfilterbc
    USE NETCDF
    implicit none
    integer :: idx, ncid, dimid, vid, nx, ny 
    integer :: istr,iend,jstr,jend,dx,dy
    real(kind=8) :: trx0,hmin
    real(kind=8), allocatable :: h(:,:),mask(:,:)
    real(kind=8), allocatable :: rx0(:,:)
    character(len=200) :: filein,fileout

    filein = 'ROMS_grid_05NM_wider.nc'
    istr = 100
    iend = 260
    jstr = 400
    jend = 500
    trx0 = 0.05

    hmin = 5.0
    
    dx = 1+iend-istr
    dy = 1+jend-jstr

    idx = index(filein,'.nc')
    write(fileout,'(a)') filein(1:idx-1)//'_filter.nc'
!    write(*,*) trim(fileout)

    CALL system('cp '//trim(filein)//' '//trim(fileout))

    CALL check(nf90_open(trim(filein),NF90_NOWRITE,ncid),10)
    CALL check(nf90_inq_dimid(ncid, "xi_rho", dimid),11)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nx),12)
    CALL check(nf90_inq_dimid(ncid, "eta_rho", dimid),13)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=ny),12)
    allocate( h(nx,ny), mask(nx,ny), rx0(nx,ny) )
    CALL check(nf90_inq_varid(ncid,"h",vid),17)
    CALL check(nf90_get_var(ncid,vid,h),18)
    CALL check(nf90_inq_varid(ncid,"mask_rho",vid),17)
    CALL check(nf90_get_var(ncid,vid,mask),18)
    CALL check(nf90_close(ncid),35)

    write(*,*) trim(fileout)

    where(h.lt.hmin) h = hmin

    CALL calc_rx0(nx,ny,istr,iend,jstr,jend,mask,h,rx0)
    CALL laplacian(nx,ny,istr,iend,jstr,jend,mask,h,trx0)
    CALL calc_rx0(nx,ny,istr,iend,jstr,jend,mask,h,rx0)

    write(*,*) maxval(rx0)

    CALL check(nf90_open(trim(fileout),NF90_WRITE,ncid),10)
    CALL check(nf90_inq_varid(ncid,"h",vid),17)
    CALL check(nf90_put_var(ncid,vid,h),18)
    CALL check(nf90_close(ncid),35)

    deallocate( h, mask, rx0 )

END PROGRAM hfilterbc

SUBROUTINE check(status,label)
    USE NETCDF
    implicit none
    INTEGER, INTENT(in) :: status , label

    IF (status .ne. nf90_noerr) then
        print*, trim(nf90_strerror(status)),' label=',label
        stop "Stopped"
    END IF
END SUBROUTINE check

SUBROUTINE calc_rx0(nx,ny,istr,iend,jstr,jend,mask,h,rx0)
    implicit none
    integer, intent(in) :: nx,ny,istr,iend,jstr,jend
    real(kind=8), intent(in) :: h(nx,ny), mask(nx,ny)
    real(kind=8), intent(out) :: rx0(nx,ny)

    integer :: i,j
    real(kind=8) :: rx0_tmp

    rx0 = 0


    do i = istr,iend
    do j = jstr,jend
      if (mask(i,j).gt.0.5) then
        if(mask(i+1,j).gt.0.5) then
          rx0_tmp = abs(h(i,j)-h(i+1,j))/(h(i,j)+h(i+1,j))
          rx0(i,j) = max(rx0_tmp,rx0(i,j))
        end if
        if(mask(i-1,j).gt.0.5) then
          rx0_tmp = abs(h(i,j)-h(i-1,j))/(h(i,j)+h(i-1,j))
          rx0(i,j) = max(rx0_tmp,rx0(i,j))
        end if
        if(mask(i,j+1).gt.0.5) then
          rx0_tmp = abs(h(i,j)-h(i,j+1))/(h(i,j)+h(i,j+1))
          rx0(i,j) = max(rx0_tmp,rx0(i,j))
        end if
        if(mask(i,j-1).gt.0.5) then
          rx0_tmp = abs(h(i,j)-h(i,j-1))/(h(i,j)+h(i,j-1))
          rx0(i,j) = max(rx0_tmp,rx0(i,j))
        end if
      end if
    end do
    end do
END SUBROUTINE calc_rx0

SUBROUTINE dc(nx,ny,mask,h,trx0)
    integer, intent(in) :: nx,ny
    real(kind=8), intent(in) :: mask(nx,ny), trx0
    real(kind=8), intent(inout) :: h(nx,ny)

    real(kind=8) :: r,h0,h1
    integer :: i,j,n
    do n = 1,1000
    do i = 2,nx-1
    do j = 2,ny-1
      if (mask(i,j).gt.0.5) then
        h0 = h(i,j)

        if (mask(i-1,j).gt.0.5) then
          h1 = h(i-1,j)
          r = (h0-h1)/(h0+h1)
          if ( r.gt.trx0) then
            h(i,j) = h1*(1+trx0)/(1-trx0)
          end if
        end if

        if (mask(i+1,j).gt.0.5) then
          h1 = h(i+1,j)
          r = (h0-h1)/(h0+h1)
          if ( r.gt.trx0) then
            h(i,j) = h1*(1+trx0)/(1-trx0)
          end if
        end if

        if (mask(i,j-1).gt.0.5) then
          h1 = h(i,j-1)
          r = (h0-h1)/(h0+h1)
          if ( r.gt.trx0) then
            h(i,j) = h1*(1+trx0)/(1-trx0)
          end if
        end if

        if (mask(i,j+1).gt.0.5) then
          h1 = h(i,j+1)
          r = (h0-h1)/(h0+h1)
          if ( r.gt.trx0) then
            h(i,j) = h1*(1+trx0)/(1-trx0)
          end if
        end if

      end if
    end do
    end do
    end do
    
     
END SUBROUTINE dc

SUBROUTINE laplacian(nx,ny,istr,iend,jstr,jend,mask,h,trx0)
    USE NETCDF
    implicit none
    integer, intent(in) :: nx,ny,istr,iend,jstr,jend
    real(kind=8), intent(in) :: mask(nx,ny), trx0
    real(kind=8), intent(inout) :: h(nx,ny)
    
    real(kind=8) :: rx0(nx,ny),rmax,dh,rmax0
    integer :: i,j,n,cnt,varid,xi_d,eta_d,vid,ncid

    CALL calc_rx0(nx,ny,istr,iend,jstr,jend,mask,h,rx0) 
    rmax = maxval(rx0)
    rmax0 = rmax

    do n = 1,2000
        write(*,*) rmax

        do i = istr,iend
        do j = jstr,jend
            if (rx0(i,j).gt.trx0) then
              cnt = 0
              dh = 0.0
              if (mask(i+1,j).gt.0.5) then
                dh = dh+(h(i+1,j)-h(i,j))
                cnt = cnt+1
              end if
              if (mask(i-1,j).gt.0.5) then
                dh = dh+(h(i-1,j)-h(i,j))
                cnt = cnt+1
              end if
              if (mask(i,j+1).gt.0.5) then
                dh = dh+(h(i,j+1)-h(i,j))
                cnt = cnt+1
              end if
              if (mask(i,j-1).gt.0.5) then
                dh = dh+(h(i,j-1)-h(i,j))
                cnt = cnt+1
              end if
              h(i,j) = h(i,j)+0.5*dh/real(cnt,8)
              CALL calc_rx0(nx,ny,istr,iend,jstr,jend,mask,h,rx0)
            end if 
        end do
        end do        

        CALL calc_rx0(nx,ny,istr,iend,jstr,jend,mask,h,rx0) 
        rmax = maxval(rx0)
        if (rmax.le.trx0) EXIT
    end do

    CALL dc(nx,ny,mask,h,rmax0)

    CALL calc_rx0(nx,ny,2,nx-1,2,ny-1,mask,h,rx0)

    CALL check( nf90_create( 'diagnostic.nc',NF90_NETCDF4,ncid ), 300 )
    CALL check( nf90_def_dim( ncid, 'xi_rho', nx, xi_d ), 301 )
    CALL check( nf90_def_dim( ncid, 'eta_rho', ny, eta_d ), 301 )
    CALL check( nf90_def_var( ncid, 'rx0', NF90_DOUBLE, (/xi_d,eta_d/), vid ), 316 )
    CALL check( nf90_enddef(ncid), 340 )
    CALL check( nf90_put_var( ncid, vid, rx0 ), 342 )
    CALL check( nf90_close(ncid), 350 )

END SUBROUTINE laplacian
