PROGRAM ROMS_bathy
    USE NETCDF
    implicit none
    integer, parameter :: isrc=250
    integer :: nx,ny,ncid,dimid,vid,overall
    integer :: i,j,k,l,cnt,istr,iend,jstr
    integer :: x_dimid,y_dimid,dimids2(2),varid(2)
    integer, allocatable :: mask(:,:)
    real(kind=8), allocatable :: h(:,:),rx0(:,:),pm(:,:),pn(:,:),A(:,:)
    real(kind=8) :: hmin,trx0
    character(len=250) :: gfile, infile
    logical, allocatable :: mk(:,:)

    hmin = 1.0
    trx0 = 0.08
    mk = .false.
    infile = '../ROMS_grids/ROMS_grid_025NMv2.nc'
    gfile = 'ROMS_grid_025NMv2.nc'

    CALL check(nf90_open(trim(infile),NF90_NOWRITE,ncid),1)

    CALL check(nf90_inq_dimid(ncid, "xi_rho", dimid),2)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nx),2)
    CALL check(nf90_inq_dimid(ncid, "eta_rho", dimid),3)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=ny),4)

    allocate( h(nx,ny), mask(nx,ny), mk(nx,ny), rx0(nx,ny) )
    allocate( pm(nx,ny), pn(nx,ny), A(nx,ny) )

    CALL check(nf90_inq_varid(ncid,"h",vid),5)
    CALL check(nf90_get_var(ncid,vid,h),6)
    CALL check(nf90_inq_varid(ncid,"mask_rho",vid),7)
    CALL check(nf90_get_var(ncid,vid,mask),8)
    CALL check(nf90_inq_varid(ncid,"pm",vid),7)
    CALL check(nf90_get_var(ncid,vid,pm),8)
    CALL check(nf90_inq_varid(ncid,"pn",vid),7)
    CALL check(nf90_get_var(ncid,vid,pn),8)

    CALL check(nf90_close(ncid),9)


    A(:,:) = 1/(pm(:,:)*pn(:,:))

    where(h.lt.hmin) h = hmin 
    where(mask.eq.0) h = hmin
    where(mask.eq.1) mk = .true.

    CALL calc_rx0(nx,ny,mk,h,rx0)
    write(*,*) 'rx0=',maxval(rx0)

    do j = ny,1,-1
      if (mask(isrc,j).eq.1) then
        jstr = j+1

        do i = isrc,1,-1
          if (mask(i,j).eq.0) then
            istr = i+1
            exit
          end if
        end do
    
        do i = isrc,nx
          if (mask(i,j).eq.0) then
            iend = i-1
            exit
          end if
        end do

        exit
      end if
    end do

    write(*,*) istr,iend,jstr

    mask(istr:iend,jstr) = 3

    CALL meo_scheme(nx,ny,mask,mk,h,A,trx0)
    do i = 1,5
      CALL gauss(nx,ny,h,mk)
    end do
    CALL calc_rx0(nx,ny,mk,h,rx0)
    write(*,*) 'rx0=',maxval(rx0)

    CALL check(nf90_create( "diag.nc", NF90_NETCDF4, ncid ), 200)
    CALL check(nf90_def_dim( ncid, 'x', nx, x_dimid ), 201)
    CALL check(nf90_def_dim( ncid, 'y', ny, y_dimid ), 202)
    dimids2 = (/ x_dimid, y_dimid /)
    CALL check(nf90_def_var( ncid, 'rx0', NF90_DOUBLE, dimids2, varid(1)), 205)
    CALL check(nf90_def_var( ncid, 'mask', NF90_INT, dimids2, varid(2)), 205)
    CALL check( nf90_enddef( ncid ), 207 )
    CALL check( nf90_put_var( ncid, varid(1), rx0 ), 208 )
    CALL check( nf90_put_var( ncid, varid(2), mask ), 208 )
    CALL check(nf90_close( ncid ), 209 ) 

    CALL check(nf90_open(trim(gfile),NF90_WRITE,ncid),10)
    CALL check(nf90_inq_varid(ncid,"h",vid),11)
    CALL check(nf90_put_var(ncid,vid,h),12)
    CALL check(nf90_close(ncid),15)
 
    deallocate( h, mask, mk, rx0 )
    deallocate( pm, pn, A )
END PROGRAM ROMS_bathy

SUBROUTINE check(status,label)
    USE NETCDF
    implicit none
    INTEGER, INTENT(in) :: status , label

    IF (status .ne. nf90_noerr) then
        print*, trim(nf90_strerror(status)),' label=',label
        stop "Stopped"
    END IF
END SUBROUTINE check

SUBROUTINE calc_rx0(nx,ny,mask,h,rx0)
    implicit none
    integer, intent(in) :: nx,ny
    logical, intent(in) :: mask(nx,ny)
    real(kind=8), intent(in) :: h(nx,ny)
    real(kind=8), intent(out) :: rx0(nx,ny)

    integer :: i,j
    real(kind=8) :: rx0_tmp

    rx0 = 0

    do i = 2,nx-1
    do j = 2,nx-1
      if (mask(i,j)) then
        if(mask(i+1,j)) then
          rx0_tmp = abs(h(i,j)-h(i+1,j))/(h(i,j)+h(i+1,j))
          rx0(i,j) = max(rx0_tmp,rx0(i,j))
        end if
        if(mask(i-1,j)) then
          rx0_tmp = abs(h(i,j)-h(i-1,j))/(h(i,j)+h(i-1,j))
          rx0(i,j) = max(rx0_tmp,rx0(i,j))
        end if
        if(mask(i,j+1)) then
          rx0_tmp = abs(h(i,j)-h(i,j+1))/(h(i,j)+h(i,j+1))
          rx0(i,j) = max(rx0_tmp,rx0(i,j))
        end if
        if(mask(i,j-1)) then
          rx0_tmp = abs(h(i,j)-h(i,j-1))/(h(i,j)+h(i,j-1))
          rx0(i,j) = max(rx0_tmp,rx0(i,j))
        end if
      end if
    end do
    end do
END SUBROUTINE calc_rx0

SUBROUTINE meo_scheme(nx,ny,mask,mk,h,A,trx0)
    implicit none
    integer, intent(in) :: nx,ny
    integer, intent(inout) :: mask(nx,ny)
    real(kind=8), intent(in) :: A(nx,ny),trx0
    logical, intent(in) :: mk(nx,ny)
    real(kind=8), intent(inout) :: h(nx,ny)

    integer :: i,j,k,l,m,dx,dy,idx(4,2),dx2,dy2
    real(kind=8) :: V, rx0(nx,ny),calcv1,calcv2
    logical :: find

    idx(:,1) = (/1,0,-1,0/)
    idx(:,2) = (/0,-1,0,1/)

    do k = 1, 1000
    do i = 2,nx-1
    do j = 2,ny-1
      if (mask(i,j).eq.1) then

        do l = 1,4
        dx = i+idx(l,1)
        dy = j+idx(l,2)

        if(mask(dx,dy).eq.1) then
          if ( (h(i,j)-h(dx,dy))/(h(i,j)+h(dx,dy)).gt.trx0 ) then
            V = calcv1(h(i,j),h(dx,dy),A(i,j),A(dx,dy),trx0)
!            write(*,*) 'i,j=',i,j,'dh=',V/A(i,j)
            h(i,j) = h(i,j)-V/A(i,j)
            h(dx,dy) = h(dx,dy)+V/A(dx,dy)
          end if
        end if

!        if(mask(dx,dy).eq.2 .or. mask(dx,dy).eq.0) then
!          if ( (h(i,j)-h(dx,dy))/(h(i,j)+h(dx,dy)).gt.trx0 ) then
!            write(*,*) 'i,j=',i,j,'dh=',V/A(i,j)
!            h(i,j) = h(dx,dy)*(1+trx0)/(1-trx0)
              
!            do m = 1,3
!              dx2=i+idx(mod(l+m-1,4)+1,1)
!              dy2=j+idx(mod(l+m-1,4)+1,2)
!              if (mask(dx2,dy2).eq.1) then
!                h(dx2,dy2) = h(dx2,dy2)+V/A(dx2,dy2)
!                EXIT
!              end if
!            end do
            
!          end if
!        end if

        end do

      end if
    end do
    end do
    CALL calc_rx0(nx,ny,mk,h,rx0)
    write(*,*) 'rx0=',maxval(rx0)
    if (maxval(rx0).lt.0.1) EXIT
    end do

END SUBROUTINE meo_scheme

FUNCTION calcv1(h1,h2,a1,a2,r) result(v)
    implicit none
    real(kind=8), intent(in) :: h1,h2,a1,a2,r
    real(kind=8) :: v

    v = a1*a2*(h1-h2-r*h1-r*h2)/(a1*r+a1-a2*r+a2)
END FUNCTION calcv1

FUNCTION calcv2(h1,h2,a1,r) result(v)
    implicit none
    real(kind=8), intent(in) :: h1,h2,a1,r
    real(kind=8) :: v

    v = a1*(h1-h2-r*h1-r*h2)/(1-r)
END FUNCTION calcv2

SUBROUTINE gauss(nx,ny,h,mk)
    implicit none
    integer, intent(in) :: nx,ny
    real(kind=8), intent(inout) :: h(nx,ny)
    logical, intent(in) :: mk(nx,ny)

    integer :: i,j,wg(5,5)
    logical :: mk2(5,5)
    real(kind=8) :: win(5,5),tmp(nx,ny) 

    wg(1,:) = (/1,4,7,4,1/)
    wg(2,:) = (/4,16,26,16,4/)
    wg(3,:) = (/7,26,41,26,7/)
    wg(4,:) = (/4,16,26,16,4/)
    wg(5,:) = (/1,4,7,4,1/)


    do i = 3, nx-2
    do j = 3, ny-2
      if (mk(i,j)) then
        mk2 = mk(i-2:i+2,j-2:j+2)
        win = h(i-2:i+2,j-2:j+2)
        tmp(i,j) = sum(win(:,:)*real(wg(:,:),8),mask=mk2)/real(sum(wg,mask=mk2),8)
      end if
    end do
    end do
    where(mk) h=tmp
    
END SUBROUTINE gauss
