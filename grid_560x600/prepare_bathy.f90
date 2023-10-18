PROGRAM prepare_bathy
    USE NETCDF
    USE omp_lib
    implicit none
    integer :: nx,ny,ncid,dimid,vid,varid(2)
    integer :: x_dimid, y_dimid, z_dimid, t_dimid
    integer :: dimids3(3), dimids4(4)
    integer :: onx,ony,i,j,k,l,x,y,reclen
    integer, allocatable :: mask(:,:), omask(:,:), idx(:,:,:,:)
    real(kind=8) :: point(2), corners(4,2), rx0
    real(kind=8), allocatable :: h(:,:), pm(:,:), pn(:,:), A(:,:)
    real(kind=8), allocatable :: lon(:,:), lat(:,:), olon(:,:), olat(:,:)
    real(kind=8), allocatable :: oh(:,:), W(:,:,:), dis_array(:,:)
    character(len=100) :: infile, outfile
    character(len=20) :: file1, file2
    logical :: score, fexit, ex1, ex2
    logical, allocatable :: omasklog(:,:)
  
    !SET STACKSIZE EQUAL 1024 MB per thread
    CALL KMP_SET_STACKSIZE_S(3221225472)
    !END SET
 
    infile = 'ROMS_grid_full_mask_025NM.nc'
    outfile = 'ROMS_grid_2_3km_560x600_NetCDF4v2.nc'
    file1 = "idx.bin"
    file2 = "W.bin"
    
!READ INPUT GRID
    CALL check(nf90_open(trim(infile),NF90_NOWRITE,ncid),1)

    CALL check(nf90_inq_dimid(ncid, "xi_rho", dimid),2)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nx),2)
    CALL check(nf90_inq_dimid(ncid, "eta_rho", dimid),3)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=ny),4)

    write(*,*) 'nx/ny', nx, ny
    allocate( h(nx,ny), mask(nx,ny), lon(nx,ny), lat(nx,ny) )
    allocate( dis_array(nx,ny) )

    CALL check(nf90_inq_varid(ncid,"hraw",vid),5)
    CALL check(nf90_get_var(ncid,vid,h),6)
    CALL check(nf90_inq_varid(ncid,"mask_rho",vid),7)
    CALL check(nf90_get_var(ncid,vid,mask),8)
    CALL check(nf90_inq_varid(ncid,"lon_rho",vid),7)
    CALL check(nf90_get_var(ncid,vid,lon),8)
    CALL check(nf90_inq_varid(ncid,"lat_rho",vid),7)
    CALL check(nf90_get_var(ncid,vid,lat),8)

    CALL check(nf90_close(ncid),9)
!END READ

!READ INPUT GRID
    CALL check(nf90_open(trim(outfile),NF90_NOWRITE,ncid),1)

    CALL check(nf90_inq_dimid(ncid, "xi_rho", dimid),2)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=onx),2)
    CALL check(nf90_inq_dimid(ncid, "eta_rho", dimid),3)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=ony),4)

    allocate( oh(onx,ony), omask(onx,ony), olon(onx,ony), olat(onx,ony) )
    allocate( idx(onx,ony,4,2), W(onx,ony,4), omasklog(onx,ony) )
    allocate( pm(onx,ony), pn(onx,ony), A(onx,ony) )
    

    CALL check(nf90_inq_varid(ncid,"mask_rho",vid),7)
    CALL check(nf90_get_var(ncid,vid,omask),8)
    CALL check(nf90_inq_varid(ncid,"lon_rho",vid),7)
    CALL check(nf90_get_var(ncid,vid,olon),8)
    CALL check(nf90_inq_varid(ncid,"lat_rho",vid),7)
    CALL check(nf90_get_var(ncid,vid,olat),8)
    CALL check(nf90_inq_varid(ncid,"pm",vid),7)
    CALL check(nf90_get_var(ncid,vid,pm),8)
    CALL check(nf90_inq_varid(ncid,"pn",vid),7)
    CALL check(nf90_get_var(ncid,vid,pn),8)

    CALL check(nf90_close(ncid),9)

    A(:,:) = 1/(pm(:,:)*pn(:,:))
!END READ

!INTERP COEF
    INQUIRE(FILE=trim(file1),EXIST=ex1)
    INQUIRE(FILE=trim(file2),EXIST=ex2)
    if (ex1 .and. ex2) then

      inquire(iolength = reclen) idx

      OPEN(10,FILE=trim(file1),ACCESS='DIRECT', recl=reclen, FORM='UNFORMATTED')
      read(10, rec=1) idx
      CLOSE(10)

      inquire(iolength = reclen) W

      OPEN(10,FILE=trim(file2),ACCESS='DIRECT', recl=reclen, FORM='UNFORMATTED')
      read(10, rec=1) W
      CLOSE(10)

    else

    idx = 0
    W = 0
!$OMP PARALLEL DO DEFAULT(PRIVATE),SHARED(lon,lat,olon,olat),&
!$OMP& SHARED(nx,ny,idx,W,omask), FIRSTPRIVATE(onx,ony) SCHEDULE(DYNAMIC)
    do j = 1, ony
    do i = 1, onx
      if (omask(i,j).eq.1) then
    
        dis_array(:,:) = sqrt( (olon(i,j)-lon(:,:))**2+(olat(i,j)-lat(:,:))**2 )
        idx(i,j,1,:) = minloc( dis_array )
        x = idx(i,j,1,1)
        y = idx(i,j,1,2)
        if (x.gt.1 .and. x.lt.nx .and. y.gt.1 .and. y.lt.ny) then
          point(1) = olon(i,j)
          point(2) = olat(i,j)
          fexit = .false.
          do k = x-1,x
            if (fexit) EXIT
          do l = y-1,y
            if (fexit) EXIT
            corners(1,:) = (/lon(k,l),lat(k,l)/)
            corners(2,:) = (/lon(k,l+1),lat(k,l+1)/)
            corners(3,:) = (/lon(k+1,l+1),lat(k+1,l+1)/)
            corners(4,:) = (/lon(k+1,l),lat(k+1,l)/)
            CALL in_convex_polygon(corners,point,score)
            if (score) then
              fexit = .true.
              idx(i,j,1,:) = (/k,l/)
              idx(i,j,2,:) = (/k,l+1/)
              idx(i,j,3,:) = (/k+1,l+1/)
              idx(i,j,4,:) = (/k+1,l/)
              call calc_w2( lon(k:k+1,l:l+1), lat(k:k+1,l:l+1), point, W(i,j,:) )
            end if
          end do
          end do
        else
          idx(i,j,1,:) = 0
        end if

      end if
    end do
    end do
!$OMP END PARALLEL DO
      inquire(iolength = reclen) idx

      OPEN(10,FILE=trim(file1),ACCESS='DIRECT', recl=reclen, FORM='UNFORMATTED')
      write(10, rec=1) idx
      CLOSE(10)

      inquire(iolength = reclen) W

      OPEN(10,FILE=trim(file2),ACCESS='DIRECT', recl=reclen, FORM='UNFORMATTED')
      write(10, rec=1) W
      CLOSE(10)

    end if !files exist
!END INTERP COEF

    CALL check(nf90_create( "interp_coef.nc", NF90_NETCDF4, ncid ), 200)
    CALL check(nf90_def_dim( ncid, 'x', onx, x_dimid ), 201)
    CALL check(nf90_def_dim( ncid, 'y', ony, y_dimid ), 202)
    CALL check(nf90_def_dim( ncid, 'corner', 4, z_dimid ), 202)
    CALL check(nf90_def_dim( ncid, 'index', 2, t_dimid ), 202)
    dimids3 = (/ x_dimid, y_dimid, z_dimid /)
    dimids4 = (/ x_dimid, y_dimid, z_dimid, t_dimid /)
    CALL check(nf90_def_var( ncid, 'idx', NF90_INT, dimids4, varid(1)), 205)
    CALL check(nf90_def_var( ncid, 'w', NF90_DOUBLE, dimids3, varid(2)), 205)
    CALL check( nf90_enddef( ncid ), 207 )
    CALL check( nf90_put_var( ncid, varid(1), idx ), 209 )
    CALL check( nf90_put_var( ncid, varid(2), W ), 209 )
    CALL check(nf90_close( ncid ), 210 )

    write(*,*) 'min W:', minval(W)
    write(*,*) 'max W:', maxval(W)

    where( omask.eq.0 ) oh = 2.0
    CALL calcme(nx,ny,onx,ony,h,oh,W,idx)
    
    omasklog = .false.
    where( omask.eq.1 ) omasklog=.true.
!    CALL calc_rx0(onx,ony,omasklog,oh,rx0)
!    write(*,*) 'rx0=',rx0
    CALL meo_scheme(onx,ony,omasklog,oh,A,0.12)

    CALL check(nf90_open( trim(outfile), nf90_write, ncid ), 200)
    CALL check(nf90_inq_varid(ncid, "h", varid(1)),201)
    CALL check(nf90_put_var(ncid, varid(1), oh),202)
    CALL check(nf90_close( ncid ), 203 )

    deallocate( h, mask, lon, lat, pm, pn, A )
    deallocate( oh, omask, olon, olat )
    deallocate( dis_array, idx, W, omasklog )
END PROGRAM prepare_bathy

SUBROUTINE check(status,label)
    USE NETCDF
    implicit none
    INTEGER, INTENT(in) :: status , label

    IF (status .ne. nf90_noerr) then
        print*, trim(nf90_strerror(status)),' label=',label
        stop "Stopped"
    END IF
END SUBROUTINE check

SUBROUTINE calcme(nxin,nyin,nxout,nyout,inarray,outarray,W,idx)
    implicit none
    integer, intent(in) :: nxin,nyin,nxout,nyout
    real(kind=8), intent(in) :: inarray(nxin,nxin),W(nxout,nyout,4)
    real(kind=8), intent(out) :: outarray(nxout,nyout)
    integer, intent(in) :: idx(nxout,nyout,4,2)

    real(kind=8) :: tmp
    integer :: i,j,k,x,y

    do i = 1,nxout
    do j = 1,nyout
      if (idx(i,j,1,1).ne.0) then
        tmp = 0
        do k = 1,4
          x = idx(i,j,k,1)
          y = idx(i,j,k,2)
          if (x.gt.0 .and. y.gt.0) then
            tmp = tmp+inarray(x,y)*W(i,j,k)
          end if
        end do
        outarray(i,j) = tmp
      end if
    end do
    end do
END SUBROUTINE calcme

SUBROUTINE calc_w2(lon,lat,point_xy,w)
    implicit none
    real(kind=8), intent(in) :: lon(2,2),lat(2,2)
    real(kind=8), intent(in) :: point_xy(2)
    real(kind=8), intent(out) :: w(4)

    integer :: idx(2)
    real(kind=8) :: x, y, dis_array(2,2)

    x = point_xy(1)
    y = point_xy(2)
    dis_array(:,:) = sqrt((x-lon(:,:))**2+(y-lat(:,:))**2)
    idx = minloc( dis_array )

    if ( dis_array(idx(1),idx(2)).eq.0 ) then
        w = 0
        if (x.eq.1 .and. y.eq.1) then
            w(1) = 1
        elseif (x.eq.1 .and. y.eq.2) then
            w(2) = 1
        elseif (x.eq.2 .and. y.eq.2) then
            w(3) = 1
        elseif (x.eq.2 .and. y.eq.1) then
            w(4) = 1
        end if
    else
        dis_array(:,:) = 1/dis_array(:,:)
        w(1) = dis_array(1,1)/sum(dis_array)
        w(2) = dis_array(1,2)/sum(dis_array)
        w(3) = dis_array(2,2)/sum(dis_array)
        w(4) = dis_array(2,1)/sum(dis_array)
    end if

    
END SUBROUTINE calc_w2

SUBROUTINE calc_w(lon,lat,point_xy,w)
    implicit none
    real(kind=8), intent(in) :: lon(2,2),lat(2,2)
    real(kind=8), intent(in) :: point_xy(2)
    real(kind=8), intent(out) :: w(4)

    integer, parameter :: long_double = SELECTED_REAL_KIND(22)
    real(kind=long_double) :: x, y
    real(kind=long_double) :: A, B, C, s, t

    x = point_xy(1)
    y = point_xy(2)

    A = ( lon(1,1)-lon(1,2) )*( lat(2,1)-lat(2,2) )-&
        ( lat(1,1)-lat(1,2) )*( lon(2,1)-lon(2,2) )

    B = y*( ( lon(2,1)-lon(2,2) )-( lon(1,1)-lon(1,2) ) )-&
        x*( ( lat(2,1)-lat(2,2) )-( lat(1,1)-lat(1,2) ) )+&
        ( lon(1,1)-lon(1,2) )*lat(2,2)-&
        ( lat(1,1)-lat(1,2) )*lon(2,2)+&
        ( lat(2,1)-lat(2,2) )*lon(1,2)-&
        ( lon(2,1)-lon(2,2) )*lat(1,2)

    C = y*( lon(2,2)-lon(1,2) )-&
        x*( lat(2,2)-lat(1,2) )+&
        lon(1,2)*lat(2,2)-&
        lat(1,2)*lon(2,2)


    t = ( -B+sqrt(B**2-4*A*C) )/(2*A)
    s = (  y-lat(1,2)-( lat(1,1)-lat(1,2) )*t  )/&
        (  lat(2,2)+( lat(2,1)-lat(2,2) )*t-&
        lat(1,2)-( lat(1,1)-lat(1,2) )*t  )

    if ( (t<0).OR.(t>1) ) then
      t = ( -B-sqrt(B**2-4*A*C) )/(2*A)
      s = (  y-lat(1,2)-( lat(1,1)-lat(1,2) )*t  )/&
          (  lat(2,2)+( lat(2,1)-lat(2,2) )*t-&
          lat(1,2)-( lat(1,1)-lat(1,2) )*t  )
    end if

    w(1) = (1-s)*t
    w(2) = (1-s)*(1-t)
    w(3) = s*(1-t)
    w(4) = s*t

!    if ( maxval(w).gt.1 .or. minval(w).lt.0 ) then
!         write(*,*) 'ERROR',lon,lat,point_xy
!        write(*,*) 'A,b,C',A,B,C
!        write(*,*) 't,s',t,s
!        write(*,*) (lat(2,2)+( lat(2,1)-lat(2,2) )*t),(lat(1,2)+( lat(1,1)-lat(1,2) )*t)
!    end if
END SUBROUTINE calc_w

SUBROUTINE in_convex_polygon(corners,point,score)
    implicit none
    real(kind=8), intent(in) :: corners(4,2)
    real(kind=8), intent(in) :: point(2)
    logical, intent(out) :: score
    logical :: state,side
    real(kind=8) :: xi,yi,xj,yj,d
    integer*1 :: i

    score=.TRUE.
    state=.TRUE.

    do i = 1,4
      xi = corners(i,1)
      yi = corners(i,2)
      xj = corners(modulo(i,4)+1,1)
      yj = corners(modulo(i,4)+1,2)
      d = (point(1)-xi)*(yj-yi)-(point(2)-yi)*(xj-xi)
      if ( d==0.0 ) then
        go to 10
      else
        if (state) then
          state=.FALSE.
          side=(d>0.0)
        else if ((d>0.0)/=side) then
          score=.FALSE.
          exit
        end if
      end if
      10 continue
    end do
END SUBROUTINE in_convex_polygon

SUBROUTINE meo_scheme(nx,ny,mask,h,A,trx0)
    implicit none
    integer, intent(in) :: nx,ny
    real, intent(in) :: trx0
    real(kind=8), intent(in) :: A(nx,ny)
    logical, intent(in) :: mask(nx,ny)
    real(kind=8), intent(inout) :: h(nx,ny)

    integer :: i,j,k,l,m,dx,dy,idx(4,2),dx2,dy2
    real(kind=8) :: V, rx0,calcv1,calcv2
    logical :: find

    write(*,*) 'trx0:',trx0

    idx(:,1) = (/1,0,-1,0/)
    idx(:,2) = (/0,-1,0,1/)

    do k = 1, 1000
    do i = 2,nx-1
    do j = 2,ny-1
      if ( mask(i,j) ) then

        do l = 1,4
        dx = i+idx(l,1)
        dy = j+idx(l,2)

        if( mask(dx,dy) ) then
          if ( (h(i,j)-h(dx,dy))/(h(i,j)+h(dx,dy)).gt.trx0 ) then
            V = calcv1(h(i,j),h(dx,dy),A(i,j),A(dx,dy),trx0)
            h(i,j) = h(i,j)-V/A(i,j)
            h(dx,dy) = h(dx,dy)+V/A(dx,dy)
          end if
        end if
        end do

      end if
    end do
    end do
    CALL calc_rx0(nx,ny,mask,h,rx0)
    write(*,*) 'rx0=',rx0
    if (rx0.lt.(trx0+0.01)) EXIT
    end do
END SUBROUTINE meo_scheme

FUNCTION calcv1(h1,h2,a1,a2,r) result(v)
    implicit none
    real(kind=8), intent(in) :: h1,h2,a1,a2,r
    real(kind=8) :: v

    v = a1*a2*(h1-h2-r*h1-r*h2)/(a1*r+a1-a2*r+a2)
END FUNCTION calcv1

SUBROUTINE calc_rx0(nx,ny,mask,h,rx0)
    implicit none
    integer, intent(in) :: nx,ny
    logical, intent(in) :: mask(nx,ny)
    real(kind=8), intent(in) :: h(nx,ny)
    real(kind=8), intent(out) :: rx0

    integer :: i,j
    real(kind=8) :: rx0_tmp,rx02d(nx,ny)

    rx02d = 0

    do i = 2,nx-1
    do j = 2,nx-1
      if (mask(i,j)) then
        if(mask(i+1,j)) then
          rx0_tmp = abs(h(i,j)-h(i+1,j))/(h(i,j)+h(i+1,j))
          rx02d(i,j) = max(rx0_tmp,rx02d(i,j))
        end if
        if(mask(i-1,j)) then
          rx0_tmp = abs(h(i,j)-h(i-1,j))/(h(i,j)+h(i-1,j))
          rx02d(i,j) = max(rx0_tmp,rx02d(i,j))
        end if
        if(mask(i,j+1)) then
          rx0_tmp = abs(h(i,j)-h(i,j+1))/(h(i,j)+h(i,j+1))
          rx02d(i,j) = max(rx0_tmp,rx02d(i,j))
        end if
        if(mask(i,j-1)) then
          rx0_tmp = abs(h(i,j)-h(i,j-1))/(h(i,j)+h(i,j-1))
          rx02d(i,j) = max(rx0_tmp,rx02d(i,j))
        end if
      end if
    end do
    end do

    rx0 = maxval(rx02d)
END SUBROUTINE calc_rx0
