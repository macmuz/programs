PROGRAM interp_salt
    USE NETCDF
    implicit none
    integer, parameter :: nxin=600,nyin=640,nz=66,nxout=1636,nyout=1280,time=12
    integer :: i,j,reclen,reclen2,z,t,x,y
    integer :: ncid,varid(4),xdim,ydim,zdim,tdim
    integer :: ncid2,varid2(3),xdim2,ydim2,zdim2,tdim2
    integer :: idx(nxout,nyout,2),intmask(nxout,nyout)
    real(kind=8) :: pi
    real(kind=8) :: tlonin(nxin,nyin),tlatin(nxin,nyin)
    real(kind=8) :: ulonout(nxout,nyout),ulatout(nxout,nyout)
    real(kind=8) :: tlonout(nxout,nyout),tlatout(nxout,nyout)
    real(kind=8) :: pt(2),alphas(nxout,nyout,4)
    real(kind=8) :: datain(nxin,nyin),dataout(nxout,nyout) 
    character(len=100) :: gridin,gridout,input,output
    logical :: mask(nxout,nyout)

    gridin = 'grid.bs2v1.ocn.20110427.ieeer8'
    gridout = 'grid_baltic_1636x1280.bin'
    input = 'salt_3d_bs2v3_monthly_20130618.ieeer8'
    output = 'salt_3d_bs2v3_monthly_20130618_1636x1280'
    
    pi = 4.D0*DATAN(1.D0)
 
    inquire(iolength = reclen) tlonin
    open(10,FILE=trim(gridin),ACCESS='DIRECT', recl=reclen, FORM='UNFORMATTED')
    read(10, rec=8) tlatin
    read(10, rec=9) tlonin
    close(10)

    write(*,*) minval(tlatin)*180/pi,maxval(tlatin)*180/pi
    write(*,*) minval(tlonin)*180/pi,maxval(tlonin)*180/pi
   
    inquire(iolength = reclen2) ulonout
    open(10,FILE=trim(gridout),ACCESS='DIRECT', recl=reclen2, FORM='UNFORMATTED')
    read(10, rec=1) ulatout
    read(10, rec=2) ulonout
    close(10)

    write(*,*) minval(ulatout)*180/pi,maxval(ulatout)*180/pi
    write(*,*) minval(ulonout)*180/pi,maxval(ulonout)*180/pi

    CALL Tlatlon(nxout,nyout,ulonout,ulatout,tlonout,tlatout)

    write(*,*) minval(tlatout)*180/pi,maxval(tlatout)*180/pi
    write(*,*) minval(tlonout)*180/pi,maxval(tlonout)*180/pi
  
    do i = 1,nxout 
    do j = 1,nyout 
!    do i = 100,120 
!    do j = 100,120
      pt = (/tlonout(i,j),tlatout(i,j)/)
      CALL find_idx(nxin,nyin,tlonin,tlatin,pt,idx(i,j,:),alphas(i,j,:))
!      write(*,*) pt*180/pi 
!      write(*,*) idx(i,j,:) 
!      write(*,*) alphas(i,j,:) 
    end do
    end do

    mask = .false.
    do i = 1,nxout
    do j = 1,nyout
      if (idx(i,j,1).ne.0) mask(i,j)=.true.
    end do
    end do
    intmask = 0
    where(mask) intmask = 1

!CREATE OUTPUT NETCDF
    CALL check( nf90_create( 'input.nc',NF90_NETCDF4,ncid2 ), 300 )
    CALL check( nf90_def_dim( ncid2, 'nx', nxin, xdim2 ), 301 )
    CALL check( nf90_def_dim( ncid2, 'ny', nyin, ydim2 ), 301 )
    CALL check( nf90_def_dim( ncid2, 'nz', nz, zdim2 ), 301 )
    CALL check( nf90_def_dim( ncid2, 'time', NF90_UNLIMITED, tdim2 ), 301 )

    CALL check(nf90_def_var( ncid2, 'lon', NF90_DOUBLE, (/xdim2,ydim2/), varid2(1)), 327)
    CALL check( nf90_put_att( ncid2, varid2(1), 'long_name',&
        'longitude of t-points' ), 328 )
    CALL check( nf90_put_att( ncid2, varid2(1), 'units',&
        'degree_east' ), 329 )

    CALL check(nf90_def_var( ncid2, 'lat', NF90_DOUBLE, (/xdim2,ydim2/), varid2(2)), 327)
    CALL check( nf90_put_att( ncid2, varid2(2), 'long_name',&
        'latitude of t-points' ), 328 )
    CALL check( nf90_put_att( ncid2, varid2(2), 'units',&
        'degree_north' ), 329 )

    CALL check(nf90_def_var( ncid2, 'salt', NF90_DOUBLE, (/xdim2,ydim2,zdim2,tdim2/), varid2(3)), 327)
    CALL check( nf90_put_att( ncid2, varid2(3), 'long_name',&
        'salinity' ), 328 )

    CALL check( nf90_enddef(ncid2), 340 )
    CALL check( nf90_put_var( ncid2, varid2(1), tlonin*180/pi ), 347 )
    CALL check( nf90_put_var( ncid2, varid2(2), tlatin*180/pi ), 347 )
!END CREATE OUTPUT NETCDF

!CREATE OUTPUT NETCDF
    CALL check( nf90_create( trim(output)//'.nc',NF90_NETCDF4,ncid ), 300 )
    CALL check( nf90_def_dim( ncid, 'nx', nxout, xdim ), 301 )
    CALL check( nf90_def_dim( ncid, 'ny', nyout, ydim ), 301 )
    CALL check( nf90_def_dim( ncid, 'nz', nz, zdim ), 301 )
    CALL check( nf90_def_dim( ncid, 'time', NF90_UNLIMITED, tdim ), 301 )
   
    CALL check(nf90_def_var( ncid, 'lon', NF90_DOUBLE, (/xdim,ydim/), varid(1)), 327)
    CALL check( nf90_put_att( ncid, varid(1), 'long_name',&
        'longitude of t-points' ), 328 )
    CALL check( nf90_put_att( ncid, varid(1), 'units',&
        'degree_east' ), 329 ) 

    CALL check(nf90_def_var( ncid, 'lat', NF90_DOUBLE, (/xdim,ydim/), varid(2)), 327)
    CALL check( nf90_put_att( ncid, varid(2), 'long_name',&
        'latitude of t-points' ), 328 )
    CALL check( nf90_put_att( ncid, varid(2), 'units',&
        'degree_north' ), 329 )
 
    CALL check(nf90_def_var( ncid, 'salt', NF90_DOUBLE, (/xdim,ydim,zdim,tdim/), varid(3)), 327)
    CALL check( nf90_put_att( ncid, varid(3), 'long_name', 'salinity' ), 328 )

    CALL check(nf90_def_var( ncid, 'mask', NF90_INT, (/xdim,ydim/), varid(4)), 327)
    CALL check( nf90_put_att( ncid, varid(4), 'long_name', 'mask' ), 328 )

    CALL check( nf90_enddef(ncid), 340 )
    CALL check( nf90_put_var( ncid, varid(1), tlonout*180/pi ), 347 )
    CALL check( nf90_put_var( ncid, varid(2), tlatout*180/pi ), 347 )
    CALL check( nf90_put_var( ncid, varid(4), intmask ), 347 )
!END CREATE OUTPUT NETCDF

    inquire(iolength = reclen) datain
    inquire(iolength = reclen2) dataout

    open(10,FILE=trim(input),ACCESS='DIRECT', recl=reclen, FORM='UNFORMATTED')
    open(11,FILE=trim(output)//'.ieeer8',ACCESS='DIRECT',&
         recl=reclen2, FORM='UNFORMATTED')

    do t = 1,time
    do z = 1,nz

      read(10, rec=(t-1)*nz+z) datain
      CALL check( nf90_put_var( ncid2, varid2(3), datain, &
        start=(/1,1,z,t/), count=(/nxin,nyin,1,1/) ), 347 )

      dataout = 0.0
      do i = 1,nxout
      do j = 1,nyout
        if (idx(i,j,1).ne.0) then
          x = idx(i,j,1)
          y = idx(i,j,2)
          dataout(i,j) = datain(x,y)*alphas(i,j,1)+datain(x+1,y)*alphas(i,j,2)+&
            datain(x+1,y+1)*alphas(i,j,3)+datain(x,y+1)*alphas(i,j,4)
        end if
      end do
      end do
!      CALL extrap(dataout,mask,nxout,nyout,100,2) 

      write(11, rec=(t-1)*nz+z) dataout
      CALL check( nf90_put_var( ncid, varid(3), dataout, &
        start=(/1,1,z,t/), count=(/nxout,nyout,1,1/) ), 347 )

      write(*,*) t,z
    end do
    end do

    close(11)
    close(10)

    CALL check(nf90_close(ncid),350)
    CALL check(nf90_close(ncid2),350)
 
END PROGRAM interp_salt

SUBROUTINE check(status,label)
    USE NETCDF
    implicit none
    INTEGER, INTENT(in) :: status , label

    IF (status .ne. nf90_noerr) then
        print*, trim(nf90_strerror(status)),'label=',label
        stop "Stopped"
    END IF
END SUBROUTINE check

subroutine Tlatlon(nx,ny,lonu,latu,tlon,tlat)
    implicit none
    integer, intent(in) :: nx,ny
    real(kind=8), intent(in) :: lonu(nx,ny),latu(nx,ny)
    real(kind=8), intent(out) :: tlon(nx,ny),tlat(nx,ny)
    
    integer :: i,j
    real(kind=8) :: ULON(0:nx,0:ny),ULAT(0:nx,0:ny)
    real(kind=8) :: z1,x1,y1,z2,x2,y2,z3,x3,y3,z4,x4,y4,tx,ty,tz,da

    ULON(1:nx,1:ny) = lonu 
    ULAT(1:nx,1:ny) = latu

    ULON(0,1:ny) = 2*ULON(1,1:ny)-ULON(2,1:ny)
    ULON(:,0) = 2*ULON(:,1)-ULON(:,2)

    ULAT(1:nx,0) = 2*ULAT(1:nx,1)-ULAT(1:nx,2)
    ULAT(0,:) = 2*ULAT(1,:)-ULAT(2,:)

    do i = 1,nx
      do j = 1,ny
        z1 = cos(ULAT(i-1,j-1))
        x1 = cos(ULON(i-1,j-1))*z1
        y1 = sin(ULON(i-1,j-1))*z1
        z1 = sin(ULAT(i-1,j-1))

        z2 = cos(ULAT(i,j-1))
        x2 = cos(ULON(i,j-1))*z2
        y2 = sin(ULON(i,j-1))*z2
        z2 = sin(ULAT(i,j-1))

        z3 = cos(ULAT(i-1,j))
        x3 = cos(ULON(i-1,j))*z3
        y3 = sin(ULON(i-1,j))*z3
        z3 = sin(ULAT(i-1,j))

        z4 = cos(ULAT(i,j))
        x4 = cos(ULON(i,j))*z4
        y4 = sin(ULON(i,j))*z4
        z4 = sin(ULAT(i,j))

        tx = (x1+x2+x3+x4)*0.25
        ty = (y1+y2+y3+y4)*0.25
        tz = (z1+z2+z3+z4)*0.25
        da = sqrt(tx**2+ty**2+tz**2)

        tz = tz/da

        ! TLON in radians East
        tlon(i,j) = 0.0
        if (tx /= 0.0 .or. ty /= 0.0) tlon(i,j) = atan2(ty,tx)

        ! TLAT in radians North
        tlat(i,j) = asin(tz)
      end do
    end do
    
endsubroutine Tlatlon

SUBROUTINE find_idx(ni,nj,lon,lat,pt,idx,alphas)
    implicit none
    integer, intent(in) :: ni,nj
    real(kind=8), intent(in) :: lon(ni,nj),lat(ni,nj),pt(2)
    integer, intent(out) :: idx(2)
    real(kind=8), intent(out) :: alphas(4)

    real(kind=8) :: disarray(ni,nj),corners(4,2)
    integer :: tmp(2),k,l,i,j,k1,k2,l1,l2
    logical :: fexit,score
    
    idx = 0
    alphas = 0.0

    do i = 1,ni
    do j = 1,nj
      disarray(i,j) = sqrt((lon(i,j)-pt(1))**2+(lat(i,j)-pt(2))**2)
    end do
    end do
    tmp = minloc(disarray)

    if (tmp(1).eq.1 .and. tmp(2).eq.1) then
      k1 = 1
      k2 = 1
      l1 = 1
      l2 = 1
    elseif (tmp(1).eq.ni .and. tmp(2).eq.1) then
      k1 = ni-1
      k2 = ni-1
      l1 = 1
      l2 = 1
    elseif (tmp(1).eq.ni .and. tmp(2).eq.nj) then
      k1 = ni-1
      k2 = ni-1
      l1 = nj-1
      l2 = nj-1
    elseif (tmp(1).eq.1 .and. tmp(2).eq.nj) then
      k1 = 1
      k2 = 1
      l1 = nj-1
      l2 = nj-1
    elseif (tmp(1).eq.1) then
      k1 = 1
      k2 = 1
      l1 = tmp(2)-1
      l2 = tmp(2)
    elseif (tmp(2).eq.1) then
      k1 = tmp(1)-1
      k2 = tmp(1)
      l1 = 1
      l2 = 1
    elseif (tmp(1).eq.ni) then
      k1 = ni-1
      k2 = ni-1
      l1 = tmp(2)-1
      l2 = tmp(2)
    elseif (tmp(2).eq.nj) then
      k1 = tmp(1)-1
      k2 = tmp(1)
      l1 = nj-1
      l2 = nj-1
    else
      k1 = tmp(1)-1
      k2 = tmp(1)
      l1 = tmp(2)-1
      l2 = tmp(2)
    end if

    fexit = .false.
    do k = tmp(1)-1,tmp(1)
        if (fexit) exit
        do l = tmp(2)-1,tmp(2)
            corners(1,:) = (/lon(k,l),lat(k,l)/)
            corners(2,:) = (/lon(k+1,l),lat(k+1,l)/)
            corners(3,:) = (/lon(k+1,l+1),lat(k+1,l+1)/)
            corners(4,:) = (/lon(k,l+1),lat(k,l+1)/)
            CALL in_convex_polygon(corners,pt,score)
            if (score) then
                idx = (/k,l/)
                CALL calc_w(corners,pt,alphas)
!                write(*,*) 'ALPHAS',alphas
                fexit = .true.
                exit
            endif
        end do
    end do

END SUBROUTINE find_idx

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

SUBROUTINE calc_w(corners,pt,alphas)
    implicit none
    real(kind=8), intent(in) :: corners(4,2),pt(2)
    real(kind=8), intent(out) :: alphas(4)

    real(kind=8) :: v1(2),v2(2),v3(2),v4(2),a(2),b(2),c(2),d(2)
    real(kind=8) :: x,y,dx(2),tol,f(2),Df(2,2),W,Wx,Wy
    real(kind=8) :: alpha1,alpha2,alpha3,alpha4
    integer :: iter

    tol = 1e-12
    iter = 0
    x = 0.5
    y = 0.5
    dx = (/0.1,0.1/)

    v1 = corners(1,:)
    v2 = corners(2,:)
    v3 = corners(3,:)
    v4 = corners(4,:)

    a = v1-pt
    b = v2-v1
    c = v4-v1
    d = v1-v2-v4+v3

    do while(sqrt(dx(1)**2+dx(2)**2).gt.tol .and. iter.lt.20)
      f = a + b*x + c*y + d*x*y
      Df(:,1)=(b + d*y)
      Df(:,2)=(c + d*x)

      W = Df(1,1)*Df(2,2)-Df(2,1)*Df(1,2)
      Wx = -f(1)*Df(2,2)+f(2)*Df(1,2)
      Wy = Df(1,1)*(-1)*f(2)+Df(2,1)*f(1)

      dx(1) = Wx/W
      dx(2) = Wy/W

      x=x+dx(1)
      y=y+dx(2)

      iter = iter+1
      if (sqrt(dx(1)**2+dx(2)**2).gt.10) then
        iter = 20
      end if
    end do

    if (iter < 20) then
        alpha1=(1-y)*(1-x);
        alpha2=x*(1-y);
        alpha3=y*x;
        alpha4=(1-x)*y;
        alphas=(/alpha1,alpha2,alpha3,alpha4/)
    else
        alphas=(/-1,-1,-1,-1/) !wrong values
    endif

END SUBROUTINE calc_w

!    val = alphas(1)*val2d(1,1)+alphas(2)*val2d(2,1)+&
!        alphas(3)*val2d(2,2)+alphas(4)*val2d(1,2)

SUBROUTINE extrap(a,mask,lon,lat,maxscn,met)
    implicit none
    integer, intent(in) :: lon,lat,maxscn,met
    real(kind=8), intent(inout) :: a(lon,lat)
    logical, intent(in) :: mask(lon,lat)

    integer :: i,j,n,cnt,overall
    real(kind=8) :: relc,ave
    real(kind=8), dimension(lon,lat) :: sor,res
    logical :: mask_tmp(lon,lat),mask_tmp2(lon,lat)

    relc=1.0
    sor = 0.0
    where(.not.mask) sor=relc

    select case(met)
    case(0)
      where(.not.mask) a=0.0
    case(1)
      cnt = 0
      ave = 0.0
      do i=1,lon
      do j=1,lat
        if (mask(i,j)) then
            ave=ave+a(i,j)
            cnt=cnt+1
        end if
      end do
      end do
      if ( cnt.GT.0 ) ave = ave/real(cnt,8)
      where(.not.mask) a=ave
    case(2)
      mask_tmp2 = mask
      mask_tmp = mask_tmp2
      do

      overall = 0
      do i = 1, lon
      do j = 1, lat

        if (.not.mask_tmp(i,j)) then

          cnt = 0
          ave = 0.0

          if ( i.gt.1 ) then
            if ( mask_tmp(i-1,j) ) then
              ave = ave+a(i-1,j)
              cnt = cnt+1
            end if
          end if

          if ( j.gt.1 ) then
            if ( mask_tmp(i,j-1) ) then
              ave = ave+a(i,j-1)
              cnt = cnt+1
            end if
          end if

          if ( i.lt.lon ) then
            if ( mask_tmp(i+1,j) ) then
              ave = ave+a(i+1,j)
              cnt = cnt+1
            end if
          end if

          if ( j.lt.lat ) then
            if ( mask_tmp(i,j+1) ) then
              ave = ave+a(i,j+1)
              cnt = cnt+1
            end if
          end if

          if ( cnt.gt.0 ) then
            a(i,j) = ave/(real(cnt,8))
            overall = overall+cnt
            mask_tmp2(i,j) = .true.
          end if

        end if

      end do
      end do

      mask_tmp = mask_tmp2
      if ( overall.eq.0 ) EXIT

      end do
    end select

    do n=1,maxscn

!!$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j),SCHEDULE(DYNAMIC)
        do i=2,lon-1
          do j=2,lat-1
            res(i,j)=0.25*(a(i-1,j)+a(i+1,j)+a(i,j-1)+a(i,j+1))-a(i,j)
          end do
        end do
!!$OMP END PARALLEL DO

        do i=2,lon-1
          res(i,1)=0.3333*(a(i-1,1)+a(i+1,1)+a(i,2))-a(i,1)
          res(i,lat)=0.3333*(a(i-1,lat)+a(i+1,lat)+a(i,lat-1))-a(i,lat)
        end do

        do j=2,lat-1
          res(1,j)=0.3333*(a(1,j-1)+a(1,j+1)+a(2,j))-a(1,j)
          res(lon,j)=0.3333*(a(lon,j-1)+a(lon,j+1)+a(lon-1,j))-a(lon,j)
        end do

        res(1,1)=0.5*(a(1,2)+a(2,1))-a(1,1)
        res(lon,1)=0.5*(a(lon,2)+a(lon-1,1))-a(lon,1)
        res(1,lat)=0.5*(a(1,lat-1)+a(2,lat))-a(1,lat)
        res(lon,lat)=0.5*(a(lon,lat-1)+a(lon-1,lat))-a(lon,lat)

        res=res*sor
        a=a+res
    end do
END SUBROUTINE extrap
