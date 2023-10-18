PROGRAM section
    USE NETCDF
    implicit none
    integer :: ncid,vid,outncid,dimid,ncid2
    integer :: xdimid,zdimid,tdimid,dimids3(3)
    integer :: invid(3),outvid(5)
    integer :: ptdeg,nlines,io,i,j,k,nx,ny,nz,nt,t
    integer :: pos,npt,tmp,cnt,narg,pos2,nw
    integer, allocatable :: near(:,:)
    real(kind=8), allocatable :: points(:,:),coord(:,:),w(:,:)
    real(kind=8) :: dx,dy,hc,fv,d
    real(kind=8), allocatable :: h(:,:),cs_r(:),s_r(:),depth(:,:,:),time(:)
    real(kind=8), allocatable :: zeta(:,:,:),temp(:,:,:,:),salt(:,:,:,:)
    real(kind=8), allocatable :: mask(:,:),dp3d(:,:,:),lon(:,:),lat(:,:)
    real(kind=8), allocatable :: otemp(:,:,:),osalt(:,:,:),odp(:,:,:)
    real(kind=8), allocatable :: distance(:)
    logical, allocatable :: m(:,:)
    character(len=30) :: inpfile
    character(len=50) :: buffer,outfile
    character(len=250) :: infile

    ptdeg = 30
    nw = 30
    inpfile = 'input2.in'
    fv = -999.0

    !############# command line arguments #################################
    narg = command_argument_count()
    if (narg.ne.1) then
        write(*,*) "Program must have 1 argument: filename"
        stop
    end if
    call get_command_argument(1,infile)
!    write(*,*) infile
    pos = index(infile, '/',.true.)
    pos2 = index(infile, '.nc',.true.)
    write(outfile,'(a,a)') infile(pos+1:pos2-1),'_section.nc'
!    write(*,*) trim(outfile)
    !######################################################################


    nlines = 0
    OPEN (1, file = trim(inpfile))
    DO
      READ(1,*,iostat=io)
      IF (io/=0) EXIT
      nlines = nlines + 1
    END DO
    CLOSE (1)

    allocate( points(nlines,2) )
    
    OPEN (1, file = trim(inpfile))
    DO i = 1, nlines
      READ(1,'(A)') buffer
      pos = scan(buffer, ' ')
      read(buffer(1:pos),*) points(i,1)
      read(buffer(pos+1:),*) points(i,2)
!      write(*,*) points(i,:)
    END DO
    CLOSE (1)

    do i = 1, nlines-1
        tmp = nint(sqrt((points(i+1,1)-points(i,1))**2+(points(i+1,2)-&
            points(i,2))**2)*ptdeg)
        if (i.eq.1) then
          npt = tmp
        else
          npt = npt+tmp-1
        end if 
    end do
    
    allocate(coord(npt,2),near(npt,2),w(npt,2),distance(npt))

    cnt = 1
    do i = 1, nlines-1
        tmp = nint(sqrt((points(i+1,1)-points(i,1))**2+(points(i+1,2)-&
            points(i,2))**2)*ptdeg)
        dy = (points(i+1,1)-points(i,1))/(tmp-1)
        dx = (points(i+1,2)-points(i,2))/(tmp-1)
        coord(cnt,1) = points(i,1)
        coord(cnt,2) = points(i,2)
        cnt = cnt+1
        if (tmp.gt.2) then
          do j = 1, tmp-2
            coord(cnt,1) = coord(cnt-1,1)+dy
            coord(cnt,2) = coord(cnt-1,2)+dx
            cnt = cnt+1
          end do
        end if
    end do 
    coord(cnt,1) = points(nlines,1)
    coord(cnt,2) = points(nlines,2)

    distance(1) = 0.0
    do i = 2, npt
      CALL dist(coord(i-1,:),coord(i,:),d) 
      distance(i) = distance(i-1)+d
    end do
   
!    CALL dist(coord(1,:),coord(2,:),d) 
!    write(*,*) coord(1,:)
!    write(*,*) coord(2,:)
!    write(*,*) d

    !READ GRID
!    write(*,*) infile
    CALL check(nf90_open(trim(infile),NF90_NOWRITE,ncid),310)

    CALL check(nf90_inq_dimid(ncid, "xi_rho", dimid),311)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nx),312)

    CALL check(nf90_inq_dimid(ncid, "eta_rho", dimid),313)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=ny),314)

    CALL check(nf90_inq_dimid(ncid, "s_rho", dimid),315)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nz),316)

    CALL check(nf90_inq_dimid(ncid, "ocean_time", dimid),315)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nt),316)

    allocate( h(nx,ny), cs_r(nz), s_r(nz), depth(npt,nz,nt), time(nt) )
    allocate( zeta(nx,ny,nt), temp(nx,ny,nz,nt), salt(nx,ny,nz,nt) )
    allocate( mask(nx,ny), dp3d(nx,ny,nz), lon(nx,ny), lat(nx,ny) )
    allocate( otemp(npt,nz,nt),osalt(npt,nz,nt),odp(npt,nz,nt), m(nx,ny) )

    CALL check(nf90_inq_varid(ncid,"hc",vid),317)
    CALL check(nf90_get_var(ncid,vid,hc),318)

    CALL check(nf90_inq_varid(ncid,"h",vid),319)
    CALL check(nf90_get_var(ncid,vid,h),320)

    CALL check(nf90_inq_varid(ncid,"Cs_r",vid),321)
    CALL check(nf90_get_var(ncid,vid,cs_r),322)

    CALL check(nf90_inq_varid(ncid,"s_rho",vid),323)
    CALL check(nf90_get_var(ncid,vid,s_r),324)

    CALL check(nf90_inq_varid(ncid,"zeta",vid),325)
    CALL check(nf90_get_var(ncid,vid,zeta),326)

    CALL check(nf90_inq_varid(ncid,"mask_rho",vid),325)
    CALL check(nf90_get_var(ncid,vid,mask),326)

    CALL check(nf90_inq_varid(ncid,"lon_rho",vid),325)
    CALL check(nf90_get_var(ncid,vid,lon),326)

    CALL check(nf90_inq_varid(ncid,"lat_rho",vid),325)
    CALL check(nf90_get_var(ncid,vid,lat),326)

    CALL check(nf90_inq_varid(ncid,"ocean_time",invid(1)),327)
    CALL check(nf90_get_var(ncid,invid(1),time),328)

    CALL check(nf90_inq_varid(ncid,"temp",invid(2)),329)
    CALL check(nf90_get_var(ncid,invid(2),temp),330)

    CALL check(nf90_inq_varid(ncid,"salt",invid(3)),331)
    CALL check(nf90_get_var(ncid,invid(3),salt),332)

    !END READ GRID

    CALL check( nf90_create( trim(outfile),NF90_NETCDF4,ncid2 ), 500 )
    CALL check( nf90_def_dim( ncid2, 'npt', npt, xdimid ), 501 )
    CALL check( nf90_def_dim( ncid2, 's_rho', nz, zdimid ), 502 )
    CALL check( nf90_def_dim( ncid2, 'ocean_time', NF90_UNLIMITED, tdimid ), 502 )
    dimids3 = (/ xdimid, zdimid, tdimid /)
    CALL check(nf90_def_var( ncid2, 'ocean_time', NF90_DOUBLE,&
        (/tdimid/), outvid(1) ), 508)
    CALL check( nf90_copy_att(ncid, invid(1), 'long_name', ncid2, outvid(1)), 516 )
    CALL check( nf90_copy_att(ncid, invid(1), 'units', ncid2, outvid(1)), 517 )
    CALL check( nf90_copy_att(ncid, invid(1), 'calendar', ncid2, outvid(1)), 518 )
    CALL check( nf90_copy_att(ncid, invid(1), 'field', ncid2, outvid(1)), 519 )


    CALL check( nf90_def_var( ncid2, 'temp', NF90_DOUBLE,&
        dimids3, outvid(2) ), 508 )
    CALL check( nf90_copy_att(ncid, invid(2), 'long_name', ncid2, outvid(2)), 520 )
    CALL check( nf90_copy_att(ncid, invid(2), 'units', ncid2, outvid(2)), 521 )
    CALL check( nf90_copy_att(ncid, invid(2), 'time', ncid2, outvid(2)), 522 )
    CALL check( nf90_put_att(ncid2, outvid(2), '_FillValue', fv), 901 )   

    CALL check( nf90_def_var( ncid2, 'salt', NF90_DOUBLE,&
        dimids3, outvid(3) ), 508 )
    CALL check( nf90_copy_att(ncid, invid(3), 'long_name', ncid2, outvid(3)), 523 )
    CALL check( nf90_copy_att(ncid, invid(3), 'time', ncid2, outvid(3)), 525 )
    CALL check( nf90_put_att(ncid2, outvid(3), '_FillValue', fv), 901 )   

    CALL check( nf90_def_var( ncid2, 'depth', NF90_DOUBLE,&
        dimids3, outvid(4) ), 508 )
    CALL check( nf90_put_att(ncid2, outvid(4), 'long_name', 'depth (m)'), 901 )   
    CALL check( nf90_put_att(ncid2, outvid(4), 'units', 'm'), 901 )  

    CALL check( nf90_def_var( ncid2, 'distance', NF90_DOUBLE,&
        (/xdimid/), outvid(5) ), 508 )
    CALL check( nf90_put_att(ncid2, outvid(5), 'long_name', 'distance (km)'), 901 )
    CALL check( nf90_put_att(ncid2, outvid(5), 'units', 'km'), 901 ) 

    CALL check( nf90_enddef(ncid2), 516 )

    do i = 1, npt
      CALL prepare(nx,ny,lon,lat,coord(i,2),coord(i,1),near(i,:),w(i,:))
!      write(*,*) coord(i,:),near(i,:),w(i,:)
    end do

    m = .true.
    where(mask.lt.0.5) m = .false.    

    do t = 1, nt
      where(mask.lt.0.5) zeta(:,:,t) = 0.0
      do i = 1, nx
        do j = 1, ny
          do k = 1, nz
            dp3d(i,j,k) = zeta(i,j,t)+(zeta(i,j,t)+h(i,j))*(hc*s_r(k)+h(i,j)*cs_r(k))/(hc+h(i,j))
          end do
        end do
      end do
      do k = 1, nz
        where(mask.lt.0.5) dp3d(:,:,k) = 0.0
        CALL interp_linear(nx,ny,npt,dp3d(:,:,k),near,W,odp(:,k,t))
        where(odp(:,k,t).gt.0) odp(:,k,t)=0.0
        CALL interp(nx,ny,npt,nw,temp(:,:,k,t),near,W,odp(:,k,t),&
            fv,m,otemp(:,k,t)) 
        CALL interp(nx,ny,npt,nw,salt(:,:,k,t),near,W,odp(:,k,t),&
            fv,m,osalt(:,k,t)) 
        write(*,*) t,k
      end do
    end do

!    write(*,*) dp3d(300,150,1), temp(300,150,1,1), salt(300,150,1,1)
    CALL check( nf90_put_var( ncid2, outvid(1), time ), 517 )
    CALL check( nf90_put_var( ncid2, outvid(2), otemp ), 517 )
    CALL check( nf90_put_var( ncid2, outvid(3), osalt ), 517 )
    CALL check( nf90_put_var( ncid2, outvid(4), odp ), 517 )
    CALL check( nf90_put_var( ncid2, outvid(5), distance ), 517 )
    CALL check(nf90_close(ncid2),359)
    CALL check(nf90_close(ncid),360)

    deallocate( points, coord, h, cs_r, s_r, depth )
    deallocate( zeta, temp, salt, time, mask, dp3d )
    deallocate( lon, lat, near, w, otemp, osalt, odp, m )
    deallocate( distance )
END PROGRAM section

SUBROUTINE check(status,label)
    USE NETCDF
    implicit none
    INTEGER, INTENT(in) :: status , label

    IF (status .ne. nf90_noerr) then
        print*, trim(nf90_strerror(status)),' label= ',label
        stop "Stopped"
    END IF
END SUBROUTINE check

SUBROUTINE prepare(nx,ny,lon,lat,olon,olat,near,W)
    implicit none
    integer, intent(in) :: nx,ny
    real(kind=8), intent(in) :: lon(nx,ny),lat(nx,ny),olon,olat
    integer, intent(out) :: near(2)
    real(kind=8), intent(out) :: W(2)

    integer :: idx(2),k,l,x,y
    real(kind=8) :: dis_array(nx,ny)
    real(kind=8) :: point(2),corners(4,2)
    logical :: score,fexit


          do k = 1,nx
            do l = 1,ny
              dis_array(k,l) = sqrt((lon(k,l)-olon)**2+(lat(k,l)-olat)**2)
            end do
          end do

          idx = minloc( dis_array )
          x = idx(1)
          y = idx(2)
          point(1) = olon
          point(2) = olat

          if (x.gt.1 .and. x.lt.nx .and. y.gt.1 .and. y.lt.ny) then
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
                near(1) = k
                near(2) = l
                CALL calc_w(lon(k:k+1,l:l+1),lat(k:k+1,l:l+1),point,W)
                fexit = .true.
              end if
            end do
            end do
          else
            write(*,*) olon,olat,'not found'
          end if

END SUBROUTINE prepare

SUBROUTINE calc_w(lon,lat,point,w)
    implicit none
    real(kind=8), intent(in) :: lon(2,2), lat(2,2), point(2)
    real(kind=8), intent(out) :: w(2)

    real(kind=8) :: a,b,c,d,e,f,g,h
    real(kind=8) :: alpha,beta

    a = -lon(1,1) + lon(2,1)
    b = -lon(1,1) + lon(1,2)
    c = lon(1,1) - lon(1,2) - lon(2,1) + lon(2,2)
    d = point(1) - lon(1,1)
    e = -lat(1,1) + lat(2,1)
    f = -lat(1,1) + lat(1,2)
    g = lat(1,1) - lat(1,2) - lat(2,1) + lat(2,2)
    h = point(2) - lat(1,1)

    alpha = -(b*e - a*f + d*g - c*h + sqrt(-4*(c*e - a*g)*(d*f - b*h) + &
        (b*e - a*f + d*g - c*h)**2))/(2*c*e - 2*a*g)
    beta  = (b*e - a*f - d*g + c*h + sqrt(-4*(c*e - a*g)*(d*f - b*h) + &
        (b*e - a*f + d*g - c*h)**2))/(2*c*f - 2*b*g)


    if (alpha .lt. 0 .or. beta .lt. 0 .or. &
        alpha .gt. 1 .or. beta .gt. 1) then
      alpha = (-b*e + a*f - d*g + c*h + sqrt(-4*(c*e - a*g)*(d*f - b*h) + &
        (b*e - a*f + d*g - c*h)**2))/(2*c*e - 2*a*g)
      beta  = -((-b*e + a*f + d*g - c*h + sqrt(-4*(c*e - a*g)*(d*f - b*h) + &
        (b*e - a*f + d*g - c*h)**2))/( 2*c*f - 2*b*g))
    end if

    w(1) = alpha
    w(2) = beta
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

SUBROUTINE interp(nx,ny,npt,nw,input,idx,W,dp,fv,mask,output)
    implicit none
    integer, intent(in) :: nx,ny,npt,nw,idx(npt,2)
    real(kind=8), intent(in) :: input(nx,ny),W(npt,2),dp(npt),fv
    logical, intent(in) :: mask(nx,ny)
    real(kind=8), intent(out) :: output(npt)

    integer :: i,j,n,x,stx,spx,sty,spy
    real(kind=8) :: xarg(nx),yarg(ny)
    real(kind=8) :: splint,xp,yp
    real(kind=8) :: y2a(nx,ny),line(nw),y2line(nw)

    do i = 1,nx
      xarg(i) = i
    end do
    do j = 1,ny
      yarg(j) = j
    end do

    CALL extrap(input,mask,nx,ny,100,2)

    CALL splie2(yarg,nx,ny,input,y2a)
    do n = 1, npt
      if (dp(n).lt.0) then
      xp = real(idx(n,1),8)+W(n,1)
      yp = real(idx(n,2),8)+W(n,2)

      stx = idx(n,1)-int(nw/2)+1
      if (stx.lt.1) stx = 1
      spx = idx(n,1)+nw-int(nw/2)
      if (spx.gt.nx) spx = nx

      sty = idx(n,2)-int(nw/2)+1
      if (sty.lt.1) sty = 1
      spy = idx(n,2)+nw-int(nw/2)
      if (spy.gt.ny) spy = ny

      do x = stx, spx
        line(x+1-stx) = splint(yarg(sty:spy),input(x,sty:spy),y2a(x,sty:spy),nw,yp)
      end do
      call spline(xarg(stx:spx),line,nw,y2line)
      output(n) = splint(xarg(stx:spx),line,y2line,nw,xp)
      else
        output(n) = fv
      end if
    end do
END SUBROUTINE interp

SUBROUTINE interp_linear(nx,ny,npt,input,idx,W,output)
    implicit none
    integer, intent(in) :: nx,ny,npt,idx(npt,2)
    real(kind=8), intent(in) :: input(nx,ny),W(npt,2)
    real(kind=8), intent(out) :: output(npt)

    integer :: n
    real(kind=8) :: p11,p12,p21,p22,w11,w12,w21,w22

    do n = 1, npt
      p11 = input(idx(n,1),idx(n,2))
      p12 = input(idx(n,1),idx(n,2)+1)
      p21 = input(idx(n,1)+1,idx(n,2))
      p22 = input(idx(n,1)+1,idx(n,2)+1)
      w11 = (1-W(n,1))*(1-W(n,2))
      w12 = (1-W(n,1))*W(n,2)
      w21 = W(n,1)*(1-W(n,2))
      w22 = W(n,1)*W(n,2)
      output(n) = p11*w11+p12*w12+p21*w21+p22*w22
    end do
    
END SUBROUTINE interp_linear

FUNCTION splint(xa,ya,y2a,n,x)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    REAL(kind=8), DIMENSION(n),  INTENT(IN) :: xa,ya,y2a
    REAL(kind=8), INTENT(IN) ::  x

    INTEGER :: khi,klo,locate,jl,jm,ju
    REAL(kind=8) :: a,b,h,splint
    LOGICAL :: ascnd

    ascnd=(xa(n)>=xa(1))
    jl=0
    ju=n+1
    do
        if (ju-jl<=1) exit
        jm=(ju+jl)/2
        if (ascnd .eqv. (x>=xa(jm))) then
            jl=jm
        else
            ju=jm
        end if
    end do

    if (x==xa(1)) then
        locate=1
    else if (x==xa(n)) then
        locate=n-1
    else
        locate=jl
    end if

    klo=max(min(locate,n-1),1)
    khi=klo+1
    h=xa(khi)-xa(klo)
    if (h == 0.0) then
        write(*,*) 'bad xa input in splint'
        STOP
    end if
    a=(xa(khi)-x)/h
    b=(x-xa(klo))/h
    splint=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0
END FUNCTION splint

SUBROUTINE splie2(x2a,m,n,ya,y2a)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: m,n
    REAL(kind=8), DIMENSION(n), INTENT(IN) :: x2a
    REAL(kind=8), DIMENSION(m,n), INTENT(IN) :: ya
    REAL(kind=8), DIMENSION(m,n), INTENT(OUT) :: y2a

    INTEGER :: j

    DO j=1,m
        call spline(x2a,ya(j,:),n,y2a(j,:))
    END DO
END SUBROUTINE splie2

SUBROUTINE spline(x,y,n,y2)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    REAL(kind=8), DIMENSION(n), INTENT(IN) :: x,y
    REAL(kind=8), DIMENSION(n), INTENT(OUT) :: y2

    REAL(kind=8), DIMENSION(n) :: a,b,c,r

    c(1:n-1)=x(2:n)-x(1:n-1)
    r(1:n-1)=6.0*((y(2:n)-y(1:n-1))/c(1:n-1))
    r(2:n-1)=r(2:n-1)-r(1:n-2)
    a(2:n-1)=c(1:n-2)
    b(2:n-1)=2.0*(c(2:n-1)+a(2:n-1))
    b(1)=1.0
    b(n)=1.0
    r(1)=0.0
    c(1)=0.0
    r(n)=0.0
    a(n)=0.0
    call tridag(n,a(2:n),b(1:n),c(1:n-1),r(1:n),y2(1:n))
END SUBROUTINE spline

SUBROUTINE tridag(n,a,b,c,r,u)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    REAL(kind=8), DIMENSION(n), INTENT(IN) :: b,r
    REAL(kind=8), DIMENSION(n-1), INTENT(IN) :: a,c
    REAL(kind=8), DIMENSION(n), INTENT(OUT) :: u

    REAL(kind=8) :: gam(n),bet
    INTEGER :: j

    bet=b(1)

    if (bet==0.0) then
        write(*,*) 'tridag: Error at code stage 1'
        STOP
    end if

    u(1)=r(1)/bet
    do j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j-1)*gam(j)
        if (bet==0.0) then
            write(*,*) 'tridag: Error at code stage 2'
            STOP
        end if
        u(j)=(r(j)-a(j-1)*u(j-1))/bet
    end do
    do j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
    end do
END SUBROUTINE tridag

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

        do i=2,lon-1
          do j=2,lat-1
            res(i,j)=0.25*(a(i-1,j)+a(i+1,j)+a(i,j-1)+a(i,j+1))-a(i,j)
          end do
        end do

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

SUBROUTINE dist(a,b,d)
    implicit none
    real(kind=8), intent(in) :: a(2),b(2)
    real(kind=8), intent(out) :: d

    real(kind=8) :: x,c,r,pi

    r = 6371.0
    pi = 4.D0*DATAN(1.D0) 
    x = (sin(0.5*(b(1)*pi/180-a(1)*pi/180)))**2+cos(a(1)*pi/180)*&
        cos(b(1)*pi/180)*(sin(0.5*(b(2)*pi/180-a(2)*pi/180)))**2
    c = 2*atan2(sqrt(x),sqrt(1-x))
    d = r*c
    
END SUBROUTINE dist
