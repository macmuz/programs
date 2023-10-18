PROGRAM interp
    USE netcdf
    implicit none
    include 'mpif.h'
    integer :: my_task, size_Of_Cluster, ierror
    integer, parameter :: master_task = 0
    integer :: ncid,varid,dimids2(2),dimids3(3),dimids4(4)
    integer :: dim1, dim2, x_dimid, y_dimid, vid(3)
    integer :: c_dimid, d_dimid, i, j, x, y, m, n
    integer :: istr,iend,jstr,jend
    integer :: indim1, indim2
    integer, allocatable :: idx(:,:,:,:)
    real(kind=8), allocatable :: W(:,:,:)
    real(kind=8), allocatable :: inlon(:,:), inlat(:,:)
    real(kind=8), allocatable :: outlon(:,:), outlat(:,:)
    real(kind=8), allocatable :: dis_array(:,:)
    character(len=250) :: ingrid, outgrid
    real(kind=8) :: corners(4,2)
    real(kind=8) :: point(2),tmp
    logical :: score

    call MPI_INIT(ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,size_Of_Cluster,ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD,my_task,ierror)


    ingrid = '/users/work/mmuzyka/CSDIR/input_560x600/ROMS_grid_125NM_may.nc'
    outgrid = '/users/work/mmuzyka/CSDIR/input_560x600/ROMS_grid_125NM_may.nc'

    CALL check(nf90_open(trim(ingrid), NF90_NOWRITE, ncid),100)
    CALL check(nf90_inq_varid(ncid, 'lon_rho', varid),101)
    CALL check(nf90_inquire_variable(ncid, varid, dimids=dimids2),102)
    CALL check(nf90_inquire_dimension(ncid, dimids2(1), len=indim1),103)
    CALL check(nf90_inquire_dimension(ncid, dimids2(2), len=indim2),104)
    allocate( inlon(indim1,indim2), inlat(indim1,indim2) )
    allocate( dis_array(indim1,indim2) )
    CALL check(nf90_get_var(ncid, varid, inlon),105)
    CALL check(nf90_inq_varid(ncid, 'lat_rho', varid),106)
    CALL check(nf90_get_var(ncid, varid, inlat),107)


    CALL check(nf90_open(trim(outgrid), NF90_NOWRITE, ncid),110)
    CALL check(nf90_inq_varid(ncid, 'lon_rho', varid),111)
    CALL check(nf90_inquire_variable(ncid, varid, dimids=dimids2),112)
    CALL check(nf90_inquire_dimension(ncid, dimids2(1), len=dim1),113)
    CALL check(nf90_inquire_dimension(ncid, dimids2(2), len=dim2),114)
    allocate( outlon(dim1,dim2), outlat(dim1,dim2) )
    allocate( idx(dim1,dim2,4,2), W(dim1,dim2,4) )
    CALL check(nf90_get_var(ncid, varid, outlon),115)
    CALL check(nf90_inq_varid(ncid, 'lat_rho', varid),116)
    CALL check(nf90_get_var(ncid, varid, outlat),117)


    idx = 0
    W = 0
    W(:,:,1) = 1.0

    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    do j = 1, dim2
    if (mod(j,size_Of_Cluster).eq.my_task) then
    do i = 1, dim1
        dis_array(:,:) = sqrt((outlon(i,j)-inlon(:,:))**2+(outlat(i,j)-inlat(:,:))**2)
        idx(i,j,1,:) = minloc( dis_array )
        x = idx(i,j,1,1)
        y = idx(i,j,1,2)
        if (x.gt.1 .and. x.lt.indim1 .and. y.gt.1 .and. y.lt.indim2) then
          istr = x-1
          iend = x
          jstr = y-1
          jend = y
        elseif ( x.eq.indim1 .and. y.gt.1 .and. y.lt.indim2 ) then
          istr = x-1
          iend = istr
          jstr = y-1
          jend = y
        elseif ( x.eq.1 .and. y.gt.1 .and. y.lt.indim2 ) then
          istr = x
          iend = istr
          jstr = y-1
          jend = y
        elseif ( x.gt.1 .and. x.lt.indim1 .and. y.eq.indim2 ) then
          istr = x-1
          iend = x
          jstr = y-1
          jend = Jstr
        elseif ( x.gt.1 .and. x.lt.indim1 .and. y.eq.1 ) then
          istr = x-1
          iend = x
          jstr = y
          jend = Jstr
        elseif ( x.eq.1 .and. y.eq.1 ) then
          istr = x
          iend = istr
          jstr = y
          jend = jstr
        elseif ( x.eq.1 .and. y.eq. indim2 ) then
          istr = x
          iend = istr
          jstr = y-1
          jend = jstr 
        elseif ( x.eq.indim1 .and. y.eq. indim2 ) then
          istr = x-1
          iend = istr
          jstr = y-1
          jend = jstr 
        elseif ( x.eq.indim1 .and. y.eq. 1 ) then
          istr = x-1
          iend = istr
          jstr = y
          jend = jstr 
        end if
       
        point(1) = outlon(i,j)
        point(2) = outlat(i,j)
        do m = istr,iend
        do n = jstr,jend
          corners(1,:) = (/inlon(m,n),inlat(m,n)/)
          corners(4,:) = (/inlon(m,n+1),inlat(m,n+1)/)
          corners(3,:) = (/inlon(m+1,n+1),inlat(m+1,n+1)/)
          corners(2,:) = (/inlon(m+1,n),inlat(m+1,n)/)
          CALL in_convex_polygon(corners,point,score)
          if (score) then
            idx(i,j,1,:) = (/m,n/)
            idx(i,j,4,:) = (/m,n+1/)
            idx(i,j,3,:) = (/m+1,n+1/)
            idx(i,j,2,:) = (/m+1,n/)
!            CALL coef2(corners,point,W(i,j,:))
            CALL calc_w(corners,point,W(i,j,:))
            go to 1001
          end if
        end do
        end do 
1001    continue
    end do
    if (my_task.eq.master_task) write(*,*) j
    end if
    end do

    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    if (my_task==master_task) then

      do j = 1,dim2
        if (mod(j,size_Of_Cluster).ne.master_task) then
          call MPI_RECV(idx(:,j,:,:), 4*2*dim1, MPI_INTEGER, &
             mod(j,size_Of_Cluster), j, MPI_COMM_WORLD, &
             MPI_STATUS_IGNORE, ierror)
        end if
      end do

    else

      do j = 1,dim2
        if (mod(j,size_Of_Cluster).eq.my_task) then
          call MPI_SEND(idx(:,j,:,:), 4*2*dim1, MPI_INTEGER, &
               master_task, j, MPI_COMM_WORLD, ierror)
        end if
      end do

    end if

    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    if (my_task==master_task) then

      do j = 1,dim2
        if (mod(j,size_Of_Cluster).ne.master_task) then
          call MPI_RECV(W(:,j,:), 4*dim1, MPI_DOUBLE_PRECISION, &
             mod(j,size_Of_Cluster), j, MPI_COMM_WORLD, &
             MPI_STATUS_IGNORE, ierror)
        end if
      end do

    else

      do j = 1,dim2
        if (mod(j,size_Of_Cluster).eq.my_task) then
          call MPI_SEND(W(:,j,:), 4*dim1, MPI_DOUBLE_PRECISION, &
               master_task, j, MPI_COMM_WORLD, ierror)
        end if
      end do

    end if


    if (my_task==master_task) then
      CALL check(nf90_create( 'out_1_25NM.nc', NF90_CLOBBER, ncid ), 200)
      CALL check(nf90_def_dim( ncid, 'nx', dim1, x_dimid ), 201)
      CALL check(nf90_def_dim( ncid, 'ny', dim2, y_dimid ), 202)
      CALL check(nf90_def_dim( ncid, 'corners', 4, c_dimid ), 203)
      CALL check(nf90_def_dim( ncid, 'dim', 2, d_dimid ), 204)
      dimids2 = (/ x_dimid, y_dimid /)
      dimids3 = (/ x_dimid, y_dimid, c_dimid /)
      dimids4 = (/ x_dimid, y_dimid, c_dimid, d_dimid /)
      CALL check(nf90_def_var( ncid, 'idx', NF90_SHORT, dimids4, vid(1)), 205)
      CALL check(nf90_def_var( ncid, 'W', NF90_FLOAT, dimids3, vid(2)), 206)
      CALL check( nf90_enddef( ncid ), 207 )
      CALL check( nf90_put_var( ncid, vid(1), int(idx,2) ), 208 )
      CALL check( nf90_put_var( ncid, vid(2), real(W,4) ), 209 )
      CALL check(nf90_close( ncid ), 210 )
    end if
    

    deallocate(dis_array)
    deallocate(idx,W)
    deallocate(outlon,outlat)
    deallocate(inlon,inlat)
    call MPI_FINALIZE(ierror)
END PROGRAM interp

SUBROUTINE check(status,label)
    USE NETCDF
    implicit none
    INTEGER, INTENT(in) :: status , label

    IF (status .ne. nf90_noerr) then
        print*, trim(nf90_strerror(status)),'label=',label
        stop "Stopped"
    END IF
END SUBROUTINE check

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

SUBROUTINE coef(corners,point,W)
    implicit none
    real(kind=8), intent(in) :: corners(4,2)
    real(kind=8), intent(in) :: point(2)
    real(kind=8), intent(out) :: W(4)

    real(kind=8) :: a,b,c,d,e,f,g,h,alpha,beta

    a = -corners(1,1)+corners(4,1)
    b = -corners(1,1)+corners(2,1)
    c = corners(1,1)-corners(2,1)-corners(4,1)+corners(3,1)
    b = point(1)-corners(1,1)
    e = -corners(1,2)+corners(4,2)
    f = -corners(1,2)+corners(2,2)
    g = corners(1,2)-corners(2,2)-corners(4,2)+corners(3,2)
    h = point(2)-corners(1,2)

    alpha = (-1)*(b*e-a*f+d*g-c*h+sqrt(-4*(c*e-a*g)*(d*f-b*h)+&
            (b*e-a*f+d*g-c*h)**2))/(2*c*e-2*a*g)
    beta = (b*e-a*f-d*g+c*h+sqrt(-4*(c*e-a*g)*(d*f-b*h)+&
            (b*e-a*f+d*g-c*h)**2))/(2*c*f-2*b*g)

    write(*,*) alpha,beta

    if (alpha.lt.0 .or. alpha.gt.1 .or. &
        beta.lt.0 .or. beta.gt.1 ) then
      alpha = ((-1)*b*e+a*f-d*g+c*h+sqrt(-4*(c*e-a*g)*(d*f-b*h)+&
              (b*e-a*f+d*g-c*h)**2))/(2*c*e-2*a*g)
      beta = (-1)*((-b*e+a*f+d*g-c*h+sqrt(-4*(c*e-a*g)*(d*f-b*h)+&
              (b*e-a*f+d*g-c*h)**2))/(2*c*f-2*b*g))
      write(*,*) alpha,beta
    end if

    W(1) = (1-alpha)*(1-beta)
    W(2) = (1-alpha)*beta
    W(3) = alpha*(1-beta)
    W(4) = alpha*beta
END SUBROUTINE coef

SUBROUTINE coef2(corners,point,W)
    implicit none
    real(kind=8), intent(in) :: corners(4,2)
    real(kind=8), intent(in) :: point(2)
    real(kind=8), intent(out) :: W(4)

    real(kind=8) :: A,B,C,t,s
    real(kind=8), dimension(2) :: p1,p2,p3,p4

    p1 = (/corners(2,1),corners(2,2)/)
    p2 = (/corners(3,1),corners(3,2)/)
    p3 = (/corners(1,1),corners(1,2)/)
    p4 = (/corners(4,1),corners(4,2)/)

    A = (p3(1)-p1(1))*(p4(2)-p2(2))-(p3(2)-p1(2))*(p4(1)-p2(1))
    B = point(2)*((p4(1)-p2(1))-(p3(1)-p1(1)))-&
        point(1)*((p4(2)-p2(2))-(p3(2)-p1(2)))+&
        (p3(1)-p1(1))*p2(2)-(p3(2)-p1(2))*p2(1)+&
        (p4(2)-p2(2))*p1(1)-(p4(1)-p2(1))*p1(2)
    C = point(2)*(p2(1)-p1(1))-point(1)*(p2(2)-p1(2))+&
        p1(1)*p2(2)-p1(2)*p2(1)

    t = (-B+sqrt(B**2-4*A*C))/(2*A)
    if ((t<0).OR.(t>1)) t=(-B-sqrt(B**2-4*A*C))/(2*A)
    s = (point(2)-p1(2)-(p3(2)-p1(2))*t)/(p2(2)+(p4(2)-p2(2))*t-&
        p1(2)-(p3(2)-p1(2))*t)

    W(1) = (1-s)*t
    W(2) = (1-s)*(1-t)
    W(3) = s*(1-t)
    W(4) = s*t
END SUBROUTINE coef2

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
