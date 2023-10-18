PROGRAM meo_scheme
    USE netcdf
    implicit none
    integer, parameter :: nx=1370,ny=1460,fh=15
    integer :: ios=0, pos, N, narg, ncid, varid, i, j, k
    integer, allocatable :: mask(:,:)
    logical, allocatable :: masklog(:,:), mask2log(:,:)
    real(kind=8) :: hc,theta_s,theta_b,ds,rx1,rx0,trx0
    real(kind=8), allocatable :: h(:,:), depth(:,:,:), s_w(:), cs_w(:)
    real(kind=8), allocatable :: pm(:,:), pn(:,:), A(:,:)
    character(len=50) :: buffer, inpfile, label 
    character(len=150) :: buffer2,filename 

!############# command line arguments #################################
    narg = command_argument_count()
    if (narg.ne.3) then
        write(*,*) "Program must have 3 arguments: &
                    &filename, inputfile and trx0"
        stop
    end if
    call get_command_argument(1,buffer2)
    write(filename, *) buffer2
    call get_command_argument(2,buffer)
    read(buffer, *) inpfile
    call get_command_argument(3,buffer)
    read(buffer, *) trx0

    filename = adjustl(filename)
    write(*,*) trim(filename), trx0
!######################################################################

!##########parameters#reading#########
    open(fh, file=trim(inpfile))

    do while (ios == 0)
     read(fh, '(A)', iostat=ios) buffer
     if (ios == 0) then

        ! Find the first instance of whitespace.  Split label and data.
        pos = scan(buffer, '    ')
        label = buffer(1:pos)
        buffer = buffer(pos+1:)

        select case (label)
        case ('N')
           read(buffer, *, iostat=ios) N
        case ('hc')
           read(buffer, *, iostat=ios) hc
        case ('theta_s')
           read(buffer, *, iostat=ios) theta_s
        case ('theta_b')
           read(buffer, *, iostat=ios) theta_b
        end select
     end if
    end do

    close(fh)
    print *, 'Read N: ', N
    print *, 'Read hc: ', hc
    print *, 'Read theta_s: ', theta_s
    print *, 'Read theta_b: ', theta_b
!#####################################
    
    allocate( mask(nx,ny), h(nx,ny), depth(nx,ny,N+1) )
    allocate( masklog(nx,ny), mask2log(nx,ny) )
    allocate( pm(nx,ny), pn(nx,ny), A(nx,ny) )
    
    CALL check(nf90_open( trim(filename), nf90_nowrite, ncid ), 200)
    CALL check(nf90_inq_varid( ncid, "mask_rho", varid ),201)
    CALL check(nf90_get_var (ncid, varid, mask ),202)
    CALL check(nf90_inq_varid( ncid, "hraw", varid ),203)
    CALL check(nf90_get_var (ncid, varid, h ),204)
    CALL check(nf90_inq_varid(ncid,"pm",varid),7)
    CALL check(nf90_get_var(ncid,varid,pm),8)
    CALL check(nf90_inq_varid(ncid,"pn",varid),7)
    CALL check(nf90_get_var(ncid,varid,pn),8)
    CALL check(nf90_close( ncid ), 205 )

    A(:,:) = 1/(pm(:,:)*pn(:,:))

!########## calc depth ################################################
    allocate( s_w(N+1), cs_w(N+1) )

    depth = 0
    s_w(1) = -1
    s_w(N+1) = 0
    ds = 1.0/real(N,8)

    do i = 2,N
      s_w(i) = s_w(i-1)+ds
    end do

    if (theta_s .gt. 0) then
      do i = 1,N+1
        cs_w(i) = (1-cosh(theta_s*s_w(i)))/(cosh(theta_s)-1)
      end do
    else
      do i = 1,N+1
        cs_w(i) = -1*(s_w(i)**2)
      end do
    end if

    if (theta_b .gt. 0) then
      do i = 1,N+1
        cs_w(i) = (exp(theta_b*cs_w(i))-1)/(1-exp(-1*theta_b))
      end do
    end if

    do i = 1,nx
    do j = 1,ny
      if ( mask(i,j).eq.1 ) then
        do k = 1,N+1
          depth(i,j,k) = h(i,j)*( hc*s_w(k)+h(i,j)*cs_w(k) )/( hc+h(i,j) )
        end do
      end if
    end do
    end do

!######################################################################

    CALL calc_rx1(nx,ny,N,mask,depth,rx1)
    write(*,*) 'rx1=',rx1

    masklog = .false.
    where( mask.eq.1 ) masklog = .true.
    CALL calc_rx0(nx,ny,masklog,h,rx0)
    write(*,*) 'rx0=',rx0

    CALL meo(nx,ny,masklog,h,A,trx0)
!    mask2log = .false.
!    where(masklog(60:75,588:608)) mask2log(60:75,588:608) = .true.
!    CALL mb_scheme(nx,ny,mask2log,h,trx0)


    do i = 1,nx
    do j = 1,ny
      if ( mask(i,j).eq.1 ) then
        do k = 1,N+1
          depth(i,j,k) = h(i,j)*( hc*s_w(k)+h(i,j)*cs_w(k) )/( hc+h(i,j) )
        end do
      end if
    end do
    end do
    CALL calc_rx1(nx,ny,N,mask,depth,rx1)
    write(*,*) 'rx1=',rx1
    
    CALL check(nf90_open(trim(filename),NF90_WRITE,ncid),10)
    CALL check(nf90_inq_varid(ncid,"h",varid),11)
    CALL check(nf90_put_var(ncid,varid,h),12)
    CALL check(nf90_close(ncid),15)

    deallocate( s_w, cs_w )
    deallocate( mask, h, depth, masklog )
    deallocate( pm, pn, A, mask2log )
END PROGRAM meo_scheme

SUBROUTINE check(status,label)
    USE NETCDF
    implicit none
    INTEGER, INTENT(in) :: status , label

    IF (status .ne. nf90_noerr) then
        print*, trim(nf90_strerror(status)),'label=',label
        stop "Stopped"
    END IF
END SUBROUTINE check

SUBROUTINE calc_rx1(nx,ny,nz,mask,depth,rx1max)
    implicit none
    integer, intent(in) :: nx,ny,nz,mask(nx,ny)
    real(kind=8), intent(in) :: depth(nx,ny,nz+1)
    real(kind=8), intent(out) :: rx1max

    integer :: i,j,k
    real(kind=8) :: rx1_tmp(nz)
    real(kind=8) :: rx1(nx,ny)

    rx1 = 0

    do i = 2,nx-1
    do j = 2,ny-1
      if (mask(i,j).eq.1) then
        if(mask(i+1,j).eq.1) then
        
          do k = 2,nz+1
            rx1_tmp(k-1) = &
    abs( depth(i,j,k)-depth(i+1,j,k)+depth(i,j,k-1)-depth(i+1,j,k-1) )/&
       ( depth(i,j,k)+depth(i+1,j,k)-depth(i,j,k-1)-depth(i+1,j,k-1) )
          end do
          rx1(i,j) = max(maxval(rx1_tmp),rx1(i,j))

        end if
        if(mask(i-1,j).eq.1) then

          do k = 2,nz+1
            rx1_tmp(k-1) = &
    abs( depth(i,j,k)-depth(i-1,j,k)+depth(i,j,k-1)-depth(i-1,j,k-1) )/&
       ( depth(i,j,k)+depth(i-1,j,k)-depth(i,j,k-1)-depth(i-1,j,k-1) )
          end do
          rx1(i,j) = max(maxval(rx1_tmp),rx1(i,j))

        end if
        if(mask(i,j+1).eq.1) then

          do k = 2,nz+1
            rx1_tmp(k-1) = &
    abs( depth(i,j,k)-depth(i,j+1,k)+depth(i,j,k-1)-depth(i,j+1,k-1) )/&
       ( depth(i,j,k)+depth(i,j+1,k)-depth(i,j,k-1)-depth(i,j+1,k-1) )
          end do
          rx1(i,j) = max(maxval(rx1_tmp),rx1(i,j))

        end if
        if(mask(i,j-1).eq.1) then

          do k = 2,nz+1
            rx1_tmp(k-1) = &
    abs( depth(i,j,k)-depth(i,j-1,k)+depth(i,j,k-1)-depth(i,j-1,k-1) )/&
       ( depth(i,j,k)+depth(i,j-1,k)-depth(i,j,k-1)-depth(i,j-1,k-1) )
          end do
          rx1(i,j) = max(maxval(rx1_tmp),rx1(i,j))

        end if
      end if
    end do
    end do
    rx1max = maxval(rx1)
ENDSUBROUTINE calc_rx1

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

SUBROUTINE meo(nx,ny,mask,h,A,trx0)
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
END SUBROUTINE meo

FUNCTION calcv1(h1,h2,a1,a2,r) result(v)
    implicit none
    real(kind=8), intent(in) :: h1,h2,a1,a2,r
    real(kind=8) :: v

    v = a1*a2*(h1-h2-r*h1-r*h2)/(a1*r+a1-a2*r+a2)
END FUNCTION calcv1

SUBROUTINE mb_scheme(nx,ny,mask,h,target_rx0)
    implicit none
    integer, intent(in) :: nx, ny
    logical, dimension(nx,ny), intent(in) :: mask
    real(kind=8), dimension(nx,ny), intent(inout) :: h
    real(kind=8), intent(in) :: target_rx0

    integer :: i,j
    real(kind=8) :: rx0

    do
      CALL calc_rx0(nx,ny,mask,h,rx0)

      write(*,*) "Next loop",rx0

      if ( rx0.le.(target_rx0+0.001) ) EXIT

      do j = 2,ny-1
      do i = 2,nx-1
      if (mask(i,j)) then
        if(mask(i+1,j)) then
          if ( (h(i,j)-h(i+1,j))/(h(i,j)+h(i+1,j)).gt.target_rx0 ) then
            h(i+1,j) = h(i,j)*(1-target_rx0)/(1+target_rx0)
          end if
        end if
        if(mask(i-1,j)) then
          if ( (h(i,j)-h(i-1,j))/(h(i,j)+h(i-1,j)).gt.target_rx0 ) then
            h(i-1,j) = h(i,j)*(1-target_rx0)/(1+target_rx0)
          end if
        end if
        if(mask(i,j+1)) then
          if ( (h(i,j)-h(i,j+1))/(h(i,j)+h(i,j+1)).gt.target_rx0 ) then
            h(i,j+1) = h(i,j)*(1-target_rx0)/(1+target_rx0)
          end if
        end if
        if(mask(i,j-1)) then
          if ( (h(i,j)-h(i,j-1))/(h(i,j)+h(i,j-1)).gt.target_rx0 ) then
            h(i,j-1) = h(i,j)*(1-target_rx0)/(1+target_rx0)
          end if
        end if
      end if
      end do
      end do
    end do
END SUBROUTINE mb_scheme
