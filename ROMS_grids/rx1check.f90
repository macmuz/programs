PROGRAM rx1check
    USE netcdf
    implicit none
    integer, parameter :: nx=2700,ny=3200,fh=15
    integer :: ios=0, pos, N, narg, ncid, varid, i, j, k
    integer, allocatable :: mask(:,:)
    real(kind=8) :: hc,theta_s,theta_b,ds,rx1
    real(kind=8), allocatable :: h(:,:), depth(:,:,:), s_w(:), cs_w(:)
    character(len=50) :: buffer, filename, inpfile, label 

!############# command line arguments #################################
    narg = command_argument_count()
    if (narg.ne.2) then
        write(*,*) "Program must have 2 arguments: &
                    &filename and inputfile"
        stop
    end if
    call get_command_argument(1,buffer)
    read(buffer, *) filename
    call get_command_argument(2,buffer)
    read(buffer, *) inpfile
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
    
    CALL check(nf90_open( trim(filename), nf90_nowrite, ncid ), 200)
    CALL check(nf90_inq_varid( ncid, "mask_rho", varid ),201)
    CALL check(nf90_get_var (ncid, varid, mask ),202)
    CALL check(nf90_inq_varid( ncid, "h", varid ),203)
    CALL check(nf90_get_var (ncid, varid, h ),204)
    CALL check(nf90_close( ncid ), 205 )

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

    deallocate( s_w, cs_w )
!######################################################################

    CALL calc_rx1(nx,ny,N,mask,depth,rx1)
    write(*,*) 'rx1=',rx1

    deallocate( mask, h, depth )
END PROGRAM rx1check

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
