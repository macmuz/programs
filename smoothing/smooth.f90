PROGRAM smooth
    USE netcdf
    implicit none 
    include 'mpif.h'

    character(len=150) :: inputfile,paramfile
    character(len=100) :: buffer,label
    real(kind=8), allocatable :: mask(:,:),hraw(:,:)
    real(kind=8), allocatable :: pm(:,:),pn(:,:),invA(:,:)
    real(kind=8) :: hc,theta_s,theta_b,target_rx0,target_rx1
    real(kind=8), allocatable :: s_w(:),cs_w(:),depth(:,:,:)
    real(kind=8) :: ds,dx,rx0_tmp
    real(kind=8), allocatable :: rx0(:,:),rx1(:,:)
    real(kind=8), allocatable :: mb_out(:,:),dc_out(:,:),meo_out(:,:)
    real(kind=8) :: tstart, tend
    integer :: ncid,varid,dim1,dim2,dimids2(2),vid(4),i,j,k
    integer :: x_dimid, y_dimid
    integer :: N,pos,ios = 0,line = 0
    integer :: my_task, size_Of_Cluster, ierror
    integer, parameter :: master_task = 0
    integer, parameter :: fh=15
    logical, allocatable :: wet(:,:)

    INTERFACE
      SUBROUTINE calc_rx0(mask,h,rx0,my_task,size_Of_Cluster,&
        ierror,master)
        logical, dimension(:,:) :: mask
        real(kind=8), dimension(:,:) :: h
        real(kind=8), dimension(:,:) :: rx0
        integer :: my_task,size_Of_Cluster,ierror,master
      END SUBROUTINE
    END INTERFACE 

    INTERFACE
      SUBROUTINE calc_rx1(mask,depths,rx1,my_task,size_Of_Cluster,&
        ierror,master)
        logical, dimension(:,:) :: mask
        real(kind=8), dimension(:,:,:) :: depths
        real(kind=8), dimension(:,:) :: rx1
        integer :: my_task,size_Of_Cluster,ierror,master
      END SUBROUTINE
    END INTERFACE 

    INTERFACE
      SUBROUTINE mb_scheme(mask,h,target_rx0,my_task,size_Of_Cluster,&
        ierror,master)
        logical, dimension(:,:) :: mask
        real(kind=8), dimension(:,:) :: h
        real(kind=8) :: target_rx0
        integer :: my_task,size_Of_Cluster,ierror,master
      END SUBROUTINE
    END INTERFACE 

    INTERFACE
      SUBROUTINE dc_scheme(mask,h,target_rx0,my_task,size_Of_Cluster,&
        ierror,master)
        logical, dimension(:,:) :: mask
        real(kind=8), dimension(:,:) :: h
        real(kind=8) :: target_rx0
        integer :: my_task,size_Of_Cluster,ierror,master
      END SUBROUTINE
    END INTERFACE 

    call MPI_INIT(ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,size_Of_Cluster,ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD,my_task,ierror)

!##################################
    inputfile = 'ROMS_linear_10m_filter.nc'
    paramfile = 'parameters.in'
    
!##################################

    CALL check(nf90_open(trim(inputfile), NF90_NOWRITE, ncid),100)
    CALL check(nf90_inq_varid(ncid, 'mask_rho', varid),101)
    CALL check(nf90_inquire_variable(ncid, varid, dimids=dimids2),102)
    CALL check(nf90_inquire_dimension(ncid, dimids2(1), len=dim1),103)
    CALL check(nf90_inquire_dimension(ncid, dimids2(2), len=dim2),104)
    allocate( mask(dim1,dim2), hraw(dim1,dim2) )
    allocate( pm(dim1,dim2), pn(dim1,dim2), invA(dim1,dim2) )
    CALL check(nf90_get_var(ncid, varid, mask),105)
    CALL check(nf90_inq_varid(ncid, 'hraw', varid),106)
    CALL check(nf90_get_var(ncid, varid, hraw),107)
    CALL check(nf90_inq_varid(ncid, 'pm', varid),108)
    CALL check(nf90_get_var(ncid, varid, pm),109)
    CALL check(nf90_inq_varid(ncid, 'pn', varid),110)
    CALL check(nf90_get_var(ncid, varid, pn),111)
    CALL check(nf90_close(ncid),112)

    allocate( wet(dim1,dim2) )
    wet = .false.

    where(mask.gt.0.5) wet = .true.
    deallocate( mask )

    FORALL (i=1:dim1,j=1:dim2) invA(i,j) = pm(i,j)*pn(i,j)
!##########parameters#reading#########
    open(fh, file=trim(paramfile))
   
    do while (ios == 0)
     read(fh, '(A)', iostat=ios) buffer
     if (ios == 0) then
        line = line + 1

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
        case ('target_rx0')
           read(buffer, *, iostat=ios) target_rx0
        case ('target_rx1')
           read(buffer, *, iostat=ios) target_rx1
        end select
     end if
    end do
 
    close(fh)
    if (my_task==master_task) then
      print *, 'Read N: ', N
      print *, 'Read hc: ', hc
      print *, 'Read theta_s: ', theta_s
      print *, 'Read theta_b: ', theta_b
      print *, 'Read target_rx0: ', target_rx0
      print *, 'Read target_rx1: ', target_rx1
    end if
!#####################################

    allocate( s_w(N+1), cs_w(N+1), depth(dim1,dim2,N+1) )
    
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

    if (my_task==master_task) then
      do i = 1,N+1
        write(*,*) s_w(i),cs_w(i)
      end do
    end if

    do i = 1,dim1
    do j = 1,dim2
      if (wet(i,j)) then
        do k = 1,N+1
          depth(i,j,k) = hraw(i,j)*( hc*s_w(k)+hraw(i,j)*cs_w(k) )/( hc+hraw(i,j) )
        end do
      end if
    end do
    end do
   
!#########################################

    allocate( rx0(dim1,dim2), rx1(dim1,dim2) )
    allocate( mb_out(dim1,dim2), dc_out(dim1,dim2), meo_out(dim1,dim2) )

    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    tstart  = MPI_Wtime()
  
    CALL calc_rx0(wet,hraw,rx0,my_task,size_Of_Cluster,ierror,master_task)
    if (my_task==master_task) then
        write(*,*) 'rx0=',maxval(rx0)
    end if

    call MPI_BARRIER(MPI_COMM_WORLD,ierror)

    CALL calc_rx1(wet,depth,rx1,my_task,size_Of_Cluster,ierror,master_task)
    if (my_task==master_task) then
        write(*,*) 'rx1=',maxval(rx1)
    end if


    mb_out = hraw
    dc_out = hraw
    meo_out = hraw


    dx = 0.01
    rx0_tmp = target_rx0
    do
        do i = 1,dim1
        do j = 1,dim2
          if (wet(i,j)) then
            do k = 1,N+1
                  depth(i,j,k) = mb_out(i,j)*( hc*s_w(k)+mb_out(i,j)*&
                    cs_w(k) )/( hc+mb_out(i,j) )
            end do
          end if
        end do
        end do
      CALL calc_rx1(wet,depth,rx1,my_task,size_Of_Cluster,ierror,master_task)

      do j = 1,size(mask,2)
          CALL MPI_BCAST(rx1(:,j), size(mask,1), MPI_DOUBLE_PRECISION,&
                 master_task, MPI_COMM_WORLD, ierror)
      end do

      if (my_task==master_task) write(*,*) "rx1=",maxval(rx1)
      if (maxval(rx1).le.(target_rx1+0.001) .or. rx0_tmp.le.dx) EXIT

      if (my_task==master_task) write(*,*) "rx0_target:",rx0_tmp
      mb_out = hraw
      CALL mb_scheme(wet,mb_out,rx0_tmp,my_task,size_Of_Cluster,&
            ierror,master_task)
      rx0_tmp = rx0_tmp-dx

    end do

    CALL meo_scheme(wet,meo_out,invA,target_rx0)
!    do
!        do i = 1,dim1
!        do j = 1,dim2
!          if (wet(i,j)) then
!            do k = 1,N+1
!                  depth(i,j,k) = dc_out(i,j)*( hc*s_w(k)+dc_out(i,j)*&
!                    cs_w(k) )/( hc+dc_out(i,j) )
!            end do
!          end if
!        end do
!        end do
!      CALL calc_rx1(wet,depth,rx1,my_task,size_Of_Cluster,ierror,master_task)
!
!      do j = 1,size(mask,2)
!          CALL MPI_BCAST(rx1(:,j), size(mask,1), MPI_DOUBLE_PRECISION,&
!                 master_task, MPI_COMM_WORLD, ierror)
!      end do

!      if (my_task==master_task) write(*,*) "rx1=",maxval(rx1)
!      if (maxval(rx1).le.(target_rx1+0.001) .or. rx0_tmp.le.dx) EXIT

!      if (my_task==master_task) write(*,*) "rx0_target:",rx0_tmp
!      dc_out = hraw
!      CALL dc_scheme(wet,dc_out,rx0_tmp,my_task,size_Of_Cluster,&
!            ierror,master_task)
!      rx0_tmp = rx0_tmp-dx

!    end do



    if (my_task==master_task) then
      CALL check(nf90_create( 'rx_out_10m.nc', NF90_CLOBBER, ncid ), 200)

      CALL check(nf90_def_dim( ncid, 'x', dim1, x_dimid ), 201)
      CALL check(nf90_def_dim( ncid, 'y', dim2, y_dimid ), 202)
      dimids2 = (/ x_dimid, y_dimid /)
      CALL check(nf90_def_var( ncid, 'rx0', NF90_DOUBLE, dimids2, vid(1)), 203)
      CALL check(nf90_def_var( ncid, 'rx1', NF90_DOUBLE, dimids2, vid(2)), 204)
      CALL check(nf90_def_var( ncid, 'mb', NF90_DOUBLE, dimids2, vid(3)), 205)
      CALL check(nf90_def_var( ncid, 'meo', NF90_DOUBLE, dimids2, vid(4)), 206)

      CALL check( nf90_enddef( ncid ), 207 )
      CALL check( nf90_put_var( ncid, vid(1), rx0 ), 208)
      CALL check( nf90_put_var( ncid, vid(2), rx1 ), 209)
      CALL check( nf90_put_var( ncid, vid(3), mb_out ), 210)
      CALL check( nf90_put_var( ncid, vid(4), meo_out ), 211)
      CALL check(nf90_close( ncid ), 212)
    end if


    deallocate( rx0, rx1 )
    deallocate( mb_out, dc_out, meo_out )
    deallocate( hraw, wet )
    deallocate( pm, pn, invA )
    deallocate( s_w, cs_w, depth)
   
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    tend  = MPI_Wtime()
 
    call MPI_FINALIZE(ierror)

    if (my_task==master_task) then
      write(*,*) "Runtime = ",tend-tstart
    end if

END PROGRAM smooth

SUBROUTINE calc_rx0(mask,h,rx0,my_task,size_Of_Cluster,ierror,master)
    implicit none
    include 'mpif.h'
    logical, dimension(:,:), intent(in) :: mask
    real(kind=8), dimension(:,:), intent(in) :: h
    real(kind=8), dimension(:,:), intent(out) :: rx0
    integer, intent(in) :: my_task,size_Of_Cluster,ierror,master
   
    integer :: i,j 
    real(kind=8) :: rx0_tmp

    rx0 = 0

    do i = 2,size(mask,1)-1
    if (mod(i,size_Of_Cluster).eq.my_task) then
    do j = 2,size(mask,2)-1
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
    end if
    end do

    if (my_task==master) then
      do i = 2,size(mask,1)-1
        if (mod(i,size_Of_Cluster).ne.master) then
            call MPI_RECV(rx0(i,:), size(mask,2), MPI_DOUBLE_PRECISION, &
             mod(i,size_Of_Cluster), i, MPI_COMM_WORLD, &
             MPI_STATUS_IGNORE, ierror)
        end if
      end do
    else
      do i = 2,size(mask,1)-1
        if (mod(i,size_Of_Cluster).eq.my_task) then
            call MPI_SEND(rx0(i,:), size(mask,2), MPI_DOUBLE_PRECISION, &
               master, i, MPI_COMM_WORLD, ierror)
        end if
      end do
    end if
END SUBROUTINE calc_rx0

SUBROUTINE calc_rx1(mask,depths,rx1,my_task,size_Of_Cluster,ierror,master)
    implicit none
    include 'mpif.h'
    logical, dimension(:,:), intent(in) :: mask
    real(kind=8), dimension(:,:,:), intent(in) :: depths
    real(kind=8), dimension(:,:), intent(out) :: rx1
    integer, intent(in) :: my_task,size_Of_Cluster,ierror,master

    integer :: i,j,k
    real(kind=8) :: rx1_tmp(size(depths,3)-1)

    rx1 = 0

    do i = 2,size(mask,1)-1
    if (mod(i,size_Of_Cluster).eq.my_task) then
    do j = 2,size(mask,2)-1
      if (mask(i,j)) then
        if(mask(i+1,j)) then
          do k = 2,nz
            rx1_tmp(k-1) = &
    abs( depths(i,j,k)-depths(i+1,j,k)+depths(i,j,k-1)-depths(i+1,j,k-1) )/&
       ( depths(i,j,k)+depths(i+1,j,k)-depths(i,j,k-1)-depths(i+1,j,k-1) )
          end do
          rx1(i,j) = max(maxval(rx1_tmp),rx1(i,j))
        end if
        if(mask(i-1,j)) then
          do k = 2,nz
            rx1_tmp(k-1) = &
    abs( depths(i,j,k)-depths(i-1,j,k)+depths(i,j,k-1)-depths(i-1,j,k-1) )/&
       ( depths(i,j,k)+depths(i-1,j,k)-depths(i,j,k-1)-depths(i-1,j,k-1) )
          end do
          rx1(i,j) = max(maxval(rx1_tmp),rx1(i,j))
        end if
        if(mask(i,j+1)) then
          do k = 2,nz
            rx1_tmp(k-1) = &
    abs( depths(i,j,k)-depths(i,j+1,k)+depths(i,j,k-1)-depths(i,j+1,k-1) )/&
       ( depths(i,j,k)+depths(i,j+1,k)-depths(i,j,k-1)-depths(i,j+1,k-1) )
          end do
          rx1(i,j) = max(maxval(rx1_tmp),rx1(i,j))
        end if
        if(mask(i,j-1)) then
          do k = 2,nz
            rx1_tmp(k-1) = &
    abs( depths(i,j,k)-depths(i,j-1,k)+depths(i,j,k-1)-depths(i,j-1,k-1) )/&
       ( depths(i,j,k)+depths(i,j-1,k)-depths(i,j,k-1)-depths(i,j-1,k-1) )
          end do
          rx1(i,j) = max(maxval(rx1_tmp),rx1(i,j))
        end if
      end if
    end do
    end if
    end do

    if (my_task==master) then
      do i = 2,size(mask,1)-1
        if (mod(i,size_Of_Cluster).ne.master) then
            call MPI_RECV(rx1(i,:), size(mask,2), MPI_DOUBLE_PRECISION, &
             mod(i,size_Of_Cluster), i, MPI_COMM_WORLD, &
             MPI_STATUS_IGNORE, ierror)
        end if
      end do
    else
      do i = 2,size(mask,1)-1
        if (mod(i,size_Of_Cluster).eq.my_task) then
            call MPI_SEND(rx1(i,:), size(mask,2), MPI_DOUBLE_PRECISION, &
               master, i, MPI_COMM_WORLD, ierror)
        end if
      end do
    end if
END SUBROUTINE calc_rx1

SUBROUTINE mb_scheme(mask,h,target_rx0,my_task,size_Of_Cluster,ierror,master)
    implicit none
    include 'mpif.h'
    logical, dimension(:,:), intent(in) :: mask
    real(kind=8), dimension(:,:), intent(inout) :: h
    real(kind=8), intent(in) :: target_rx0
    integer, intent(in) :: my_task,size_Of_Cluster,ierror,master

    integer :: i,j
    real(kind=8), dimension(:,:), allocatable :: rx0

    INTERFACE
      SUBROUTINE calc_rx0(mask,h,rx0,my_task,size_Of_Cluster,ierror,master)
        logical, dimension(:,:) :: mask
        real(kind=8), dimension(:,:) :: h
        real(kind=8), dimension(:,:) :: rx0
        integer :: my_task,size_Of_Cluster,ierror,master
      END SUBROUTINE
    END INTERFACE

    allocate( rx0(size(mask,1),size(mask,2)) )

    do
      CALL calc_rx0(mask,h,rx0,my_task,size_Of_Cluster,ierror,master)

      do j = 1,size(mask,2)
          CALL MPI_BCAST(rx0(:,j), size(mask,1), MPI_DOUBLE_PRECISION,&
                 master, MPI_COMM_WORLD, ierror)
      end do

      if (my_task==master) write(*,*) "Next loop",maxval(rx0)

      if ( maxval(rx0).le.(target_rx0+0.001) ) EXIT

      do j = 2,size(mask,2)-1
      if (mod(j,size_Of_Cluster).eq.my_task) then
      do i = 2,size(mask,1)-1
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
      end if
      end do
      end if
      end do

    if (my_task==master) then
      do j = 2,size(mask,2)-1
        if (mod(j,size_Of_Cluster).ne.master) then
            call MPI_RECV(h(:,j), size(mask,1), MPI_DOUBLE_PRECISION, &
             mod(j,size_Of_Cluster), j, MPI_COMM_WORLD, &
             MPI_STATUS_IGNORE, ierror)
        end if
      end do
    else
      do j = 2,size(mask,2)-1
        if (mod(j,size_Of_Cluster).eq.my_task) then
            call MPI_SEND(h(:,j), size(mask,1), MPI_DOUBLE_PRECISION, &
               master, j, MPI_COMM_WORLD, ierror)
        end if
      end do
    end if

      do j = 2,size(mask,2)-1
          CALL MPI_BCAST(h(:,j), size(mask,1), MPI_DOUBLE_PRECISION,&
                 master, MPI_COMM_WORLD, ierror)
      end do

      do i = 2,size(mask,1)-1
      if (mod(i,size_Of_Cluster).eq.my_task) then
      do j = 2,size(mask,2)-1
      if (mask(i,j)) then
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
      end if
      end do

    if (my_task==master) then
      do i = 2,size(mask,1)-1
        if (mod(i,size_Of_Cluster).ne.master) then
            call MPI_RECV(h(i,:), size(mask,2), MPI_DOUBLE_PRECISION, &
             mod(i,size_Of_Cluster), i, MPI_COMM_WORLD, &
             MPI_STATUS_IGNORE, ierror)
        end if
      end do
    else
      do i = 2,size(mask,1)-1
        if (mod(i,size_Of_Cluster).eq.my_task) then
            call MPI_SEND(h(i,:), size(mask,2), MPI_DOUBLE_PRECISION, &
               master, i, MPI_COMM_WORLD, ierror)
        end if
      end do
    end if

      do i = 2,size(mask,1)-1
          CALL MPI_BCAST(h(i,:), size(mask,2), MPI_DOUBLE_PRECISION,&
                 master, MPI_COMM_WORLD, ierror)
      end do

    end do

    deallocate( rx0 )
END SUBROUTINE mb_scheme


SUBROUTINE dc_scheme(mask,h,target_rx0,my_task,size_Of_Cluster,ierror,master)
    implicit none
    include 'mpif.h'
    logical, dimension(:,:), intent(in) :: mask
    real(kind=8), dimension(:,:), intent(inout) :: h
    real(kind=8), intent(in) :: target_rx0
    integer, intent(in) :: my_task,size_Of_Cluster,ierror,master

    integer :: i,j
    real(kind=8), dimension(:,:), allocatable :: rx0

    INTERFACE
      SUBROUTINE calc_rx0(mask,h,rx0,my_task,size_Of_Cluster,ierror,master)
        logical, dimension(:,:) :: mask
        real(kind=8), dimension(:,:) :: h
        real(kind=8), dimension(:,:) :: rx0
        integer :: my_task,size_Of_Cluster,ierror,master
      END SUBROUTINE
    END INTERFACE

    allocate( rx0(size(mask,1),size(mask,2)) )

    do
      CALL calc_rx0(mask,h,rx0,my_task,size_Of_Cluster,ierror,master)

      do j = 1,size(mask,2)
          CALL MPI_BCAST(rx0(:,j), size(mask,1), MPI_DOUBLE_PRECISION,&
                 master, MPI_COMM_WORLD, ierror)
      end do

      if (my_task==master) write(*,*) "Next loop",maxval(rx0)

      if ( maxval(rx0).le.(target_rx0+0.001) ) EXIT

      do j = 2,size(mask,2)-1
      if (mod(j,size_Of_Cluster).eq.my_task) then
      do i = 2,size(mask,1)-1
      if (mask(i,j)) then
        if(mask(i+1,j)) then
          if ( (h(i,j)-h(i+1,j))/(h(i,j)+h(i+1,j)).gt.target_rx0 ) then
            h(i,j) = h(i+1,j)*(1+target_rx0)/(1-target_rx0)
          end if
        end if
        if(mask(i-1,j)) then
          if ( (h(i,j)-h(i-1,j))/(h(i,j)+h(i-1,j)).gt.target_rx0 ) then
            h(i,j) = h(i-1,j)*(1+target_rx0)/(1-target_rx0)
          end if
        end if
        if(mask(i,j+1)) then
          if ( (h(i,j)-h(i,j+1))/(h(i,j)+h(i,j+1)).gt.target_rx0 ) then
            h(i,j) = h(i,j+1)*(1+target_rx0)/(1-target_rx0)
          end if
        end if
        if(mask(i,j-1)) then
          if ( (h(i,j)-h(i,j-1))/(h(i,j)+h(i,j-1)).gt.target_rx0 ) then
            h(i,j) = h(i,j-1)*(1+target_rx0)/(1-target_rx0)
          end if
        end if
      end if
      end do
      end if
      end do

    if (my_task==master) then
      do j = 2,size(mask,2)-1
        if (mod(j,size_Of_Cluster).ne.master) then
            call MPI_RECV(h(:,j), size(mask,1), MPI_DOUBLE_PRECISION, &
             mod(j,size_Of_Cluster), j, MPI_COMM_WORLD, &
             MPI_STATUS_IGNORE, ierror)
        end if
      end do
    else
      do j = 2,size(mask,2)-1
        if (mod(j,size_Of_Cluster).eq.my_task) then
            call MPI_SEND(h(:,j), size(mask,1), MPI_DOUBLE_PRECISION, &
               master, j, MPI_COMM_WORLD, ierror)
        end if
      end do
    end if

      do j = 2,size(mask,2)-1
          CALL MPI_BCAST(h(:,j), size(mask,1), MPI_DOUBLE_PRECISION,&
                 master, MPI_COMM_WORLD, ierror)
      end do

    end do

    deallocate( rx0 )
END SUBROUTINE dc_scheme


SUBROUTINE meo_scheme(mask,h,invA,target_rx0)
    implicit none
    include 'mpif.h'
    logical, dimension(:,:), intent(in) :: mask
    real(kind=8), dimension(:,:), intent(in) ::invA
    real(kind=8), dimension(:,:), intent(inout) :: h
    real(kind=8), intent(in) :: target_rx0

    integer :: i,j
    real(kind=8) :: rf,V
    real(kind=8), dimension(:,:), allocatable :: rx0
    integer :: my_task, size_Of_Cluster, ierror
    integer, parameter :: master = 0

    INTERFACE
      SUBROUTINE calc_rx0(mask,h,rx0,my_task,size_Of_Cluster,ierror,master)
        logical, dimension(:,:) :: mask
        real(kind=8), dimension(:,:) :: h
        real(kind=8), dimension(:,:) :: rx0
        integer :: my_task,size_Of_Cluster,ierror,master
      END SUBROUTINE
    END INTERFACE

    call MPI_COMM_SIZE(MPI_COMM_WORLD,size_Of_Cluster,ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD,my_task,ierror)

    rf = (1+target_rx0)/(1-target_rx0)

    allocate( rx0(size(mask,1),size(mask,2)) )

    do
      CALL calc_rx0(mask,h,rx0,my_task,size_Of_Cluster,ierror,master)

      do j = 1,size(mask,2)
          CALL MPI_BCAST(rx0(:,j), size(mask,1), MPI_DOUBLE_PRECISION,&
                 master, MPI_COMM_WORLD, ierror)
      end do

      if (my_task==master) write(*,*) "Next loop",maxval(rx0)

      if ( maxval(rx0).le.(target_rx0+0.001) ) EXIT

      do j = 2,size(mask,2)-1
      if (mod(j,size_Of_Cluster).eq.my_task) then
      do i = 2,size(mask,1)-1
      if (mask(i,j)) then
        if(mask(i+1,j)) then
          if ( (h(i,j)-h(i+1,j))/(h(i,j)+h(i+1,j)).gt.target_rx0 ) then
            V = (h(i,j)-rf*h(i+1,j))/(rf*(invA(i,j)+invA(i+1,j)))
            h(i,j) = h(i,j)-V/invA(i,j)
            h(i+1,j) = h(i+1,j)+V/invA(i+1,j)
          end if
        end if
        if(mask(i-1,j)) then
          if ( (h(i,j)-h(i-1,j))/(h(i,j)+h(i-1,j)).gt.target_rx0 ) then
            V = (h(i,j)-rf*h(i-1,j))/(rf*(invA(i,j)+invA(i-1,j)))
            h(i,j) = h(i,j)-V/invA(i,j)
            h(i-1,j) = h(i-1,j)+V/invA(i-1,j)
          end if
        end if
      end if
      end do
      end if
      end do

    if (my_task==master) then
      do j = 2,size(mask,2)-1
        if (mod(j,size_Of_Cluster).ne.master) then
            call MPI_RECV(h(:,j), size(mask,1), MPI_DOUBLE_PRECISION, &
             mod(j,size_Of_Cluster), j, MPI_COMM_WORLD, &
             MPI_STATUS_IGNORE, ierror)
        end if
      end do
    else
      do j = 2,size(mask,2)-1
        if (mod(j,size_Of_Cluster).eq.my_task) then
            call MPI_SEND(h(:,j), size(mask,1), MPI_DOUBLE_PRECISION, &
               master, j, MPI_COMM_WORLD, ierror)
        end if
      end do
    end if

      do j = 2,size(mask,2)-1
          CALL MPI_BCAST(h(:,j), size(mask,1), MPI_DOUBLE_PRECISION,&
                 master, MPI_COMM_WORLD, ierror)
      end do

      do i = 2,size(mask,1)-1
      if (mod(i,size_Of_Cluster).eq.my_task) then
      do j = 2,size(mask,2)-1
      if (mask(i,j)) then
        if(mask(i,j+1)) then
          if ( (h(i,j)-h(i,j+1))/(h(i,j)+h(i,j+1)).gt.target_rx0 ) then
            V = (h(i,j)-rf*h(i,j+1))/(rf*(invA(i,j)+invA(i,j+1)))
            h(i,j) = h(i,j)-V/invA(i,j)
            h(i,j+1) = h(i,j+1)+V/invA(i,j+1)
          end if
        end if
        if(mask(i,j-1)) then
          if ( (h(i,j)-h(i,j-1))/(h(i,j)+h(i,j-1)).gt.target_rx0 ) then
            V = (h(i,j)-rf*h(i,j-1))/(rf*(invA(i,j)+invA(i,j-1)))
            h(i,j) = h(i,j)-V/invA(i,j)
            h(i,j-1) = h(i,j-1)+V/invA(i,j-1)
          end if
        end if
      end if
      end do
      end if
      end do

    if (my_task==master) then
      do i = 2,size(mask,1)-1
        if (mod(i,size_Of_Cluster).ne.master) then
            call MPI_RECV(h(i,:), size(mask,2), MPI_DOUBLE_PRECISION, &
             mod(i,size_Of_Cluster), i, MPI_COMM_WORLD, &
             MPI_STATUS_IGNORE, ierror)
        end if
      end do
    else
      do i = 2,size(mask,1)-1
        if (mod(i,size_Of_Cluster).eq.my_task) then
            call MPI_SEND(h(i,:), size(mask,2), MPI_DOUBLE_PRECISION, &
               master, i, MPI_COMM_WORLD, ierror)
        end if
      end do
    end if

      do i = 2,size(mask,1)-1
          CALL MPI_BCAST(h(i,:), size(mask,2), MPI_DOUBLE_PRECISION,&
                 master, MPI_COMM_WORLD, ierror)
      end do

    end do

    deallocate( rx0 )
END SUBROUTINE meo_scheme

SUBROUTINE check(status,label)
    USE NETCDF
    implicit none
    INTEGER, INTENT(in) :: status , label

    IF (status .ne. nf90_noerr) then
        print*, trim(nf90_strerror(status)),' label=',label
        stop "Stopped"
    END IF
END SUBROUTINE check

