PROGRAM ROMS_grid
    USE NETCDF
    implicit none
    integer, parameter :: itr=3,wc(2)=(/1400,500/),norm_no=3,diag_no=2
    integer :: nx,ny,ncid,dimid,vid,overall
    integer :: i,j,k,l,cnt,idx(itr*4,2)
    integer, allocatable :: mask_int(:,:)
    real(kind=8), allocatable :: h(:,:),mask(:,:)
    character(len=250) :: infile,outfile
    logical :: nb(4),diag

    diag = .true.
    infile = 'ROMS_grid_025NM_minus3.nc'
    outfile = 'ROMS_grid_025NM_minus3_fixed.nc'

    CALL check(nf90_open(trim(infile),NF90_NOWRITE,ncid),1)

    CALL check(nf90_inq_dimid(ncid, "xi_rho", dimid),2)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nx),2)
    CALL check(nf90_inq_dimid(ncid, "eta_rho", dimid),3)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=ny),4)

    allocate( h(nx,ny), mask(nx,ny), mask_int(nx,ny) )

    CALL check(nf90_inq_varid(ncid,"h",vid),5)
    CALL check(nf90_get_var(ncid,vid,h),6)
    CALL check(nf90_inq_varid(ncid,"mask_rho",vid),7)
    CALL check(nf90_get_var(ncid,vid,mask),8)

    CALL check(nf90_close(ncid),9)


!ISLAND REMOVE   
    mask_int = 0
    where( mask.gt.0.5 ) mask_int = 1
goto 1002

    do i = 1, nx
    do j = 1, ny
!    do i = 246, 246
!    do j = 1058, 1058
      if ( mask_int(i,j).eq.0 ) then

        cnt = 1
        idx = 0
        idx(1,:) = (/i,j/)
        do k = 1, itr*4
!          write(*,*) 'k=',k,'cnt=',cnt,'idx=',idx(k,:)
          if ( cnt.gt.itr .or. idx(k,1).eq.0 ) then
            do l = 1,k-1
!              write(*,*) 'refill: ',idx(l,:)
              mask_int(idx(l,1),idx(l,2)) = 0 
            end do
            EXIT
          end if
          CALL find_nb(nx,ny,mask_int,idx(k,1),idx(k,2),nb)
           
          if ( nb(1) ) then
            cnt = cnt+1
            idx(cnt,:) = (/idx(k,1),idx(k,2)+1/)
!            write(*,*) idx(cnt,:)
          end if
          if ( nb(2) ) then
            cnt = cnt+1
            idx(cnt,:) = (/idx(k,1)+1,idx(k,2)/)
!            write(*,*) idx(cnt,:)
          end if
          if ( nb(3) ) then
            cnt = cnt+1
            idx(cnt,:) = (/idx(k,1),idx(k,2)-1/)
!            write(*,*) idx(cnt,:)
          end if
          if ( nb(4) ) then
            cnt = cnt+1
            idx(cnt,:) = (/idx(k,1)-1,idx(k,2)/)
!            write(*,*) idx(cnt,:)
          end if
        end do

        if ( cnt.le.itr ) then
          do k = 1, cnt
            mask_int(idx(k,1),idx(k,2)) = 1
          end do
        end if

      end if
    end do
    end do

    where( mask_int.eq.2 ) mask_int = 0
    where( mask_int.eq.1 ) mask = 1.0
   
!    goto 1001
!CHANNEL REMOVE

    do
      cnt = 0
      if (diag) then

        do i = 1,nx-1-diag_no
        do j = 1,ny-1-diag_no
        if (mask_int(i,j).eq.0 .and. mask_int(i+1,j+1).eq.1) then
          do k = 1,diag_no
            if ( mask_int(i+1+k,j+1+k).eq.0) then
              cnt = cnt+1
              do l = k,1,-1
                mask_int(i+l,j+l) = 0
              end do
              exit
            end if
          end do
        end if
        end do
        end do

        do i = 2+diag_no,nx
        do j = 1,ny-1-diag_no
        if (mask_int(i,j).eq.0 .and. mask_int(i-1,j+1).eq.1) then
          do k = 1,diag_no
            if ( mask_int(i-1-k,j+1+k).eq.0) then
              cnt = cnt+1
              do l = k,1,-1
                mask_int(i-l,j+l) = 0
              end do
              exit
            end if
          end do
        end if
        end do
        end do
      end if
    
      do i = 1,nx
      do j = 1,ny-1-norm_no
      if (mask_int(i,j).eq.0 .and. mask_int(i,j+1).eq.1) then
        do k = 1,norm_no
          if ( mask_int(i,j+1+k).eq.0) then
            cnt = cnt+1
            do l = k,1,-1
              mask_int(i,j+l) = 0
            end do
            exit
          end if
        end do
      end if
      end do
      end do

      do j = 1,ny
      do i = 1,nx-1-norm_no
      if (mask_int(i,j).eq.0 .and. mask_int(i+1,j).eq.1) then
        do k = 1,norm_no
          if ( mask_int(i+1+k,j).eq.0) then
            cnt = cnt+1
            do l = k,1,-1
              mask_int(i+l,j) = 0
            end do
            exit
          end if
        end do
      end if
      end do
      end do
 
      write(*,*) 'cnt=',cnt
      if (cnt.eq.0) exit
    end do 

    mask = 0.0
    where(mask_int.eq.1) mask = 1.0

1002 continue
!RECONNECTION
    overall = 0
    mask_int(wc(1),wc(2)) = 2
    do
      cnt = 0
      do i = 2,nx-1
      do j = 2,ny-1
        if (mask_int(i,j).eq.2 .and. mask_int(i+1,j).eq.1) then
           mask_int(i+1,j) = 2
           cnt = cnt+1
        end if
        if (mask_int(i,j).eq.2 .and. mask_int(i-1,j).eq.1) then
           mask_int(i-1,j) = 2
           cnt = cnt+1
        end if
        if (mask_int(i,j).eq.2 .and. mask_int(i,j+1).eq.1) then
           mask_int(i,j+1) = 2
           cnt = cnt+1
        end if
        if (mask_int(i,j).eq.2 .and. mask_int(i,j-1).eq.1) then
           mask_int(i,j-1) = 2
           cnt = cnt+1
        end if
      end do
      end do

      overall = overall+cnt
      write(*,*) 'cnt=',cnt
      write(*,*) 'overall=',overall
      if (cnt.eq.0) exit
    end do

    where(mask_int.ne.2) mask_int = 0
    where(mask_int.eq.2) mask_int = 1

!1001 continue
    mask = 0.0
    where(mask_int.eq.1) mask = 1.0
    where( h.lt.1.0 ) h = 1.0
    where(mask_int.eq.0) h = 1.0
 

    CALL check(nf90_open(trim(outfile),NF90_WRITE,ncid),10)
!    CALL check(nf90_inq_varid(ncid,"h",vid),11)
!    CALL check(nf90_put_var(ncid,vid,h),12)
    CALL check(nf90_inq_varid(ncid,"mask_rho",vid),13)
    CALL check(nf90_put_var(ncid,vid,mask),14)
    CALL check(nf90_close(ncid),15)
 
    deallocate( h, mask, mask_int )
END PROGRAM ROMS_grid

SUBROUTINE check(status,label)
    USE NETCDF
    implicit none
    INTEGER, INTENT(in) :: status , label

    IF (status .ne. nf90_noerr) then
        print*, trim(nf90_strerror(status)),' label=',label
        stop "Stopped"
    END IF
END SUBROUTINE check

SUBROUTINE find_nb(nx,ny,mask,pi,pj,nb)
    implicit none
    integer, intent(in) :: nx,ny,pi,pj
    integer, intent(inout) :: mask(nx,ny)
    logical, intent(out) :: nb(4)

    integer :: i,j

    nb = .false.

!    write(*,*) 'pi,pj=',pi,pj
    if ( pj.lt.ny ) then
      if ( mask(pi,pj+1).eq.0 ) then 
        nb(1)=.true.
!        write(*,*) 'land:',pi,pj+1 
      end if
    end if
    if ( pi.lt.nx ) then
      if ( mask(pi+1,pj).eq.0 ) then 
        nb(2)=.true. 
!        write(*,*) 'land:',pi+1,pj 
      end if
    end if
    if ( pj.gt.1 ) then
      if ( mask(pi,pj-1).eq.0 ) then
        nb(3)=.true. 
!        write(*,*) 'land:',pi,pj-1
      end if
    end if
    if ( pi.gt.1 ) then
      if ( mask(pi-1,pj).eq.0 ) then
        nb(4)=.true. 
!        write(*,*) 'land:',pi-1,pj 
      end if
    end if
    
    mask(pi,pj) = 2
END SUBROUTINE find_nb
