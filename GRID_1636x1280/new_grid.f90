PROGRAM new_grid
    USE netcdf
    implicit none
    integer, parameter :: imt=1636,jmt=1280
    integer :: ncid,varid(6),xdim,ydim,dimids2(2)
    integer :: reclen,nx,ny,dimid,vid,i,j,d
    integer :: bath(imt,jmt),mask(imt,jmt)
    integer :: idx(1),jdx(1),tmp,mask2(imt,jmt),mask3(imt,jmt)
    integer, allocatable :: hi_mask(:,:,:)
    real(kind=8) :: lon(imt,jmt),lat(imt,jmt),trlon,rad2deg
    real(kind=8), allocatable :: elev(:,:,:),hi_lon(:,:),hi_lat(:,:)
    character(len=50) :: f_lon,f_lat,f_bath,f_elev(2)

    f_lon = 'grid_ro_rec_02.bin'
    f_lat = 'grid_ro_rec_01.bin'
    f_bath = 'kmt.bs0575v4.ocn.20221106.ieeei4'
    f_elev(1) = 'D5_2020.nc'
    f_elev(2) = 'D6_2020.nc'

    rad2deg = 180/(4.D0*DATAN(1.D0))

    inquire(iolength = reclen) bath
    OPEN(10,FILE=trim(f_bath),ACCESS='DIRECT', recl=reclen, &
        FORM='UNFORMATTED')
    read(10, rec=1) bath
    CLOSE(10)

    inquire(iolength = reclen) lon
    OPEN(10,FILE=trim(f_lon),ACCESS='DIRECT', recl=reclen, &
        FORM='UNFORMATTED')
    read(10, rec=1) lon
    CLOSE(10)
    lon = lon*rad2deg

    inquire(iolength = reclen) lat
    OPEN(10,FILE=trim(f_lat),ACCESS='DIRECT', recl=reclen, &
        FORM='UNFORMATTED')
    read(10, rec=1) lat
    CLOSE(10)
    lat = lat*rad2deg

    CALL check(nf90_open(trim(f_elev(1)),NF90_NOWRITE,ncid),10)
    CALL check(nf90_inq_dimid(ncid, "lon", dimid),11)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nx),12)
    CALL check(nf90_inq_dimid(ncid, "lat", dimid),13)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=ny),14)
    allocate( elev(nx,ny,2), hi_lon(nx,2), hi_lat(ny,2), hi_mask(nx,ny,2) )
    CALL check(nf90_inq_varid(ncid,"elevation",vid),17)
    CALL check(nf90_get_var(ncid,vid,elev(:,:,1)),18)
    CALL check(nf90_inq_varid(ncid,"lon",vid),19)
    CALL check(nf90_get_var(ncid,vid,hi_lon(:,1)),20)
    CALL check(nf90_inq_varid(ncid,"lat",vid),21)
    CALL check(nf90_get_var(ncid,vid,hi_lat(:,1)),22)
    CALL check(nf90_close(ncid),35)

!    write(*,*) minval(hi_lon(:,1)),minval(hi_lat(:,1))
!    write(*,*) maxval(hi_lon(:,1)),maxval(hi_lat(:,1))
!    write(*,*) minval(elev(:,:,1)),maxval(elev(:,:,1))

    CALL check(nf90_open(trim(f_elev(2)),NF90_NOWRITE,ncid),10)
    CALL check(nf90_inq_varid(ncid,"elevation",vid),17)
    CALL check(nf90_get_var(ncid,vid,elev(:,:,2)),18)
    CALL check(nf90_inq_varid(ncid,"lon",vid),19)
    CALL check(nf90_get_var(ncid,vid,hi_lon(:,2)),20)
    CALL check(nf90_inq_varid(ncid,"lat",vid),21)
    CALL check(nf90_get_var(ncid,vid,hi_lat(:,2)),22)
    CALL check(nf90_close(ncid),35)

    hi_mask = 1
    where(elev.ne.elev) hi_mask = 0

    trlon = maxval(hi_lon(:,1))
    mask = 1
 
    do i = 1,imt
    do j = 1,jmt
      d = 2
      if (lon(i,j).lt.trlon) d=1
      idx = minloc(abs(hi_lon(:,d)-lon(i,j))) 
      jdx = minloc(abs(hi_lat(:,d)-lat(i,j)))
!      write(*,*) idx,jdx,hi_lon(idx,d),hi_lat(jdx,d),lon(i,j),lat(i,j)
!      write(*,*) hi_mask(idx,jdx,d)
      if (hi_mask(idx(1),jdx(1),d).eq.0) mask(i,j) = 0
    end do
    if (mod(i,100).eq.0) write(*,*) i
    end do
!    write(*,*) minval(hi_lon(:,2)),minval(hi_lat(:,2))
!    write(*,*) maxval(hi_lon(:,2)),maxval(hi_lat(:,2))
!    write(*,*) minval(elev(:,:,2)),maxval(elev(:,:,2))

    CALL wc(imt,jmt,mask,1200,600)
    mask2 = mask
    mask3 = mask
    CALL connect(imt,jmt,mask2,3)
    CALL wc(imt,jmt,mask2,1200,600)
    CALL connect(imt,jmt,mask3,5)
    CALL wc(imt,jmt,mask3,1200,600)
    

    CALL check( nf90_create( 'bathy.nc',NF90_NETCDF4,ncid ), 300 )
    CALL check( nf90_def_dim( ncid, 'imt', imt, xdim ), 301 )
    CALL check( nf90_def_dim( ncid, 'jmt', jmt, ydim ), 301 )
    dimids2 = (/ xdim, ydim /) 
    CALL check(nf90_def_var( ncid, 'h', NF90_INT,&
        dimids2, varid(1) ), 304)
    CALL check( nf90_put_att( ncid, varid(1), 'units',&
        'meter' ), 314 )
    CALL check(nf90_def_var( ncid, 'lon', NF90_DOUBLE,&
        dimids2, varid(2) ), 304)
    CALL check( nf90_put_att( ncid, varid(2), 'long_name',&
        'longitude' ), 314 )
    CALL check( nf90_put_att( ncid, varid(2), 'units',&
        'degree_east' ), 314 )
    CALL check(nf90_def_var( ncid, 'lat', NF90_DOUBLE,&
        dimids2, varid(3) ), 304)
    CALL check( nf90_put_att( ncid, varid(3), 'long_name',&
        'latitude' ), 314 )
    CALL check( nf90_put_att( ncid, varid(3), 'units',&
        'degree_north' ), 314 )
    CALL check(nf90_def_var( ncid, 'mask', NF90_INT,&
        dimids2, varid(4) ), 304)
    CALL check(nf90_def_var( ncid, 'mask2', NF90_INT,&
        dimids2, varid(5) ), 304)
    CALL check(nf90_def_var( ncid, 'mask3', NF90_INT,&
        dimids2, varid(6) ), 304)
    CALL check( nf90_enddef(ncid), 341 )
    CALL check( nf90_put_var( ncid, varid(1), bath ), 352 )
    CALL check( nf90_put_var( ncid, varid(2), lon ), 352 )
    CALL check( nf90_put_var( ncid, varid(3), lat ), 352 )
    CALL check( nf90_put_var( ncid, varid(4), mask ), 352 )
    CALL check( nf90_put_var( ncid, varid(5), mask2 ), 352 )
    CALL check( nf90_put_var( ncid, varid(6), mask3 ), 352 )
    CALL check( nf90_close(ncid), 353 )
    
    deallocate( elev, hi_lon, hi_lat, hi_mask )

END PROGRAM new_grid

SUBROUTINE check(status,label)
    USE NETCDF
    implicit none
    INTEGER, INTENT(in) :: status , label

    IF (status .ne. nf90_noerr) then
        print*, trim(nf90_strerror(status)),' label=',label
        stop "Stopped"
    END IF
END SUBROUTINE check

SUBROUTINE wc(nx,ny,mask,wi,wj)
    integer, intent(in) :: nx,ny,wi,wj
    integer, intent(inout) :: mask(nx,ny)
   
    integer :: i,j,cnt,tmp_mask(nx,ny)
    tmp_mask = 0
    tmp_mask(wi,wj) = 1 
 
    do
      cnt = 0
      do i = 2,nx-1
      do j = 2,ny-1
        if (tmp_mask(i,j).eq.0.and.mask(i,j).eq.1) then
            if (tmp_mask(i,j-1).eq.1) then
              tmp_mask(i,j) = 1
              cnt = cnt+1
            elseif (tmp_mask(i-1,j).eq.1) then
              tmp_mask(i,j) = 1
              cnt = cnt+1
            elseif (tmp_mask(i,j+1).eq.1) then
              tmp_mask(i,j) = 1
              cnt = cnt+1
            elseif (tmp_mask(i+1,j).eq.1) then
              tmp_mask(i,j) = 1
              cnt = cnt+1
            end if
        end if  
      end do
      end do
      write(*,*) cnt
      if (cnt.eq.0) EXIT
    end do 

    do
      cnt = 0 
      do i = 2,nx-1
        if (tmp_mask(i,1).eq.0.and.mask(i,1).eq.1) then
            if (tmp_mask(i-1,1).eq.1) then
              tmp_mask(i,1) = 1
              cnt = cnt+1
            elseif (tmp_mask(i,2).eq.1) then
              tmp_mask(i,1) = 1
              cnt = cnt+1
            elseif (tmp_mask(i+1,1).eq.1) then
              tmp_mask(i,1) = 1
              cnt = cnt+1
            end if 
        end if
        if (tmp_mask(i,ny).eq.0.and.mask(i,ny).eq.1) then
            if (tmp_mask(i-1,ny).eq.1) then
              tmp_mask(i,ny) = 1
              cnt = cnt+1
            elseif (tmp_mask(i,ny-1).eq.1) then
              tmp_mask(i,ny) = 1
              cnt = cnt+1
            elseif (tmp_mask(i+1,ny).eq.1) then
              tmp_mask(i,ny) = 1
              cnt = cnt+1
            end if
        end if
      end do
      if (cnt.eq.0) EXIT
    end do

    do
      cnt = 0
      do j = 2,ny-1
        if (tmp_mask(1,j).eq.0.and.mask(1,j).eq.1) then
            if (tmp_mask(1,j-1).eq.1) then
              tmp_mask(1,j) = 1
              cnt = cnt+1
            elseif (tmp_mask(2,j).eq.1) then
              tmp_mask(1,j) = 1
              cnt = cnt+1
            elseif (tmp_mask(1,j+1).eq.1) then
              tmp_mask(1,j) = 1
              cnt = cnt+1
            end if
        end if
        if (tmp_mask(nx,j).eq.0.and.mask(nx,j).eq.1) then
            if (tmp_mask(nx,j-1).eq.1) then
              tmp_mask(nx,j) = 1
              cnt = cnt+1
            elseif (tmp_mask(nx-1,j).eq.1) then
              tmp_mask(nx,j) = 1
              cnt = cnt+1
            elseif (tmp_mask(nx,j+1).eq.1) then
              tmp_mask(nx,j) = 1
              cnt = cnt+1
            end if
        end if
      end do
      if (cnt.eq.0) EXIT
    end do

    if (tmp_mask(1,1).eq.0.and.mask(1,1).eq.1) then
      if (tmp_mask(1,2).eq.1) then
        tmp_mask(1,1) = 1
      elseif (tmp_mask(2,1).eq.1) then
        tmp_mask(1,1) = 1
      end if
    end if
    if (tmp_mask(nx,1).eq.0.and.mask(nx,1).eq.1) then
      if (tmp_mask(nx-1,1).eq.1) then
        tmp_mask(nx,1) = 1
      elseif (tmp_mask(nx,2).eq.1) then
        tmp_mask(nx,1) = 1
      end if
    end if
    if (tmp_mask(1,ny).eq.0.and.mask(1,ny).eq.1) then
      if (tmp_mask(1,ny-1).eq.1) then
        tmp_mask(1,ny) = 1
      elseif (tmp_mask(2,ny).eq.1) then
        tmp_mask(1,ny) = 1
      end if
    end if
    if (tmp_mask(nx,ny).eq.0.and.mask(nx,ny).eq.1) then
      if (tmp_mask(nx-1,ny).eq.1) then
        tmp_mask(nx,ny) = 1
      elseif (tmp_mask(nx,ny-1).eq.1) then
        tmp_mask(nx,ny) = 1
      end if
    end if

    mask = tmp_mask
END SUBROUTINE wc

SUBROUTINE connect(nx,ny,mask,dx)
    integer, intent(in) :: nx,ny,dx
    integer, intent(inout) :: mask(nx,ny)

    integer :: i,j,k,cnt,tmp_mask(nx,ny)
    logical :: fill
    tmp_mask = 1

    do i = 1,nx
    do j = 1,ny
      if (mask(i,j).eq.0) then
        if (i.eq.1 .or. i.eq.nx .or. j.eq.1 .or. j.eq.ny) then
          tmp_mask(i,j) = 0
        elseif (mask(i+1,j).eq.0 .and. mask(i-1,j).eq.0 .and. &
                mask(i,j+1).eq.0 .and. mask(i,j-1).eq.0) then
          tmp_mask(i,j) = 0
        end if
      end if
    end do
    end do

    do
      cnt = 0
      do i = 2+dx,nx-1-dx
      do j = 2+dx,ny-1-dx
        if (tmp_mask(i,j).eq.1 .and. mask(i,j).eq.0) then

          if (mask(i+1,j).eq.1) then
            fill = .false.
            do k = dx+1,2,-1
              if ( mask(i+k,j).eq.0 ) fill = .true.
              if ( mask(i+k-1,j).eq.1 .and. fill ) then
                mask(i+k-1,j) = 0
                cnt = cnt+1
              end if
            end do
          end if

          if (mask(i-1,j).eq.1) then
            fill = .false.
            do k = dx+1,2,-1
              if ( mask(i-k,j).eq.0 ) fill = .true.
              if ( mask(i-k+1,j).eq.1 .and. fill ) then
                mask(i-k+1,j) = 0
                cnt = cnt+1
              end if
            end do
          end if

          if (mask(i,j+1).eq.1) then
            fill = .false.
            do k = dx+1,2,-1
              if ( mask(i,j+k).eq.0 ) fill = .true.
              if ( mask(i,j+k-1).eq.1 .and. fill ) then
                mask(i,j+k-1) = 0
                cnt = cnt+1
              end if
            end do
          end if

          if (mask(i,j-1).eq.1) then
            fill = .false.
            do k = dx+1,2,-1
              if ( mask(i,j-k).eq.0 ) fill = .true.
              if ( mask(i,j-k+1).eq.1 .and. fill ) then
                mask(i,j-k+1) = 0
                cnt = cnt+1
              end if
            end do
          end if
    
        end if
      end do
      end do
      write(*,*) cnt
      if (cnt.eq.0) EXIT
    end do

END SUBROUTINE
