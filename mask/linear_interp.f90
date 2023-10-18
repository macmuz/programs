PROGRAM interpolation
    USE NETCDF
    implicit none
    include 'mpif.h'
    integer :: my_task, size_Of_Cluster, ierror
    integer, parameter :: master_task = 0
    character(len=150) :: filename,outfile(2),gridfile
    character(len=15) :: varname(2),dim1(2),dim2(2)
    character(len=15) :: lonname(3),latname(3)
    integer :: ncid,varid,dimid,nlons(2),nlats(2),vid(2)
    integer :: nlonso,nlatso
    integer :: x_dimid,y_dimid,dimids2(2)
    integer :: i,j,k,idx(2),nx
    integer :: istr1,iend1,jstr1,jend1
    integer :: istr2,iend2,jstr2,jend2
    real(kind=8), allocatable :: lon1(:,:),lon2(:,:),lon(:,:)
    real(kind=8), allocatable :: lat1(:,:),lat2(:,:),lat(:,:)
    real(kind=8), allocatable :: bathy1(:,:),bathy2(:,:)
    real(kind=8), allocatable :: bathy(:,:),mask_rho(:,:)
    real(kind=8), allocatable :: dis_array1(:,:),dis_array2(:,:)
!    real(kind=8), allocatable :: y2der1(:,:),yy_tmp1(:),y2_tmp1(:)
!    real(kind=8), allocatable :: y2der2(:,:),yy_tmp2(:),y2_tmp2(:)
!    real(kind=8), allocatable :: y2der3(:,:),yy_tmp3(:),y2_tmp3(:)
    real(kind=8) :: fv,tmp,splint,hmin
    real(kind=8) :: corner_x(2),corner_y(2),f(2,2)
    logical, allocatable :: mymask1(:,:),mymask2(:,:),mymaskout(:,:)

    call MPI_INIT(ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,size_Of_Cluster,ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD,my_task,ierror)

!##################################
    filename = 'myncfile.nc'
    gridfile = 'ROMS_grid_without_bathy_025NM.nc'
    varname(1) = 'data1'
    varname(2) = 'data2'
    lonname(1) = 'tlon1'
    lonname(2) = 'tlon2'
    lonname(3) = 'lon_rho'
    latname(1) = 'tlat1'
    latname(2) = 'tlat2'
    latname(3) = 'lat_rho'
    fv = -999.0
    hmin = 2.0
    outfile(1) = 'out_extrap_025NM.nc'
    outfile(2) = 'ROMS_bathy_linear_2m_025NM.nc'
!##################################

    CALL check(nf90_open(trim(filename), NF90_NOWRITE, ncid),100)

    CALL check(nf90_inq_varid(ncid, trim(varname(1)), varid),101)
    CALL check(nf90_inquire_variable(ncid, varid, dimids=dimids2),102)
    CALL check(nf90_inquire_dimension(ncid, dimids2(1), len=nlons(1), &
                                                     name=dim1(1)),103)
    CALL check(nf90_inquire_dimension(ncid, dimids2(2), len=nlats(1), &
                                                     name=dim2(1)),104)
    allocate(bathy1(nlons(1),nlats(1)),mymask1(nlons(1),nlats(1)))
    allocate(lon1(nlons(1),nlats(1)),lat1(nlons(1),nlats(1)))
    allocate(dis_array1(nlons(1),nlats(1)))
!    allocate(y2der1(nlons(1),nlats(1)),yy_tmp1(nlons(1)),y2_tmp1(nlons(1)))
    CALL check(nf90_get_var(ncid, varid, bathy1),105)

    CALL check(nf90_inq_varid(ncid, trim(lonname(1)), varid),106)
    CALL check(nf90_get_var(ncid, varid, lon1),107)
    CALL check(nf90_inq_varid(ncid, trim(latname(1)), varid),108)
    CALL check(nf90_get_var(ncid, varid, lat1),109)


    CALL check(nf90_inq_varid(ncid, trim(varname(2)), varid),110)
    CALL check(nf90_inquire_variable(ncid, varid, dimids=dimids2),111)
    CALL check(nf90_inquire_dimension(ncid, dimids2(1), len=nlons(2), &
                                                     name=dim1(2)),112)
    CALL check(nf90_inquire_dimension(ncid, dimids2(2), len=nlats(2), &
                                                     name=dim2(2)),113)
    allocate(bathy2(nlons(2),nlats(2)),mymask2(nlons(2),nlats(2)))
    allocate(lon2(nlons(2),nlats(2)),lat2(nlons(2),nlats(2)))
    allocate(dis_array2(nlons(2),nlats(2)))
!    allocate(y2der2(nlons(2),nlats(2)),yy_tmp2(nlons(2)),y2_tmp2(nlons(2)))
    CALL check(nf90_get_var(ncid, varid, bathy2),114)

    CALL check(nf90_inq_varid(ncid, trim(lonname(2)), varid),115)
    CALL check(nf90_get_var(ncid, varid, lon2),116)
    CALL check(nf90_inq_varid(ncid, trim(latname(2)), varid),117)
    CALL check(nf90_get_var(ncid, varid, lat2),118)

    CALL check(nf90_close(ncid),119)

    
    CALL check(nf90_open(trim(gridfile), NF90_NOWRITE, ncid),120)
    CALL check(nf90_inq_varid(ncid, trim(lonname(3)), varid),121)
    CALL check(nf90_inquire_variable(ncid, varid, dimids=dimids2),122)
    CALL check(nf90_inquire_dimension(ncid, dimids2(1), len=nlonso),123)
    CALL check(nf90_inquire_dimension(ncid, dimids2(2), len=nlatso),124)
    allocate(lon(nlonso,nlatso),lat(nlonso,nlatso))
    allocate(bathy(nlonso,nlatso),mymaskout(nlonso,nlatso))
    allocate(mask_rho(nlonso,nlatso))
    CALL check(nf90_get_var(ncid, varid, lon),125)
    CALL check(nf90_inq_varid(ncid, trim(latname(3)), varid),126)
    CALL check(nf90_get_var(ncid, varid, lat),127)
    CALL check(nf90_inq_varid(ncid, 'mask_rho', varid),126)
    CALL check(nf90_get_var(ncid, varid, mask_rho),127)
    CALL check(nf90_close(ncid),128)


    mymask1=.TRUE.
    where(bathy1.NE.bathy1) mymask1=.FALSE.
    mymask2=.TRUE.
    where(bathy2.NE.bathy2) mymask2=.FALSE.

    mymaskout=.TRUE.
    where(mask_rho.lt.0.5) mymaskout=.FALSE.


    do i=3,nlons(1)
    do j=1,nlats(1)
      if (mymask1(i,j)) then
        tmp=lon1(i-2,1)
        goto 100
      end if 
    end do
    end do
!100 write(*,*) tmp
100 continue
    do i=1,nlonso
      if (lon(i,1).gt.tmp) then
        istr1=i
        EXIT
      end if
    end do
!    write(*,*) istr1

    do i=nlons(1)-2,1,-1
    do j=1,nlats(1)
      if (mymask1(i,j)) then
        tmp=lon1(i+2,1)
        goto 101
      end if
    end do
    end do
!101 write(*,*) tmp
101 continue
    do i=nlonso,1,-1
      if (lon(i,1).lt.tmp) then
        iend1=i
        EXIT
      end if
    end do
!    write(*,*) iend1

    do j=nlats(1)-2,1,-1
    do i=1,nlons(1)
      if (mymask1(i,j)) then
        tmp=lat1(1,j+2)
        goto 102
      end if
    end do
    end do
!102 write(*,*) tmp
102 continue
    do j=nlatso,1,-1
      if (lat(1,j).lt.tmp) then
        jend1=j
        EXIT
      end if
    end do
!    write(*,*) jend1

    tmp=lat1(1,1)
!    write(*,*) tmp
    do j=1,nlatso
      if (lat(1,j).gt.tmp) then
        jstr1=j
        EXIT
      end if
    end do
!    write(*,*) jstr1

    do i=nlons(2)-2,1,-1
    do j=1,nlats(2)
      if (mymask2(i,j)) then
        tmp=lon2(i+2,1)
        goto 103
      end if
    end do
    end do
!103 write(*,*) tmp
103 continue
    do i=nlonso,1,-1
      if (lon(i,1).lt.tmp) then
        iend2=i
        EXIT
      end if
    end do
!    write(*,*) iend2

    do j=3,nlats(2)
    do i=1,nlons(2)
      if (mymask2(i,j)) then
        tmp=lat2(1,j-2)
        goto 104
      end if
    end do
    end do
!104 write(*,*) tmp
104 continue
    do j=1,nlatso
      if (lat(1,j).gt.tmp) then
        jstr2=j
        EXIT
      end if
    end do
!    write(*,*) jstr2

    istr2=1
    jend2=jstr1-1
    if (my_task==master_task) then
      write(*,*) istr1,iend1,jstr1,jend1
      write(*,*) istr2,iend2,jstr2,jend2
    end if
    
 
    where(.not.mymask1) bathy1=hmin
    where(bathy1.lt.hmin) bathy1=hmin
    CALL extrap(bathy1,mymask1,nlons(1),nlats(1),1000)
    where(.not.mymask2) bathy2=hmin
    where(bathy2.lt.hmin) bathy2=hmin
    CALL extrap(bathy2,mymask2,nlons(2),nlats(2),1000)

    if (my_task==master_task) then
      CALL check(nf90_create( trim(outfile(1)), NF90_CLOBBER, ncid ), 200)

      CALL check(nf90_def_dim( ncid, trim(dim1(1)), nlons(1), x_dimid ), 201)
      CALL check(nf90_def_dim( ncid, trim(dim2(1)), nlats(1), y_dimid ), 202)
      dimids2 = (/ x_dimid, y_dimid /)
      CALL check(nf90_def_var( ncid, trim(varname(1)), NF90_DOUBLE, &
                            dimids2, vid(1)), 203)
      CALL check( nf90_put_att( ncid, vid(1), '_FillValue', fv), 204)

      CALL check(nf90_def_dim( ncid, trim(dim1(2)), nlons(2), x_dimid ), 201)
      CALL check(nf90_def_dim( ncid, trim(dim2(2)), nlats(2), y_dimid ), 202)
      dimids2 = (/ x_dimid, y_dimid /)
      CALL check(nf90_def_var( ncid, trim(varname(2)), NF90_DOUBLE, &
                            dimids2, vid(2)), 203)
      CALL check( nf90_put_att( ncid, vid(2), '_FillValue', fv), 204)

      CALL check( nf90_enddef( ncid ), 205 )
      CALL check( nf90_put_var( ncid, vid(1), bathy1 ), 206)
      CALL check( nf90_put_var( ncid, vid(2), bathy2 ), 206)
      CALL check(nf90_close( ncid ), 207)
    end if


    bathy = hmin

    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    do j=jstr1,jend1
      if (mod(j,size_Of_Cluster).eq.my_task) then
        do i=istr1,iend1
          if (mymaskout(i,j)) then

            do k = 1, nlons(1)-1
              if ( lon(i,j).ge.lon1(k,1) .and. &
                   lon(i,j).le.lon1(k+1,1) ) then
                corner_x(1) = lon1(k,1)
                corner_x(2) = lon1(k+1,1)
                nx = k
                EXIT
              end if
            end do

            do k = 1, nlats(1)-1
              if ( lat(i,j).ge.lat1(1,k) .and. &
                   lat(i,j).le.lat1(1,k+1) ) then
                corner_y(1) = lat1(1,k)
                corner_y(2) = lat1(1,k+1)
                f(1,1) = bathy1(nx,k) 
                f(1,2) = bathy1(nx,k+1) 
                f(2,1) = bathy1(nx+1,k) 
                f(2,2) = bathy1(nx+1,k+1) 
                EXIT
              end if
            end do

            CALL linear_interp(corner_x,corner_y,f,&
                    (/lon(i,j),lat(i,j)/),bathy(i,j))            

          end if
        end do
      end if
    end do

    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    do j=jstr2,jend2
      if (mod(j,size_Of_Cluster).eq.my_task) then
        do i=istr2,iend2
          if (mymaskout(i,j)) then

            do k = 1, nlons(2)-1
              if ( lon(i,j).ge.lon2(k,1) .and. &
                   lon(i,j).le.lon2(k+1,1) ) then
                corner_x(1) = lon2(k,1)
                corner_x(2) = lon2(k+1,1)
                nx = k
                EXIT
              end if
            end do

            do k = 1, nlats(2)-1
              if ( lat(i,j).ge.lat2(1,k) .and. &
                   lat(i,j).le.lat2(1,k+1) ) then
                corner_y(1) = lat2(1,k)
                corner_y(2) = lat2(1,k+1)
                f(1,1) = bathy2(nx,k)
                f(1,2) = bathy2(nx,k+1)
                f(2,1) = bathy2(nx+1,k)
                f(2,2) = bathy2(nx+1,k+1)
                EXIT
              end if
            end do

            CALL linear_interp(corner_x,corner_y,f,&
                    (/lon(i,j),lat(i,j)/),bathy(i,j))

          end if
        end do
      end if
    end do

!    call splie2(topo_lat,nlons(3),nlats(3),topo,y2der3)
!    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
!    do j=1,nlatso
!    if (mod(j,size_Of_Cluster).eq.my_task) then
!      write(*,*) "Process ", my_task, " is processing j=:",j
!      do i=1,nlons(3)
!        yy_tmp3(i)=splint(topo_lat,topo(i,:),y2der3(i,:),nlats(3),lat(1,j))
!      end do
!      call spline(topo_lon,yy_tmp3,nlons(3),y2_tmp3)
!      do i=1,nlonso
!        bathy(i,j)=splint(topo_lon,yy_tmp3,y2_tmp3,nlons(3),lon(i,j))
!      end do
!    end if
!    end do


!    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
!    nx = nlonso
!    if (my_task==master_task) then

!      do j = 1,nlatso
!        if (mod(j,size_Of_Cluster).ne.master_task) then
!          call MPI_RECV(bathy(:,j), nx, MPI_DOUBLE_PRECISION, &
!             mod(j,size_Of_Cluster), j, MPI_COMM_WORLD, &
!             MPI_STATUS_IGNORE, ierror)
!        end if
!      end do

!    else

!      do j = 1,nlatso
!        if (mod(j,size_Of_Cluster).eq.my_task) then
!          call MPI_SEND(bathy(:,j), nx, MPI_DOUBLE_PRECISION, &
!               master_task, j, MPI_COMM_WORLD, ierror)
!        end if
!      end do

!    end if
!    call MPI_BARRIER(MPI_COMM_WORLD,ierror)



!    if (my_task==master_task) then
!      write(*,*) 'before interp1'
!    end if

!    call splie2(lat1(1,:),nlons(1),nlats(1),bathy1,y2der1)
!    if (my_task==master_task) then
!      write(*,*) 'after 1st splie2'
!    end if

!    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
!    do j=jstr1,jend1
!    if (mod(j,size_Of_Cluster).eq.my_task) then
!      write(*,*) "Process ", my_task, " is processing j=:",j
!      do i=1,nlons(1)
!        yy_tmp1(i)=splint(lat1(1,:),bathy1(i,:),y2der1(i,:),nlats(1),lat(1,j))
!      end do
!      call spline(lon1(:,1),yy_tmp1,nlons(1),y2_tmp1)
!      do i=istr1,iend1
!        if (mymaskout(i,j)) then
!          bathy(i,j)=splint(lon1(:,1),yy_tmp1,y2_tmp1,nlons(1),lon(i,j))
!        end if
!      end do
!    end if
!    end do


!    if (my_task==master_task) then
!      write(*,*) 'before interp2'
!    end if

!    call splie2(lat2(1,:),nlons(2),nlats(2),bathy2,y2der2)
!    if (my_task==master_task) then
!      write(*,*) 'after 2nd splie2'
!    end if

!    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
!    do j=jstr2,jend2
!    if (mod(j,size_Of_Cluster).eq.my_task) then
!      write(*,*) "Process ", my_task, " is processing j=:",j
!      do i=1,nlons(2)
!        yy_tmp2(i)=splint(lat2(1,:),bathy2(i,:),y2der2(i,:),nlats(2),lat(1,j))
!      end do
!      call spline(lon2(:,1),yy_tmp2,nlons(2),y2_tmp2)
!      do i=istr2,iend2
!        if (mymaskout(i,j)) then
!          bathy(i,j)=splint(lon2(:,1),yy_tmp2,y2_tmp2,nlons(2),lon(i,j))
!        end if
!      end do
!    end if
!    end do
!    call MPI_BARRIER(MPI_COMM_WORLD,ierror)

!    do j=jstr1,jend1
!      if (mod(j,size_Of_Cluster).eq.my_task) then
!        do i=istr1,iend1
!          dis_array1(:,:)=sqrt((lon(i,j)-lon1(:,:))**2+(lat(i,j)-lat1(:,:))**2)
!          idx=minloc(dis_array1)
!          if(.not.mymask1(idx(1),idx(2))) bathy(i,j)=fv
!        end do  
!      end if
!    end do

!    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
!    do j=jstr2,jend2
!      if (mod(j,size_Of_Cluster).eq.my_task) then
!        do i=istr2,iend2
!          dis_array2(:,:)=sqrt((lon(i,j)-lon2(:,:))**2+(lat(i,j)-lat2(:,:))**2)
!          idx=minloc(dis_array2)
!          if(.not.mymask2(idx(1),idx(2))) bathy(i,j)=fv
!        end do
!      end if 
!    end do

    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    nx = 1+iend1-istr1
    if (my_task==master_task) then

      do j = jstr1,jend1
        if (mod(j,size_Of_Cluster).ne.master_task) then
          call MPI_RECV(bathy(istr1:iend1,j), nx, MPI_DOUBLE_PRECISION, &
             mod(j,size_Of_Cluster), j, MPI_COMM_WORLD, &
             MPI_STATUS_IGNORE, ierror)
        end if
      end do
    
    else

      do j = jstr1,jend1
        if (mod(j,size_Of_Cluster).eq.my_task) then
          call MPI_SEND(bathy(istr1:iend1,j), nx, MPI_DOUBLE_PRECISION, &
               master_task, j, MPI_COMM_WORLD, ierror)
        end if
      end do  

    end if
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    nx = 1+iend2-istr2
    if (my_task==master_task) then

      do j = jstr2,jend2
        if (mod(j,size_Of_Cluster).ne.master_task) then
          call MPI_RECV(bathy(istr2:iend2,j), nx, MPI_DOUBLE_PRECISION, &
             mod(j,size_Of_Cluster), j, MPI_COMM_WORLD, &
             MPI_STATUS_IGNORE, ierror)
        end if
      end do

    else

      do j = jstr2,jend2
        if (mod(j,size_Of_Cluster).eq.my_task) then
          call MPI_SEND(bathy(istr2:iend2,j), nx, MPI_DOUBLE_PRECISION, &
               master_task, j, MPI_COMM_WORLD, ierror)
        end if
      end do

    end if

    if (my_task==master_task) then
      CALL check(nf90_create( trim(outfile(2)), NF90_CLOBBER, ncid ), 200)
      CALL check(nf90_def_dim( ncid, 'xi_rho', nlonso, x_dimid ), 201)
      CALL check(nf90_def_dim( ncid, 'eta_rho', nlatso, y_dimid ), 202)
      dimids2 = (/ x_dimid, y_dimid /)
      CALL check(nf90_def_var( ncid, 'h', NF90_DOUBLE, &
                            dimids2, vid(1)), 203)
      CALL check( nf90_put_att( ncid, vid(1), 'long_names', &
                'bathymetry at RHO-points'), 204)
      CALL check( nf90_put_att( ncid, vid(1), 'units', 'meter'), 205)
      CALL check( nf90_put_att( ncid, vid(1), '_FillValue', fv), 205)
      CALL check(nf90_def_var( ncid, 'hraw', NF90_DOUBLE, &
                            dimids2, vid(2)), 206)
      CALL check( nf90_put_att( ncid, vid(2), 'long_names', &
                'Working bathymetry at RHO-points'), 207)
      CALL check( nf90_put_att( ncid, vid(2), 'units', 'meter'), 208)
      CALL check( nf90_put_att( ncid, vid(2), '_FillValue', fv), 208)
      CALL check( nf90_enddef( ncid ), 209)
      CALL check( nf90_put_var( ncid, vid(1), bathy ), 210)
      CALL check( nf90_put_var( ncid, vid(2), bathy ), 211)
      CALL check(nf90_close( ncid ), 212)
    end if

    call MPI_FINALIZE(ierror)

    deallocate(bathy1,mymask1)
    deallocate(bathy2,mymask2)
    deallocate(lon,lon1,lon2)
    deallocate(lat,lat1,lat2)
    deallocate(bathy,mymaskout,mask_rho)
    deallocate(dis_array1,dis_array2)
!    deallocate(y2der1,yy_tmp1,y2_tmp1)
!    deallocate(y2der2,yy_tmp2,y2_tmp2)
!    deallocate(y2der3,yy_tmp3,y2_tmp3)
END PROGRAM interpolation

SUBROUTINE check(status,label)
    USE NETCDF
    implicit none
    INTEGER, INTENT(in) :: status , label

    IF (status .ne. nf90_noerr) then
        print*, trim(nf90_strerror(status)),'label=',label
        stop "Stopped"
    END IF
END SUBROUTINE check

SUBROUTINE extrap(a,seamask,lon,lat,maxscn)
    implicit none
    integer, intent(in) :: lon,lat,maxscn
    real(kind=8), intent(inout) :: a(lon,lat)
    logical, intent(in) :: seamask(lon,lat)

    integer :: i,j,n
    real(kind=8) :: relc
    real(kind=8), dimension(lon,lat) :: sor,res

    relc=1.0
    do i=1,lon
        do j=1,lat
            if (.not.seamask(i,j)) then
                sor(i,j)=relc
            else
                sor(i,j)=0.0
            endif
        end do
    end do

    do n=1,maxscn
        do i=2,lon-1
            do j=2,lat-1
                res(i,j)=0.25*(a(i-1,j)+a(i+1,j)+a(i,j-1)+a(i,j+1))-a(i,j)
            end do
        end do
        do i=2,lon-1
            do j=2,lat-1
                res(i,j)=res(i,j)*sor(i,j)
                a(i,j)=a(i,j)+res(i,j)
            end do
        end do
        do j=1,lat
            a(1,j)=a(2,j)
            a(lon,j)=a(lon-1,j)
        end do
        do i=1,lon
            a(i,1)=a(i,2)
            a(i,lat)=a(i,lat-1)
        end do
    end do
END SUBROUTINE extrap

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
END SUBROUTINE

SUBROUTINE linear_interp(corner_x,corner_y,f,point_xy,val)
    implicit none
    real(kind=8), intent(in) :: corner_x(2),corner_y(2)
    real(kind=8), intent(in) :: f(2,2),point_xy(2)
    real(kind=8), intent(out) :: val

    real(kind=8) :: x,x1,x2,y,y1,y2
    real(kind=8) :: fr1,fr2

    x = point_xy(1)
    y = point_xy(2)
    x1 = corner_x(1) 
    x2 = corner_x(2) 
    y1 = corner_y(1) 
    y2 = corner_y(2) 

    fr1 = f(1,1)*(x2-x)/(x2-x1)+f(2,1)*(x-x1)/(x2-x1)
    fr2 = f(1,2)*(x2-x)/(x2-x1)+f(2,2)*(x-x1)/(x2-x1)
    
    val = fr1*(y2-y)/(y2-y1)+fr2*(y-y1)/(y2-y1)
END SUBROUTINE linear_interp
