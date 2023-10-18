PROGRAM ROMS_ini_ebaltic
    USE NETCDF
    USE omp_lib
    implicit none
    integer, parameter :: xin=600,yin=640,zin=66
    integer :: t_dimid,date(3),julian
    integer :: xu_dimid,yu_dimid,xv_dimid,yv_dimid
    integer :: dimids4(4),varid(12),dimids2(2)
    integer :: dimids3u(3),dimids4u(4),dimids3v(3),dimids4v(4)
    integer :: nx,ny,Nlvl,i,j,k,l,reclen,x,y
    integer :: narg,YEAR,MON,DAY,idx,ncid,vid,dimid
    integer :: x_dimid,y_dimid,z_dimid,dimids3(3)
    integer, allocatable :: idx_rho(:,:,:,:),mask_rho(:,:)
    integer, allocatable :: idx_u(:,:,:,:),idx_v(:,:,:,:)
    integer, allocatable :: masku(:,:),maskv(:,:)
    real(kind=8), allocatable :: W_u(:,:,:),W_v(:,:,:)
    real(kind=8), allocatable :: W_rho(:,:,:),uvelin(:,:,:),vvelin(:,:,:)
    real(kind=8), allocatable :: loninu(:,:),latinu(:,:),zetain(:,:)
    real(kind=8), allocatable :: lonin(:,:),latin(:,:),depth(:)
    real(kind=8), allocatable :: tempin(:,:,:),saltin(:,:,:)
    real(kind=8), allocatable :: tmpzin(:),tmpxyzin(:,:,:),dis_array(:,:)
    real(kind=8) :: hc,point(2),corners(4,2),time(1),line,tmp_rho,fv,fvu,zw(2)
    real(kind=8), allocatable :: Cs_r(:),s_rho(:),h(:,:),rho(:,:,:,:)
    real(kind=8), allocatable :: Cs_w(:),s_w(:),hz(:,:,:)
    real(kind=8), allocatable :: z_u(:,:,:),z_v(:,:,:),hz_u(:,:,:),hz_v(:,:,:)
    real(kind=8), allocatable :: lon(:,:),lat(:,:),z_rho(:,:,:)
    real(kind=8), allocatable :: lonu(:,:),lonv(:,:),latu(:,:),latv(:,:)
    real(kind=8), allocatable :: temp(:,:,:,:),salt(:,:,:,:),tmp(:,:,:)
    real(kind=8), allocatable :: zeta(:,:,:),uvel(:,:,:,:),vvel(:,:,:,:)
    real(kind=8), allocatable :: ubar(:,:,:),vbar(:,:,:)
    real(kind=8) :: tmp_sum,tmp_h
    real, allocatable :: zero(:,:,:,:)
    character(len=250) :: grid,outgrid,infile
    character(len=20) :: buffer,file1(3),file2(3)
    logical :: ex1,ex2,score,fexit,maskin(xin,yin),maskin3d(xin,yin,zin)
    integer :: mask(xin,yin)
    logical, allocatable :: maskout(:,:)

    !SET STACKSIZE EQUAL 1024 MB per thread
    CALL KMP_SET_STACKSIZE_S(3221225472)
    !END SET

    allocate(lonin(xin,yin),latin(xin,yin),depth(zin+1))
    allocate(loninu(xin,yin),latinu(xin,yin),zetain(xin,yin))
    allocate(uvelin(xin,yin,zin),vvelin(xin,yin,zin))
    allocate(tempin(xin,yin,zin),saltin(xin,yin,zin))
    allocate(tmpzin(zin),tmpxyzin(xin,yin,zin),dis_array(xin,yin))

!READ INPUT FILENAME AND DATE
    narg = command_argument_count()
    if (narg.ne.4) then
      write(*,*) "Program must have 4 arguments: gridfile,Y,M,D"
      stop
    end if
    call get_command_argument(1,grid)
    call get_command_argument(2,buffer)
    read(buffer, *) YEAR
    call get_command_argument(3,buffer)
    read(buffer, *) MON
    call get_command_argument(4,buffer)
    read(buffer, *) DAY

    date=(/YEAR,MON,DAY/)

    idx = index(grid,'.nc')
    write(outgrid,'(A,A)') trim(grid(1:idx-1)),'_2003v2_initial_full.nc'
    write(infile,'(A,i4.4,"-",i2.2,"-",i2.2,A)') &
        'PLGNG001.pop.h.',YEAR,MON,DAY,'-03600.nc'
    file1(1) = "idx_rho2.bin"
    file2(1) = "W_rho2.bin"
    file1(2) = "idx_u2.bin"
    file2(2) = "W_u2.bin"
    file1(3) = "idx_v2.bin"
    file2(3) = "W_v2.bin"

!    write(*,*) trim(outgrid)
    write(*,*) trim(infile)
!END READ COMMAN LINE


    !READ grid
    CALL check(nf90_open(trim(grid),NF90_NOWRITE,ncid),10)
    CALL check(nf90_inq_dimid(ncid, "xi_rho", dimid),11)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nx),12)
    CALL check(nf90_inq_dimid(ncid, "eta_rho", dimid),13)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=ny),14)
    CALL check(nf90_inq_dimid(ncid, "s_rho", dimid),15)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=Nlvl),16)
    allocate(Cs_r(Nlvl),s_rho(Nlvl),h(nx,ny))
    allocate(Cs_w(Nlvl+1),s_w(Nlvl+1),hz(nx,ny,Nlvl))
    allocate( z_u(nx-1,ny,Nlvl),z_v(nx,ny-1,Nlvl) )
    allocate( hz_u(nx-1,ny,Nlvl),hz_v(nx,ny-1,Nlvl) )
    allocate(lon(nx,ny),lat(nx,ny),z_rho(nx,ny,Nlvl))
    allocate(lonu(nx-1,ny),latu(nx-1,ny),masku(nx-1,ny))
    allocate(lonv(nx,ny-1),latv(nx,ny-1),maskv(nx,ny-1))
    allocate(idx_rho(nx,ny,4,2), W_rho(nx,ny,4))
    allocate(idx_u(nx-1,ny,4,2), W_u(nx-1,ny,4))
    allocate(idx_v(nx,ny-1,4,2), W_v(nx,ny-1,4))
    allocate(temp(nx,ny,Nlvl,1),salt(nx,ny,Nlvl,1))
    allocate(tmp(nx,ny,zin+1),maskout(nx,ny),mask_rho(nx,ny))
    allocate(zero(nx,ny,Nlvl,1),rho(nx,ny,Nlvl,1))
    allocate(zeta(nx,ny,1),uvel(nx-1,ny,Nlvl,1),vvel(nx,ny-1,Nlvl,1))
    allocate(ubar(nx-1,ny,1),vbar(nx,ny-1,1))
    CALL check(nf90_inq_varid(ncid,"h",vid),17)
    CALL check(nf90_get_var(ncid,vid,h),18)
    CALL check(nf90_inq_varid(ncid,"hc",vid),19)
    CALL check(nf90_get_var(ncid,vid,hc),20)
    CALL check(nf90_inq_varid(ncid,"lon_rho",vid),21)
    CALL check(nf90_get_var(ncid,vid,lon),22)
    CALL check(nf90_inq_varid(ncid,"lat_rho",vid),27)
    CALL check(nf90_get_var(ncid,vid,lat),28)
    CALL check(nf90_inq_varid(ncid,"lon_u",vid),21)
    CALL check(nf90_get_var(ncid,vid,lonu),22)
    CALL check(nf90_inq_varid(ncid,"lat_u",vid),27)
    CALL check(nf90_get_var(ncid,vid,latu),28)
    CALL check(nf90_inq_varid(ncid,"lon_v",vid),21)
    CALL check(nf90_get_var(ncid,vid,lonv),22)
    CALL check(nf90_inq_varid(ncid,"lat_v",vid),27)
    CALL check(nf90_get_var(ncid,vid,latv),28)
    CALL check(nf90_inq_varid(ncid,"Cs_r",vid),31)
    CALL check(nf90_get_var(ncid,vid,Cs_r),32)
    CALL check(nf90_inq_varid(ncid,"s_rho",vid),33)
    CALL check(nf90_get_var(ncid,vid,s_rho),34)
    CALL check(nf90_inq_varid(ncid,"Cs_w",vid),31)
    CALL check(nf90_get_var(ncid,vid,Cs_w),32)
    CALL check(nf90_inq_varid(ncid,"s_w",vid),33)
    CALL check(nf90_get_var(ncid,vid,s_w),34)
    CALL check(nf90_inq_varid(ncid,"mask_rho",vid),33)
    CALL check(nf90_get_var(ncid,vid,mask_rho),34)
    CALL check(nf90_inq_varid(ncid,"mask_u",vid),33)
    CALL check(nf90_get_var(ncid,vid,masku),34)
    CALL check(nf90_inq_varid(ncid,"mask_v",vid),33)
    CALL check(nf90_get_var(ncid,vid,maskv),34)
    CALL check(nf90_close(ncid),35)
    !END READ grid

!CALCULATIONS OF SIGMA COORD VARIABLE
    do i = 1,nx 
      do j = 1,ny
        do k = 1,Nlvl
          z_rho(i,j,k) = h(i,j)*(hc*s_rho(k)+h(i,j)*Cs_r(k))/(hc+h(i,j))
          zw(1) = h(i,j)*(hc*s_w(k+1)+h(i,j)*Cs_w(k+1))/(hc+h(i,j))
          zw(2) = h(i,j)*(hc*s_w(k)+h(i,j)*Cs_w(k))/(hc+h(i,j))
          hz(i,j,k)= zw(1)-zw(2)
        end do
      end do
    end do
    
    do i = 1,nx-1 
      do j = 1,ny
        do k = 1,Nlvl
          z_u(i,j,k) = 0.5*(z_rho(i,j,k)+z_rho(i+1,j,k))
          hz_u(i,j,k) = 0.5*(hz(i,j,k)+hz(i+1,j,k))
        end do
      end do
    end do
    
    do i = 1,nx 
      do j = 1,ny-1
        do k = 1,Nlvl
          z_v(i,j,k) = 0.5*(z_rho(i,j,k)+z_rho(i,j+1,k))
          hz_v(i,j,k) = 0.5*(hz(i,j,k)+hz(i,j+1,k))
        end do
      end do
    end do
!END CALCULATIONS OF SIGMA COORD VARIABLES

    do k = 1,Nlvl
      write(*,*) 'k=',k,'depth=',z_rho(1400,500,k) 
    end do

    write(*,*) 'nx,ny,N',nx,ny,Nlvl

    !READ initial from ebaltic
    CALL check(nf90_open(trim(infile),NF90_NOWRITE,ncid),40)

    CALL check(nf90_inq_varid(ncid,"TLONG",vid),41)
    CALL check(nf90_get_var(ncid,vid,lonin),42)
    CALL check(nf90_inq_varid(ncid,"TLAT",vid),43)
    CALL check(nf90_get_var(ncid,vid,latin),44)
    CALL check(nf90_inq_varid(ncid,"ULONG",vid),41)
    CALL check(nf90_get_var(ncid,vid,loninu),42)
    CALL check(nf90_inq_varid(ncid,"ULAT",vid),43)
    CALL check(nf90_get_var(ncid,vid,latinu),44)
    CALL check(nf90_inq_varid(ncid,"z_t",vid),45)
    CALL check(nf90_get_var(ncid,vid,tmpzin),46)
    do k = 1,zin
      depth(k) = (-1)*tmpzin(1+zin-k)/100.0
    end do
    depth(zin+1) = 0.0
    
    CALL check(nf90_inq_varid(ncid,"SSH",vid),47)
    CALL check(nf90_get_var(ncid,vid,zetain,&
        start=(/1,1,1/),count=(/xin,yin,1/)),46)

    CALL check(nf90_inq_varid(ncid,"TEMP",vid),47)
    CALL check(nf90_get_att(ncid,vid,"_FillValue",fv),471)
    CALL check(nf90_get_var(ncid,vid,tmpxyzin,&
        start=(/1,1,1,1/),count=(/xin,yin,zin,1/)),48)    
    do k = 1,zin
      tempin(:,:,k) = tmpxyzin(:,:,1+zin-k)
    end do
    CALL check(nf90_inq_varid(ncid,"SALT",vid),49)
    CALL check(nf90_get_var(ncid,vid,tmpxyzin,&
        start=(/1,1,1,1/),count=(/xin,yin,zin,1/)),50)    
    do k = 1,zin
      saltin(:,:,k) = tmpxyzin(:,:,1+zin-k)
    end do

    CALL check(nf90_inq_varid(ncid,"UVEL",vid),47)
    CALL check(nf90_get_att(ncid,vid,"_FillValue",fvu),471)
    CALL check(nf90_get_var(ncid,vid,tmpxyzin,&
        start=(/1,1,1,1/),count=(/xin,yin,zin,1/)),48)
    do k = 1,zin
      uvelin(:,:,k) = tmpxyzin(:,:,1+zin-k)
    end do
    CALL check(nf90_inq_varid(ncid,"VVEL",vid),49)
    CALL check(nf90_get_var(ncid,vid,tmpxyzin,&
        start=(/1,1,1,1/),count=(/xin,yin,zin,1/)),50)
    do k = 1,zin
      vvelin(:,:,k) = tmpxyzin(:,:,1+zin-k)
    end do

    CALL check(nf90_close(ncid),51)
    !END READ initial from ebaltic



    do k = 1,zin+1
      write(*,*) 'k=',k,'depth=',depth(k) 
    end do


    !PREPARE bilinear interpolation coefficient
!___RHO___
    INQUIRE(FILE=trim(file1(1)),EXIST=ex1)
    INQUIRE(FILE=trim(file2(1)),EXIST=ex2)

    if (ex1 .and. ex2) then

      inquire(iolength = reclen) idx_rho

      OPEN(10,FILE=trim(file1(1)),ACCESS='DIRECT', recl=reclen, FORM='UNFORMATTED')
      read(10, rec=1) idx_rho
      CLOSE(10)

      inquire(iolength = reclen) W_rho

      OPEN(10,FILE=trim(file2(1)),ACCESS='DIRECT', recl=reclen, FORM='UNFORMATTED')
      read(10, rec=1) W_rho
      CLOSE(10)

    else
      !CALC coeff
      idx_rho = 0
      W_rho = 0
!$OMP PARALLEL DO DEFAULT(PRIVATE),SHARED(lon,lat,lonin,latin),&
!$OMP& SHARED(idx_rho,W_rho), FIRSTPRIVATE(nx,ny) SCHEDULE(DYNAMIC)
      do j = 1, ny
      do i = 1, nx
        dis_array(:,:) = sqrt( (lon(i,j)-lonin(:,:))**2+(lat(i,j)-latin(:,:))**2 )
        idx_rho(i,j,1,:) = minloc( dis_array )
        x = idx_rho(i,j,1,1)
        y = idx_rho(i,j,1,2)
        if (x.gt.1 .and. x.lt.xin .and. y.gt.1 .and. y.lt.yin) then
          point(1) = lon(i,j)
          point(2) = lat(i,j)
          fexit = .false.
          do k = x-1,x
            if (fexit) EXIT
          do l = y-1,y
            if (fexit) EXIT
            corners(1,:) = (/lonin(k,l),latin(k,l)/)
            corners(2,:) = (/lonin(k,l+1),latin(k,l+1)/)
            corners(3,:) = (/lonin(k+1,l+1),latin(k+1,l+1)/)
            corners(4,:) = (/lonin(k+1,l),latin(k+1,l)/)
            CALL in_convex_polygon(corners,point,score)
            if (score) then
              fexit = .true.
              idx_rho(i,j,1,:) = (/k,l/)
              idx_rho(i,j,2,:) = (/k,l+1/)
              idx_rho(i,j,3,:) = (/k+1,l+1/)
              idx_rho(i,j,4,:) = (/k+1,l/)
              call calc_w( lonin(k:k+1,l:l+1), latin(k:k+1,l:l+1), point, W_rho(i,j,:) )
            end if
          end do
          end do
        else
          idx_rho(i,j,1,:) = 0
        end if
      end do
        write(*,*) 'j=',j
      end do
!$OMP END PARALLEL DO
      !END CALC

      inquire(iolength = reclen) idx_rho

      OPEN(10,FILE=trim(file1(1)),ACCESS='DIRECT', recl=reclen, FORM='UNFORMATTED')
      write(10, rec=1) idx_rho
      CLOSE(10)

      inquire(iolength = reclen) W_rho

      OPEN(10,FILE=trim(file2(1)),ACCESS='DIRECT', recl=reclen, FORM='UNFORMATTED')
      write(10, rec=1) W_rho
      CLOSE(10)

    end if !files exist
!___RHO___

!___U___
    INQUIRE(FILE=trim(file1(2)),EXIST=ex1)
    INQUIRE(FILE=trim(file2(2)),EXIST=ex2)

    if (ex1 .and. ex2) then

      inquire(iolength = reclen) idx_u

      OPEN(10,FILE=trim(file1(2)),ACCESS='DIRECT', recl=reclen, FORM='UNFORMATTED')
      read(10, rec=1) idx_u
      CLOSE(10)

      inquire(iolength = reclen) W_u

      OPEN(10,FILE=trim(file2(2)),ACCESS='DIRECT', recl=reclen, FORM='UNFORMATTED')
      read(10, rec=1) W_u
      CLOSE(10)

    else
      !CALC coeff
      idx_u = 0
      W_u = 0
!$OMP PARALLEL DO DEFAULT(PRIVATE),SHARED(lonu,latu,loninu,latinu),&
!$OMP& SHARED(idx_u,W_u), FIRSTPRIVATE(nx,ny) SCHEDULE(DYNAMIC)
      do j = 1, ny
      do i = 1, nx-1
        dis_array(:,:) = sqrt( (lonu(i,j)-loninu(:,:))**2+(latu(i,j)-latinu(:,:))**2 )
        idx_u(i,j,1,:) = minloc( dis_array )
        x = idx_u(i,j,1,1)
        y = idx_u(i,j,1,2)
        if (x.gt.1 .and. x.lt.xin .and. y.gt.1 .and. y.lt.yin) then
          point(1) = lonu(i,j)
          point(2) = latu(i,j)
          fexit = .false.
          do k = x-1,x
            if (fexit) EXIT
          do l = y-1,y
            if (fexit) EXIT
            corners(1,:) = (/loninu(k,l),latinu(k,l)/)
            corners(2,:) = (/loninu(k,l+1),latinu(k,l+1)/)
            corners(3,:) = (/loninu(k+1,l+1),latinu(k+1,l+1)/)
            corners(4,:) = (/loninu(k+1,l),latinu(k+1,l)/)
            CALL in_convex_polygon(corners,point,score)
            if (score) then
              fexit = .true.
              idx_u(i,j,1,:) = (/k,l/)
              idx_u(i,j,2,:) = (/k,l+1/)
              idx_u(i,j,3,:) = (/k+1,l+1/)
              idx_u(i,j,4,:) = (/k+1,l/)
              call calc_w( loninu(k:k+1,l:l+1), latinu(k:k+1,l:l+1), point, W_u(i,j,:) )
            end if
          end do
          end do
        else
          idx_u(i,j,1,:) = 0
        end if
      end do
        write(*,*) 'j=',j
      end do
!$OMP END PARALLEL DO
      !END CALC

      inquire(iolength = reclen) idx_u

      OPEN(10,FILE=trim(file1(2)),ACCESS='DIRECT', recl=reclen, FORM='UNFORMATTED')
      write(10, rec=1) idx_u
      CLOSE(10)

      inquire(iolength = reclen) W_u

      OPEN(10,FILE=trim(file2(2)),ACCESS='DIRECT', recl=reclen, FORM='UNFORMATTED')
      write(10, rec=1) W_u
      CLOSE(10)

    end if !files exist
!___U___

!___V___
    INQUIRE(FILE=trim(file1(3)),EXIST=ex1)
    INQUIRE(FILE=trim(file2(3)),EXIST=ex2)

    if (ex1 .and. ex2) then

      inquire(iolength = reclen) idx_v

      OPEN(10,FILE=trim(file1(3)),ACCESS='DIRECT', recl=reclen, FORM='UNFORMATTED')
      read(10, rec=1) idx_v
      CLOSE(10)

      inquire(iolength = reclen) W_v

      OPEN(10,FILE=trim(file2(3)),ACCESS='DIRECT', recl=reclen, FORM='UNFORMATTED')
      read(10, rec=1) W_v
      CLOSE(10)

    else
      !CALC coeff
      idx_v = 0
      W_v = 0
!$OMP PARALLEL DO DEFAULT(PRIVATE),SHARED(lonv,latv,loninu,latinu),&
!$OMP& SHARED(idx_v,W_v), FIRSTPRIVATE(nx,ny) SCHEDULE(DYNAMIC)
      do j = 1, ny-1
      do i = 1, nx
        dis_array(:,:) = sqrt( (lonv(i,j)-loninu(:,:))**2+(latv(i,j)-latinu(:,:))**2 )
        idx_v(i,j,1,:) = minloc( dis_array )
        x = idx_v(i,j,1,1)
        y = idx_v(i,j,1,2)
        if (x.gt.1 .and. x.lt.xin .and. y.gt.1 .and. y.lt.yin) then
          point(1) = lonv(i,j)
          point(2) = latv(i,j)
          fexit = .false.
          do k = x-1,x
            if (fexit) EXIT
          do l = y-1,y
            if (fexit) EXIT
            corners(1,:) = (/loninu(k,l),latinu(k,l)/)
            corners(2,:) = (/loninu(k,l+1),latinu(k,l+1)/)
            corners(3,:) = (/loninu(k+1,l+1),latinu(k+1,l+1)/)
            corners(4,:) = (/loninu(k+1,l),latinu(k+1,l)/)
            CALL in_convex_polygon(corners,point,score)
            if (score) then
              fexit = .true.
              idx_v(i,j,1,:) = (/k,l/)
              idx_v(i,j,2,:) = (/k,l+1/)
              idx_v(i,j,3,:) = (/k+1,l+1/)
              idx_v(i,j,4,:) = (/k+1,l/)
              call calc_w( loninu(k:k+1,l:l+1), latinu(k:k+1,l:l+1), point, W_v(i,j,:) )
            end if
          end do
          end do
        else
          idx_v(i,j,1,:) = 0
        end if
      end do
        write(*,*) 'j=',j
      end do
!$OMP END PARALLEL DO
      !END CALC

      inquire(iolength = reclen) idx_v

      OPEN(10,FILE=trim(file1(3)),ACCESS='DIRECT', recl=reclen, FORM='UNFORMATTED')
      write(10, rec=1) idx_v
      CLOSE(10)

      inquire(iolength = reclen) W_v

      OPEN(10,FILE=trim(file2(3)),ACCESS='DIRECT', recl=reclen, FORM='UNFORMATTED')
      write(10, rec=1) W_v
      CLOSE(10)

    end if !files exist
!___V___

    !END PREPARE

    write(*,*) minval(idx_rho),maxval(idx_rho)
    write(*,*) minval(W_rho),maxval(W_rho)
    write(*,*) tempin(200,100,1),tempin(200,100,zin)


    CALL check(nf90_create( "interp_coef.nc", NF90_NETCDF4, ncid ), 200)
    CALL check(nf90_def_dim( ncid, 'x', nx, x_dimid ), 201)
    CALL check(nf90_def_dim( ncid, 'y', ny, y_dimid ), 202)
    CALL check(nf90_def_dim( ncid, 'xu', nx-1, xu_dimid ), 201)
    CALL check(nf90_def_dim( ncid, 'yv', ny-1, yv_dimid ), 202)
    CALL check(nf90_def_dim( ncid, 'corner', 4, z_dimid ), 202)
    CALL check(nf90_def_dim( ncid, 'index', 2, t_dimid ), 202)
    dimids3 = (/ x_dimid, y_dimid, z_dimid /)
    dimids4 = (/ x_dimid, y_dimid, z_dimid, t_dimid /)
    dimids3u = (/ xu_dimid, y_dimid, z_dimid /)
    dimids4u = (/ xu_dimid, y_dimid, z_dimid, t_dimid /)
    dimids3v = (/ x_dimid, yv_dimid, z_dimid /)
    dimids4v = (/ x_dimid, yv_dimid, z_dimid, t_dimid /)
    CALL check(nf90_def_var( ncid, 'idx', NF90_INT, dimids4, varid(1)), 205)
    CALL check(nf90_def_var( ncid, 'w', NF90_DOUBLE, dimids3, varid(2)), 205)
    CALL check(nf90_def_var( ncid, 'idxu', NF90_INT, dimids4u, varid(3)), 205)
    CALL check(nf90_def_var( ncid, 'wu', NF90_DOUBLE, dimids3u, varid(4)), 205)
    CALL check(nf90_def_var( ncid, 'idxv', NF90_INT, dimids4v, varid(5)), 205)
    CALL check(nf90_def_var( ncid, 'wv', NF90_DOUBLE, dimids3v, varid(6)), 205)
    CALL check( nf90_enddef( ncid ), 207 )
    CALL check( nf90_put_var( ncid, varid(1), idx_rho ), 209 )
    CALL check( nf90_put_var( ncid, varid(2), W_rho ), 209 )
    CALL check( nf90_put_var( ncid, varid(3), idx_u ), 209 )
    CALL check( nf90_put_var( ncid, varid(4), W_u ), 209 )
    CALL check( nf90_put_var( ncid, varid(5), idx_v ), 209 )
    CALL check( nf90_put_var( ncid, varid(6), W_v ), 209 )
    CALL check(nf90_close( ncid ), 210 )

!    goto 1001

    
    maskin3d = .true.
    where(tempin.eq.fv) maskin3d = .false.
    where(maskin3d) saltin = saltin*1000
    CALL extrap3d(tempin,maskin3d,xin,yin,zin,5000,0)
    CALL extrap3d(saltin,maskin3d,xin,yin,zin,1000,2)

    maskin = maskin3d(:,:,zin)
    CALL extrap(zetain,maskin,xin,yin,1000,2)
    zetain = zetain*0.01
    CALL calcme(xin,yin,nx,ny,zetain,zeta(:,:,1),W_rho,idx_rho)
    maskout = .true.
    where(idx_rho(:,:,1,1).eq.0) maskout = .false.
    CALL extrap(zeta(:,:,1),maskout,nx,ny,500,2)

    maskin3d = .true.
    where(uvelin.eq.fvu) maskin3d = .false.
    where(maskin3d) uvelin = uvelin*0.01
    where(maskin3d) vvelin = vvelin*0.01
    CALL extrap3d(uvelin,maskin3d,xin,yin,zin,5000,0)
    CALL extrap3d(vvelin,maskin3d,xin,yin,zin,5000,0)

    CALL check(nf90_create( "input_extrap.nc", NF90_CLOBBER, ncid ), 200)
    CALL check(nf90_def_dim( ncid, 'x', xin, x_dimid ), 201)
    CALL check(nf90_def_dim( ncid, 'y', yin, y_dimid ), 202)
    CALL check(nf90_def_dim( ncid, 'z', zin, z_dimid ), 202)
    dimids3 = (/ x_dimid, y_dimid, z_dimid /)
    CALL check(nf90_def_var( ncid, 'temp', NF90_DOUBLE, dimids3, varid(1)), 205)
    CALL check(nf90_def_var( ncid, 'salt', NF90_DOUBLE, dimids3, varid(2)), 205)
    CALL check( nf90_enddef( ncid ), 207 )
    CALL check( nf90_put_var( ncid, varid(1), tempin ), 209 )
    CALL check( nf90_put_var( ncid, varid(2), saltin ), 209 )
    CALL check(nf90_close( ncid ), 210 )

   
    CALL check(nf90_create( "output_horizontal.nc", NF90_NETCDF4, ncid ), 200)
    CALL check(nf90_def_dim( ncid, 'x', nx, x_dimid ), 201)
    CALL check(nf90_def_dim( ncid, 'y', ny, y_dimid ), 202)
    CALL check(nf90_def_dim( ncid, 'z', zin+1, z_dimid ), 202)
    dimids3 = (/ x_dimid, y_dimid, z_dimid /)
    CALL check(nf90_def_var( ncid, 'temp', NF90_DOUBLE, dimids3, varid(1)), 205)
    CALL check(nf90_def_var( ncid, 'salt', NF90_DOUBLE, dimids3, varid(2)), 205)
    CALL check( nf90_enddef( ncid ), 207 )

    write(*,*) 'DIAGNOSTIC 1'
!$OMP PARALLEL DO DEFAULT(PRIVATE),SHARED(tempin,tmp,W_rho,idx_rho,maskout),&
!$OMP& FIRSTPRIVATE(nx,ny) SCHEDULE(DYNAMIC)
    do k = 1,zin
      CALL calcme(xin,yin,nx,ny,tempin(:,:,k),tmp(:,:,k),W_rho,idx_rho)
      CALL extrap(tmp(:,:,k),maskout,nx,ny,500,2)
      write(*,*) 'temp k=',k
    end do
!$OMP END PARALLEL DO
    tmp(:,:,zin+1)=tmp(:,:,zin)
    CALL vertical_interp(nx,ny,Nlvl,zin+1,depth,z_rho,tmp,temp(:,:,:,1))

    CALL check( nf90_put_var( ncid, varid(1), tmp ), 209 )

    write(*,*) 'DIAGNOSTIC 2'
!$OMP PARALLEL DO DEFAULT(PRIVATE),SHARED(saltin,tmp,W_rho,idx_rho,maskout),&
!$OMP& FIRSTPRIVATE(nx,ny) SCHEDULE(DYNAMIC)
    do k = 1,zin
      CALL calcme(xin,yin,nx,ny,saltin(:,:,k),tmp(:,:,k),W_rho,idx_rho)
      CALL extrap(tmp(:,:,k),maskout,nx,ny,500,2)
      write(*,*) 'salt k=',k
    end do
!$OMP END PARALLEL DO
    tmp(:,:,zin+1)=tmp(:,:,zin)
    CALL vertical_interp(nx,ny,Nlvl,zin+1,depth,z_rho,tmp,salt(:,:,:,1))
    
 
    CALL check( nf90_put_var( ncid, varid(2), tmp ), 209 )
    CALL check(nf90_close( ncid ), 211 )


!INTERPOLATE VELOCITES
    maskout = .true.
    where(idx_u(:,:,1,1).eq.0) maskout(1:nx-1,:) = .false.
!$OMP PARALLEL DO DEFAULT(PRIVATE),SHARED(uvelin,tmp,W_u,idx_u,maskout),&
!$OMP& FIRSTPRIVATE(nx,ny) SCHEDULE(DYNAMIC)
    do k = 1,zin
      CALL calcme(xin,yin,nx-1,ny,uvelin(:,:,k),tmp(1:nx-1,:,k),W_u,idx_u)
      CALL extrap(tmp(1:nx-1,:,k),maskout(1:nx-1,:),nx-1,ny,500,2)
      write(*,*) 'uvel k=',k
    end do
!$OMP END PARALLEL DO
    tmp(1:nx-1,:,zin+1)=tmp(1:nx-1,:,zin)
    CALL vertical_interp(nx-1,ny,Nlvl,zin+1,depth,z_u,tmp(1:nx-1,:,:),uvel(:,:,:,1))


    maskout = .true.
    where(idx_v(:,:,1,1).eq.0) maskout(:,1:ny-1) = .false.
!$OMP PARALLEL DO DEFAULT(PRIVATE),SHARED(vvelin,tmp,W_v,idx_v,maskout),&
!$OMP& FIRSTPRIVATE(nx,ny) SCHEDULE(DYNAMIC)
    do k = 1,zin
      CALL calcme(xin,yin,nx,ny-1,vvelin(:,:,k),tmp(:,1:ny-1,k),W_v,idx_v)
      CALL extrap(tmp(:,1:ny-1,k),maskout(:,1:ny-1),nx,ny-1,500,2)
      write(*,*) 'vvel k=',k
    end do
!$OMP END PARALLEL DO
    tmp(:,1:ny-1,zin+1)=tmp(:,1:ny-1,zin)
    CALL vertical_interp(nx,ny-1,Nlvl,zin+1,depth,z_v,tmp(:,1:ny-1,:),vvel(:,:,:,1))
!END VELOCITIES

!CREATE COMMON OCEAN_TIME
    CALL tjd(date,julian)
    time = real(julian,8)
!END CREATE
    zero = 0.0


    write(*,*) 'temp: ',minval(temp),maxval(temp)
    write(*,*) 'salt: ',minval(salt),maxval(salt)

    CALL dens(nx,ny,Nlvl,temp(:,:,:,1),salt(:,:,:,1),z_rho,rho(:,:,:,1))

    do j = 1, ny
      do i = 1, nx
        if (mask_rho(i,j).eq.1) then
          tmp_rho = rho(i,j,Nlvl,1)
          do k = Nlvl-1,1,-1
            if (rho(i,j,k,1).lt.tmp_rho) then
              CALL bubble_sort(Nlvl,rho(i,j,:,1),temp(i,j,:,1),salt(i,j,:,1))
              write(*,*) 'rho increases with depth at (i,j)=',i,j
              EXIT
            end if
            tmp_rho = rho(i,j,k,1) 
          end do
        end if
      end do 
    end do


    ubar = 0.0
    vbar = 0.0
    do i = 1,nx-1
      do j = 1,ny
        tmp_sum = 0.0
        tmp_h = 0.0
        do k = 1,Nlvl
          tmp_sum = tmp_sum+hz_u(i,j,k)*uvel(i,j,k,1)
          tmp_h = tmp_h+hz_u(i,j,k)
        end do
        if (tmp_h .GT. 0) ubar(i,j,1) = tmp_sum/tmp_h
      end do
    end do
    do i = 1,nx
      do j = 1,ny-1
        tmp_sum = 0.0
        tmp_h = 0.0
        do k = 1,Nlvl
          tmp_sum = tmp_sum+hz_v(i,j,k)*vvel(i,j,k,1)
          tmp_h = tmp_h+hz_v(i,j,k)
        end do
        if (tmp_h .GT. 0) vbar(i,j,1) = tmp_sum/tmp_h
      end do
    end do

!WRITE INI
    CALL check( nf90_create( trim(outgrid),NF90_NETCDF4,ncid ), 300 )
    CALL check( nf90_def_dim( ncid, 'xi_rho', nx, x_dimid ), 301 )
    CALL check( nf90_def_dim( ncid, 'xi_u', nx-1, xu_dimid ), 301 )
    CALL check( nf90_def_dim( ncid, 'xi_v', nx, xv_dimid ), 301 )
    CALL check( nf90_def_dim( ncid, 'eta_rho', ny, y_dimid ), 302 )
    CALL check( nf90_def_dim( ncid, 'eta_u', ny, yu_dimid ), 302 )
    CALL check( nf90_def_dim( ncid, 'eta_v', ny-1, yv_dimid ), 302 )
    CALL check( nf90_def_dim( ncid, 's_rho', Nlvl, z_dimid ), 302 )
    CALL check( nf90_def_dim( ncid, 'ocean_time', NF90_UNLIMITED, t_dimid ), 303)
    dimids2 = (/ x_dimid, y_dimid /)
    dimids3 = (/ x_dimid, y_dimid, t_dimid /)
    dimids4 = (/ x_dimid, y_dimid, z_dimid, t_dimid /)
    dimids3u = (/ xu_dimid, yu_dimid, t_dimid /)
    dimids4u = (/ xu_dimid, yu_dimid, z_dimid, t_dimid /)
    dimids3v = (/ xv_dimid, yv_dimid, t_dimid /)
    dimids4v = (/ xv_dimid, yv_dimid, z_dimid, t_dimid /)

    CALL check(nf90_def_var( ncid, 'ocean_time', NF90_DOUBLE,&
        (/t_dimid/), varid(1) ), 304)
    CALL check( nf90_put_att( ncid, varid(1), 'long_name',&
        'time since initialization' ), 305 )
    CALL check( nf90_put_att( ncid, varid(1), 'units',&
        'days since 1968-05-23 00:00:00 GMT' ), 306 )
    CALL check( nf90_put_att( ncid, varid(1), 'calendar',&
        'gregorian' ), 307 )

    CALL check(nf90_def_var( ncid, 'zeta', NF90_FLOAT,&
        dimids3, varid(2) ), 308)
    CALL check( nf90_put_att( ncid, varid(2), 'long_name',&
        'free-surface' ), 309 )
    CALL check( nf90_put_att( ncid, varid(2), 'units',&
        'meter' ), 310 )
    CALL check( nf90_put_att( ncid, varid(2), 'field',&
        'free-surface, scalar, series' ), 311 )
    CALL check( nf90_put_att( ncid, varid(2), 'time',&
        'ocean_time' ), 312 )

    CALL check(nf90_def_var( ncid, 'salt', NF90_FLOAT,&
        dimids4, varid(3) ), 313)
    CALL check( nf90_put_att( ncid, varid(3), 'long_name',&
        'salinity' ), 314 )
    CALL check( nf90_put_att( ncid, varid(3), 'units',&
        'PSU' ), 315 )
    CALL check( nf90_put_att( ncid, varid(3), 'field',&
        'salinity, scalar, series' ), 316 )
    CALL check( nf90_put_att( ncid, varid(3), 'time',&
        'ocean_time' ), 317 ) 

    CALL check(nf90_def_var( ncid, 'depth', NF90_DOUBLE,&
        (/ x_dimid, y_dimid, z_dimid /), varid(4) ), 313)
    CALL check( nf90_put_att( ncid, varid(4), 'long_name',&
        'depth' ), 314 )
    CALL check( nf90_put_att( ncid, varid(4), 'units',&
        'meter' ), 315 )

    CALL check(nf90_def_var( ncid, 'temp', NF90_FLOAT,&
        dimids4, varid(5) ), 318)
    CALL check( nf90_put_att( ncid, varid(5), 'long_name',&
        'potential temperature' ), 319 )
    CALL check( nf90_put_att( ncid, varid(5), 'units',&
        'Celsius' ), 320 )
    CALL check( nf90_put_att( ncid, varid(5), 'field',&
        'temperature, scalar, series' ), 321 )
    CALL check( nf90_put_att( ncid, varid(5), 'time',&
        'ocean_time' ), 322 )

    CALL check(nf90_def_var( ncid, 'u', NF90_FLOAT,&
        dimids4u, varid(6) ), 323)
    CALL check( nf90_put_att( ncid, varid(6), 'long_name',&
        'u-momentum component' ), 324 )
    CALL check( nf90_put_att( ncid, varid(6), 'units',&
        'meter second-1' ), 325 )
    CALL check( nf90_put_att( ncid, varid(6), 'field',&
        'u-velocity, scalar, series' ), 326 )
    CALL check( nf90_put_att( ncid, varid(6), 'time',&
        'ocean_time' ), 327 )

    CALL check(nf90_def_var( ncid, 'ubar', NF90_FLOAT,&
        dimids3u, varid(7) ), 328)
    CALL check( nf90_put_att( ncid, varid(7), 'long_name',&
        'vertically integrated u-momentum component' ), 329 )
    CALL check( nf90_put_att( ncid, varid(7), 'units',&
        'meter second-1' ), 330 )
    CALL check( nf90_put_att( ncid, varid(7), 'field',&
        'ubar-velocity, scalar, series' ), 331 )
    CALL check( nf90_put_att( ncid, varid(7), 'time',&
        'ocean_time' ), 332 )

    CALL check(nf90_def_var( ncid, 'v', NF90_FLOAT,&
        dimids4v, varid(8) ), 333)
    CALL check( nf90_put_att( ncid, varid(8), 'long_name',&
        'v-momentum component' ), 334 )
    CALL check( nf90_put_att( ncid, varid(8), 'units',&
        'meter second-1' ), 335 )
    CALL check( nf90_put_att( ncid, varid(8), 'field',&
        'v-velocity, scalar, series' ), 336 )
    CALL check( nf90_put_att( ncid, varid(8), 'time',&
        'ocean_time' ), 337 )

    CALL check(nf90_def_var( ncid, 'vbar', NF90_FLOAT,&
        dimids3v, varid(9) ), 338)
    CALL check( nf90_put_att( ncid, varid(9), 'long_name',&
        'vertically integrated v-momentum component' ), 339 )
    CALL check( nf90_put_att( ncid, varid(9), 'units',&
        'meter second-1' ), 340 )
    CALL check( nf90_put_att( ncid, varid(9), 'field',&
        'vbar-velocity, scalar, series' ), 341 )
    CALL check( nf90_put_att( ncid, varid(9), 'time',&
        'ocean_time' ), 342 )

    CALL check(nf90_def_var( ncid, 'lon_rho', NF90_DOUBLE,&
        dimids2, varid(10) ), 313)
    CALL check( nf90_put_att( ncid, varid(10), 'long_name',&
        'longitude of RHO-points' ), 314 )
    CALL check( nf90_put_att( ncid, varid(10), 'units',&
        'degree_east' ), 315 )

    CALL check(nf90_def_var( ncid, 'lat_rho', NF90_DOUBLE,&
        dimids2, varid(11) ), 313)
    CALL check( nf90_put_att( ncid, varid(11), 'long_name',&
        'latitude of RHO-points' ), 314 )
    CALL check( nf90_put_att( ncid, varid(11), 'units',&
        'degree_north' ), 315 )

    CALL check(nf90_def_var( ncid, 'rho', NF90_FLOAT,&
        dimids4, varid(12) ), 313)
    CALL check( nf90_put_att( ncid, varid(12), 'long_name',&
        'density anomaly' ), 314 )
    CALL check( nf90_put_att( ncid, varid(12), 'units',&
        'kilogram meter-3' ), 315 )
    CALL check( nf90_put_att( ncid, varid(12), 'field',&
        'density, scalar, series' ), 341 )
    CALL check( nf90_put_att( ncid, varid(12), 'time',&
        'ocean_time' ), 337 )

    CALL check( nf90_enddef(ncid), 341 )
    CALL check( nf90_put_var( ncid, varid(1), time ), 342 )
    CALL check( nf90_put_var( ncid, varid(2), real(zeta,4) ), 343 )
    CALL check( nf90_put_var( ncid, varid(3), real(salt,4) ), 344 )
    CALL check( nf90_put_var( ncid, varid(4), z_rho ), 345 )
    CALL check( nf90_put_var( ncid, varid(5), real(temp,4) ), 346 )
    CALL check( nf90_put_var( ncid, varid(6), real(uvel,4) ), 347 )
    CALL check( nf90_put_var( ncid, varid(7), real(ubar,4) ), 348 )
    CALL check( nf90_put_var( ncid, varid(8), real(vvel,4) ), 349 )
    CALL check( nf90_put_var( ncid, varid(9), real(vbar,4) ), 350 )
    CALL check( nf90_put_var( ncid, varid(10), lon ), 351 )
    CALL check( nf90_put_var( ncid, varid(11), lat ), 352 )
    CALL check( nf90_put_var( ncid, varid(12), real(rho,4) ), 352 )
    CALL check( nf90_close(ncid), 353 )
!END WRITE INI

1001 continue
    deallocate(idx_u,idx_v,W_u,W_v)
    deallocate(loninu,latinu,zetain,uvelin,vvelin)
    deallocate(Cs_r,s_rho,h,lon,lat,z_rho)
    deallocate(Cs_w,s_w,hz,masku,maskv)
    deallocate(lonu,latu,lonv,latv)
    deallocate(hz_u,hz_v,z_u,z_v)
    deallocate(idx_rho,W_rho,rho,mask_rho)
    deallocate(lonin,latin,depth,tempin,saltin)
    deallocate(tmpzin,tmpxyzin,dis_array)
    deallocate(temp,salt,tmp,maskout,zero)
    deallocate(zeta,uvel,vvel,ubar,vbar)
END PROGRAM ROMS_ini_ebaltic

SUBROUTINE check(status,label)
    USE NETCDF
    implicit none
    INTEGER, INTENT(in) :: status , label

    IF (status .ne. nf90_noerr) then
        print*, trim(nf90_strerror(status)),' label=',label
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

SUBROUTINE calc_w(lon,lat,point_xy,w)
    implicit none
    real(kind=8), intent(in) :: lon(2,2),lat(2,2)
    real(kind=8), intent(in) :: point_xy(2)
    real(kind=8), intent(out) :: w(4)

    real(kind=8) :: x, y
    real(kind=8) :: A, B, C, s, t

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
END SUBROUTINE calc_w

SUBROUTINE extrap3d(a,mask,lon,lat,n,maxscn,met)
    USE omp_lib
    implicit none
    integer, intent(in) :: lon,lat,n,maxscn,met
    real(kind=8), intent(inout) :: a(lon,lat,n)
    logical, intent(in) :: mask(lon,lat,n)

    integer :: i,j,k,l,cnt,overall
    real(kind=8) :: relc,ave
    real(kind=8), dimension(lon,lat,n) :: sor,res
    logical :: mask_tmp(lon,lat,n),mask_tmp2(lon,lat,n)

    relc = 0.6
    sor = 0.0
    where(.not.mask) sor=relc

    !FILLING LAND WITH LAYER AVERAGE
    select case(met)
    case(0)
      where(.not.mask) a=0.0
    case(1)
      do k = 1,n

      cnt = 0
      ave = 0.0        

      do i=1,lon
      do j=1,lat
        if (mask(i,j,k)) then
            ave=ave+a(i,j,k)
            cnt=cnt+1
        end if
      end do
      end do
      if ( cnt.GT.0 ) ave = ave/real(cnt,8)
      where(.not.mask(:,:,k)) a(:,:,k)=ave

      end do
    case(2)
      mask_tmp2 = mask
      mask_tmp = mask_tmp2
      do

      overall = 0
      do i = 1, lon
      do j = 1, lat
      do k = 1, n
        
        if (.not.mask_tmp(i,j,k)) then

          cnt = 0
          ave = 0.0

          if ( i.gt.1 ) then
            if ( mask_tmp(i-1,j,k) ) then
              ave = ave+a(i-1,j,k)
              cnt = cnt+1 
            end if
          end if

          if ( j.gt.1 ) then
            if ( mask_tmp(i,j-1,k) ) then
              ave = ave+a(i,j-1,k)
              cnt = cnt+1 
            end if
          end if

          if ( k.gt.1 ) then
            if ( mask_tmp(i,j,k-1) ) then
              ave = ave+a(i,j,k-1)
              cnt = cnt+1 
            end if
          end if

          if ( i.lt.lon ) then
            if ( mask_tmp(i+1,j,k) ) then
              ave = ave+a(i+1,j,k)
              cnt = cnt+1 
            end if
          end if

          if ( j.lt.lat ) then
            if ( mask_tmp(i,j+1,k) ) then
              ave = ave+a(i,j+1,k)
              cnt = cnt+1 
            end if
          end if

          if ( k.lt.n ) then
            if ( mask_tmp(i,j,k+1) ) then
              ave = ave+a(i,j,k+1)
              cnt = cnt+1 
            end if
          end if

          if ( cnt.gt.0 ) then
            a(i,j,k) = ave/(real(cnt,8))
            overall = overall+cnt
            mask_tmp2(i,j,k) = .true.
          end if

        end if 

      end do
      end do
      end do
      mask_tmp = mask_tmp2
      write(*,*) 'extrap filling: ',overall
      if ( overall.eq.0 ) EXIT

      end do
    end select
    !END FILLING

    do l = 1,maxscn

!$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k),SCHEDULE(DYNAMIC)
      do i=2,lon-1
        do j=2,lat-1
          res(i,j,n)=0.2*(a(i-1,j,n)+a(i+1,j,n)+a(i,j-1,n)+a(i,j+1,n)+&
            a(i,j,n-1))-a(i,j,n)
          res(i,j,1)=0.2*(a(i-1,j,1)+a(i+1,j,1)+a(i,j-1,1)+a(i,j+1,1)+&
            a(i,j,2))-a(i,j,1)
        do k = 2,n-1
          res(i,j,k)=0.1667*(a(i-1,j,k)+a(i+1,j,k)+a(i,j-1,k)+a(i,j+1,k)+&
            a(i,j,k+1)+a(i,j,k-1))-a(i,j,k)
        end do 
        end do
      end do
!$OMP END PARALLEL DO
 
      do k = 2,n-1
      do i=2,lon-1
        res(i,1,k)=0.2*(a(i-1,1,k)+a(i+1,1,k)+a(i,1,k+1)+a(i,1,k-1)+a(i,2,k))-a(i,1,k)
        res(i,lat,k)=0.2*(a(i-1,lat,k)+a(i+1,lat,k)+a(i,lat,k+1)+a(i,lat,k-1)+a(i,lat-1,k))-a(i,lat,k)
      end do
      do j=2,lat-1
        res(1,j,k)=0.2*(a(1,j-1,k)+a(1,j+1,k)+a(1,j,k+1)+a(1,j,k-1)+a(2,j,k))-a(1,j,k)
        res(lon,j,k)=0.2*(a(lon,j-1,k)+a(lon,j+1,k)+a(lon,j,k+1)+a(lon,j,k-1)+a(lon-1,j,k))-a(lon,j,k)
      end do
      end do

      do i=2,lon-1
        res(i,1,1)=0.25*(a(i-1,1,1)+a(i+1,1,1)+a(i,1,2)+a(i,2,1))-a(i,1,1)
        res(i,lat,1)=0.25*(a(i-1,lat,1)+a(i+1,lat,1)+a(i,lat,2)+a(i,lat-1,1))-a(i,lat,1)
        res(i,1,n)=0.25*(a(i-1,1,n)+a(i+1,1,n)+a(i,1,n-1)+a(i,2,n))-a(i,1,n)
        res(i,lat,n)=0.25*(a(i-1,lat,n)+a(i+1,lat,n)+a(i,lat-1,n)+a(i,lat,n-1))-a(i,lat,n)
      end do

      do j=2,lat-1
        res(1,j,1)=0.25*(a(1,j-1,1)+a(1,j+1,1)+a(1,j,2)+a(2,j,1))-a(1,j,1)
        res(lon,j,1)=0.25*(a(lon,j-1,1)+a(lon,j+1,1)+a(lon,j,2)+a(lon-1,j,1))-a(lon,j,1)
        res(1,j,n)=0.25*(a(1,j-1,n)+a(1,j+1,n)+a(1,j,n-1)+a(2,j,n))-a(1,j,n)
        res(lon,j,n)=0.25*(a(lon,j-1,n)+a(lon,j+1,n)+a(lon,j,n-1)+a(lon-1,j,n))-a(lon,j,n)
      end do

      do k = 2,n-1
       res(1,1,k)=0.25*(a(1,1,k-1)+a(1,1,k+1)+a(2,1,k)+a(1,2,k))-a(1,1,k)
       res(lon,1,k)=0.25*(a(lon,1,k-1)+a(lon,1,k+1)+a(lon-1,1,k)+a(lon,2,k))-a(lon,1,k)
       res(1,lat,k)=0.25*(a(1,lat,k-1)+a(1,lat,k+1)+a(2,lat,k)+a(1,lat-1,k))-a(1,lat,k)
       res(lon,lat,k)=0.25*(a(lon,lat,k-1)+a(lon,lat,k+1)+a(lon-1,lat,k)+a(lon,lat-1,k))-a(lon,lat,k)
      end do

      res(1,1,1)=0.3333*(a(2,1,1)+a(1,2,1)+a(1,1,2))-a(1,1,1)
      res(lon,1,1)=0.3333*(a(lon-1,1,1)+a(lon,2,1)+a(lon,1,2))-a(lon,1,1)
      res(1,lat,1)=0.3333*(a(2,lat,1)+a(1,lat-1,1)+a(1,lat,2))-a(1,lat,1)
      res(lon,lat,1)=0.3333*(a(lon-1,lat,1)+a(lon,lat-1,1)+a(lon,lat,2))-a(lon,lat,1)
      res(1,1,n)=0.3333*(a(2,1,n)+a(1,2,n)+a(1,1,n-1))-a(1,1,n)
      res(lon,1,n)=0.3333*(a(lon-1,1,n)+a(lon,2,n)+a(lon,1,n-1))-a(lon,1,n)
      res(1,lat,n)=0.3333*(a(2,lat,n)+a(1,lat-1,n)+a(1,lat,n-1))-a(1,lat,n)
      res(lon,lat,n)=0.3333*(a(lon-1,lat,n)+a(lon,lat-1,n)+a(lon,lat,n-1))-a(lon,lat,n)
      
       
      res=res*sor
      a=a+res

      if (mod(l,100).eq.0) write(*,*) 'extrap3d, n=',l
    end do    

END SUBROUTINE extrap3d

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

!$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j),SCHEDULE(DYNAMIC)
        do i=2,lon-1
          do j=2,lat-1
            res(i,j)=0.25*(a(i-1,j)+a(i+1,j)+a(i,j-1)+a(i,j+1))-a(i,j)
          end do
        end do
!$OMP END PARALLEL DO

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

        if (mod(n,100).eq.0) write(*,*) 'extrap2d, n=',n
    end do
END SUBROUTINE extrap

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

SUBROUTINE vertical_interp(nx,ny,Nlvl,zin,depth,z_rho,input,output)
    implicit none
    integer, intent(in) :: nx,ny,Nlvl,zin
    real(kind=8), allocatable :: temp(:,:,:),salt(:,:,:),tmp(:,:,:)
    real(kind=8), intent(in) :: depth(zin),z_rho(nx,ny,Nlvl)
    real(kind=8), intent(in) :: input(nx,ny,zin) 
    real(kind=8), intent(out) :: output(nx,ny,Nlvl)

    integer :: i,j,k,m
    real(kind=8) :: y0,y1,x0,x1,x

!$OMP PARALLEL DO DEFAULT(PRIVATE),SHARED(output,z_rho,depth,input),&
!$OMP& FIRSTPRIVATE(nx,ny,Nlvl,zin) SCHEDULE(DYNAMIC)
    do i = 1, nx
    do j = 1, ny
    do k = Nlvl,1,-1
      do m = zin,2,-1
        if (z_rho(i,j,k).le.depth(m) .and. z_rho(i,j,k).gt.depth(m-1)) then
          x = z_rho(i,j,k)
          x0 = depth(m-1)
          x1 = depth(m)
          y0 = input(i,j,m-1) 
          y1 = input(i,j,m) 
          output(i,j,k) = y0+(y1-y0)*(x-x0)/(x1-x0)
          EXIT
        else
!          write(*,*) 'ERROR, depth not found'
           continue
        end if
      end do
    end do
    end do
    write(*,*) 'vert_interp, i=',i
    end do
!$OMP END PARALLEL DO
END SUBROUTINE vertical_interp

SUBROUTINE tjd(date,julian)
    implicit none
    integer, intent(in) :: date(3)
    integer, intent(out) :: julian

    integer :: date_0(3),date_tmp(3)

    date_0 = (/1968,5,23/)
    date_tmp(1) = date_0(1)
    date_tmp(2) = date_0(2)
    date_tmp(3) = date_0(3)

    julian = 0

    do
        if( date_tmp(1).eq.date(1) .and.&
            date_tmp(2).eq.date(2) .and.&
            date_tmp(3).eq.date(3) ) EXIT
        julian = julian+1
        CALL add_day(date_tmp,.true.)
    end do
END SUBROUTINE tjd

SUBROUTINE add_day(date,use_leap)
    implicit none
    integer, intent(inout) :: date(3)
    logical, intent(in) :: use_leap

    logical :: leap
    integer :: days_in_month(12)

    days_in_month = (/31,28,31,30,31,30,31,31,30,31,30,31/)

    if ( date(2).eq.2 .and. date(3).eq.28 .and. use_leap) then
        leap = .false.
        if ( mod(date(1),4).eq.0 ) leap = .true.
        if ( mod(date(1),100).eq.0 ) leap = .false.
        if ( mod(date(1),400).eq.0 ) leap = .true.

        if ( leap ) days_in_month(2) = 29
    end if

    date(3) = date(3)+1
    if ( date(3).gt.days_in_month(date(2)) ) then
        date(3) = 1
        date(2) = date(2)+1
        if ( date(2).gt.12 ) then
            date(2) = 1
            date(1) = date(1)+1
        end if
    end if

END SUBROUTINE add_day

SUBROUTINE gauss(nx,ny,h)
    implicit none
    integer, intent(in) :: nx,ny
    real(kind=8), intent(inout) :: h(nx,ny)

    integer :: i,j,wg(5,5)
    real(kind=8) :: win(5,5),tmp(nx,ny)

    wg(1,:) = (/1,4,7,4,1/)
    wg(2,:) = (/4,16,26,16,4/)
    wg(3,:) = (/7,26,41,26,7/)
    wg(4,:) = (/4,16,26,16,4/)
    wg(5,:) = (/1,4,7,4,1/)


    tmp=h
    do i = 3, nx-2
    do j = 3, ny-2
        win = h(i-2:i+2,j-2:j+2)
        tmp(i,j) = sum(win(:,:)*real(wg(:,:)))/273.0
    end do
    end do
    h=tmp

END SUBROUTINE gauss

FUNCTION line(x1,y1,x2,y2,x) result(y)
    implicit none
    real(kind=8), intent(in) :: x1,y1,x2,y2,x
    real(kind=8) :: y

    y = x*(y1-y2)/(x1-x2)+(y1-x1*(y1-y2)/(x1-x2))
END FUNCTION line

SUBROUTINE dens(nx,ny,nz,temp,salt,depth,rho)
    implicit none
    integer, intent(in) :: nx,ny,nz
    real(kind=8), intent(in) :: temp(nx,ny,nz),salt(nx,ny,nz),depth(nx,ny,nz)
    real(kind=8), intent(out) :: rho(nx,ny,nz)

    integer :: i,j,k
    real(kind=8) :: Tp, Tpr10, Ts, Tt, sqrtTs, cff
    real(kind=8), dimension(0:9) :: C
    real(kind=8) :: den,den1,bulk,bulk0,bulk1,bulk2

    real(kind=8), parameter :: A00 = +1.909256e+04
    real(kind=8), parameter :: A01 = +2.098925e+02
    real(kind=8), parameter :: A02 = -3.041638e+00
    real(kind=8), parameter :: A03 = -1.852732e-03
    real(kind=8), parameter :: A04 = -1.361629e-05
    real(kind=8), parameter :: B00 = +1.044077e+02
    real(kind=8), parameter :: B01 = -6.500517e+00
    real(kind=8), parameter :: B02 = +1.553190e-01
    real(kind=8), parameter :: B03 = +2.326469e-04
    real(kind=8), parameter :: D00 = -5.587545e+00
    real(kind=8), parameter :: D01 = +7.390729e-01
    real(kind=8), parameter :: D02 = -1.909078e-02
    real(kind=8), parameter :: E00 = +4.721788e-01
    real(kind=8), parameter :: E01 = +1.028859e-02
    real(kind=8), parameter :: E02 = -2.512549e-04
    real(kind=8), parameter :: E03 = -5.939910e-07
    real(kind=8), parameter :: F00 = -1.571896e-02
    real(kind=8), parameter :: F01 = -2.598241e-04
    real(kind=8), parameter :: F02 = +7.267926e-06
    real(kind=8), parameter :: G00 = +2.042967e-03
    real(kind=8), parameter :: G01 = +1.045941e-05
    real(kind=8), parameter :: G02 = -5.782165e-10
    real(kind=8), parameter :: G03 = +1.296821e-07
    real(kind=8), parameter :: H00 = -2.595994e-07
    real(kind=8), parameter :: H01 = -1.248266e-09
    real(kind=8), parameter :: H02 = -3.508914e-09
    real(kind=8), parameter :: Q00 = +9.99842594e+02
    real(kind=8), parameter :: Q01 = +6.793952e-02
    real(kind=8), parameter :: Q02 = -9.095290e-03
    real(kind=8), parameter :: Q03 = +1.001685e-04
    real(kind=8), parameter :: Q04 = -1.120083e-06
    real(kind=8), parameter :: Q05 = +6.536332e-09
    real(kind=8), parameter :: U00 = +8.24493e-01
    real(kind=8), parameter :: U01 = -4.08990e-03
    real(kind=8), parameter :: U02 = +7.64380e-05
    real(kind=8), parameter :: U03 = -8.24670e-07
    real(kind=8), parameter :: U04 = +5.38750e-09
    real(kind=8), parameter :: V00 = -5.72466e-03
    real(kind=8), parameter :: V01 = +1.02270e-04
    real(kind=8), parameter :: V02 = -1.65460e-06
    real(kind=8), parameter :: W00 = +4.8314e-04


    DO j=1,ny
      DO k=1,nz
        DO i=1,nx

          Tt=MAX(-2.0,temp(i,j,k))
          Ts=MAX(0.0,salt(i,j,k))
          sqrtTs=SQRT(Ts)
          Tp=depth(i,j,k)
          Tpr10=0.1*Tp


!
!-----------------------------------------------------------------------
!  Compute density (kg/m3) at standard one atmosphere pressure.
!-----------------------------------------------------------------------
!
          C(0)=Q00+Tt*(Q01+Tt*(Q02+Tt*(Q03+Tt*(Q04+Tt*Q05))))
          C(1)=U00+Tt*(U01+Tt*(U02+Tt*(U03+Tt*U04)))
          C(2)=V00+Tt*(V01+Tt*V02)
          den1=C(0)+Ts*(C(1)+sqrtTs*C(2)+Ts*W00)

!
!-----------------------------------------------------------------------
!  Compute secant bulk modulus.
!-----------------------------------------------------------------------
!
          C(3)=A00+Tt*(A01+Tt*(A02+Tt*(A03+Tt*A04)))
          C(4)=B00+Tt*(B01+Tt*(B02+Tt*B03))
          C(5)=D00+Tt*(D01+Tt*D02)
          C(6)=E00+Tt*(E01+Tt*(E02+Tt*E03))
          C(7)=F00+Tt*(F01+Tt*F02)
          C(8)=G01+Tt*(G02+Tt*G03)
          C(9)=H00+Tt*(H01+Tt*H02)
    
          bulk0=C(3)+Ts*(C(4)+sqrtTs*C(5))
          bulk1=C(6)+Ts*(C(7)+sqrtTs*G00)
          bulk2=C(8)+Ts*C(9)
          bulk =bulk0-Tp*(bulk1-Tp*bulk2)
!
!-----------------------------------------------------------------------
!  Compute local "in situ" density anomaly (kg/m3 - 1000).
!-----------------------------------------------------------------------
!   
          cff=1.0/(bulk+Tpr10)
          den=den1*bulk*cff
          den=den-1000.0

          rho(i,j,k)=den
        END DO
      END DO
    END DO

END SUBROUTINE dens

SUBROUTINE bubble_sort(N,rho,temp,salt)
    implicit none
    integer, intent(in) :: N
    real(kind=8), intent(inout) :: rho(N),temp(N),salt(N)

    integer :: i,cnt
    real(kind=8) :: tmp

    do
      cnt = 0
      do i = 1,N-1
        if ( rho(i).lt.rho(i+1) ) then
          cnt = cnt+1
          
          tmp = rho(i)
          rho(i) = rho(i+1)
          rho(i+1) = tmp

          tmp = temp(i)
          temp(i) = temp(i+1)
          temp(i+1) = tmp

          tmp = salt(i)
          salt(i) = salt(i+1)
          salt(i+1) = tmp
        end if
      end do
      if (cnt.eq.0) EXIT
    end do
END SUBROUTINE bubble_sort
