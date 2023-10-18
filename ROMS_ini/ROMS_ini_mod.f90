PROGRAM ROMS_ini
    USE NETCDF
    USE omp_lib
    implicit none
    integer, parameter :: xin=763,yin=774,zin=56
    integer :: t_dimid,date(3),julian
    integer :: xu_dimid,yu_dimid,xv_dimid,yv_dimid
    integer :: dimids4(4),varid(11),dimids2(2)
    integer :: dimids3u(3),dimids4u(4),dimids3v(3),dimids4v(4)
    integer :: nx,ny,Nlvl,i,j,k,l,reclen,x,y
    integer :: narg,YEAR,MON,DAY,idx,ncid,vid,dimid
    integer :: x_dimid,y_dimid,z_dimid,dimids3(3)
    integer, allocatable :: idx_rho(:,:,:,:)
    real(kind=8), allocatable :: W_rho(:,:,:)
    real(kind=8), allocatable :: lonin(:),latin(:),depth(:)
    real(kind=8), allocatable :: tempin(:,:,:),saltin(:,:,:)
    real(kind=8), allocatable :: uoin(:,:,:),voin(:,:,:)
    real(kind=8), allocatable :: tmpzin(:),tmpxyzin(:,:,:),dis_array(:,:)
    real(kind=8) :: hc,point(2),corners(4,2),time(1),sst,sss,line,fv
    real(kind=8), allocatable :: Cs_r(:),s_rho(:),h(:,:)
    real(kind=8), allocatable :: lon(:,:),lat(:,:),z_rho(:,:,:)
    real(kind=8), allocatable :: temp(:,:,:,:),salt(:,:,:,:),tmp(:,:,:)
    real, allocatable :: zero(:,:,:,:)
    character(len=250) :: grid,outgrid,nemofile
    character(len=20) :: buffer,file1,file2
    logical :: ex1,ex2,score,fexit,maskin(xin,yin)
    integer :: mask(xin,yin)
    logical, allocatable :: maskout(:,:)

    !SET STACKSIZE EQUAL 1024 MB per thread
    CALL KMP_SET_STACKSIZE_S(3221225472)
    !END SET

    allocate(lonin(xin),latin(yin),depth(zin+1))
    allocate(tempin(xin,yin,zin),saltin(xin,yin,zin))
    allocate(uoin(xin,yin,zin),voin(xin,yin,zin))
    allocate(tmpzin(zin),tmpxyzin(xin,yin,zin),dis_array(xin,yin))

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
    write(outgrid,'(A,A)') trim(grid(1:idx-1)),'_initial.nc'
    write(nemofile,'(A,i4.4,i2.2,i2.2,A)') &
        'BAL-MYP-NEMO_PHY-DailyMeans-',YEAR,MON,DAY,'.nc'
    file1 = "idx_rho.bin"
    file2 = "W_rho.bin"

!    write(*,*) trim(outgrid)
!    write(*,*) trim(nemofile)

    !READ grid
    CALL check(nf90_open(trim(grid),NF90_NOWRITE,ncid),10)
    CALL check(nf90_inq_dimid(ncid, "xi_rho", dimid),11)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nx),12)
    CALL check(nf90_inq_dimid(ncid, "eta_rho", dimid),13)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=ny),14)
    CALL check(nf90_inq_dimid(ncid, "s_rho", dimid),15)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=Nlvl),16)
    allocate(Cs_r(Nlvl),s_rho(Nlvl),h(nx,ny))
    allocate(lon(nx,ny),lat(nx,ny),z_rho(nx,ny,Nlvl))
    allocate(idx_rho(nx,ny,4,2), W_rho(nx,ny,4))
    allocate(temp(nx,ny,Nlvl,1),salt(nx,ny,Nlvl,1))
    allocate(tmp(nx,ny,zin+1),maskout(nx,ny))
    allocate(zero(nx,ny,Nlvl,1))
    CALL check(nf90_inq_varid(ncid,"h",vid),17)
    CALL check(nf90_get_var(ncid,vid,h),18)
    CALL check(nf90_inq_varid(ncid,"hc",vid),19)
    CALL check(nf90_get_var(ncid,vid,hc),20)
    CALL check(nf90_inq_varid(ncid,"lon_rho",vid),21)
    CALL check(nf90_get_var(ncid,vid,lon),22)
    CALL check(nf90_inq_varid(ncid,"lat_rho",vid),27)
    CALL check(nf90_get_var(ncid,vid,lat),28)
    CALL check(nf90_inq_varid(ncid,"Cs_r",vid),31)
    CALL check(nf90_get_var(ncid,vid,Cs_r),32)
    CALL check(nf90_inq_varid(ncid,"s_rho",vid),33)
    CALL check(nf90_get_var(ncid,vid,s_rho),34)
    CALL check(nf90_close(ncid),35)
    !END READ grid

    do i = 1,nx 
    do j = 1,ny 
    do k = 1,Nlvl
      z_rho(i,j,k) = h(i,j)*(hc*s_rho(k)+h(i,j)*Cs_r(k))/(hc+h(i,j))
    end do
    end do
    end do

    do k = 1,Nlvl
    write(*,*) 'k=',k,'depth=',z_rho(200,200,k) 
    end do

    write(*,*) 'nx,ny,N',nx,ny,Nlvl

    !READ initial from NEMO
    CALL check(nf90_open(trim(nemofile),NF90_NOWRITE,ncid),40)

    CALL check(nf90_inq_varid(ncid,"lon",vid),41)
    CALL check(nf90_get_var(ncid,vid,lonin),42)
    CALL check(nf90_inq_varid(ncid,"lat",vid),43)
    CALL check(nf90_get_var(ncid,vid,latin),44)
    CALL check(nf90_inq_varid(ncid,"depth",vid),45)
    CALL check(nf90_get_var(ncid,vid,tmpzin),46)
    do k = 1,zin
      depth(k) = (-1)*tmpzin(1+zin-k)
    end do
    depth(zin+1) = 0.0
    
    CALL check(nf90_inq_varid(ncid,"thetao",vid),47)
    CALL check(nf90_get_att(ncid,vid,'_FillValue',fv),48)
    CALL check(nf90_get_var(ncid,vid,tmpxyzin),48)    
    do k = 1,zin
      tempin(:,:,k) = tmpxyzin(:,:,1+zin-k)
    end do
    CALL check(nf90_inq_varid(ncid,"so",vid),49)
    CALL check(nf90_get_var(ncid,vid,tmpxyzin),50)    
    do k = 1,zin
      saltin(:,:,k) = tmpxyzin(:,:,1+zin-k)
    end do
    CALL check(nf90_inq_varid(ncid,"uo",vid),49)
    CALL check(nf90_get_var(ncid,vid,tmpxyzin),50)    
    do k = 1,zin
      uoin(:,:,k) = tmpxyzin(:,:,1+zin-k)
    end do
    CALL check(nf90_inq_varid(ncid,"vo",vid),49)
    CALL check(nf90_get_var(ncid,vid,tmpxyzin),50)    
    do k = 1,zin
      voin(:,:,k) = tmpxyzin(:,:,1+zin-k)
    end do

    CALL check(nf90_close(ncid),51)
    !END READ initial from NEMO


    do k = 1,zin+1
    write(*,*) 'k=',k,'depth=',depth(k) 
    end do

    !PREPARE bilinear interpolation coefficient
    INQUIRE(FILE=trim(file1),EXIST=ex1)
    INQUIRE(FILE=trim(file2),EXIST=ex2)

    if (ex1 .and. ex2) then

      inquire(iolength = reclen) idx_rho

      OPEN(10,FILE=trim(file1),ACCESS='DIRECT', recl=reclen, FORM='UNFORMATTED')
      read(10, rec=1) idx_rho
      CLOSE(10)

      inquire(iolength = reclen) W_rho

      OPEN(10,FILE=trim(file2),ACCESS='DIRECT', recl=reclen, FORM='UNFORMATTED')
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
        do k = 1, xin
        do l = 1, yin
          dis_array(k,l) = sqrt((lon(i,j)-lonin(k))**2+(lat(i,j)-latin(l))**2)
        end do
        end do
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
            corners(1,:) = (/lonin(k),latin(l)/)
            corners(2,:) = (/lonin(k),latin(l+1)/)
            corners(3,:) = (/lonin(k+1),latin(l+1)/)
            corners(4,:) = (/lonin(k+1),latin(l)/)
            CALL in_convex_polygon(corners,point,score)
            if (score) then
              fexit = .true.
              idx_rho(i,j,1,:) = (/k,l/)
              idx_rho(i,j,2,:) = (/k,l+1/)
              idx_rho(i,j,3,:) = (/k+1,l+1/)
              idx_rho(i,j,4,:) = (/k+1,l/)
              call calc_w( (/lonin(k),lonin(k+1)/),(/latin(l),latin(l+1)/), point, W_rho(i,j,:) )
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

      OPEN(10,FILE=trim(file1),ACCESS='DIRECT', recl=reclen, FORM='UNFORMATTED')
      write(10, rec=1) idx_rho
      CLOSE(10)

      inquire(iolength = reclen) W_rho

      OPEN(10,FILE=trim(file2),ACCESS='DIRECT', recl=reclen, FORM='UNFORMATTED')
      write(10, rec=1) W_rho
      CLOSE(10)

    end if !files exist
    !END PREPARE

    write(*,*) minval(idx_rho),maxval(idx_rho)
    write(*,*) minval(W_rho),maxval(W_rho)
    write(*,*) tempin(200,100,1),tempin(200,100,zin)


!$OMP PARALLEL DO DEFAULT(PRIVATE),SHARED(tempin,saltin),&
!$OMP& SCHEDULE(DYNAMIC)
    do k= zin,1,-1
      maskin = .false.
      mask = 0.0
      where(tempin(:,:,k).gt.-800.0) maskin = .true.
      where(maskin) mask = 1.0
      write(*,*) k,sum(mask)
      write(*,*) 'k=',k,'sum=',sum(mask),'(1,1)=',tempin(1,1,k)
      CALL extrap(tempin(:,:,k),maskin,xin,yin,100,0)
      CALL extrap(saltin(:,:,k),maskin,xin,yin,100,1)
    end do
!$OMP END PARALLEL DO

!    do k= zin,1,-1
!      CALL gauss(xin,yin,tempin(:,:,k))
!      CALL gauss(xin,yin,saltin(:,:,k))
!    end do
!    tempin(:,:,1) = tempin(:,:,2)
!    saltin(:,:,1) = saltin(:,:,2)

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

!    goto 1001
goto 111
   
    CALL check(nf90_create( "output_horizontal.nc", NF90_NETCDF4, ncid ), 200)
    CALL check(nf90_def_dim( ncid, 'x', nx, x_dimid ), 201)
    CALL check(nf90_def_dim( ncid, 'y', ny, y_dimid ), 202)
    CALL check(nf90_def_dim( ncid, 'z', zin+1, z_dimid ), 202)
    dimids3 = (/ x_dimid, y_dimid, z_dimid /)
    CALL check(nf90_def_var( ncid, 'temp', NF90_DOUBLE, dimids3, varid(1)), 205)
    CALL check(nf90_def_var( ncid, 'salt', NF90_DOUBLE, dimids3, varid(2)), 205)
    CALL check( nf90_enddef( ncid ), 207 )

    maskout = .true.
    where(idx_rho(:,:,1,1).eq.0) maskout = .false.
    write(*,*) 'DIAGNOSTIC 1'
!$OMP PARALLEL DO DEFAULT(PRIVATE),SHARED(tempin,tmp,W_rho,idx_rho,maskout),&
!$OMP& FIRSTPRIVATE(nx,ny) SCHEDULE(DYNAMIC)
    do k = 1,zin
      CALL calcme(xin,yin,nx,ny,tempin(:,:,k),tmp(:,:,k),W_rho,idx_rho)
      CALL extrap(tmp(:,:,k),maskout,nx,ny,100,.true.)
      write(*,*) 'k=',k
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
      CALL extrap(tmp(:,:,k),maskout,nx,ny,100,.true.)
      write(*,*) 'k=',k
    end do
!$OMP END PARALLEL DO
    tmp(:,:,zin+1)=tmp(:,:,zin)
    CALL vertical_interp(nx,ny,Nlvl,zin+1,depth,z_rho,tmp,salt(:,:,:,1))
    
 
    CALL check( nf90_put_var( ncid, varid(2), tmp ), 209 )
    CALL check(nf90_close( ncid ), 211 )

!CREATE COMMON OCEAN_TIME
    CALL tjd(date,julian)
    time = real(julian,8)
!END CREATE
    zero = 0.0


    write(*,*) 'temp: ',minval(temp),maxval(temp)
    write(*,*) 'salt: ',minval(salt),maxval(salt)

!    do j = 1, ny
!    do i = 1, nx
!    sst = temp(i,j,Nlvl,1)
!    sss = salt(i,j,Nlvl,1)

!    if (temp(i,j,Nlvl-1,1).gt.sst) then
!      temp(i,j,Nlvl-1,1) = sst
!    end if
!    sst = temp(i,j,Nlvl-1,1)

!    if (salt(i,j,Nlvl-1,1).gt.sss) then
!      salt(i,j,Nlvl-1,1) = sss
!    end if
!    sss = salt(i,j,Nlvl-1,1)

!    do k = Nlvl-2,1,-1

!      if (temp(i,j,k,1).gt.sst) then
!        temp(i,j,k,1) = line(z_rho(i,j,k+2),temp(i,j,k+2,1),&
!            z_rho(i,j,k+1),temp(i,j,k+1,1),z_rho(i,j,k))
!      end if
!      sst = temp(i,j,k,1)

!      if (salt(i,j,k,1).lt.sss) then
!        salt(i,j,k,1) = line(z_rho(i,j,k+2),salt(i,j,k+2,1),&
!            z_rho(i,j,k+1),salt(i,j,k+1,1),z_rho(i,j,k))
!      end if
!      sss = salt(i,j,k,1)

!    end do
!    end do
!    write(*,*) 'j=',j
!    end do
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
    CALL check( nf90_enddef(ncid), 341 )
    CALL check( nf90_put_var( ncid, varid(1), time ), 342 )
    CALL check( nf90_put_var( ncid, varid(2), zero(:,:,1,1) ), 343 )
    CALL check( nf90_put_var( ncid, varid(3), real(salt,4) ), 344 )
    CALL check( nf90_put_var( ncid, varid(4), z_rho ), 345 )
    CALL check( nf90_put_var( ncid, varid(5), real(temp,4) ), 346 )
    CALL check( nf90_put_var( ncid, varid(6), zero(1:nx-1,:,:,:) ), 347 )
    CALL check( nf90_put_var( ncid, varid(7), zero(1:nx-1,:,1,1) ), 348 )
    CALL check( nf90_put_var( ncid, varid(8), zero(:,1:ny-1,:,:) ), 349 )
    CALL check( nf90_put_var( ncid, varid(9), zero(:,1:ny-1,1,1) ), 350 )
    CALL check( nf90_put_var( ncid, varid(10), lon ), 351 )
    CALL check( nf90_put_var( ncid, varid(11), lat ), 352 )
    CALL check( nf90_close(ncid), 353 )
!END WRITE INI

111 continue
    deallocate(Cs_r,s_rho,h,lon,lat,z_rho)
    deallocate(idx_rho,W_rho)
    deallocate(lonin,latin,depth,tempin,saltin)
    deallocate(tmpzin,tmpxyzin,dis_array)
    deallocate(uoin,voin)
    deallocate(temp,salt,tmp,maskout,zero)
END PROGRAM ROMS_ini

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

SUBROUTINE calc_w(corner_x,corner_y,point_xy,w)
    implicit none
    real(kind=8), intent(in) :: corner_x(2),corner_y(2)
    real(kind=8), intent(in) :: point_xy(2)
    real(kind=8), intent(out) :: w(4)

    real(kind=8) :: x,x1,x2,y,y1,y2

    x = point_xy(1)
    y = point_xy(2)
    x1 = corner_x(1)
    x2 = corner_x(2)
    y1 = corner_y(1)
    y2 = corner_y(2)

    w(1) = (x2-x)/(x2-x1)*(y2-y)/(y2-y1)
    w(2) = (x2-x)/(x2-x1)*(y-y1)/(y2-y1)
    w(3) = (x-x1)/(x2-x1)*(y-y1)/(y2-y1)
    w(4) = (x-x1)/(x2-x1)*(y2-y)/(y2-y1)
END SUBROUTINE calc_w


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

!        if (mod(n,100).eq.0) write(*,*) 'extrap2d, n=',n
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
