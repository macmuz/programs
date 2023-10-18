PROGRAM forcing
    USE netcdf
    use omp_lib
    implicit none
    integer :: nx,ny,cnt
    integer :: onx,ony,narg,year,nt
    integer :: ncid,varid,dimid,i,j,m,n,l
    integer :: ncids(8),varids(9),nrec,ntarray(12)
    integer :: xdimid,ydimid,dimids2(2)
    integer, allocatable :: near(:,:,:),idx(:,:)
    real(kind=8), allocatable :: olon(:,:),olat(:,:),output(:,:,:)
    real(kind=8), allocatable :: w(:,:,:)
    real(kind=8), allocatable :: lon(:,:),lat(:,:),input(:,:,:)
    character(len=100) :: grid,igrid
    character(len=10) :: buffer
    character(len=20) :: files(2)
    logical :: ex,leap

    !CALL KMP_SET_STACKSIZE_S(3221225472)
    CALL KMP_SET_STACKSIZE_S(5368709120)
    !END SET

    files(1) = "idx_rho.bin"
    files(2) = "W_rho.bin"

!READ command line arg: YEAR
    narg = command_argument_count()
    if (narg.ne.1) then
      write(*,*) "Program must have YEAR as an argument"
      stop
    end if
    call get_command_argument(1,buffer)
    read(buffer, *) year
!END READ

    grid = "ROMS_grid_125NM_may.nc"

    leap = .false.
    if ( mod(year,4).eq.0 ) leap = .true.
    if ( mod(year,100).eq.0 ) leap = .false.
    if ( mod(year,400).eq.0 ) leap = .true.

    if (leap) then
      nrec = 366*24
    else
      nrec = 365*24
    end if

    allocate( idx(nrec,3) ) 
 
!READ GRID
    CALL check(nf90_open(trim(grid),NF90_NOWRITE,ncid),310)

    CALL check(nf90_inq_dimid(ncid, "xi_rho", dimid),311)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=onx),312)

    CALL check(nf90_inq_dimid(ncid, "eta_rho", dimid),313)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=ony),314)

    allocate( olon(onx,ony), olat(onx,ony), output(onx,ony,9) )
    allocate( near(onx,ony,2), w(onx,ony,4) )

    CALL check(nf90_inq_varid(ncid,"lon_rho",varid),317)
    CALL check(nf90_get_var(ncid,varid,olon),318)

    CALL check(nf90_inq_varid(ncid,"lat_rho",varid),319)
    CALL check(nf90_get_var(ncid,varid,olat),320)

    CALL check(nf90_close(ncid),360)
!END READ GRID

!READ lon/lat
    write(igrid,"('Tair_',i4.4,i2.2,'.nc')") year,1

    CALL check(nf90_open(trim(igrid),NF90_NOWRITE,ncid),410)
    CALL check(nf90_inq_dimid(ncid, "x", dimid),311)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nx),312)

    CALL check(nf90_inq_dimid(ncid, "y", dimid),313)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=ny),314)

    allocate( lon(nx,ny), lat(nx,ny), input(nx,ny,9) )

    CALL check(nf90_inq_varid(ncid,"longitude",varid),417)
    CALL check(nf90_get_var(ncid,varid,lon),418)
    CALL check(nf90_inq_varid(ncid,"latitude",varid),419)
    CALL check(nf90_get_var(ncid,varid,lat),420)
    CALL check(nf90_close(ncid),460)

    where(lon.gt.100.0) lon = lon-360.0
!END READ lon/lat

    CALL prepare(nx,ny,onx,ony,lon,lat,olon,olat,near,w,files)
    CALL create_output(year,onx,ony,ncids,varids)

    cnt = 0
    do m = 1,12
      write(igrid,"('Tair_',i4.4,i2.2,'.nc')") year,m
      CALL check(nf90_open(trim(igrid),NF90_NOWRITE,ncid),410)
      CALL check(nf90_inq_dimid(ncid, "time", dimid),311)
      CALL check(nf90_inquire_dimension(ncid, dimid, len=nt),312)
      CALL check(nf90_close(ncid),460)
      ntarray(m) = nt
      do n = 1,nt
        do l = 1,6
          cnt = cnt+1
          idx(cnt,:) = (/m,n,l/) 
        end do
      end do
    end do 

!$OMP PARALLEL DO DEFAULT(PRIVATE),FIRSTPRIVATE(nrec,year,near,w),&
!$OMP& FIRSTPRIVATE(ntarray,idx,nx,ny,onx,ony,ncids,varids), SCHEDULE(DYNAMIC)
    do i = 1,nrec
      nt = ntarray(idx(i,1))

      !$OMP CRITICAL
      CALL read_data(nx,ny,year,idx(i,1),idx(i,2),idx(i,3),nt,input)
      write(*,*) i
      !$OMP END CRITICAL
      
      where(input(:,:,1:4).lt.0.0) input(:,:,1:4) = 0.0
      input(:,:,5) = input(:,:,5)-273.15
      input(:,:,6) = input(:,:,6)/100.0
      where(input(:,:,7).eq.0.0) input(:,:,7) = 100.0
      where(input(:,:,7).lt.0.0) input(:,:,7) = 0.0
      where(input(:,:,7).gt.100.0) input(:,:,7) = 100.0

      CALL humid(nx,ny,input(:,:,5:7))
      CALL wind(nx,ny,input(:,:,8:9))

      CALL interp2(nx,ny,onx,ony,input,near,w,output)

      !$OMP CRITICAL
      CALL write_output(onx,ony,i,output,ncids,varids)
      !$OMP END CRITICAL

    end do
!$OMP END PARALLEL DO
 
    CALL close_output(ncids) 

    deallocate(idx,olon,olat,lon,lat,near,w,input,output)
END PROGRAM forcing

SUBROUTINE check(status,label)
    USE NETCDF
    implicit none
    INTEGER, INTENT(in) :: status , label

    IF (status .ne. nf90_noerr) then
        print*, trim(nf90_strerror(status)),' label= ',label
        stop "Stopped"
    END IF
END SUBROUTINE check

SUBROUTINE wind(nx,ny,mydata)
    implicit none
    integer, intent(in) :: nx,ny
    real(kind=8), intent(inout) :: mydata(nx,ny,2) 

    real(kind=8) :: wdir(nx,ny),wspeed(nx,ny),pi

    pi = 4.0*ATAN(1.0)

    wdir = (270.0-mydata(:,:,1))*pi/180.0
    wspeed = mydata(:,:,2)

    mydata(:,:,1) = wspeed(:,:)*cos(wdir(:,:))
    mydata(:,:,2) = wspeed(:,:)*sin(wdir(:,:))
END SUBROUTINE wind

SUBROUTINE humid(nx,ny,mydata)
    implicit none
    integer, intent(in) :: nx,ny
    real(kind=8), intent(inout) :: mydata(nx,ny,3)

    real(kind=8) :: cff(nx,ny)

    cff(:,:)=(1.0007+3.46E-6*mydata(:,:,2))*6.1121*                   &
     &        EXP(17.502*mydata(:,:,1)/(240.97+mydata(:,:,1)))
    cff(:,:)=cff(:,:)*mydata(:,:,3)/100.0
    mydata(:,:,3)=0.62197*(cff(:,:)/(mydata(:,:,2)-0.378*cff(:,:)))
    where(mydata(:,:,3).lt.0.0) mydata(:,:,3)=0.0
END SUBROUTINE humid

SUBROUTINE create_output(year,nx,ny,ncids,varids)
    USE NETCDF
    implicit none
    integer, intent(in) :: year,nx,ny
    integer, intent(out) :: ncids(8),varids(9)  

    integer :: i,varid,nt
    integer :: date(3),date0(3),julian
    integer :: xdimid,ydimid,tdimid,dimids3(3)
    real(kind=8), allocatable :: time(:)
    character(len=10) :: names(8),vnam(9),tdim(8)
    character(len=30) :: outfile,units(9)
    character(len=100) :: longnames(9),tnames(8)
    logical :: leap

    leap = .false.
    if ( mod(year,4).eq.0 ) leap = .true.
    if ( mod(year,100).eq.0 ) leap = .false.
    if ( mod(year,400).eq.0 ) leap = .true.

    if (leap) then
        nt = 366*24
    else
        nt = 365*24
    end if
    allocate( time(nt) )

    names(1) = "lwrad_down"
    names(2) = "swrad_down"
    names(3) = "swrad"
    names(4) = "rain"
    names(5) = "Tair"
    names(6) = "Pair"
    names(7) = "Qair"
    names(8) = "wind"

    tdim(1) = "lrf_time"
    tdim(2) = "srf_time"
    tdim(3) = "srf_time"
    tdim(4) = "rain_time"
    tdim(5) = "tair_time"
    tdim(6) = "pair_time"
    tdim(7) = "qair_time"
    tdim(8) = "wind_time"

    tnames(1) = "time for heat flux"
    tnames(2) = "time for heat flux"
    tnames(3) = "time for heat flux"
    tnames(4) = "time for surface air temperature"
    tnames(5) = "time for surface air temperature"
    tnames(6) = "time for surface air pressure"
    tnames(7) = "time for surface air specific humidity"
    tnames(8) = "10-meter wind time"

    vnam(1) = "lwrad_down"
    vnam(2) = "swrad_down"
    vnam(3) = "swrad"
    vnam(4) = "rain"
    vnam(5) = "Tair"
    vnam(6) = "Pair"
    vnam(7) = "Qair"
    vnam(8) = "Uwind"
    vnam(9) = "Vwind"

    longnames(1) = "longwave downward radiation"
    longnames(2) = "shortwave downward radiation"
    longnames(3) = "solar shortwave radiation"
    longnames(4) = "rain fall rate"
    longnames(5) = "surface air temperature"
    longnames(6) = "surface air pressure"
    longnames(7) = "surface air specific humidity"
    longnames(8) = "10-meter u-wind component"
    longnames(9) = "10-meter v-wind component"

    units(1) = "Watts meter-2"
    units(2) = "Watts meter-2"
    units(3) = "Watts meter-2"
    units(4) = "kilogram meter-2 sec"
    units(5) = "Celsius"
    units(6) = "millibar"
    units(7) = "kg/kg"
    units(8) = "meter second-1"
    units(9) = "meter second-1"
  
    date0=(/1968,5,23/)
    date=(/year,1,1/)
    CALL elapsed(date,date0,julian)
    time(1) = real(julian,8)+1.0/24.0
    do i = 2, nt
      time(i) = time(i-1)+1.0/24.0
    end do
 
    102 FORMAT('baltic_',a,'_',i4.4,'.nc')

    do i = 1,8
      write(outfile,102) trim(names(i)),year

      CALL check( nf90_create( trim(outfile),NF90_NETCDF4,ncids(i) ), 500 )

      CALL check( nf90_def_dim( ncids(i), 'xi_rho', nx, xdimid ), 501 )
      CALL check( nf90_def_dim( ncids(i), 'eta_rho', ny, ydimid ), 502 )
      CALL check( nf90_def_dim( ncids(i), trim(tdim(i)), NF90_UNLIMITED, tdimid ),503)

      dimids3 = (/ xdimid, ydimid, tdimid /)

      CALL check(nf90_def_var( ncids(i), trim(tdim(i)), NF90_DOUBLE,&
        (/tdimid/), varid ), 504)
      CALL check( nf90_put_att( ncids(i), varid, 'long_name',&
        trim(tnames(i)) ), 505 )
      CALL check( nf90_put_att( ncids(i), varid, 'units',&
        "days since 1968-05-23 00:00:00 GMT" ), 506 )
      CALL check( nf90_put_att( ncids(i), varid, 'calendar',&
        "gregorian" ), 507 )

      CALL check(nf90_def_var( ncids(i), trim(vnam(i)), NF90_FLOAT,&
        dimids3, varids(i) ), 508)
      CALL check( nf90_put_att( ncids(i), varids(i), 'long_name',&
        trim(longnames(i)) ), 509 )
      CALL check( nf90_put_att( ncids(i), varids(i), 'units',&
        trim(units(i)) ), 510 )
      CALL check( nf90_put_att( ncids(i), varids(i), 'time',&
        trim(tdim(i)) ), 511 )

      if (i.eq.8) then
        CALL check(nf90_def_var( ncids(i), trim(vnam(9)), NF90_FLOAT,&
          dimids3, varids(9) ), 512)
        CALL check( nf90_put_att( ncids(i), varids(9), 'long_name',&
          trim(longnames(9)) ), 513 )
        CALL check( nf90_put_att( ncids(i), varids(9), 'units',&
          trim(units(9)) ), 514 )
        CALL check( nf90_put_att( ncids(i), varids(9), 'time',&
          trim(tdim(i)) ), 515 )
      end if

      CALL check( nf90_enddef(ncids(i)), 516 )

      CALL check( nf90_put_var( ncids(i), varid, time ), 517 )

    end do

    deallocate(time)
END SUBROUTINE create_output

SUBROUTINE write_output(nx,ny,t,output,ncids,varids)
    USE NETCDF
    implicit none
    integer, intent(in) :: nx,ny,t,ncids(8),varids(9)
    real(kind=8) :: output(nx,ny,9)

    integer :: i

    do i = 1,8
      CALL check( nf90_put_var( ncids(i), varids(i), output(:,:,i),&
        start=(/1,1,t/), count=(/nx,ny,1/) ), 600 )
    end do
    CALL check( nf90_put_var( ncids(8), varids(9), output(:,:,9),&
      start=(/1,1,t/), count=(/nx,ny,1/) ), 600 )
 
END SUBROUTINE write_output

SUBROUTINE write_rain(nx,ny,output,idx,nt,ncid,varid)
    USE NETCDF
    implicit none
    integer, intent(in) :: nx,ny,idx,nt,ncid,varid
    real(kind=8) :: output(nx,ny)

    integer :: i

    if (idx.eq.1) then
      do i = 1,6
        CALL check( nf90_put_var( ncid, varid, output/86400.0,&
          start=(/1,1,i/) ), 600 )
      end do
    elseif (idx.eq.nt+1) then
      do i = 1,18
        CALL check( nf90_put_var( ncid, varid, output/86400.0,&
          start=(/1,1,(idx-2)*24+6+i/) ), 600 )
      end do
    else
      do i = 1,24
        CALL check( nf90_put_var( ncid, varid, output/86400.0,&
          start=(/1,1,(idx-2)*24+6+i/) ), 600 )
      end do
    end if

END SUBROUTINE write_rain

SUBROUTINE longwave(nx,ny,tair,cloud,output)
    implicit none
    integer, intent(in) :: nx,ny
    real(kind=8), intent(in) :: tair(nx,ny),cloud(nx,ny)
    real(kind=8), intent(out) :: output(nx,ny)

    real(kind=8) :: stefan_boltzmann,Tffresh

    stefan_boltzmann = 567.0e-10
    Tffresh = 273.15

    output(:,:) = stefan_boltzmann*(tair(:,:)+Tffresh)**4 &
             * (1.0 - 0.261 &
             * exp(-7.77e-4*tair(:,:)**2)) &
             * (1.0 + 0.275*(cloud(:,:)/100.0))
END SUBROUTINE longwave

SUBROUTINE shortwave(nx,ny,year,n,h,lon,lat,qair,cloud,output)
    implicit none
    integer, intent(in) :: nx,ny,year,n,h
    real(kind=8), intent(in) :: lon(nx,ny),lat(nx,ny)
    real(kind=8), intent(in) :: qair(nx,ny),cloud(nx,ny)
    real(kind=8), intent(out) :: output(nx,ny)

    integer :: i,j
    integer :: mhr,yday,dayyr,ydayUTC
    real(kind=8) :: pi,deg2rad
    real(kind=8) :: hour_angle,solar_time,declin,cosZ,e,d
    logical :: isleap

    pi = 4.0*ATAN(1.0) 
    deg2rad = pi/180.0

    dayyr = 365
    CALL leap(year,isleap)
    if(isleap) dayyr = 366

    mhr = modulo(modulo((n-1)*6,24)+h,24)
    ydayUTC = ceiling(real(n,8)/4.0)
    if (mhr.eq.0) ydayUTC = ydayUTC+1
    if (ydayUTC.gt.dayyr) then
      ydayUTC = 1
      dayyr = 365
      CALL leap(year+1,isleap)
      if(isleap) dayyr = 366
    end if

    do i = 1,nx
    do j = 1,ny
      yday = ydayUTC
      solar_time = real(mhr,8)+lon(i,j)/15.0
      if (solar_time .ge. 24.0) then
        yday = ydayUTC+1
        if (yday.gt.dayyr) then
          yday = 1
          dayyr = 365
          CALL leap(year+1,isleap)
          if(isleap) dayyr = 366
        end if
        solar_time = solar_time - 24.0
      end if
      hour_angle = (12.0 - solar_time)*pi/12.0

      declin = 23.44*cos((172.0-yday)*2.0*pi/dayyr)*deg2rad 
      cosZ = sin(lat(i,j)*deg2rad)*sin(declin) &
             + cos(lat(i,j)*deg2rad)*cos(declin)*cos(hour_angle)
      cosZ = max(cosZ,0.0)
      e = 1.e5*qair(i,j)/(0.622 + 0.378*qair(i,j))
      d = (cosZ+2.7)*e*1.e-5+1.085*cosZ+0.1
      output(i,j) = 1353.0*cosZ**2/d
      output(i,j) = max(output(i,j),0.0)
      output(i,j) = output(i,j)*(1.0-0.6*(cloud(i,j)/100.0)**3)
    end do
    end do
END SUBROUTINE

SUBROUTINE temporal_interp(nx,ny,input,output)
    USE NETCDF
    implicit none
    integer, intent(in) :: nx,ny
    real(kind=8), intent(in) :: input(nx,ny,2) 
    real(kind=8), intent(out) :: output(nx,ny,6)

    integer :: i
    real(kind=8) :: dx(nx,ny)    

    output(:,:,6) = input(:,:,2)
    dx(:,:) = (input(:,:,2)-input(:,:,1))/6

    output(:,:,1) = input(:,:,1)+dx(:,:)
    do i = 2,5
      output(:,:,i) = output(:,:,i-1)+dx(:,:) 
    end do 
END SUBROUTINE temporal_interp

SUBROUTINE close_output(ncids)
    USE NETCDF
    implicit none
    integer, intent(in) :: ncids(8)

    integer :: i

    do i = 1,8
      CALL check( nf90_close(ncids(i)), 600+i )
    end do
END SUBROUTINE close_output

SUBROUTINE read_data(nx,ny,year,mon,n,l,maxn,input)
    USE NETCDF
    implicit none
    integer, intent(in) :: nx,ny,year,mon,n,l,maxn
    real(kind=8), intent(out) :: input(nx,ny,9)

    integer :: ncid,varid,i,year2,mon2,ncid2,varid2
    real(kind=8) :: tmp2d(nx,ny,3),ave2d(nx,ny,2)
    character(len=10) :: var(9),varname(9)
    character(len=50) :: filename
    logical :: ex,mask(nx,ny)

    var(1) = 'lwrad_down'
    var(2) = 'swrad_down'
    var(3) = 'swrad'
    var(4) = 'rain'
    var(5) = 'Tair'
    var(6) = 'Pair'
    var(7) = 'Qair'
    var(8) = 'wdir'
    var(9) = 'wspeed'

    varname(1) = 'strd'
    varname(2) = 'ssrd'
    varname(3) = 'ssr'
    varname(4) = 'tp'
    varname(5) = 't2m'
    varname(6) = 'msl'
    varname(7) = 'r2'
    varname(8) = 'wdir10'
    varname(9) = 'si10'

    101 format (a,'_',i4.4,i2.2,'.nc')
    do i = 1,4

      write(filename,101) trim(var(i)),year,mon
      CALL check(nf90_open(trim(filename),NF90_NOWRITE,ncid),710)
      CALL check(nf90_inq_varid(ncid,trim(varname(i)),varid),419)
      CALL check(nf90_get_var(ncid,varid,tmp2d(:,:,2),&
            start=(/1,1,l,n/),count=(/nx,ny,1,1/)),420)

      if ( l.lt.6) then
        CALL check(nf90_get_var(ncid,varid,tmp2d(:,:,3),&
            start=(/1,1,l+1,n/),count=(/nx,ny,1,1/)),420)
        ave2d(:,:,2) = (tmp2d(:,:,3)-tmp2d(:,:,2))/3600.0

        if (l.eq.1) then
          ave2d(:,:,1) = tmp2d(:,:,2)/3600.0 
        else
          CALL check(nf90_get_var(ncid,varid,tmp2d(:,:,1),&
              start=(/1,1,l-1,n/),count=(/nx,ny,1,1/)),420)
          ave2d(:,:,1) = (tmp2d(:,:,2)-tmp2d(:,:,1))/3600.0
        end if

      else
        CALL check(nf90_get_var(ncid,varid,tmp2d(:,:,1),&
               start=(/1,1,l-1,n/),count=(/nx,ny,1,1/)),420)
        ave2d(:,:,1) = (tmp2d(:,:,2)-tmp2d(:,:,1))/3600.0

        if (n.lt.maxn) then
          CALL check(nf90_get_var(ncid,varid,tmp2d(:,:,3),&
            start=(/1,1,1,n+1/),count=(/nx,ny,1,1/)),420)
          ave2d(:,:,2) = tmp2d(:,:,3)/3600.0 
        else
          mon2 = mon+1
          if (mon2.gt.12) then
            mon2 = 1
            year2 = year+1
          else
            year2 = year
          end if

          write(filename,101) trim(var(i)),year2,mon2
          INQUIRE(FILE=trim(filename),EXIST=ex)
          if (ex) then
            CALL check(nf90_open(trim(filename),NF90_NOWRITE,ncid2),711)
            CALL check(nf90_inq_varid(ncid2,trim(varname(i)),varid2),419)
            CALL check(nf90_get_var(ncid2,varid2,tmp2d(:,:,3),&
                start=(/1,1,1,1/),count=(/nx,ny,1,1/)),420)
            CALL check(nf90_close(ncid2),460)
            ave2d(:,:,2) = tmp2d(:,:,3)/3600.0 
          else
            ave2d(:,:,2) = ave2d(:,:,1)
          end if 
        end if
        
      end if

      input(:,:,i) = (ave2d(:,:,1)+ave2d(:,:,2))*0.5
 
      CALL check(nf90_close(ncid),460)
    end do

    do i = 5,9

      write(filename,101) trim(var(i)),year,mon

      CALL check(nf90_open(trim(filename),NF90_NOWRITE,ncid),712)
      CALL check(nf90_inq_varid(ncid,trim(varname(i)),varid),419)
      CALL check(nf90_get_var(ncid,varid,input(:,:,i),&
            start=(/1,1,l,n/),count=(/nx,ny,1,1/)),420)
      CALL check(nf90_close(ncid),460)
    end do

END SUBROUTINE read_data


SUBROUTINE read_rain(nx,ny,year,idx,input)
    USE NETCDF
    implicit none
    integer, intent(in) :: nx,ny,year,idx
    real(kind=8), intent(out) :: input(nx,ny)

    integer :: ncid,varid
    character(len=50) :: filename

    101 format ('UERRA_rain_',i4.4,'.nc')
    write(filename,101) year
    CALL check(nf90_open(trim(filename),NF90_NOWRITE,ncid),713)
    CALL check(nf90_inq_varid(ncid,'tp',varid),419)
    CALL check(nf90_get_var(ncid,varid,input,&
        start=(/500,550,idx/),count=(/nx,ny,1/)),420)
    CALL check(nf90_close(ncid),460)
END SUBROUTINE read_rain

SUBROUTINE leap(year,answer)
    implicit none
    integer, intent(in) :: year
    logical, intent(out) :: answer

    answer = .false.
    if( mod(year,4).eq.0 ) answer = .true.
    if( mod(year,4).eq.100 ) answer = .false.
    if( mod(year,4).eq.400 ) answer = .true.
END SUBROUTINE leap

SUBROUTINE elapsed(date,date_0,days)
    implicit none
    integer, intent(in) :: date(3),date_0(3)
    integer, intent(out) :: days

    integer :: date_tmp(3)

    date_tmp(1) = date_0(1)
    date_tmp(2) = date_0(2)
    date_tmp(3) = date_0(3)

    days = 0

    do
        if( date_tmp(1).eq.date(1) .and.&
            date_tmp(2).eq.date(2) .and.&
            date_tmp(3).eq.date(3) ) EXIT
        days = days+1
        CALL add_day(date_tmp,.true.)
    end do
END SUBROUTINE elapsed

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

SUBROUTINE prepare(nx,ny,onx,ony,lon,lat,olon,olat,near,W,files)
    implicit none
    integer, intent(in) :: nx,ny,onx,ony
    real(kind=8), intent(in) :: lon(nx,ny),lat(nx,ny),olon(onx,ony),olat(onx,ony)
    integer, intent(out) :: near(onx,ony,2)
    character(len=20), intent(in) :: files(2) 
    real(kind=8), intent(out) :: W(onx,ony,4)
 
    integer :: i,j,reclen
    real(kind=8) :: point(2)
    logical :: ex1,ex2

    INQUIRE(FILE=trim(files(1)),EXIST=ex1)
    INQUIRE(FILE=trim(files(2)),EXIST=ex2)

    if (ex1 .and. ex2) then

      inquire(iolength = reclen) near

      OPEN(10,FILE=trim(files(1)),ACCESS='DIRECT', recl=reclen, FORM='UNFORMATTED')
      read(10, rec=1) near
      CLOSE(10)

      inquire(iolength = reclen) W

      OPEN(10,FILE=trim(files(2)),ACCESS='DIRECT', recl=reclen, FORM='UNFORMATTED')
      read(10, rec=1) W
      CLOSE(10)

    else

      do i = 1,onx
        do j = 1,ony
          point = (/olon(i,j),olat(i,j)/)
          CALL find_idx(nx,ny,lon,lat,point,near(i,j,:),W(i,j,:))
        end do
      end do

      inquire(iolength = reclen) near

      OPEN(10,FILE=trim(files(1)),ACCESS='DIRECT', recl=reclen, FORM='UNFORMATTED')
      write(10, rec=1) near
      CLOSE(10)

      inquire(iolength = reclen) W

      OPEN(10,FILE=trim(files(2)),ACCESS='DIRECT', recl=reclen, FORM='UNFORMATTED')
      write(10, rec=1) W
      CLOSE(10)

    end if 
 
END SUBROUTINE prepare
    
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

SUBROUTINE interp(nx,ny,onx,ony,input,idx,W,output)
    implicit none
    integer, intent(in) :: nx,ny,onx,ony,idx(onx,ony,2)
    real(kind=8), intent(in) :: input(nx,ny),W(onx,ony,2)
    real(kind=8), intent(out) :: output(onx,ony)

    integer :: i,j,x
    real(kind=8) :: xarg(nx),yarg(ny)
    real(kind=8) :: splint,xp,yp
    real(kind=8) :: y2a(nx,ny),line(nx),y2line(nx)

    do i = 1,nx
      xarg(i) = i
    end do
    do j = 1,ny
      yarg(j) = j
    end do


    CALL splie2(yarg,nx,ny,input,y2a)
    do i = 1, onx
      do j = 1, ony
        xp = real(idx(i,j,1),8)+W(i,j,1)
        yp = real(idx(i,j,2),8)+W(i,j,2)
        do x = 1, nx
          line(x) = splint(yarg,input(x,:),y2a(x,:),ny,yp)
        end do
        call spline(xarg,line,nx,y2line)
        output(i,j) = splint(xarg,line,y2line,nx,xp)
      end do
    end do
END SUBROUTINE interp

SUBROUTINE interp2(nx,ny,onx,ony,input,idx,W,output)
    implicit none
    integer, intent(in) :: nx,ny,onx,ony,idx(onx,ony,2)
    real(kind=8), intent(in) :: input(nx,ny,9),W(onx,ony,4)
    real(kind=8), intent(out) :: output(onx,ony,9)
    
    integer :: i,j,v,x,y
    real(kind=8) :: a(4)

    do i = 1,onx
    do j = 1,ony
      x = idx(i,j,1)
      y = idx(i,j,2)
      a(1) = W(i,j,1)
      a(2) = W(i,j,2)
      a(3) = W(i,j,3)
      a(4) = W(i,j,4)
      do v = 1,9
        output(i,j,v) = a(1)*input(x,y,v)+a(2)*input(x+1,y,v)+&
          a(3)*input(x+1,y+1,v)+a(4)*input(x,y+1,v)
      end do
    end do
    end do

END SUBROUTINE interp2


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

SUBROUTINE find_idx(ni,nj,lon,lat,pt,idx,alphas)
    implicit none
    integer, intent(in) :: ni,nj
    real(kind=8), intent(in) :: lon(ni,nj),lat(ni,nj),pt(2)
    integer, intent(out) :: idx(2)
    real(kind=8), intent(out) :: alphas(4)

    real(kind=8) :: disarray(ni,nj),corners(4,2)
    integer :: tmp(2),k,l,i,j
    logical :: fexit,score

    do i = 1,ni
    do j = 1,nj
      disarray(i,j) = sqrt((lon(i,j)-pt(1))**2+(lat(i,j)-pt(2))**2)
    end do
    end do
    tmp = minloc(disarray)

    fexit = .false.
    do k = tmp(1)-1,tmp(1)
        if (fexit) exit
        do l = tmp(2)-1,tmp(2)
            corners(1,:) = (/lon(k,l),lat(k,l)/)
            corners(2,:) = (/lon(k+1,l),lat(k+1,l)/)
            corners(3,:) = (/lon(k+1,l+1),lat(k+1,l+1)/)
            corners(4,:) = (/lon(k,l+1),lat(k,l+1)/)
            CALL in_convex_polygon(corners,pt,score)
            if (score) then
                idx = (/k,l/)
                CALL calc_w(corners,pt,alphas)
                fexit = .true.
                exit
            endif
        end do
    end do

END SUBROUTINE find_idx

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
