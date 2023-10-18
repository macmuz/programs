PROGRAM NEMO
    USE NETCDF
    implicit none
    integer :: ncid,dimid,vid,ni,nj,nk,ncid2
    integer :: nt,np,nz,varid(6)
    integer :: tdimid,xdimid,zdimid
    integer :: k,narg,id(2),fvi
    integer, allocatable :: tmp2d(:,:)
    real(kind=8), allocatable :: depth(:),lon(:),lat(:)
    real, allocatable :: pos_lon(:),pos_lat(:),time(:)
    real, allocatable :: temp(:,:),salt(:,:),deph(:,:)
    real ::  fvd, af, sf
    character(len=200) :: path,by5file,by15file,filen,fileout

    path = '/users/work/mmuzyka/programs/ROMS_bc/new_files/1993/01/&
            &BAL-MYP-NEMO_PHY-DailyMeans-19930101.nc'

!READ command line arg: YEAR
    narg = command_argument_count()
    if (narg.ne.1) then
      write(*,*) "Program must have filename as argument"
      stop
    end if
    call get_command_argument(1,filen)
    id(1) = index(filen,'/')
    id(2) = index(filen,'.nc')
    write(fileout,'(A,A)') trim(filen(id(1)+1:id(2)-1)),'_NEMO.nc'
!READ command line    
    write(*,*) trim(filen)
    write(*,*) trim(fileout)

    CALL check(nf90_open(trim(path), NF90_NOWRITE, ncid),310)

    CALL check(nf90_inq_dimid(ncid, "depth", dimid),313)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nk),314)

    CALL check(nf90_inq_dimid(ncid, "lon", dimid),313)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=ni),314)

    CALL check(nf90_inq_dimid(ncid, "lat", dimid),313)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nj),314)

    allocate( depth(nk), lon(ni), lat(nj) )

    CALL check(nf90_inq_varid(ncid,"depth",vid),31)
    CALL check(nf90_get_var(ncid,vid,depth),32)

    CALL check(nf90_inq_varid(ncid,"lon",vid),31)
    CALL check(nf90_get_var(ncid,vid,lon),32)

    CALL check(nf90_inq_varid(ncid,"lat",vid),31)
    CALL check(nf90_get_var(ncid,vid,lat),32)

    CALL check(nf90_close(ncid),360)

!####read file
    CALL check(nf90_open(trim(filen), NF90_NOWRITE, ncid),310)

    CALL check(nf90_inq_dimid(ncid, "TIME", dimid),314)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nt),314)

    CALL check(nf90_inq_dimid(ncid, "POSITION", dimid),315)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=np),314)

    CALL check(nf90_inq_dimid(ncid, "DEPTH", dimid),316)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nz),314)

    allocate( pos_lon(np), pos_lat(np), time(nt), tmp2d(nz,np) )
    allocate( temp(nk,np), salt(nk,np), deph(nz,np) )
    temp = -999.0
    salt = -999.0

    CALL check(nf90_inq_varid(ncid,"TIME",vid),33)
    CALL check(nf90_get_var(ncid,vid,time),34)

    CALL check(nf90_inq_varid(ncid,"LONGITUDE",vid),33)
    CALL check(nf90_get_var(ncid,vid,pos_lon),34)

    CALL check(nf90_inq_varid(ncid,"LATITUDE",vid),33)
    CALL check(nf90_get_var(ncid,vid,pos_lat),34)

    CALL check(nf90_inq_varid(ncid,"DEPH",vid),33)
    CALL check(nf90_get_var(ncid,vid,deph),34)
    CALL check(nf90_get_att(ncid,vid,'_FillValue',fvd),34)

    CALL check(nf90_inq_varid(ncid,"TEMP",vid),33)
    CALL check(nf90_get_var(ncid,vid,tmp2d),34)
    CALL check(nf90_get_att(ncid,vid,'_FillValue',fvi),34)
    CALL check(nf90_get_att(ncid,vid,'add_offset',af),34)
    CALL check(nf90_get_att(ncid,vid,'scale_factor',sf),34)

    CALL move(nz,np,deph,fvd,tmp2d,fvi,af,sf,nk,depth,temp)

    CALL check(nf90_inq_varid(ncid,"PSAL",vid),33)
    CALL check(nf90_get_var(ncid,vid,tmp2d),34)
    CALL check(nf90_get_att(ncid,vid,'_FillValue',fvi),34)
    CALL check(nf90_get_att(ncid,vid,'add_offset',af),34)
    CALL check(nf90_get_att(ncid,vid,'scale_factor',sf),34)

    CALL move(nz,np,deph,fvd,tmp2d,fvi,af,sf,nk,depth,salt)
!############create file###########
    CALL check(nf90_create(trim(fileout),NF90_NETCDF4,ncid2),310)

    CALL check( nf90_def_dim( ncid2, 'time', NF90_UNLIMITED, tdimid ),503)
    CALL check( nf90_def_dim( ncid2, 'position', np, xdimid ),503)
    CALL check( nf90_def_dim( ncid2, 'depth', nk, zdimid ),503)

    CALL check( nf90_def_var( ncid2, 'time', NF90_DOUBLE, (/tdimid/), varid(1)), 504)
    CALL check(nf90_inq_varid(ncid,'TIME',vid),120)
    CALL check( nf90_copy_att(ncid, vid, 'standard_name', ncid2, varid(1)), 504)
    CALL check( nf90_copy_att(ncid, vid, 'long_name', ncid2, varid(1)), 504)
    CALL check( nf90_copy_att(ncid, vid, 'units', ncid2, varid(1)), 504)
    CALL check( nf90_copy_att(ncid, vid, 'calendar', ncid2, varid(1)),504)
    CALL check( nf90_copy_att(ncid, vid, 'axis', ncid2, varid(1)), 504)
 
    CALL check( nf90_def_var( ncid2, 'pos_lon', NF90_FLOAT, (/xdimid/), varid(2)), 504)
    CALL check(nf90_inq_varid(ncid,'LONGITUDE',vid),120)
    CALL check( nf90_copy_att(ncid, vid, 'standard_name', ncid2, varid(2)), 504)
    CALL check( nf90_copy_att(ncid, vid, 'long_name', ncid2, varid(2)), 504)
    CALL check( nf90_copy_att(ncid, vid, 'units', ncid2, varid(2)), 504)

    CALL check( nf90_def_var( ncid2, 'pos_lat', NF90_FLOAT, (/xdimid/), varid(3)), 504)
    CALL check(nf90_inq_varid(ncid,'LATITUDE',vid),120)
    CALL check( nf90_copy_att(ncid, vid, 'standard_name', ncid2, varid(3)), 504)
    CALL check( nf90_copy_att(ncid, vid, 'long_name', ncid2, varid(3)), 504)
    CALL check( nf90_copy_att(ncid, vid, 'units', ncid2, varid(3)), 504)

    CALL check( nf90_def_var( ncid2, 'depth', NF90_FLOAT, (/zdimid/), varid(4)), 504)
    CALL check( nf90_put_att( ncid2, varid(4), 'standard_name', 'depth'), 504)
    CALL check( nf90_put_att( ncid2, varid(4), 'long_name',&
        'Depth of each measurement'), 504)
    CALL check( nf90_put_att( ncid2, varid(4), 'units', 'meters'), 504)
    CALL check( nf90_put_att( ncid2, varid(4), 'positive', 'down'), 504)
    CALL check( nf90_put_att( ncid2, varid(4), 'axis', 'Z'), 504)

    CALL check( nf90_def_var( ncid2, 'temp', NF90_FLOAT, (/zdimid,tdimid/), varid(5)), 504)
    CALL check( nf90_put_att( ncid2, varid(5), 'standard_name',&
        'sea_water_temperature'), 504)
    CALL check( nf90_put_att( ncid2, varid(5), 'long_name',&
        'Sea temperature'), 504)
    CALL check( nf90_put_att( ncid2, varid(5), 'units', 'degrees_C'), 504)
    CALL check( nf90_put_att( ncid2, varid(5), '_FillValue', -999.0), 504)

    CALL check( nf90_def_var( ncid2, 'salt', NF90_FLOAT, (/zdimid,tdimid/), varid(6)), 504)
    CALL check( nf90_put_att( ncid2, varid(6), 'standard_name',&
        'sea_water_practical_salinity'), 504)
    CALL check( nf90_put_att( ncid2, varid(6), 'long_name',&
        'Practical salinity'), 504)
    CALL check( nf90_put_att( ncid2, varid(6), 'units', '0.001'), 504)
    CALL check( nf90_put_att( ncid2, varid(6), '_FillValue', -999.0), 504)

    CALL check( nf90_enddef(ncid2), 516 )
    
    CALL check(nf90_put_var(ncid2,varid(1),time),18)
    CALL check(nf90_put_var(ncid2,varid(2),pos_lon),18)
    CALL check(nf90_put_var(ncid2,varid(3),pos_lat),18)
    CALL check(nf90_put_var(ncid2,varid(4),real(depth,4)),18)
    CALL check(nf90_put_var(ncid2,varid(5),temp),18)
    CALL check(nf90_put_var(ncid2,varid(6),salt),18)

    CALL check(nf90_close(ncid2),360)
!##################################
    CALL check(nf90_close(ncid),360)
   
!    do k = 1,nk 
!      write(*,*) depth(k)
!    end do
 

    deallocate( depth, lon, lat, deph, tmp2d ) 
    deallocate( pos_lon, pos_lat, time, temp, salt )
     
END PROGRAM NEMO

SUBROUTINE check(status,label)
    USE NETCDF
    implicit none
    INTEGER, INTENT(in) :: status , label

    IF (status .ne. nf90_noerr) then
        print*, trim(nf90_strerror(status)),' label= ',label
        stop "Stopped"
    END IF
END SUBROUTINE check

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

SUBROUTINE move(nz,np,deph,fvd,tmp2d,fvi,af,sf,nk,depth,output)
    implicit none
    integer, intent(in) :: nz,np,fvi,nk
    integer, intent(in) :: tmp2d(nz,np)
    real, intent(in) :: deph(nz,np),fvd,af,sf
    real(kind=8), intent(in) :: depth(nk)
    real, intent(inout) :: output(nk,np)

    integer :: i,k,cnt,kstr,kend
    real(kind=8) :: splint, input(nz), y2(nz)

    do i = 1,np

    cnt = 0
    do k = 1,nz
        if (tmp2d(k,i).eq.fvi .and. cnt.eq.0) cycle
        if (tmp2d(k,i).eq.fvi) then
            exit
        end if
        cnt = cnt+1
    end do
    if (cnt.gt.3) then
      cnt = 0
      do k = 1,nz
        if (tmp2d(k,i).eq.fvi .and. cnt.eq.0) cycle
        if (tmp2d(k,i).eq.fvi) then
            kend = k-1
            exit
        end if
        cnt = cnt+1
        if (cnt.eq.1) kstr = k 
      end do

      do k = kstr,kend
        input(k) = real(tmp2d(k,i),8)*sf+af
      end do

      CALL spline(real(deph(kstr:kend,i),8),input(kstr:kend),&
            1+kend-kstr,y2(kstr:kend))
      do k = 1,nk
        if (depth(k).lt.deph(kstr,i).or.depth(k).gt.deph(kend,i)) cycle
        output(k,i) = splint(real(deph(kstr:kend,i),8),input(kstr:kend),&
            y2(kstr:kend),1+kend-kstr,depth(k))
      end do
    else
      write(*,*) 'MISSING',i
    end if
!    write(*,*) i,kstr,kend
!    write(*,*) real(tmp2d(kstr,i))*sf+af
!    write(*,*) real(tmp2d(kend,i))*sf+af

    end do    

 !   write(*,*) 'sf/af',sf,af,fvi,tmp2d(i,kstr)==fvi
END SUBROUTINE
