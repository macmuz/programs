PROGRAM flowrate
    USE netcdf
    implicit none
    integer :: ncid, varid, dimid, nx, ny, nz,nc
    integer :: i,j,k,v_j,v_istr,v_iend,vid(3),n,cnt
    integer :: xdimid,ydimid,zdimid,tdimid,dimids2(2)
    integer :: date(3),datestp(3),nrec,pts(2,2) 
    real(kind=8) :: hc, dx, myfx, a, b, S, time, dif
    real(kind=8), allocatable :: h(:,:),cs_w(:),zeta(:,:),depth(:,:,:) 
    real(kind=8), allocatable :: mask(:,:),s_w(:),depth_psi(:,:,:),mask_v(:,:)
    real(kind=8), allocatable :: v(:,:,:), pm(:,:), myflux(:,:)
    real(kind=8), allocatable :: zdif(:), frate(:)
    character(len=200) :: filename,path
    character(len=50) :: tunits,tcal

    100 format(a,'/ocean_avg_',i4,'-',i2.2,'-',i2.2,'.nc')
    !path = '/users/work/mmuzyka/CSDIR/metro_560x600_era2004v5/run/baltic'
    path = '/users/work/mmuzyka/CSDIR/metro_560x600_uerraGLS/run/baltic'
!    path = '/users/magazyn/mmuzyka/ROMS_CASES/metro_560x600_era5MYtest'

    v_j = 99
    v_istr = 104
    v_iend = 111
!    v_j = 221
!    v_istr = 255
!    v_iend = 271

    nrec = 4

    pts(1,:) = (/98,86/)
    pts(2,:) = (/101,125/)
!    pts(1,:) = (/238,187/)
!    pts(2,:) = (/244,285/)

    date = (/1993,1,2/)
    datestp = (/1993,11,15/)

    write(filename,100) trim(path),date(1),date(2),date(3)
!READ GRID
    CALL check(nf90_open(trim(filename),NF90_NOWRITE,ncid),310)

    CALL check(nf90_inq_dimid(ncid, "xi_rho", dimid),311)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nx),312)

    CALL check(nf90_inq_dimid(ncid, "eta_rho", dimid),313)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=ny),314)

    CALL check(nf90_inq_dimid(ncid, "s_rho", dimid),315)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nz),316)

    allocate( h(nx,ny), cs_w(nz+1), zeta(nx,ny), depth(nx,ny,nz+1) )
    allocate( mask(nx,ny), s_w(nz+1), depth_psi(nx-1,ny-1,nz+1) )
    allocate( v(nx,ny-1,nz), pm(nx,ny), myflux(nx-2,ny-1), mask_v(nx,ny-1) )

    CALL check(nf90_inq_varid(ncid,"hc",varid),317)
    CALL check(nf90_get_var(ncid,varid,hc),318)

    CALL check(nf90_inq_varid(ncid,"Cs_w",varid),319)
    CALL check(nf90_get_var(ncid,varid,cs_w),320)

    CALL check(nf90_inq_varid(ncid,"s_w",varid),319)
    CALL check(nf90_get_var(ncid,varid,s_w),320)

    CALL check(nf90_inq_varid(ncid,"h",varid),319)
    CALL check(nf90_get_var(ncid,varid,h),320)

    CALL check(nf90_inq_varid(ncid,"pm",varid),319)
    CALL check(nf90_get_var(ncid,varid,pm),320)

    CALL check(nf90_inq_varid(ncid,"zeta",varid),319)
    CALL check(nf90_get_var(ncid,varid,zeta,start=(/1,1,2/)),320)

    CALL check(nf90_inq_varid(ncid,"v",varid),319)
    CALL check(nf90_get_var(ncid,varid,v,start=(/1,1,1,2/)),320)

    CALL check(nf90_inq_varid(ncid,"mask_rho",varid),319)
    CALL check(nf90_get_var(ncid,varid,mask),320)

    CALL check(nf90_inq_varid(ncid,"mask_v",varid),319)
    CALL check(nf90_get_var(ncid,varid,mask_v),320)

    CALL check(nf90_inq_varid(ncid,"ocean_time",varid),319)
    CALL check(nf90_get_att(ncid,varid,'units',tunits),320)
    CALL check(nf90_get_att(ncid,varid,'calendar',tcal),320)

    CALL check(nf90_close(ncid),360)
!END READ GRID

    CALL check( nf90_create( 'flowrate_uerraGLS_new.nc',NF90_NETCDF4,ncid ), 500 )

    CALL check( nf90_def_dim( ncid, 'time', NF90_UNLIMITED, tdimid ), 501 )

    CALL check(nf90_def_var( ncid, 'time', NF90_DOUBLE,&
        (/tdimid/), vid(1) ), 508)
    CALL check( nf90_put_att( ncid, vid(1), 'units',&
        trim(tunits) ), 305 )
    CALL check( nf90_put_att( ncid, vid(1), 'calendar',&
        trim(tcal) ), 306 )

    CALL check(nf90_def_var( ncid, 'zetadif', NF90_DOUBLE,&
        (/tdimid/), vid(2) ), 508)
    CALL check( nf90_put_att( ncid, vid(2), 'longname',&
        'SSH diference' ), 305 )
    CALL check( nf90_put_att( ncid, vid(2), 'units',&
        'm' ), 305 )

    CALL check(nf90_def_var( ncid, 'flowrate', NF90_DOUBLE,&
        (/tdimid/), vid(3) ), 508)
    CALL check( nf90_put_att( ncid, vid(3), 'longname',&
        'volume flow rate' ), 305 )
    CALL check( nf90_put_att( ncid, vid(3), 'units',&
        'm**3/s' ), 305 )

    CALL check( nf90_enddef(ncid), 516 )

    cnt = 1
    do
        write(filename,100) trim(path),date(1),date(2),date(3)
        write(*,*) trim(filename)

        CALL check(nf90_open(trim(filename),NF90_NOWRITE,nc),310)

        do n = 1,nrec

        CALL check(nf90_inq_varid(nc,"ocean_time",varid),319)
        CALL check(nf90_get_var(nc,varid,time,start=(/n/)),320)

        CALL check(nf90_inq_varid(nc,"v",varid),319)
        CALL check(nf90_get_var(nc,varid,v,start=(/1,1,1,n/)),320)

        CALL check(nf90_inq_varid(nc,"zeta",varid),319)
        CALL check(nf90_get_var(nc,varid,zeta,start=(/1,1,n/)),320)

    do k = 1,nz
      where(mask_v.eq.0) v(:,:,k)=0.0
    end do
    
    where(mask.lt.0.5) zeta = 0.0
    do i = 1, nx   
      do j = 1, ny   
        do k = 1, nz+1   
          depth(i,j,k) = zeta(i,j)+(zeta(i,j)+h(i,j))*(hc*s_w(k)+h(i,j)*cs_w(k))/(hc+h(i,j))
        end do
      end do
    end do

    do i = 1, nx-1
      do j = 1, ny-1
        depth_psi(i,j,:) = 0.25*(depth(i,j,:)+depth(i+1,j,:)+&
            depth(i,j+1,:)+depth(i+1,j+1,:))
      end do
    end do

    myfx = 0.0
    do i = v_istr, v_iend
      dx = 0.5*(1/pm(i,v_j)+1/pm(i,v_j+1))
!      write(*,*) dx 
      do k = 1, nz
        a = depth_psi(i-1,v_j,k+1)-depth_psi(i-1,v_j,k) 
        b = depth_psi(i,v_j,k+1)-depth_psi(i,v_j,k) 
!        write(*,*) a,b
        S = 0.5*(a+b)*dx
!        write(*,*) S
        myfx = myfx + S*v(i,v_j,k)
      end do
    end do

    dif = zeta(pts(1,1),pts(1,2))-zeta(pts(2,1),pts(2,2)) 

    write(*,"(a,i3,a,i3,a,i3,a,f9.2,a)") "flux for i=",v_istr,":",v_iend,&
        " and j=",v_j," towards the north: ",myfx," m**3/s"

!    write(*,*) v(v_istr-1,v_j,1),v(v_istr,v_j,1),v(v_iend,v_j,1),v(v_iend+1,v_j,1)


!    myflux = 0.0
!    do i = 1, nx-2
!      do j = 1, ny-1
!        dx = 0.5*(1/pm(i+1,j)+1/pm(i+1,j+1))
!        do k = 1, nz
!          a = depth_psi(i,j,k+1)-depth_psi(i,j,k)
!          b = depth_psi(i+1,j,k+1)-depth_psi(i+1,j,k)
!          S = 0.5*(a+b)*dx
!          myflux(i,j) = myflux(i,j) + S*v(i+1,j,k)
!        end do
!      end do
!    end do

!    CALL check( nf90_create( 'flux.nc',NF90_NETCDF4,ncid ), 500 )
!    CALL check( nf90_def_dim( ncid, 'nx', nx-2, xdimid ), 501 )
!    CALL check( nf90_def_dim( ncid, 'ny', ny-1, ydimid ), 502 )
!    dimids2 = (/ xdimid, ydimid /)
!    CALL check(nf90_def_var( ncid, 'nflux', NF90_DOUBLE,&
!        dimids2, varid ), 508)
!    CALL check( nf90_put_att( ncid, varid, 'long_name',&
!        'flux towards the north' ), 305 )
!    CALL check( nf90_put_att( ncid, varid, 'units',&
!        'm**3/s' ), 306 )
!    CALL check( nf90_enddef(ncid), 516 )
!    CALL check( nf90_put_var( ncid, varid, myflux ), 517 )
!    CALL check( nf90_close(ncid), 600 )



        CALL check( nf90_put_var( ncid, vid(1), time, start=(/cnt/)), 517 )
        CALL check( nf90_put_var( ncid, vid(2), dif, start=(/cnt/)), 517 )
        CALL check( nf90_put_var( ncid, vid(3), myfx, start=(/cnt/)), 517 )
        
        cnt = cnt+1

        end do
        CALL check(nf90_close(nc),360)

        if (date(1).ge.datestp(1).and.&
            date(2).ge.datestp(2).and.&
            date(3).ge.datestp(3)) EXIT
        CALL add_day(date,.true.)
    end do
    CALL check( nf90_close(ncid), 600 ) 


    deallocate( h, cs_w, zeta, depth, mask, s_w, depth_psi, v, pm, myflux )

    allocate( zdif(cnt-1), frate(cnt-1) )

    CALL check(nf90_open('flowrate.nc',NF90_NOWRITE,ncid),310)

    CALL check(nf90_inq_varid(ncid,"zetadif",varid),319)
    CALL check(nf90_get_var(ncid,varid,zdif),320)

    CALL check(nf90_inq_varid(ncid,"flowrate",varid),319)
    CALL check(nf90_get_var(ncid,varid,frate),320)

    CALL check( nf90_close(ncid), 600 )
    
    zdif = zdif - sum(zdif)/(cnt-1)
    frate = frate - sum(frate)/(cnt-1)
 
    open (unit = 10, file = "fr.txt", form="FORMATTED")
    do i = 1,cnt-1
        write(10,"(f7.4,a,f)") zdif(i),achar(9),frate(i)
    end do
    close (10)

    deallocate( zdif, frate )
END PROGRAM flowrate

SUBROUTINE check(status,label)
    USE NETCDF
    implicit none
    INTEGER, INTENT(in) :: status , label

    IF (status .ne. nf90_noerr) then
        print*, trim(nf90_strerror(status)),' label= ',label
        stop "Stopped"
    END IF
END SUBROUTINE check

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
