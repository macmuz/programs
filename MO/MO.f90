PROGRAM MO
    USE netcdf
    implicit none
    integer, parameter :: npt4d=1,nvar4d=2
    integer :: date(3),datestp(3),cnt,i,j,k,t,n,nold
    integer :: date0(3),day,daystp,daytmp,dateold(3)
    integer :: ncid,dimid,ncid2,tdimid,zdimid
    integer :: ni,nj,ns,nk,nt,vid,varid(6)
    character(len=200) :: filename,path
    character(len=30) :: vname4d(nvar4d),varn,cfile
    real(kind=8) :: dplim(2),r,hc,pt4d(npt4d,2),npt
    real(kind=8) :: AOt,AOs,SFt,SFs,FV,mytime
    integer :: FVt,FVs,idx(2)
    real(kind=8) :: h(2,2),zeta(2,2),alphas(4)
    real(kind=8), allocatable :: time(:),Cs_r(:),s_r(:),lon(:,:),lat(:,:)
    real(kind=8), allocatable :: mt(:,:,:),ms(:,:,:),dp(:,:,:)
    real(kind=8), allocatable :: avet(:),aves(:),depth(:)
    real(kind=4) :: ptlon(1),ptlat(1)
    integer, allocatable :: temp(:,:),salt(:,:),cntt(:),cnts(:)
    integer(kind=1), allocatable :: temp_qc(:,:),salt_qc(:,:) 
    logical :: first
    logical, allocatable :: mask3d(:,:,:),hmask(:)

    100 format(a,'/ocean_avg_',i4,'-',i2.2,'-',i2.2,'.nc')
!    path = '/users/work/mmuzyka/CSDIR/metro_05NM_era5/run/baltic'
!    path = '/users/magazyn/mmuzyka/ROMS_CASES/metro_560x600_era5v2'
!    path = '/users/magazyn/mmuzyka/ROMS_CASES/metro_560x600_era5v2'
    path = '/users/work/mmuzyka/CSDIR/metro_560x600_uerraGLS4/run/baltic'
!    path = '/users/magazyn/mmuzyka/ROMS_CASES/metro_560x600_era5test4'
!    path = '/users/magazyn/mmuzyka/ROMS_CASES/metro_560x600_era5MYtest'
!    path = '/users/magazyn/mmuzyka/ROMS_CASES/metro_560x600_uerraMY'
!    path = '/users/work/mmuzyka/CSDIR/metro_560x600_era2004v14/run/baltic'
    cfile = 'NO_TS_MO_Arkona.nc'

    FV = -999.0
    
    date = (/2021,1,2/)
    datestp = (/2022,6,4/)
    date0 = (/1950,1,1/)

    CALL check(nf90_open(trim(cfile),NF90_NOWRITE,ncid),310)

    CALL check(nf90_inq_dimid(ncid, "TIME", dimid),311)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nt),312)

    CALL check(nf90_inq_dimid(ncid, "DEPTH", dimid),311)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nk),312)

    allocate( time(nt), depth(nk), cntt(nk), cnts(nk) )
    allocate( temp(nk,nt), temp_qc(nk,nt), avet(nk) )
    allocate( salt(nk,nt), salt_qc(nk,nt), aves(nk), hmask(nk) )

    CALL check(nf90_inq_varid(ncid,"TIME",vid),27)
    CALL check(nf90_get_var(ncid,vid,time),28)

    CALL check(nf90_inq_varid(ncid,"DEPH",vid),27)
    CALL check(nf90_get_var(ncid,vid,depth,start=(/1,1/),&
        count=(/nk,1/)),28)

    CALL check(nf90_inq_varid(ncid,"LONGITUDE",vid),27)
    CALL check(nf90_get_var(ncid,vid,ptlon,start=(/1/),&
        count=(/1/)),28)

    CALL check(nf90_inq_varid(ncid,"LATITUDE",vid),27)
    CALL check(nf90_get_var(ncid,vid,ptlat,start=(/1/),&
        count=(/1/)),28)

    CALL check(nf90_inq_varid(ncid,"TEMP",vid),27)
    CALL check(nf90_get_var(ncid,vid,temp),28)
    CALL check(nf90_get_att(ncid,vid,'_FillValue',FVt),28)
    CALL check(nf90_get_att(ncid,vid,'add_offset',AOt),28)
    CALL check(nf90_get_att(ncid,vid,'scale_factor',SFt),28)
    CALL check(nf90_inq_varid(ncid,"TEMP_QC",vid),27)
    CALL check(nf90_get_var(ncid,vid,temp_qc),28)

    CALL check(nf90_inq_varid(ncid,"PSAL",vid),27)
    CALL check(nf90_get_var(ncid,vid,salt),28)
    CALL check(nf90_get_att(ncid,vid,'_FillValue',FVs),28)
    CALL check(nf90_get_att(ncid,vid,'add_offset',AOs),28)
    CALL check(nf90_get_att(ncid,vid,'scale_factor',SFs),28)
    CALL check(nf90_inq_varid(ncid,"PSAL_QC",vid),27)
    CALL check(nf90_get_var(ncid,vid,salt_qc),28)


    write(*,*) ptlon,ptlat
    write(*,*) depth

    write(*,*) FVt,AOt,SFt,FVs,AOs,SFs

    CALL check(nf90_close(ncid),360)

   
    write(filename,100) trim(path),date(1),date(2),date(3)
    write(*,*) trim(filename) 
    CALL check(nf90_open(trim(filename),NF90_NOWRITE,ncid),310)
    
    CALL check(nf90_inq_dimid(ncid, "xi_rho", dimid),311)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=ni),312)

    CALL check(nf90_inq_dimid(ncid, "eta_rho", dimid),311)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nj),312)

    CALL check(nf90_inq_dimid(ncid, "s_rho", dimid),311)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=ns),312)

    allocate( Cs_r(ns),s_r(ns),lon(ni,nj),lat(ni,nj) )
    allocate( mt(2,2,ns),ms(2,2,ns),dp(2,2,ns) )

    CALL check(nf90_inq_varid(ncid,"lon_rho",vid),27)
    CALL check(nf90_get_var(ncid,vid,lon),28)
     
    CALL check(nf90_inq_varid(ncid,"lat_rho",vid),27)
    CALL check(nf90_get_var(ncid,vid,lat),28)
     
    CALL check(nf90_inq_varid(ncid,"Cs_r",vid),27)
    CALL check(nf90_get_var(ncid,vid,Cs_r),28)
     
    CALL check(nf90_inq_varid(ncid,"s_rho",vid),27)
    CALL check(nf90_get_var(ncid,vid,s_r),28)
     
    CALL check(nf90_inq_varid(ncid,"hc",vid),27)
    CALL check(nf90_get_var(ncid,vid,hc),28)
     

    CALL find_idx(ni,nj,lon,lat,real((/ptlon(1),ptlat(1)/),8),idx,alphas)
    write(*,*) (/ptlon(1),ptlat(1)/),idx
    write(*,*) alphas
    write(*,*) (/lon(idx(1),idx(2)),lat(idx(1),idx(2))/)
    write(*,*) (/lon(idx(1)+1,idx(2)),lat(idx(1)+1,idx(2))/)
    write(*,*) (/lon(idx(1)+1,idx(2)+1),lat(idx(1)+1,idx(2)+1)/)
    write(*,*) (/lon(idx(1),idx(2)+1),lat(idx(1),idx(2)+1)/)

    CALL check(nf90_inq_varid(ncid,"h",vid),27)
    CALL check(nf90_get_var(ncid,vid,h,start=(/idx(1),idx(2)/),&
        count=(/2,2/)),28)

!    write(*,*) h
    
    CALL check(nf90_close(ncid),360)

    zeta = 0.0
    CALL calc_dp(2,2,ns,zeta,h,hc,s_r,Cs_r,dp)


    CALL elapsed(date,date0,day) 
    CALL elapsed(datestp,date0,daystp)

    CALL check(nf90_create('CTD_Arkona_uerraGLS4.nc',NF90_NETCDF4,ncid2),310)
    CALL check( nf90_def_dim( ncid2, 'depth', nk, zdimid ),503)
    CALL check( nf90_def_dim( ncid2, 'time', NF90_UNLIMITED, tdimid ),503)

    CALL check( nf90_def_var( ncid2, 'depth', NF90_DOUBLE, (/zdimid/), varid(1) ), 504)
    CALL check( nf90_put_att( ncid2, varid(1), 'long_name', 'depth' ), 504)
    CALL check( nf90_put_att( ncid2, varid(1), 'units', 'm' ), 504)
    CALL check( nf90_put_att( ncid2, varid(1), 'axis', 'Z' ), 504)
    
    CALL check( nf90_def_var( ncid2, 'time', NF90_DOUBLE, (/tdimid/), varid(2) ), 504)
    CALL check( nf90_put_att( ncid2, varid(2), 'long_name', 'time' ), 504)
    CALL check( nf90_put_att( ncid2, varid(2), 'units',&
        'days since 1950-01-01T00:00:00Z' ), 504)
    CALL check( nf90_put_att( ncid2, varid(2), 'axis', 'T' ), 504)
    CALL check( nf90_put_att( ncid2, varid(2), 'calendar', 'standard' ), 504)
    
    CALL check( nf90_def_var( ncid2, 'temp_mo', NF90_DOUBLE, (/zdimid,tdimid/), varid(3) ), 504)
    CALL check( nf90_put_att( ncid2, varid(3), 'long_name','Sea temperature' ), 504)
    CALL check( nf90_put_att( ncid2, varid(3), 'units','degrees_C' ), 504)
    CALL check( nf90_put_att( ncid2, varid(3), '_FillValue',FV ), 504)

    CALL check( nf90_def_var( ncid2, 'salt_mo', NF90_DOUBLE, (/zdimid,tdimid/), varid(4) ), 504)
    CALL check( nf90_put_att( ncid2, varid(4), 'long_name','Practical salinity' ), 504)
    CALL check( nf90_put_att( ncid2, varid(4), 'units','0.001' ), 504)
    CALL check( nf90_put_att( ncid2, varid(4), '_FillValue',FV ), 504)

    CALL check( nf90_def_var( ncid2, 'temp', NF90_DOUBLE, (/zdimid,tdimid/), varid(5) ), 504)
    CALL check( nf90_put_att( ncid2, varid(5), 'long_name','Sea temperature' ), 504)
    CALL check( nf90_put_att( ncid2, varid(5), 'units','degrees_C' ), 504)
    CALL check( nf90_put_att( ncid2, varid(5), '_FillValue',FV ), 504)

    CALL check( nf90_def_var( ncid2, 'salt', NF90_DOUBLE, (/zdimid,tdimid/), varid(6) ), 504)
    CALL check( nf90_put_att( ncid2, varid(6), 'long_name','Practical salinity' ), 504)
    CALL check( nf90_put_att( ncid2, varid(6), 'units','0.001' ), 504)
    CALL check( nf90_put_att( ncid2, varid(6), '_FillValue',FV ), 504)

    CALL check( nf90_enddef(ncid2), 516 )
   
    CALL check( nf90_put_var( ncid2, varid(1), -1*depth ), 517 ) 

    date = date0
    cntt = 0
    cnts = 0
    cnt = 0
    n = 0
    nold = 0
    do t = 1,nt
      if (time(t).lt.day) cycle
      if (time(t).gt.daystp) exit 

      do 
        CALL elapsed(date,date0,daytmp)
        if (time(t)-real(daytmp).lt.1.0) then
          if (time(t)-real(daytmp).lt.0.25) then
            n = 1
          elseif (time(t)-real(daytmp).lt.0.5) then
            n = 2
          elseif (time(t)-real(daytmp).lt.0.75) then
            n = 3
          else
            n = 4
          end if
          if (n.ne.nold .and. nold.ne.0) then
            cnt = cnt+1
            CALL check( nf90_put_var( ncid2, varid(2), mytime, start = (/cnt/) ), 517 ) 
            do k = 1,nk
              if (cntt(k).ne.0) then
                avet(k) = (avet(k)/cntt(k))*SFt+AOt
!                write(*,*)  't', k, cntt(k), avet(k)
              else
                avet(k) = FV
              endif
              if (cnts(k).ne.0) then
                aves(k) = (aves(k)/cnts(k))*SFs+AOs
!                write(*,*)  's', k, cnts(k), aves(k)
              else
                aves(k) = FV
              endif
            end do
            !save real(aves(k))*SFs+AOs to file
            CALL check( nf90_put_var( ncid2, varid(3), avet, &
                start = (/1,cnt/), count=(/nk,1/) ), 517 ) 
            CALL check( nf90_put_var( ncid2, varid(4), aves, &
                start = (/1,cnt/), count=(/nk,1/) ), 517 ) 
            cntt=0
            cnts=0

!read dateold file record nold
            write(*,*) dateold,nold
            
            write(filename,100) trim(path),dateold(1),dateold(2),dateold(3)
            CALL check(nf90_open(trim(filename),NF90_NOWRITE,ncid),310)

            CALL check(nf90_inq_varid(ncid,"zeta",vid),27)
            CALL check(nf90_get_var(ncid,vid,zeta,start=(/idx(1),idx(2),nold/),&
                count=(/2,2,1/)),28)

            CALL check(nf90_inq_varid(ncid,"temp",vid),27)
            CALL check(nf90_get_var(ncid,vid,mt,&
                start=(/idx(1),idx(2),1,nold/),count=(/2,2,ns,1/)),28)

            CALL check(nf90_inq_varid(ncid,"salt",vid),27)
            CALL check(nf90_get_var(ncid,vid,ms,&
                start=(/idx(1),idx(2),1,nold/),count=(/2,2,ns,1/)),28)

            CALL check(nf90_close(ncid),360)

            CALL calc_dp(2,2,ns,zeta,h,hc,s_r,Cs_r,dp)

            !dp,mt/ms,hmask,alphas
            hmask = .false.
            where(avet.ne.FV) hmask = .true.
            do k = 3,nk
!                if (hmask(k)) write(*,*) 'TEMP dp=',depth(k)
!                if (hmask(k)) then
                    CALL calc_val(ns,mt,dp,alphas,depth(k),avet(k))
!                endif
            end do
            hmask = .false.
            where(aves.ne.FV) hmask = .true.
            do k = 3,nk
!                if (hmask(k)) write(*,*) 'SALT dp=',depth(k)
!                if (hmask(k)) then
                    CALL calc_val(ns,ms,dp,alphas,depth(k),aves(k))
!                endif
            end do
           

            CALL check( nf90_put_var( ncid2, varid(5), avet, &
                start = (/1,cnt/), count=(/nk,1/) ), 517 )
            CALL check( nf90_put_var( ncid2, varid(6), aves, &
                start = (/1,cnt/), count=(/nk,1/) ), 517 ) 
!end read    
            avet = 0
            aves = 0

          endif
          nold = n
          dateold = date
          mytime = real(daytmp)+(n-1)*0.25+0.125
          do k = 1,nk
            if (temp_qc(k,t).eq.1 .and. temp(k,t).ne.FVt) then
                cntt(k) = cntt(k)+1 
                avet(k) = avet(k)+temp(k,t) 
            endif
            if (salt_qc(k,t).eq.1 .and. salt(k,t).ne.FVs) then
                cnts(k) = cnts(k)+1 
                aves(k) = aves(k)+salt(k,t) 
            endif
          end do         
!          cnt = cnt+1
!          write(*,*) t,date,nint((time(t)-real(daytmp))*24),n
          exit
        endif
 
        CALL add_day(date,.true.)
      end do
    end do 
    CALL check(nf90_close(ncid2),360)

    deallocate( time, depth, cntt, cnts, avet, aves )
    deallocate( temp, temp_qc, salt, salt_qc )
    deallocate( Cs_r,s_r,lon,lat,mt,ms,dp,hmask )

END PROGRAM MO

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

SUBROUTINE calc_dp(nx,ny,nz,zeta,h,hc,s_r,Cs_r,dp)
    implicit none
    integer, intent(in) :: nx,ny,nz
    real(kind=8), intent(in) :: zeta(nx,ny),h(nx,ny),hc,s_r(nz),Cs_r(nz)
    real(kind=8), intent(out) :: dp(nx,ny,nz)

    integer :: i,j,k

    dp = 0.0

    do i = 1,nx
    do j = 1,ny
      do k = 1,nz
        dp(i,j,k) = -1*(zeta(i,j)+(zeta(i,j)+h(i,j))*(hc*s_r(k)+h(i,j)*Cs_r(k))/&
            (hc+h(i,j)))
      end do
    end do
    end do

!    do k = 1,nz
!        write(*,*) dp(1,1,k)
!    end do

END SUBROUTINE calc_dp

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
!                write(*,*) 'ALPHAS',alphas
                fexit = .true.
                exit
            endif
        end do
    end do

END SUBROUTINE find_idx

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

SUBROUTINE calc_val(nk,m,dp,alphas,h,val)
    implicit none
    integer, intent(in) :: nk
    real(kind=8), intent(in) :: m(2,2,nk),dp(2,2,nk),alphas(4),h
    real(kind=8), intent(out) :: val

    real(kind=8) :: val2d(2,2),x0,x1,y0,y1
    integer :: i,j,k
    val = 0
!    write(*,*) "DEPTH",h,dp(1,1,1),dp(1,1,nk)

    do i = 1,2
    do j = 1,2
      do k = 1,nk-1
        if (h.le.dp(i,j,k).and.h.gt.dp(i,j,k+1)) then
!        write(*,*) "DEPTH",h,dp(1,1,k),dp(1,1,k+1)
          x0 = dp(i,j,k)
          x1 = dp(i,j,k+1)
          y0 = m(i,j,k)
          y1 = m(i,j,k+1)
          val2d(i,j) = y0*(x1-h)/(x1-x0)+y1*(h-x0)/(x1-x0)
          exit
        endif
      end do
    end do
    end do

    val = alphas(1)*val2d(1,1)+alphas(2)*val2d(2,1)+&
        alphas(3)*val2d(2,2)+alphas(4)*val2d(1,2)
END SUBROUTINE calc_val
