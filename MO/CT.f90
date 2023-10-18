PROGRAM CT
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
    real(kind=8), allocatable :: mt(:,:,:),ms(:,:,:),dp(:,:,:),dp1d(:)
    real(kind=8), allocatable :: avet(:),aves(:),depth(:,:)
    real(kind=8), allocatable :: ptlon(:),ptlat(:),temp(:,:),salt(:,:)
    integer, allocatable :: cntt(:),cnts(:),tmp2d(:,:)
    integer(kind=1), allocatable :: temp_qc(:,:),salt_qc(:,:) 
    logical :: first
    logical, allocatable :: mask3d(:,:,:),hmask(:)

    100 format(a,'/ocean_avg_',i4,'-',i2.2,'-',i2.2,'.nc')
!    path = '/users/work/mmuzyka/CSDIR/metro_05NM_era5/run/baltic'
!    path = '/users/magazyn/mmuzyka/ROMS_CASES/metro_560x600_era5v2'
!    path = '/users/magazyn/mmuzyka/ROMS_CASES/metro_560x600_era5v2'
    path = '/users/work/mmuzyka/CSDIR/metro_560x600_uerraGLS4/run/baltic'
!    path = '/users/magazyn/mmuzyka/ROMS_CASES/metro_560x600_uerraMY'
!    path = '/users/magazyn/mmuzyka/ROMS_CASES/metro_560x600_era5MYtest'
!    path = '/users/work/mmuzyka/CSDIR/metro_560x600_era5MYtest/run/baltic'
    cfile = 'BO_PR_CT_SMHIBY5.nc'

    FV = -999.0
    
    date = (/2021,1,2/)
    datestp = (/2022,6,4/)
    date0 = (/1950,1,1/)

    CALL check(nf90_open(trim(cfile),NF90_NOWRITE,ncid),310)

    CALL check(nf90_inq_dimid(ncid, "TIME", dimid),311)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nt),312)

    CALL check(nf90_inq_dimid(ncid, "DEPTH", dimid),311)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nk),312)

    allocate( time(nt), depth(nk,nt), cntt(nk), cnts(nk) )
    allocate( temp(nk,nt), temp_qc(nk,nt), tmp2d(nk,nt) )
    allocate( salt(nk,nt), salt_qc(nk,nt), hmask(nk) )
    allocate( ptlon(nt), ptlat(nt) )

    CALL check(nf90_inq_varid(ncid,"TIME",vid),27)
    CALL check(nf90_get_var(ncid,vid,time),28)

    CALL check(nf90_inq_varid(ncid,"DEPH",vid),27)
    CALL check(nf90_get_var(ncid,vid,depth),28)

    CALL check(nf90_inq_varid(ncid,"LONGITUDE",vid),27)
    CALL check(nf90_get_var(ncid,vid,ptlon),28)

    CALL check(nf90_inq_varid(ncid,"LATITUDE",vid),27)
    CALL check(nf90_get_var(ncid,vid,ptlat),28)

    CALL check(nf90_inq_varid(ncid,"TEMP",vid),27)
    CALL check(nf90_get_var(ncid,vid,tmp2d),28)
    CALL check(nf90_get_att(ncid,vid,'_FillValue',FVt),28)
    CALL check(nf90_get_att(ncid,vid,'add_offset',AOt),28)
    CALL check(nf90_get_att(ncid,vid,'scale_factor',SFt),28)
    CALL check(nf90_inq_varid(ncid,"TEMP_QC",vid),27)
    CALL check(nf90_get_var(ncid,vid,temp_qc),28)
    temp = real(tmp2d,8)*SFt+AOt
    write(*,*) 'SF AO',SFt,AOt,temp(24,1),tmp2d(24,1)


    CALL check(nf90_inq_varid(ncid,"PSAL",vid),27)
    CALL check(nf90_get_var(ncid,vid,tmp2d),28)
    CALL check(nf90_get_att(ncid,vid,'_FillValue',FVs),28)
    CALL check(nf90_get_att(ncid,vid,'add_offset',AOs),28)
    CALL check(nf90_get_att(ncid,vid,'scale_factor',SFs),28)
    CALL check(nf90_inq_varid(ncid,"PSAL_QC",vid),27)
    CALL check(nf90_get_var(ncid,vid,salt_qc),28)
    salt = real(tmp2d,8)*SFs+AOs


!    write(*,*) ptlon,ptlat
!    write(*,*) depth

!    write(*,*) FVt,AOt,SFt,FVs,AOs,SFs

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
    allocate( mt(2,2,ns),ms(2,2,ns),dp(2,2,ns),dp1d(ns) )
    allocate( avet(ns),aves(ns) )

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
     

!    CALL find_idx(ni,nj,lon,lat,real((/ptlon(1),ptlat(1)/),8),idx,alphas)
!    write(*,*) (/ptlon(1),ptlat(1)/),idx
!    write(*,*) alphas
!    write(*,*) (/lon(idx(1),idx(2)),lat(idx(1),idx(2))/)
!    write(*,*) (/lon(idx(1)+1,idx(2)),lat(idx(1)+1,idx(2))/)
!    write(*,*) (/lon(idx(1)+1,idx(2)+1),lat(idx(1)+1,idx(2)+1)/)
!    write(*,*) (/lon(idx(1),idx(2)+1),lat(idx(1),idx(2)+1)/)

!    CALL check(nf90_inq_varid(ncid,"h",vid),27)
!    CALL check(nf90_get_var(ncid,vid,h,start=(/idx(1),idx(2)/),&
!        count=(/2,2/)),28)

!    write(*,*) h
    
    CALL check(nf90_close(ncid),360)

!    zeta = 0.0
!    CALL calc_dp(2,2,ns,zeta,h,hc,s_r,Cs_r,dp)


    CALL elapsed(date,date0,day) 
    CALL elapsed(datestp,date0,daystp)

    CALL check(nf90_create('CT_PR_uerraGLS4_BY5.nc',NF90_NETCDF4,ncid2),310)
    CALL check( nf90_def_dim( ncid2, 'depth', ns, zdimid ),503)
    CALL check( nf90_def_dim( ncid2, 'time', NF90_UNLIMITED, tdimid ),503)

    CALL check( nf90_def_var( ncid2, 'depth', NF90_DOUBLE, (/zdimid,tdimid/), varid(1) ), 504)
    CALL check( nf90_put_att( ncid2, varid(1), 'long_name', 'depth' ), 504)
    CALL check( nf90_put_att( ncid2, varid(1), 'units', 'm' ), 504)
    CALL check( nf90_put_att( ncid2, varid(1), 'axis', 'Z' ), 504)
    
    CALL check( nf90_def_var( ncid2, 'time', NF90_DOUBLE, (/tdimid/), varid(2) ), 504)
    CALL check( nf90_put_att( ncid2, varid(2), 'long_name', 'time' ), 504)
    CALL check( nf90_put_att( ncid2, varid(2), 'units',&
        'days since 1950-01-01T00:00:00Z' ), 504)
    CALL check( nf90_put_att( ncid2, varid(2), 'axis', 'T' ), 504)
    CALL check( nf90_put_att( ncid2, varid(2), 'calendar', 'standard' ), 504)
    
    CALL check( nf90_def_var( ncid2, 'temp', NF90_DOUBLE, (/zdimid,tdimid/), varid(3) ), 504)
    CALL check( nf90_put_att( ncid2, varid(3), 'long_name','Sea temperature' ), 504)
    CALL check( nf90_put_att( ncid2, varid(3), 'units','degrees_C' ), 504)
    CALL check( nf90_put_att( ncid2, varid(3), '_FillValue',FV ), 504)

    CALL check( nf90_def_var( ncid2, 'salt', NF90_DOUBLE, (/zdimid,tdimid/), varid(4) ), 504)
    CALL check( nf90_put_att( ncid2, varid(4), 'long_name','Practical salinity' ), 504)
    CALL check( nf90_put_att( ncid2, varid(4), 'units','0.001' ), 504)
    CALL check( nf90_put_att( ncid2, varid(4), '_FillValue',FV ), 504)

    CALL check( nf90_def_var( ncid2, 'ctd_temp', NF90_DOUBLE, (/zdimid,tdimid/), varid(5) ), 504)
    CALL check( nf90_put_att( ncid2, varid(5), 'long_name','Sea temperature' ), 504)
    CALL check( nf90_put_att( ncid2, varid(5), 'units','degrees_C' ), 504)
    CALL check( nf90_put_att( ncid2, varid(5), '_FillValue',FV ), 504)

    CALL check( nf90_def_var( ncid2, 'ctd_salt', NF90_DOUBLE, (/zdimid,tdimid/), varid(6) ), 504)
    CALL check( nf90_put_att( ncid2, varid(6), 'long_name','Practical salinity' ), 504)
    CALL check( nf90_put_att( ncid2, varid(6), 'units','0.001' ), 504)
    CALL check( nf90_put_att( ncid2, varid(6), '_FillValue',FV ), 504)

    CALL check( nf90_enddef(ncid2), 516 )
   
!   CALL check( nf90_put_var( ncid2, varid(1), -1*depth ), 517 ) 

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
!          if (n.ne.nold .and. nold.ne.0) then
          cnt = cnt+1
        
          CALL check( nf90_put_var( ncid2, varid(2), &
            real(daytmp)+((n-1)*6.0+3.0)/24.0, start = (/cnt/) ), 517 )

!            do k = 1,nk
!              if (cntt(k).ne.0) then
!                avet(k) = (avet(k)/cntt(k))*SFt+AOt
!              else
!                avet(k) = FV
!              endif
!              if (cnts(k).ne.0) then
!                aves(k) = (aves(k)/cnts(k))*SFs+AOs
!              else
!                aves(k) = FV
!              endif
!            end do
            !save real(aves(k))*SFs+AOs to file
!            CALL check( nf90_put_var( ncid2, varid(3), avet, &
!                start = (/1,cnt/), count=(/nk,1/) ), 517 ) 
!            CALL check( nf90_put_var( ncid2, varid(4), aves, &
!                start = (/1,cnt/), count=(/nk,1/) ), 517 ) 
!            cntt=0
!            cnts=0

!read dateold file record nold
 !           write(*,*) dateold,nold
           
            CALL find_idx(ni,nj,lon,lat,real((/ptlon(t),ptlat(t)/),8),idx,alphas)
            write(*,*) idx,ptlon(t),ptlat(t)
 
            write(filename,100) trim(path),date(1),date(2),date(3)
            write(*,*) t,date(1),date(2),date(3),n
            CALL check(nf90_open(trim(filename),NF90_NOWRITE,ncid),310)

            CALL check(nf90_inq_varid(ncid,"zeta",vid),27)
            CALL check(nf90_get_var(ncid,vid,zeta,start=(/idx(1),idx(2),n/),&
                count=(/2,2,1/)),28)

            CALL check(nf90_inq_varid(ncid,"temp",vid),27)
            CALL check(nf90_get_var(ncid,vid,mt,&
                start=(/idx(1),idx(2),1,n/),count=(/2,2,ns,1/)),28)

            CALL check(nf90_inq_varid(ncid,"salt",vid),27)
            CALL check(nf90_get_var(ncid,vid,ms,&
                start=(/idx(1),idx(2),1,n/),count=(/2,2,ns,1/)),28)

            CALL check(nf90_inq_varid(ncid,"h",vid),27)
            CALL check(nf90_get_var(ncid,vid,h,start=(/idx(1),idx(2)/),&
                count=(/2,2/)),28)

            CALL check(nf90_close(ncid),360)

            CALL calc_dp(2,2,ns,zeta,h,hc,s_r,Cs_r,dp)

!            do k = 1,ns
!              if (dp(1,1,k).gt.1000.0) EXIT
!              write(*,*) k, dp(1,1,k)
!            end do

            hmask = .false.
            where(depth(:,t).lt.1000.0) hmask = .true.

            write(*,*) minval(dp(1,1,:)),minval(depth(:,t),hmask)
            write(*,*) maxval(dp(1,1,:)),maxval(depth(:,t),hmask)

            avet = FV
            aves = FV

!            do k = 1,nk
!              if (depth(k,t).gt.1000.0) EXIT
              CALL calc_val(ns,mt,dp,alphas,avet,dp1d)
              CALL calc_val(ns,ms,dp,alphas,aves,dp1d)
!            enddo

            !dp,mt/ms,hmask,alphas
!            hmask = .false.
!            where(avet.ne.FV) hmask = .true.
!            do k = 3,nk
!                if (hmask(k)) write(*,*) 'TEMP dp=',depth(k)
!                if (hmask(k)) then
!                    CALL calc_val(ns,mt,dp,alphas,depth(k),avet(k))
!                endif
!            end do
!            hmask = .false.
!            where(aves.ne.FV) hmask = .true.
!            do k = 3,nk
!                if (hmask(k)) write(*,*) 'SALT dp=',depth(k)
!                if (hmask(k)) then
!                    CALL calc_val(ns,ms,dp,alphas,depth(k),aves(k))
!                endif
!            end do
           
          CALL check( nf90_put_var( ncid2, varid(1), dp1d, &
            start = (/1,cnt/), count = (/ns,1/) ), 518 )


            CALL check( nf90_put_var( ncid2, varid(3), avet, &
                start = (/1,cnt/), count=(/ns,1/) ), 519 )
            CALL check( nf90_put_var( ncid2, varid(4), aves, &
                start = (/1,cnt/), count=(/ns,1/) ), 520 )

            write(*,*) 'T+++',t
            CALL move(ns,nk,dp1d,depth(:,t),temp(:,t),avet) 
            CALL move(ns,nk,dp1d,depth(:,t),salt(:,t),aves) 

            CALL check( nf90_put_var( ncid2, varid(5), avet, &
                start = (/1,cnt/), count=(/ns,1/) ), 519 )
            CALL check( nf90_put_var( ncid2, varid(6), aves, &
                start = (/1,cnt/), count=(/ns,1/) ), 520 )

!end read    
!            avet = 0
!            aves = 0

!          endif
!          nold = n
!          dateold = date
!          mytime = real(daytmp)+(n-1)*0.25+0.125
!          do k = 1,nk
!            if (temp_qc(k,t).eq.1 .and. temp(k,t).ne.FVt) then
!                cntt(k) = cntt(k)+1 
!                avet(k) = avet(k)+temp(k,t) 
!            endif
!            if (salt_qc(k,t).eq.1 .and. salt(k,t).ne.FVs) then
!                cnts(k) = cnts(k)+1 
!                aves(k) = aves(k)+salt(k,t) 
!            endif
!          end do         
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
    deallocate( Cs_r,s_r,lon,lat,mt,ms,dp,hmask,dp1d )
    deallocate( ptlon, ptlat, tmp2d )

END PROGRAM CT

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

SUBROUTINE calc_val(nk,m,dp,alphas,val,dp1d)
    implicit none
    integer, intent(in) :: nk
    real(kind=8), intent(in) :: m(2,2,nk),dp(2,2,nk),alphas(4)
    real(kind=8), intent(out) :: val(nk),dp1d(nk)

    integer :: i,j,k
    val = 0
!    write(*,*) "DEPTH",h,dp(1,1,1),dp(1,1,nk)

    do k = 1,nk
      val(k) = alphas(1)*m(1,1,k)+alphas(2)*m(2,1,k)+&
        alphas(3)*m(2,2,k)+alphas(4)*m(1,2,k)
      dp1d(k) = alphas(1)*dp(1,1,k)+alphas(2)*dp(2,1,k)+&
        alphas(3)*dp(2,2,k)+alphas(4)*dp(1,2,k)
    end do

END SUBROUTINE calc_val

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

SUBROUTINE move(ns,nk,dp1d,depth,input,output)
    implicit none
    integer, intent(in) :: ns,nk
    real(kind=8), intent(in) :: dp1d(ns),depth(nk),input(nk)
    real(kind=8), intent(out) :: output(ns)

    integer :: k,kmax
    real(kind=8) :: splint,y2a(nk)

    output = 0

    do k=1,nk
!        write(*,*) k,input(k)
        if (input(k).lt.-1000.0) EXIT
    end do
    kmax=k-1

    CALL spline(depth(1:kmax),input(1:kmax),kmax,y2a) 
!    write(*,*) 'inp1',depth(1)
!    write(*,*) 'inpend',depth(kmax)
!    write(*,*) 'out1',dp1d(1)
!    write(*,*) 'outend',dp1d(ns)

    do k = 1,ns
      if (dp1d(k).lt.depth(1)) then
        output(k) = input(1)
      elseif (dp1d(k).gt.depth(kmax)) then
        output(k) = input(kmax)
      else
        output(k) = splint(depth(1:kmax),input(1:kmax),y2a,kmax,dp1d(k))
      end if
      
    end do 

END SUBROUTINE move
