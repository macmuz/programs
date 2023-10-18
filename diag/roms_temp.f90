PROGRAM roms_temp
    USE netcdf
    implicit none
    integer :: date(3),datestp(3),cnt,i,j,t,k,l,x,y
    integer :: ncid,dimid,ncid2,tdimid,tmp(2)
    integer :: ni,nj,nk,nt,vid,varid(3)
    character(len=200) :: filename,path
    character(len=30) :: vname(2),varn
    real(kind=8) :: hc,pt(2),npt,W(4),corners(4,2)
    real(kind=8), allocatable :: time(:),Cs_w(:),s_w(:),lon(:,:),lat(:,:)
    real(kind=8), allocatable :: h(:,:),zeta(:,:,:),dataout(:),dist(:,:),dp(:,:)
    real(kind=8), allocatable :: datatmp(:,:,:,:)
    logical :: first,score

    100 format(a,'/ocean_avg_',i4,'-',i2.2,'-',i2.2,'.nc')
!    path = '/users/work/mmuzyka/CSDIR/metro_05NM_era5/run/baltic'
!    path = '/users/magazyn/mmuzyka/ROMS_CASES/metro_560x600_era5v2'
!    path = '/users/magazyn/mmuzyka/ROMS_CASES/metro_560x600_era5v2'
!    path = '/users/work/mmuzyka/CSDIR/metro_560x600_era5test4/run/baltic'
!    path = '/users/magazyn/mmuzyka/ROMS_CASES/metro_560x600_era5test4'
!    path = '/users/magazyn/mmuzyka/ROMS_CASES/metro_560x600_era5MY'
    path = '/users/work/mmuzyka/CSDIR/metro_560x600_uerraMY/run/baltic'
    
    date = (/1993,1,2/)
    datestp = (/1994,12,31/)

    pt = (/15.8,55.32/)

    vname(1) = 'salt'
    vname(2) = 'temp'

    first = .true.
    cnt = 0
    CALL check(nf90_create('roms_1pt_stats_uerraMY.nc',NF90_NETCDF4,ncid2),310)

    do
        write(filename,100) trim(path),date(1),date(2),date(3)
        write(*,*) trim(filename)        

        CALL check(nf90_open(trim(filename),NF90_NOWRITE,ncid),310)
        if (first) then
            CALL check(nf90_inq_dimid(ncid, "xi_rho", dimid),311)
            CALL check(nf90_inquire_dimension(ncid, dimid, len=ni),312)

            CALL check(nf90_inq_dimid(ncid, "eta_rho", dimid),313)
            CALL check(nf90_inquire_dimension(ncid, dimid, len=nj),314)
        
            CALL check(nf90_inq_dimid(ncid, "s_rho", dimid),313)
            CALL check(nf90_inquire_dimension(ncid, dimid, len=nk),314)

            CALL check(nf90_inq_dimid(ncid, "ocean_time", dimid),315)
            CALL check(nf90_inquire_dimension(ncid, dimid, len=nt),316)

            allocate( time(nt), datatmp(2,2,nk,nt), dataout(nt) )
            allocate( Cs_w(nk+1), s_w(nk+1), dp(nt,nk), lon(ni,nj), lat(ni,nj) )
            allocate( h(2,2), zeta(2,2,nt), dist(ni,nj) )


            CALL check(nf90_inq_varid(ncid,"lon_rho",vid),29)
            CALL check(nf90_get_var(ncid,vid,lon),30)
            CALL check(nf90_inq_varid(ncid,"lat_rho",vid),31)
            CALL check(nf90_get_var(ncid,vid,lat),32)
            CALL check(nf90_inq_varid(ncid,"hc",vid),35)
            CALL check(nf90_get_var(ncid,vid,hc),36)
            CALL check(nf90_inq_varid(ncid,"Cs_w",vid),41)
            CALL check(nf90_get_var(ncid,vid,Cs_w),42)
            CALL check(nf90_inq_varid(ncid,"s_w",vid),43)
            CALL check(nf90_get_var(ncid,vid,s_w),44)

            dist = sqrt((lon-pt(1))**2+(lat-pt(2))**2)
            tmp = minloc(dist)
            x = tmp(1)
            y = tmp(2)
            do k = x-1,x
            do l = y-1,y
              corners(1,:) = (/lon(k,l),lat(k,l)/)
              corners(2,:) = (/lon(k,l+1),lat(k,l+1)/)
              corners(3,:) = (/lon(k+1,l+1),lat(k+1,l+1)/)
              corners(4,:) = (/lon(k+1,l),lat(k+1,l)/)
              CALL in_convex_polygon(corners,pt,score)
              if (score) then
                CALL coef2(corners,pt,W)
                write(*,*) k,l
!                write(*,*) pt
!                write(*,*) corners(1,:),W(1) 
!                write(*,*) corners(2,:),W(2) 
!                write(*,*) corners(3,:),W(3) 
!                write(*,*) corners(4,:),W(4) 
                goto 1001
              endif
            end do
            end do
1001        continue

            CALL check(nf90_inq_varid(ncid,"h",vid),33)
            CALL check(nf90_get_var(ncid,vid,h,start=(/k,l/),&
                count=(/2,2/)),34)

            CALL check(nf90_inq_varid(ncid,'ocean_time',vid),120)
            CALL check(nf90_get_var(ncid,vid,time),18)

            CALL check(nf90_inq_varid(ncid,'ocean_time',vid),120)

            CALL check( nf90_def_dim( ncid2, 'time', NF90_UNLIMITED, tdimid ),503)
            CALL check( nf90_def_var( ncid2, 'time', NF90_DOUBLE, (/tdimid/), varid(1) ), 504)
            CALL check( nf90_copy_att(ncid, vid, 'long_name', ncid2, varid(1)), 504)
            CALL check( nf90_copy_att(ncid, vid, 'units', ncid2, varid(1)), 504)
            CALL check( nf90_copy_att(ncid, vid, 'calendar', ncid2, varid(1)), 504)
            CALL check( nf90_copy_att(ncid, vid, 'field', ncid2, varid(1)), 504)

            do i = 1, 2
            CALL check( nf90_def_var( ncid2, trim(vname(i)), NF90_DOUBLE,&
              (/tdimid/),varid(1+i) ), 504)
            end do

            CALL check( nf90_enddef(ncid2), 516 )
            first = .false.
        end if

        CALL check(nf90_inq_varid(ncid,'ocean_time',vid),120)
        CALL check(nf90_get_var(ncid,vid,time),18)
        CALL check( nf90_put_var( ncid2, varid(1), time, start=(/cnt+1/) ), 517 )

        CALL check(nf90_inq_varid(ncid,'zeta',vid),120)
        CALL check(nf90_get_var(ncid,vid,zeta,start=(/k,l,1/),&
            count=(/2,2,nt/)),18)

        CALL calc_dp(nk,nt,zeta,h,hc,s_w,Cs_w,dp)

!        write(*,*) 'after dp'
        do i = 1,2
          CALL check(nf90_inq_varid(ncid,vname(i),vid),120)
          CALL check(nf90_get_var(ncid,vid,datatmp,start=(/k,l,1,1/),&
            count=(/2,2,nk,nt/)),18)
          CALL ave(nk,nt,datatmp,W,dp,dataout)
!          write(*,*) 'after ave'
          CALL check( nf90_put_var( ncid2, varid(1+i),&
            dataout, start=(/cnt+1/) ), 518 )
        end do
        
        CALL check(nf90_close(ncid),360)

        cnt = cnt+nt

        if (date(1).ge.datestp(1).and.&
            date(2).ge.datestp(2).and.&
            date(3).ge.datestp(3)) EXIT
        CALL add_day(date,.true.)
    end do

    CALL check(nf90_close(ncid2),360)
   
    deallocate(time,datatmp,Cs_w,s_w,dp,lon,lat,dataout)
    deallocate(h,zeta,dist) 
    
END PROGRAM roms_temp

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

SUBROUTINE calc_dp(nz,nt,zeta,h,hc,s_w,Cs_w,dp)
    implicit none
    integer, intent(in) :: nz,nt
    real(kind=8), intent(in) :: zeta(2,2,nt),h(2,2),hc,s_w(nz+1),Cs_w(nz+1)
    real(kind=8), intent(out) :: dp(nt,nz)

    integer :: i,j,t,z
    real(kind=8) :: tmp(2,2),lvl(nz+1)

    dp = 0.0

    do t = 1,nt
    do z = 1,nz+1
      do i = 1,2
      do j = 1,2
        tmp(i,j) = (zeta(i,j,t)+(zeta(i,j,t)+h(i,j))*(hc*s_w(z)+h(i,j)*Cs_w(z))/&
           (hc+h(i,j)))
      end do
      end do
      lvl(z) = 0.25*sum(tmp)
    end do
    do z = 1,nz
      dp(t,z) = lvl(z+1)-lvl(z)
!      write(*,*) dp(t,z)
    end do
    end do
END SUBROUTINE calc_dp

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

SUBROUTINE coef(corners,point,W)
    implicit none
    real(kind=8), intent(in) :: corners(4,2)
    real(kind=8), intent(in) :: point(2)
    real(kind=8), intent(out) :: W(4)

    real(kind=8) :: a,b,c,d,e,f,g,h,alpha,beta

    a = -corners(1,1)+corners(4,1)
    b = -corners(1,1)+corners(2,1)
    c = corners(1,1)-corners(2,1)-corners(4,1)+corners(3,1)
    b = point(1)-corners(1,1)
    e = -corners(1,2)+corners(4,2)
    f = -corners(1,2)+corners(2,2)
    g = corners(1,2)-corners(2,2)-corners(4,2)+corners(3,2)
    h = point(2)-corners(1,2)

    alpha = (-1)*(b*e-a*f+d*g-c*h+sqrt(-4*(c*e-a*g)*(d*f-b*h)+&
            (b*e-a*f+d*g-c*h)**2))/(2*c*e-2*a*g)
    beta = (b*e-a*f-d*g+c*h+sqrt(-4*(c*e-a*g)*(d*f-b*h)+&
            (b*e-a*f+d*g-c*h)**2))/(2*c*f-2*b*g)

!    write(*,*) alpha,beta

    if (alpha.lt.0 .or. alpha.gt.1 .or. &
        beta.lt.0 .or. beta.gt.1 ) then
      alpha = ((-1)*b*e+a*f-d*g+c*h+sqrt(-4*(c*e-a*g)*(d*f-b*h)+&
              (b*e-a*f+d*g-c*h)**2))/(2*c*e-2*a*g)
      beta = (-1)*((-b*e+a*f+d*g-c*h+sqrt(-4*(c*e-a*g)*(d*f-b*h)+&
              (b*e-a*f+d*g-c*h)**2))/(2*c*f-2*b*g))
!      write(*,*) alpha,beta
    end if

    W(1) = (1-alpha)*(1-beta)
    W(2) = (1-alpha)*beta
    W(3) = alpha*(1-beta)
    W(4) = alpha*beta
END SUBROUTINE coef

SUBROUTINE coef2(corners,point,W)
    implicit none
    real(kind=8), intent(in) :: corners(4,2)
    real(kind=8), intent(in) :: point(2)
    real(kind=8), intent(out) :: W(4)

    real(kind=8) :: A,B,C,t,s
    real(kind=8), dimension(2) :: p1,p2,p3,p4

    p1 = (/corners(2,1),corners(2,2)/)
    p2 = (/corners(3,1),corners(3,2)/)
    p3 = (/corners(1,1),corners(1,2)/)
    p4 = (/corners(4,1),corners(4,2)/)

    A = (p3(1)-p1(1))*(p4(2)-p2(2))-(p3(2)-p1(2))*(p4(1)-p2(1))
    B = point(2)*((p4(1)-p2(1))-(p3(1)-p1(1)))-&
        point(1)*((p4(2)-p2(2))-(p3(2)-p1(2)))+&
        (p3(1)-p1(1))*p2(2)-(p3(2)-p1(2))*p2(1)+&
        (p4(2)-p2(2))*p1(1)-(p4(1)-p2(1))*p1(2)
    C = point(2)*(p2(1)-p1(1))-point(1)*(p2(2)-p1(2))+&
        p1(1)*p2(2)-p1(2)*p2(1)

    t = (-B+sqrt(B**2-4*A*C))/(2*A)
    if ((t<0).OR.(t>1)) t=(-B-sqrt(B**2-4*A*C))/(2*A)
    s = (point(2)-p1(2)-(p3(2)-p1(2))*t)/(p2(2)+(p4(2)-p2(2))*t-&
        p1(2)-(p3(2)-p1(2))*t)

    W(1) = (1-s)*t
    W(2) = (1-s)*(1-t)
    W(3) = s*(1-t)
    W(4) = s*t
END SUBROUTINE coef2

SUBROUTINE ave(nk,nt,tmp,W,dp,output)
    implicit none
    integer, intent(in) :: nk,nt
    real(kind=8), intent(in) :: tmp(2,2,nk,nt), W(4), dp(nt,nk)
    real(kind=8), intent(out) :: output(nt)

    integer :: z,t
    
    output = 0
    do t = 1,nt
    do z = 1,nk
      output(t) = output(t)+dp(t,z)*(tmp(1,1,z,t)*W(1)+tmp(1,2,z,t)*W(2)+&
            &   tmp(2,2,z,t)*W(3)+tmp(2,1,z,t)*W(4))
    end do
    output(t) = output(t)/real(sum(dp(t,:)))
!    write(*,*) output(t)
    end do
END SUBROUTINE ave
