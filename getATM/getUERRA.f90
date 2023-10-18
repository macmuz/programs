PROGRAM getUERRA
    USE netcdf
    implicit none
    integer, parameter :: nx=150,ny=150
    integer :: cnt,idx(2),x,y
    integer :: ncid,dimid,ncid2,tdimid
    integer :: ni,nj,nt,vid,varid(7),i,n,p,t
    character(len=200) :: filename,path
    character(len=30) :: coord
    real(kind=8) :: W(2),pt(2),ao,sf
    integer, allocatable :: time(:)
    real(kind=8)  :: lon(nx,ny),lat(nx,ny)
    real(kind=8), allocatable :: data3d(:,:,:),uvel(:),vvel(:)
    real(kind=8), allocatable :: dew(:),tair(:),pair(:),qair(:)
    real(kind=8), allocatable :: tmp3d(:,:,:),RH(:)
    logical :: first,score,fexit

    100 format(a,'/UERRA_',a,'_',i4.4,'.nc')
    101 format(f7.4,'N, ',f8.4,'E')
    102 format('UERRA_mescan_qair_',i4.4,'.nc')
    103 format('UERRA_mescan_tair_',i4.4,'.nc')
    path = '/scratch/lustre/plgmacmuz/ROMS/UERRA'
    first = .true.
    cnt = 0
    
    pt = (/16.5,55.05/)
    write(coord,101) pt(2),pt(1)

!READ GRID
    write(filename,100) trim(path),'temp',1990
    CALL check(nf90_open(trim(filename),NF90_NOWRITE,ncid),310)
    CALL check(nf90_inq_varid(ncid,'longitude',vid),119)
    CALL check(nf90_get_var(ncid,vid,lon,&
      start=(/270,320/),count=(/nx,ny/) ),18)
    CALL check(nf90_inq_varid(ncid,'latitude',vid),119)
    CALL check(nf90_get_var(ncid,vid,lat,&
      start=(/270,320/),count=(/nx,ny/) ),18)
    CALL check(nf90_close(ncid),460) 
!END READ GRID

    CALL prepare(nx,ny,lon,lat,pt(1),pt(2),idx,W)
    x = idx(1)
    y = idx(2)
    write(*,*) pt
    write(*,*) x,y
    write(*,*) lon(x,y),lat(x,y),(1-W(1))*(1-W(2))
    write(*,*) lon(x+1,y),lat(x+1,y),W(1)*(1-W(2))
    write(*,*) lon(x,y+1),lat(x,y+1),(1-W(1))*W(2)
    write(*,*) lon(x+1,y+1),lat(x+1,y+1),W(1)*W(2)


    CALL check(nf90_create('atm_UERRA.nc',NF90_NETCDF4,ncid2),310)

    do i = 1990,2019
        write(filename,100) trim(path),'temp',i
!        write(filename,103) i
        write(*,*) trim(filename)
        CALL check(nf90_open(trim(filename),NF90_NOWRITE,ncid),310)

        CALL check(nf90_inq_dimid(ncid, "time", dimid),315)
        CALL check(nf90_inquire_dimension(ncid, dimid, len=nt),316)

        if (first) then

            CALL check(nf90_inq_varid(ncid,'time',vid),122)

            CALL check( nf90_def_dim( ncid2, 'time', NF90_UNLIMITED, tdimid ),503)
            CALL check( nf90_def_var( ncid2, 'time', NF90_INT, (/tdimid/), varid(1) ), 704)
            CALL check( nf90_copy_att(ncid, vid, 'units', ncid2, varid(1)), 705)
            CALL check( nf90_copy_att(ncid, vid, 'calendar', ncid2, varid(1)), 706)

            CALL check( nf90_def_var( ncid2, 'u10', NF90_FLOAT, (/tdimid/),varid(2) ), 707)
            CALL check( nf90_put_att( ncid2, varid(2), 'long_name', &
                '10 metre U wind component'), 708)
            CALL check( nf90_put_att( ncid2, varid(2), 'units', &
                'meter second-1'), 709)
            CALL check( nf90_put_att( ncid2, varid(2), 'coordinates', &
                trim(coord) ), 800)
            CALL check( nf90_def_var( ncid2, 'v10', NF90_FLOAT, (/tdimid/),varid(3) ), 801)
            CALL check( nf90_put_att( ncid2, varid(3), 'long_name', &
                '10 metre V wind component'), 504)
            CALL check( nf90_put_att( ncid2, varid(3), 'units', &
                'meter second-1'), 504)
            CALL check( nf90_put_att( ncid2, varid(3), 'coordinates', &
                trim(coord) ), 504)
            CALL check( nf90_def_var( ncid2, 'msl', NF90_FLOAT, (/tdimid/),varid(4) ), 505)
            CALL check( nf90_put_att( ncid2, varid(4), 'long_name', &
                'Mean sea level pressure'), 504)
            CALL check( nf90_put_att( ncid2, varid(4), 'units', &
                'hPa'), 504)
            CALL check( nf90_put_att( ncid2, varid(4), 'coordinates', &
                trim(coord) ), 504)
            CALL check( nf90_def_var( ncid2, 't2m', NF90_FLOAT, (/tdimid/),varid(5) ), 505)
            CALL check( nf90_put_att( ncid2, varid(5), 'long_name', &
                '2 metre temperature'), 504)
            CALL check( nf90_put_att( ncid2, varid(5), 'units', &
                'K'), 504)
            CALL check( nf90_put_att( ncid2, varid(5), 'coordinates', &
                trim(coord) ), 504)
            CALL check( nf90_def_var( ncid2, 'qair', NF90_FLOAT, (/tdimid/),varid(6) ), 505)
            CALL check( nf90_put_att( ncid2, varid(6), 'long_name', &
                'surface air specific humidity'), 504)
            CALL check( nf90_put_att( ncid2, varid(6), 'units', &
                'g/kg'), 504)
            CALL check( nf90_put_att( ncid2, varid(6), 'coordinates', &
                trim(coord) ), 504)
            CALL check( nf90_def_var( ncid2, 'RH', NF90_FLOAT, (/tdimid/),varid(7) ), 505)
            CALL check( nf90_put_att( ncid2, varid(7), 'long_name', &
                'surface air relative humidity'), 604)
            CALL check( nf90_put_att( ncid2, varid(7), 'units', &
                '%'), 605)
            CALL check( nf90_put_att( ncid2, varid(7), 'coordinates', &
                trim(coord) ), 606)

            CALL check( nf90_enddef(ncid2), 516 )
            first = .false.
        end if

        allocate( time(nt), data3d(nx,ny,nt), uvel(nt), vvel(nt), RH(nt) )
        allocate( dew(nt), tair(nt), pair(nt), qair(nt), tmp3d(nx,ny,nt) ) 

        CALL check(nf90_inq_varid(ncid,'time',vid),123)
        CALL check(nf90_get_var(ncid,vid,time),18)

!TEMP
        CALL check(nf90_inq_varid(ncid,'t2m',vid),124)
        CALL check(nf90_get_var(ncid,vid,data3d,&
          start=(/270,320,1/),count=(/nx,ny,nt/) ),18)

        CALL interp_linear(nx,ny,nt,data3d,idx,W,tair)

        CALL check(nf90_close(ncid),360)

!WDIR
        write(filename,100) trim(path),'wdir',i
        write(*,*) trim(filename)
        CALL check(nf90_open(trim(filename),NF90_NOWRITE,ncid),310) 
        CALL check(nf90_inq_varid(ncid,'wdir10',vid),124)
        CALL check(nf90_get_var(ncid,vid,data3d,&
          start=(/270,320,1/),count=(/nx,ny,nt/) ),18)
        CALL check(nf90_close(ncid),360)

!WSPEED
        write(filename,100) trim(path),'wspeed',i
        write(*,*) trim(filename)
        CALL check(nf90_open(trim(filename),NF90_NOWRITE,ncid),310)
        CALL check(nf90_inq_varid(ncid,'si10',vid),124)
        CALL check(nf90_get_var(ncid,vid,tmp3d,&
          start=(/270,320,1/),count=(/nx,ny,nt/) ),18)
        CALL check(nf90_close(ncid),360)

        CALL wind(nx,ny,nt,data3d,tmp3d) 
        CALL interp_linear(nx,ny,nt,data3d,idx,W,uvel)
        CALL interp_linear(nx,ny,nt,tmp3d,idx,W,vvel)

!PAIR
        write(filename,100) trim(path),'pair',i
        write(*,*) trim(filename)
        CALL check(nf90_open(trim(filename),NF90_NOWRITE,ncid),310)
        CALL check(nf90_inq_varid(ncid,'msl',vid),124)
        CALL check(nf90_get_var(ncid,vid,data3d,&
          start=(/270,320,1/),count=(/nx,ny,nt/) ),18)
        CALL check(nf90_close(ncid),360)

        CALL interp_linear(nx,ny,nt,data3d,idx,W,pair)

!QAIR
        write(filename,100) trim(path),'qair',i
!        write(filename,102) i
        write(*,*) trim(filename)
        CALL check(nf90_open(trim(filename),NF90_NOWRITE,ncid),310)
        CALL check(nf90_inq_varid(ncid,'r2',vid),124)
        CALL check(nf90_get_var(ncid,vid,data3d,&
          start=(/270,320,1/),count=(/nx,ny,nt/) ),18)
        CALL check(nf90_close(ncid),360)

        CALL qfill(nx,ny,nt,data3d)
        CALL interp_linear(nx,ny,nt,data3d,idx,W,RH)
        qair = RH
        CALL RH_to_spec(nt,tair,pair,qair)


        CALL check( nf90_put_var( ncid2, varid(1), time, start=(/cnt+1/) ), 517 )
        CALL check( nf90_put_var( ncid2, varid(2), uvel,&
            start=(/cnt+1/) ), 517 )
        CALL check( nf90_put_var( ncid2, varid(3), vvel,&
            start=(/cnt+1/) ), 517 )
        CALL check( nf90_put_var( ncid2, varid(4), pair/100,&
            start=(/cnt+1/) ), 517 )
        CALL check( nf90_put_var( ncid2, varid(5), tair,&
            start=(/cnt+1/) ), 517 )
        CALL check( nf90_put_var( ncid2, varid(6), qair*1000,&
            start=(/cnt+1/) ), 517 )
        CALL check( nf90_put_var( ncid2, varid(7), RH,&
            start=(/cnt+1/) ), 517 )
        
        cnt = cnt+nt
    
        deallocate(time,data3d,uvel,vvel,RH)
        deallocate(dew,tair,pair,qair,tmp3d)

    end do

    CALL check(nf90_close(ncid2),360)
   
END PROGRAM getUERRA

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

SUBROUTINE calc_w(lon,lat,point,w)
    implicit none
    real(kind=8), intent(in) :: lon(2,2), lat(2,2), point(2)
    real(kind=8), intent(out) :: w(2)

    real(kind=8) :: a,b,c,d,e,f,g,h
    real(kind=8) :: alpha,beta

    a = -lon(1,1) + lon(2,1)
    b = -lon(1,1) + lon(1,2)
    c = lon(1,1) - lon(1,2) - lon(2,1) + lon(2,2)
    d = point(1) - lon(1,1)
    e = -lat(1,1) + lat(2,1)
    f = -lat(1,1) + lat(1,2)
    g = lat(1,1) - lat(1,2) - lat(2,1) + lat(2,2)
    h = point(2) - lat(1,1)

    alpha = -(b*e - a*f + d*g - c*h + sqrt(-4*(c*e - a*g)*(d*f - b*h) + &
        (b*e - a*f + d*g - c*h)**2))/(2*c*e - 2*a*g)
    beta  = (b*e - a*f - d*g + c*h + sqrt(-4*(c*e - a*g)*(d*f - b*h) + &
        (b*e - a*f + d*g - c*h)**2))/(2*c*f - 2*b*g)


    if (alpha .lt. 0 .or. beta .lt. 0 .or. &
        alpha .gt. 1 .or. beta .gt. 1) then
      alpha = (-b*e + a*f - d*g + c*h + sqrt(-4*(c*e - a*g)*(d*f - b*h) + &
        (b*e - a*f + d*g - c*h)**2))/(2*c*e - 2*a*g)
      beta  = -((-b*e + a*f + d*g - c*h + sqrt(-4*(c*e - a*g)*(d*f - b*h) + &
        (b*e - a*f + d*g - c*h)**2))/( 2*c*f - 2*b*g))
    end if

    w(1) = alpha
    w(2) = beta
END SUBROUTINE calc_w

SUBROUTINE spec_humid(nt,dew,pair,output)
    implicit none
    integer, intent(in) :: nt
    real(kind=8), intent(in) :: dew(nt),pair(nt)
    real(kind=8), intent(out) :: output(nt)

    integer :: t
    real(kind=8) :: Rdry,Rvap,a1,a3,a4,T0
    real(kind=8) :: esat(nt)

    Rdry = 287.0597
    Rvap = 461.5250
    a1 = 611.21
    a3 = 17.502
    a4 = 32.19
    T0 = 273.16

    do t = 1, nt
      esat(t) = a1*exp(a3*((dew(t)-T0)/(dew(t)-a4)))
    end do

    do t = 1, nt
      output(t) = (Rdry*esat(t)/Rvap)/(pair(t)-(1-Rdry/Rvap)*esat(t))
    end do
END SUBROUTINE spec_humid

SUBROUTINE interp_linear(nx,ny,npt,input,idx,W,output)
    implicit none
    integer, intent(in) :: nx,ny,npt,idx(2)
    real(kind=8), intent(in) :: input(nx,ny,npt),W(2)
    real(kind=8), intent(out) :: output(npt)

    integer :: n
    real(kind=8) :: p11,p12,p21,p22,w11,w12,w21,w22

    do n = 1, npt
      p11 = input(idx(1),idx(2),n)
      p12 = input(idx(1),idx(2)+1,n)
      p21 = input(idx(1)+1,idx(2),n)
      p22 = input(idx(1)+1,idx(2)+1,n)
      w11 = (1-W(1))*(1-W(2))
      w12 = (1-W(1))*W(2)
      w21 = W(1)*(1-W(2))
      w22 = W(1)*W(2)
      output(n) = p11*w11+p12*w12+p21*w21+p22*w22
    end do

END SUBROUTINE interp_linear

SUBROUTINE prepare(nx,ny,lon,lat,olon,olat,near,W)
    implicit none
    integer, intent(in) :: nx,ny
    real(kind=8), intent(in) :: lon(nx,ny),lat(nx,ny),olon,olat
    integer, intent(out) :: near(2)
    real(kind=8), intent(out) :: W(2)

    integer :: idx(2),k,l,x,y
    real(kind=8) :: dis_array(nx,ny)
    real(kind=8) :: point(2),corners(4,2)
    logical :: score,fexit


          do k = 1,nx
            do l = 1,ny
              dis_array(k,l) = sqrt((lon(k,l)-olon)**2+(lat(k,l)-olat)**2)
            end do
          end do

          idx = minloc( dis_array )
          x = idx(1)
          y = idx(2)
          point(1) = olon
          point(2) = olat

          if (x.gt.1 .and. x.lt.nx .and. y.gt.1 .and. y.lt.ny) then
            fexit = .false.
            do k = x-1,x
              if (fexit) EXIT
            do l = y-1,y
              if (fexit) EXIT
              corners(1,:) = (/lon(k,l),lat(k,l)/)
              corners(2,:) = (/lon(k,l+1),lat(k,l+1)/)
              corners(3,:) = (/lon(k+1,l+1),lat(k+1,l+1)/)
              corners(4,:) = (/lon(k+1,l),lat(k+1,l)/)
              CALL in_convex_polygon(corners,point,score)
              if (score) then
                near(1) = k
                near(2) = l
                CALL calc_w(lon(k:k+1,l:l+1),lat(k:k+1,l:l+1),point,W)
                fexit = .true.
              end if
            end do
            end do
          else
            write(*,*) olon,olat,'not found'
          end if

END SUBROUTINE prepare

SUBROUTINE wind(nx,ny,nt,data1,data2)
    implicit none
    integer, intent(in) :: nx,ny,nt
    real(kind=8), intent(inout) :: data1(nx,ny,nt),data2(nx,ny,nt)

    integer :: t
    real(kind=8) :: wdir(nx,ny),wspeed(nx,ny),pi

    pi = 4.0*ATAN(1.0)

    do t = 1, nt
      wdir = (270.0-data1(:,:,t))*pi/180.0
      wspeed = data2(:,:,t)

      data1(:,:,t) = wspeed(:,:)*cos(wdir(:,:))
      data2(:,:,t) = wspeed(:,:)*sin(wdir(:,:))
    end do
END SUBROUTINE wind

SUBROUTINE RH_to_spec(nt,tair,pair,humid)
    implicit none
    integer, intent(in) :: nt
    real(kind=8), intent(in) :: tair(nt),pair(nt)
    real(kind=8), intent(inout) :: humid(nt)

    real(kind=8) :: cff(nt),tempC(nt),press(nt)

    tempC = tair-273.15
    press = pair/100

    cff(:)=(1.0007+3.46E-6*press(:))*6.1121*                   &
     &        EXP(17.502*tempC(:)/(240.97+tempC(:)))
    cff(:)=cff(:)*humid(:)/100.0
    humid(:)=0.62197*(cff(:)/(press(:)-0.378*cff(:)))
    where(humid.lt.0.0) humid=0.0
END SUBROUTINE RH_to_spec

SUBROUTINE qfill(nx,ny,nt,humid)
    implicit none
    integer, intent(in) :: nx,ny,nt
    real(kind=8), intent(inout) :: humid(nx,ny,nt)

    integer :: i,j,t,total,cnt
    real(kind=8) :: suma

    do t = 1, nt
      do
        total =0
        do i = 2,nx-1
        do j = 2,ny-1
          cnt = 0
          suma = 0.0
          if (humid(i,j,t).eq.0) then
            if (humid(i-1,j,t).ne.0) then
              cnt = cnt+1
              suma = suma+humid(i-1,j,t)
            end if
            if (humid(i+1,j,t).ne.0) then
              cnt = cnt+1
              suma = suma+humid(i+1,j,t)
            end if
            if (humid(i,j-1,t).ne.0) then
              cnt = cnt+1
              suma = suma+humid(i,j-1,t)
            end if
            if (humid(i,j+1,t).ne.0) then
              cnt = cnt+1
              suma = suma+humid(i,j+1,t)
            end if
            if (cnt.ge.2) then
              humid(i,j,t) = suma/real(cnt)
              total = total+1
            end if
          end if
        end do
        end do
        if (total.eq.0) EXIT
      end do
    end do

END SUBROUTINE qfill
