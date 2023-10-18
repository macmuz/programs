PROGRAM flow_ROMS
    USE NETCDF
    implicit none
    integer :: date0(3),date(3),date1(3)
    integer :: ncid,dimid,varid,nx,ny,nz,nt,vid(8)
    integer :: xd,yd,zd,td,cnt,myvarid(4)
    integer :: oncid,i,j,k,t,vi(2),vj,ui,uj(2)
    real(kind=8), allocatable :: lon(:,:),lat(:,:),f(:,:)
    real(kind=8), allocatable :: Cs_w(:),s_w(:),h(:,:),uzeta(:),vzeta(:)
    real(kind=8), allocatable :: Cs_r(:),s_r(:),zeta(:,:,:)
    real(kind=8), allocatable :: u(:,:,:,:),v(:,:,:,:),uslice(:,:),vslice(:,:)
    real(kind=8), allocatable :: lon_u(:,:),lat_u(:,:),lon_v(:,:),lat_v(:,:)
    real(kind=8), allocatable :: dx(:,:),dy(:,:),time(:)
    real(kind=8) :: hc,tref,Cp,uflow_e,uflow_w,vflow_n,vflow_s,div,myh,vinflowb
    real(kind=8) :: top,bottom,dz
    character(len=200) :: pref,path
    integer, allocatable :: mask(:,:),mask_u(:,:),mask_v(:,:)

    100 format(a,'/ocean_avg_',i4.4,'-',i2.2,'-',i2.2,'.nc')
    102 format('ROMS_rho_',i4.4,i2.2,i2.2,'.nc')

    pref = '/users/work/mmuzyka/CSDIR/metro_560x600_uerraGLS2/run/baltic'
!    pref = '/users/work/mmuzyka/CSDIR/metro_560x600_uerraGLS/run/baltic'
    date0 = (/1993,1,2/)
!    date1 = (/1999,1,6/)
    date1 = (/1994,2,6/)
    date = date0
    tref = -2.0
    Cp = 3985.0

    vi = (/48,108/)
    vj = 161
    ui = 120
    uj = (/53,88/)

    
    write(path,100) trim(pref),date(1),date(2),date(3)

    CALL check(nf90_open(trim(path),NF90_NOWRITE,ncid),10)
    CALL check(nf90_inq_dimid(ncid, "xi_rho", dimid),11)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nx),12)
    CALL check(nf90_inq_dimid(ncid, "eta_rho", dimid),13)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=ny),14)
    CALL check(nf90_inq_dimid(ncid, "s_rho", dimid),15)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nz),16)
    CALL check(nf90_inq_dimid(ncid, "ocean_time", dimid),17)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nt),18)

    allocate( lon(nx-1,ny-1), lat(nx-1,ny-1), f(nx,ny) )
    allocate( Cs_w(nz+1), s_w(nz+1), Cs_r(nz), s_r(nz), h(nx,ny) )
    allocate( mask(nx,ny),zeta(nx,ny,nt), time(nt) )
    allocate( u(nx-1,ny,nz,nt),v(nx,ny-1,nz,nt) )
    allocate( lon_u(nx-1,ny), lat_u(nx-1,ny), lon_v(nx,ny-1), lat_v(nx,ny-1))
    allocate( mask_u(nx-1,ny), mask_v(nx,ny-1), dx(nx,ny), dy(nx,ny) )
    allocate( uzeta( 1+uj(2)-uj(1) ),vzeta( 1+vi(2)-vi(1) ) )
    allocate( uslice( 1+uj(2)-uj(1), nz ),vslice( 1+vi(2)-vi(1), nz ) )

    div = 1.0/real(nt,8)

    CALL check(nf90_inq_varid(ncid,"mask_rho",varid),27)
    CALL check(nf90_get_var(ncid,varid,mask),28)
    CALL check(nf90_inq_varid(ncid,"mask_u",varid),27)
    CALL check(nf90_get_var(ncid,varid,mask_u),28)
    CALL check(nf90_inq_varid(ncid,"mask_v",varid),27)
    CALL check(nf90_get_var(ncid,varid,mask_v),28)
    CALL check(nf90_inq_varid(ncid,"lon_psi",varid),29)
    CALL check(nf90_get_var(ncid,varid,lon),30)
    CALL check(nf90_inq_varid(ncid,"lat_psi",varid),31)
    CALL check(nf90_get_var(ncid,varid,lat),32)
    CALL check(nf90_inq_varid(ncid,"lon_u",varid),29)
    CALL check(nf90_get_var(ncid,varid,lon_u),30)
    CALL check(nf90_inq_varid(ncid,"lat_u",varid),31)
    CALL check(nf90_get_var(ncid,varid,lat_u),32)
    CALL check(nf90_inq_varid(ncid,"lon_v",varid),29)
    CALL check(nf90_get_var(ncid,varid,lon_v),30)
    CALL check(nf90_inq_varid(ncid,"lat_v",varid),31)
    CALL check(nf90_get_var(ncid,varid,lat_v),32)
    CALL check(nf90_inq_varid(ncid,"h",varid),33)
    CALL check(nf90_get_var(ncid,varid,h),34)
    CALL check(nf90_inq_varid(ncid,"hc",varid),35)
    CALL check(nf90_get_var(ncid,varid,hc),36)
    CALL check(nf90_inq_varid(ncid,"Cs_w",varid),37)
    CALL check(nf90_get_var(ncid,varid,Cs_w),38)
    CALL check(nf90_inq_varid(ncid,"s_w",varid),39)
    CALL check(nf90_get_var(ncid,varid,s_w),40)
    CALL check(nf90_inq_varid(ncid,"Cs_r",varid),41)
    CALL check(nf90_get_var(ncid,varid,Cs_r),42)
    CALL check(nf90_inq_varid(ncid,"s_rho",varid),43)
    CALL check(nf90_get_var(ncid,varid,s_r),44)

    CALL check(nf90_close(ncid),35)

    CALL fieldROMS(nx,ny,lon,lat,f)
    CALL calc_dx(nx,ny,lon_u,lat_u,dx)
    CALL calc_dy(nx,ny,lon_v,lat_v,dy)

    cnt = 0
    CALL check(nf90_create( "ROMS_flow_uerraGLS2.nc", NF90_NETCDF4, oncid ), 200)
    CALL check(nf90_def_dim( oncid, 'time', NF90_UNLIMITED, td ), 201)

    CALL check(nf90_def_var( oncid, 'time', NF90_DOUBLE, (/td/), vid(1)), 217)
    CALL check(nf90_put_att( oncid, vid(1), 'standard_name', 'time'),220)
    CALL check(nf90_put_att( oncid, vid(1), 'long_name', 'Validity time'),220)
    CALL check(nf90_put_att( oncid, vid(1), 'units', 'seconds since 1968-05-23 00:00:00'),220)
    CALL check(nf90_put_att( oncid, vid(1), 'calendar', 'proleptic_gregorian'),220)
    CALL check(nf90_put_att( oncid, vid(1), 'axis', 'T'),220)

!    CALL check(nf90_def_var( oncid, 'tots', NF90_DOUBLE, (/td/), vid(2)), 217)
!    CALL check(nf90_put_att( oncid, vid(2), 'standard_name', 'total salinity'),220)
!    CALL check(nf90_put_att( oncid, vid(2), 'units', 'ton'),220)

!    CALL check(nf90_def_var( oncid, 'toth', NF90_DOUBLE, (/td/), vid(3)), 221)
!    CALL check(nf90_put_att( oncid, vid(3), 'standard_name', 'total heat for Tref=-2.0 celsius'),222)
!    CALL check(nf90_put_att( oncid, vid(3), 'units', 'GJ'),223)

    CALL check(nf90_def_var( oncid, 'uflow_e', NF90_DOUBLE, (/td/), vid(4)), 217)
    CALL check(nf90_put_att( oncid, vid(4), 'standard_name', &
        'u-section water mass in east dir'),220)
    CALL check(nf90_put_att( oncid, vid(4), 'units', 'm3/s'),220) 

    CALL check(nf90_def_var( oncid, 'uflow_w', NF90_DOUBLE, (/td/), vid(5)), 217)
    CALL check(nf90_put_att( oncid, vid(5), 'standard_name', &
        'u-section water mass in west dir'),220)
    CALL check(nf90_put_att( oncid, vid(5), 'units', 'm3/s'),220) 

    CALL check(nf90_def_var( oncid, 'vflow_n', NF90_DOUBLE, (/td/), vid(6)), 221)
    CALL check(nf90_put_att( oncid, vid(6), 'standard_name', &
        'v-section water mass in north dir'),222)
    CALL check(nf90_put_att( oncid, vid(6), 'units', 'm3/s'),223) 

    CALL check(nf90_def_var( oncid, 'vflow_s', NF90_DOUBLE, (/td/), vid(7)), 221)
    CALL check(nf90_put_att( oncid, vid(7), 'standard_name', &
        'v-section water mass in south dir'),222)  
    CALL check(nf90_put_att( oncid, vid(7), 'units', 'm3/s'),223)    

    CALL check(nf90_def_var( oncid, 'vinflowb', NF90_DOUBLE, (/td/), vid(8)), 221)
    CALL check(nf90_put_att( oncid, vid(8), 'standard_name', &
        'v-section bottom water mass in south dir'),222)  
    CALL check(nf90_put_att( oncid, vid(8), 'units', 'm3/s'),223)    

    CALL check( nf90_enddef( oncid ), 208 )
    do
        write(*,*) trim(path)
        CALL check(nf90_open(trim(path),NF90_NOWRITE,ncid),10)
        CALL check(nf90_inq_varid(ncid,"u",myvarid(1)),127)
        CALL check(nf90_inq_varid(ncid,"v",myvarid(2)),129)
        CALL check(nf90_inq_varid(ncid,"zeta",myvarid(3)),131)
        CALL check(nf90_inq_varid(ncid,"ocean_time",myvarid(4)),133)

        CALL check(nf90_get_var(ncid,myvarid(1),u),128)
        CALL check(nf90_get_var(ncid,myvarid(2),v),130)
        CALL check(nf90_get_var(ncid,myvarid(3),zeta),132)
        CALL check(nf90_get_var(ncid,myvarid(4),time),134)

        CALL check(nf90_close(ncid),35)
        do t = 1,nt
          where(mask.eq.0) zeta(:,:,t)=0.0
          do k = 1,nz
            where(mask_u.eq.0) u(:,:,k,t)=0.0
            where(mask_v.eq.0) v(:,:,k,t)=0.0
          end do
        end do

        uflow_e = 0
        uflow_w = 0
        vflow_n = 0
        vflow_s = 0
        vinflowb = 0

        uzeta = div*sum(zeta(ui,uj(1):uj(2),:),2)
        vzeta = div*sum(zeta(vi(1):vi(2),vj,:),2)

        uslice = 0.5*div*(sum(u(ui-1,uj(1):uj(2),:,:),3)+&
            sum(u(ui,uj(1):uj(2),:,:),3))
        vslice = 0.5*div*(sum(v(vi(1):vi(2),vj-1,:,:),3)+&
            sum(v(vi(1):vi(2),vj,:,:),3))
!          CALL calc_dp(nx,ny,nz,zeta,h,hc,s_r,Cs_r,mask,depth)
!          CALL dens(nx,ny,nz,temp,salt,depth,den)
!          CALL calc_v(nx,ny,nz,zeta,h,hc,s_w,Cs_w,mask,f,v)
!          mass = salt*0.001*den*v
!          heat = (temp-tref)*Cp*den*v          

        do j = 1,1+uj(2)-uj(1)
          if (mask(ui,uj(1)+j-1)) then
            myh = h(ui,uj(1)+j-1)
            top = -1*(uzeta(j)+(uzeta(j)+myh)*(hc*s_w(1)+myh*Cs_w(1))/(hc+myh))
            do k = 1,nz
              bottom = top
              top = -1*(uzeta(j)+(uzeta(j)+myh)*(hc*s_w(k+1)+myh*Cs_w(k+1))/(hc+myh))
              dz = abs(top-bottom)
              if (uslice(j,k).gt.0) then
                uflow_e = uflow_e+uslice(j,k)*dy(ui,uj(1)+j-1)*dz
              else
                uflow_w = uflow_w+uslice(j,k)*dy(ui,uj(1)+j-1)*dz
              end if
            end do
          end if
        end do
        do i = 1,1+vi(2)-vi(1)
          if (mask(vi(1)+i-1,vj)) then
            myh = h(vi(1)+i-1,vj)
            top = -1*(uzeta(j)+(uzeta(j)+myh)*(hc*s_w(1)+myh*Cs_w(1))/(hc+myh))
            do k = 1,nz
              bottom = top
              top = -1*(uzeta(j)+(uzeta(j)+myh)*(hc*s_w(k+1)+myh*Cs_w(k+1))/(hc+myh))
              dz = abs(top-bottom)
              if (vslice(i,k).gt.0) then
                vflow_n = vflow_n+vslice(i,k)*dx(vi(1)+i-1,vj)*dz
              else
                vflow_s = vflow_s+vslice(i,k)*dx(vi(1)+i-1,vj)*dz
                if (0.5*(top+bottom).gt.25.0) then
                  vinflowb = vinflowb+vslice(i,k)*dx(vi(1)+i-1,vj)*dz
                end if 
              end if
            end do
          end if
        end do
 
        cnt = cnt+1
        CALL check( nf90_put_var( oncid, vid(1), 0.5*sum(time(int(nt/2):int(nt/2)+1)),&
            start = (/cnt/) ), 209 )
!          CALL check( nf90_put_var( oncid, vid(2), sum(mass), start = (/cnt/) ), 209 )
!          CALL check( nf90_put_var( oncid, vid(3), sum(heat)*1e-9, start = (/cnt/) ), 209 )
          CALL check( nf90_put_var( oncid, vid(4), uflow_e, start = (/cnt/) ), 209 )
          CALL check( nf90_put_var( oncid, vid(5), uflow_w, start = (/cnt/) ), 209 )
          CALL check( nf90_put_var( oncid, vid(6), vflow_n, start = (/cnt/) ), 209 )
          CALL check( nf90_put_var( oncid, vid(7), vflow_s, start = (/cnt/) ), 209 )
          CALL check( nf90_put_var( oncid, vid(8), vinflowb, start = (/cnt/) ), 209 )
 
        if (all(date==date1)) EXIT
        CALL add_day(date,.true.)
        write(path,100) trim(pref),date(1),date(2),date(3)
    end do
    CALL check(nf90_close( oncid ), 213 )

!    write(*,*) salt(200,200,1),den(200,200,1),v(200,200,1)
    
    CALL check(nf90_create( "field2d.nc", NF90_NETCDF4, ncid ), 200)
    CALL check(nf90_def_dim( ncid, 'x', nx, xd ), 201)
    CALL check(nf90_def_dim( ncid, 'y', ny, yd ), 202)
    CALL check(nf90_def_dim( ncid, 'z', nz, zd ), 203)
    CALL check(nf90_def_var( ncid, 'f', NF90_DOUBLE, (/xd,yd/), vid(1)), 204) 
!    CALL check(nf90_def_var( ncid, 'salt', NF90_DOUBLE, (/xd,yd,zd/), vid(2)), 204) 
!    CALL check(nf90_def_var( ncid, 'den', NF90_DOUBLE, (/xd,yd,zd/), vid(3)), 205) 
!    CALL check(nf90_def_var( ncid, 'v', NF90_DOUBLE, (/xd,yd,zd/), vid(4)), 206) 
!    CALL check(nf90_def_var( ncid, 'depth', NF90_DOUBLE, (/xd,yd,zd/), vid(5)), 207)
    CALL check(nf90_def_var( ncid, 'dx', NF90_DOUBLE, (/xd,yd/), vid(6)), 204)
    CALL check(nf90_def_var( ncid, 'dy', NF90_DOUBLE, (/xd,yd/), vid(7)), 204) 
!    CALL check(nf90_def_var( ncid, 'x', NF90_DOUBLE, (/xd/), vid(5)), 217)
!    CALL check(nf90_put_att( ncid, vid(5), 'standard_name', 'temperature' ),220)
!    CALL check(nf90_put_att( ncid, vid(5), 'units', 'degrees_C' ),220)
!    CALL check(nf90_put_att( ncid, vid(5), 'axis', 'X' ),220)
!    CALL check(nf90_def_var( ncid, 'y', NF90_DOUBLE, (/yd/), vid(6)), 227)
!    CALL check(nf90_put_att( ncid, vid(6), 'standard_name', 'salinity' ),220)
!    CALL check(nf90_put_att( ncid, vid(6), 'units', '0.001' ),220)
!    CALL check(nf90_put_att( ncid, vid(6), 'axis', 'Y' ),220)
!    CALL check(nf90_def_var( ncid, 'z', NF90_DOUBLE, (/zd/), vid(7)), 237)
!    CALL check(nf90_put_att( ncid, vid(7), 'standard_name', 'depth' ),220)
!    CALL check(nf90_put_att( ncid, vid(7), 'units', 'm' ),220)
!    CALL check(nf90_put_att( ncid, vid(7), 'axis', 'Z' ),220)
!    CALL check(nf90_put_att( ncid, vid(7), 'positive', 'down' ),220)
    CALL check( nf90_enddef( ncid ), 208 )
!    CALL dens(43,43,nz,t,s,d,pden,den)
!    CALL check( nf90_put_var( ncid, vid(2), salt ), 209 )
!    CALL check( nf90_put_var( ncid, vid(3), den ), 210 )
!    CALL dens_POP(43,43,nz,t,s,d,pden,den)
    CALL check( nf90_put_var( ncid, vid(1), mask*f*1e-6 ), 211 )
!    CALL check( nf90_put_var( ncid, vid(4), v ), 211 )
!    CALL check( nf90_put_var( ncid, vid(5), depth ), 212 )
    CALL check( nf90_put_var( ncid, vid(6), mask*dx*1e-3 ), 211 )
    CALL check( nf90_put_var( ncid, vid(7), mask*dy*1e-3 ), 211 )
!    CALL check( nf90_put_var( ncid, vid(5), t(:,1,1) ), 212 )
!    CALL check( nf90_put_var( ncid, vid(6), s(1,:,1) ), 212 )
!    CALL check( nf90_put_var( ncid, vid(7), d(1,1,:) ), 212 )
    CALL check(nf90_close( ncid ), 213 )
     

    deallocate( lon, lat, f, uzeta, vzeta)
    deallocate( h, Cs_w, s_w, Cs_r, s_r )
    deallocate( mask, u, v, zeta ) 
    deallocate( dx,dy,time, uslice, vslice )
    deallocate( lon_u,lat_u,lon_v,lat_v,mask_u,mask_v)

END PROGRAM flow_ROMS

SUBROUTINE check(status,label)
    USE NETCDF
    implicit none
    INTEGER, INTENT(in) :: status , label

    IF (status .ne. nf90_noerr) then
        print*, trim(nf90_strerror(status)),' label=',label
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


SUBROUTINE dens(nx,ny,nz,temp,salt,depth,den)
    implicit none
    integer, intent(in) :: nx,ny,nz
    real(kind=8), intent(in) :: temp(nx,ny,nz),salt(nx,ny,nz),depth(nx,ny,nz)
    real(kind=8), intent(out) :: den(nx,ny,nz)

    integer :: i,j,k
    real(kind=8) :: Tp, Tpr10, Ts, Tt, sqrtTs, cff
    real(kind=8), dimension(0:9) :: C
    real(kind=8) :: den1,bulk,bulk0,bulk1,bulk2

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
          den(i,j,k)=den1*bulk*cff
        END DO
      END DO
    END DO

END SUBROUTINE dens

SUBROUTINE dens_POP(nx,ny,nz,temp,salt,depth,pden,den)

    implicit none
    integer, intent(in) :: nx,ny,nz
    real(kind=8), intent(in) :: temp(nx,ny,nz),salt(nx,ny,nz),depth(nx,ny,nz)
    real(kind=8), intent(out) :: pden(nx,ny,nz),den(nx,ny,nz)

    integer :: i,j,k
    real(kind=8) :: p,TQ,SQ,SQR
    real(kind=8) :: WORK1,WORK2,DENOMK

    real(kind=8) ::                                                        &
      mwjfnums0t0, mwjfnums0t1, mwjfnums0t2, mwjfnums0t3,              &
      mwjfnums1t0, mwjfnums1t1, mwjfnums2t0,                           &
      mwjfdens0t0, mwjfdens0t1, mwjfdens0t2, mwjfdens0t3, mwjfdens0t4, &
      mwjfdens1t0, mwjfdens1t1, mwjfdens1t3,                           &
      mwjfdensqt0, mwjfdensqt2

    real(kind=8), parameter ::                     &
      mwjfnp0s0t0 =   9.99843699e+2 * 0.001, &
      mwjfnp0s0t1 =   7.35212840e+0 * 0.001, &
      mwjfnp0s0t2 =  -5.45928211e-2 * 0.001, &
      mwjfnp0s0t3 =   3.98476704e-4 * 0.001, &
      mwjfnp0s1t0 =   2.96938239e+0 * 0.001, &
      mwjfnp0s1t1 =  -7.23268813e-3 * 0.001, &
      mwjfnp0s2t0 =   2.12382341e-3 * 0.001, &
      mwjfnp1s0t0 =   1.04004591e-2 * 0.001, &
      mwjfnp1s0t2 =   1.03970529e-7 * 0.001, &
      mwjfnp1s1t0 =   5.18761880e-6 * 0.001, &
      mwjfnp2s0t0 =  -3.24041825e-8 * 0.001, &
      mwjfnp2s0t2 =  -1.23869360e-11* 0.001

    real(kind=8), parameter ::          &
      mwjfdp0s0t0 =   1.0e+0,         &
      mwjfdp0s0t1 =   7.28606739e-3,  &
      mwjfdp0s0t2 =  -4.60835542e-5,  &
      mwjfdp0s0t3 =   3.68390573e-7,  &
      mwjfdp0s0t4 =   1.80809186e-10, &
      mwjfdp0s1t0 =   2.14691708e-3,  &
      mwjfdp0s1t1 =  -9.27062484e-6,  &
      mwjfdp0s1t3 =  -1.78343643e-10, &
      mwjfdp0sqt0 =   4.76534122e-6,  &
      mwjfdp0sqt2 =   1.63410736e-9,  &
      mwjfdp1s0t0 =   5.30848875e-6,  &
      mwjfdp2s0t3 =  -3.03175128e-16, &
      mwjfdp3s0t1 =  -1.27934137e-17

    DO i=1,nx
      DO j=1,ny
        DO k=1,nz
            p = 10.0*(0.059808*(exp(-0.025*depth(i,j,k)) - 1.0)+&
                0.100766*depth(i,j,k) + 2.28405e-7*depth(i,j,k)**2)
            TQ = temp(i,j,k)
            SQ = salt(i,j,k)
            SQR = sqrt(SQ)

            !***
      !*** first calculate numerator of MWJF density [P_1(S,T,p)]
      !***

      mwjfnums0t0 = mwjfnp0s0t0 + p*(mwjfnp1s0t0 + p*mwjfnp2s0t0)
      mwjfnums0t1 = mwjfnp0s0t1
      mwjfnums0t2 = mwjfnp0s0t2 + p*(mwjfnp1s0t2 + p*mwjfnp2s0t2)
      mwjfnums0t3 = mwjfnp0s0t3
      mwjfnums1t0 = mwjfnp0s1t0 + p*mwjfnp1s1t0
      mwjfnums1t1 = mwjfnp0s1t1
      mwjfnums2t0 = mwjfnp0s2t0

      WORK1 = mwjfnums0t0 + TQ * (mwjfnums0t1 + TQ * (mwjfnums0t2 + &
              mwjfnums0t3 * TQ)) + SQ * (mwjfnums1t0 +              &
              mwjfnums1t1 * TQ + mwjfnums2t0 * SQ)

      !***
      !*** now calculate denominator of MWJF density [P_2(S,T,p)]
      !***

      mwjfdens0t0 = mwjfdp0s0t0 + p*mwjfdp1s0t0
      mwjfdens0t1 = mwjfdp0s0t1 + p**3 * mwjfdp3s0t1
      mwjfdens0t2 = mwjfdp0s0t2
      mwjfdens0t3 = mwjfdp0s0t3 + p**2 * mwjfdp2s0t3
      mwjfdens0t4 = mwjfdp0s0t4
      mwjfdens1t0 = mwjfdp0s1t0
      mwjfdens1t1 = mwjfdp0s1t1
      mwjfdens1t3 = mwjfdp0s1t3
      mwjfdensqt0 = mwjfdp0sqt0
      mwjfdensqt2 = mwjfdp0sqt2

      WORK2 = mwjfdens0t0 + TQ * (mwjfdens0t1 + TQ * (mwjfdens0t2 +    &
           TQ * (mwjfdens0t3 + mwjfdens0t4 * TQ))) +                   &
           SQ * (mwjfdens1t0 + TQ * (mwjfdens1t1 + TQ*TQ*mwjfdens1t3)+ &
           SQR * (mwjfdensqt0 + TQ*TQ*mwjfdensqt2))

      DENOMK = 1.0/WORK2

            pden(i,j,k) = 1000*WORK1*DENOMK
            den(i,j,k) = 0
        END DO
      END DO
    END DO

!pressure = 0.059808_r8*(exp(-0.025_r8*depth) - c1)     &
!            + 0.100766_r8*depth + 2.28405e-7_r8*depth**2

!/users/work/jakacki/CESM/eBaltic-cesm/cesm1_0_1_io01/models/ocn/pop2/source/state_mod.F90
END SUBROUTINE dens_POP 

SUBROUTINE dislonlat(lon1,lat1,lon2,lat2,d)
    implicit none
    real(kind=8), intent(in) :: lon1,lat1,lon2,lat2
    real(kind=8), intent(out) :: d
 
    real(kind=8) :: d2r,R,dlon,dlat,a,c

    d2r = 4.D0*DATAN(1.D0)/180.D0
    R = 6371.D0
    dlon = abs(lon2-lon1)*d2r
    dlat = abs(lat2-lat1)*d2r

    a = sin(dlat*0.5)**2+cos(lat1*d2r)*cos(lat2*d2r)*& 
        sin(dlon*0.5)**2
    c = 2*atan2(sqrt(a), sqrt(1-a))
    d = R*c*1000.D0
    
END SUBROUTINE dislonlat

SUBROUTINE fieldROMS(nx,ny,lon,lat,f)
    implicit none
    integer, intent(in) :: nx,ny
    real(kind=8), intent(in) :: lon(nx-1,ny-1),lat(nx-1,ny-1)
    real(kind=8), intent(out) :: f(nx,ny)

    real(kind=8) :: a(2),b(2),c,s(2),p
    real(kind=8), dimension(0:nx,0:ny) :: mylon
    real(kind=8), dimension(0:nx,0:ny) :: mylat

    integer :: i,j

    mylon(1:nx-1,1:ny-1) = lon
    mylat(1:nx-1,1:ny-1) = lat

    mylon(1:nx-1,0) = 2*lon(1:nx-1,1)-lon(1:nx-1,2)
    mylon(1:nx-1,ny) = 2*lon(1:nx-1,ny-1)-lon(1:nx-1,ny-2)
    mylon(0,:) = 2*lon(1,:)-lon(2,:)
    mylon(nx,:) = 2*lon(nx-1,:)-lon(nx-2,:)

    mylat(0,1:ny-1) = 2*lat(1,1:ny-1)-lat(2,1:ny-1)
    mylat(nx,1:ny-1) = 2*lat(nx-1,1:ny-1)-lat(nx-2,1:ny-1)
    mylat(:,0) = 2*lat(:,1)-lat(:,2)
    mylat(:,ny) = 2*lat(:,ny-1)-lat(:,ny-2)

    do i = 1,nx
      do j = 1,ny
        CALL dislonlat(mylon(i-1,j),mylat(i-1,j),mylon(i,j),mylat(i,j),a(1)) 
        CALL dislonlat(mylon(i-1,j-1),mylat(i-1,j-1),mylon(i,j-1),mylat(i,j-1),a(2)) 
        CALL dislonlat(mylon(i-1,j-1),mylat(i-1,j-1),mylon(i-1,j),mylat(i-1,j),b(1)) 
        CALL dislonlat(mylon(i,j-1),mylat(i,j-1),mylon(i,j),mylat(i,j),b(2)) 
        CALL dislonlat(mylon(i-1,j-1),mylat(i-1,j-1),mylon(i,j),mylat(i,j),c)
        p = 0.5*(a(1)+b(1)+c)
        s(1) = sqrt(p*(p-a(1))*(p-b(1))*(p-c))
        p = 0.5*(a(2)+b(2)+c)
        s(2) = sqrt(p*(p-a(2))*(p-b(2))*(p-c))
        f(i,j) = sum(s)
      end do
    end do
END SUBROUTINE fieldROMS

SUBROUTINE calc_dx(nx,ny,lon_u,lat_u,dx)
    implicit none
    integer, intent(in) :: nx,ny
    real(kind=8), intent(in) :: lon_u(nx-1,ny),lat_u(nx-1,ny)
    real(kind=8), intent(out) :: dx(nx,ny)

    real(kind=8), dimension(0:nx,ny) :: mylon
    real(kind=8), dimension(0:nx,ny) :: mylat

    integer :: i,j

    mylon(1:nx-1,1:ny) = lon_u
    mylat(1:nx-1,1:ny) = lat_u

    mylon(0,:) = 2*lon_u(1,:)-lon_u(2,:)
    mylon(nx,:) = 2*lon_u(nx-1,:)-lon_u(nx-2,:)

    mylat(0,:) = 2*lat_u(1,:)-lat_u(2,:)
    mylat(nx,:) = 2*lat_u(nx-1,:)-lat_u(nx-2,:)

    do i = 1,nx
      do j = 1,ny
        CALL dislonlat(mylon(i-1,j),mylat(i-1,j),mylon(i,j),mylat(i,j),dx(i,j))
      end do
    end do
END SUBROUTINE calc_dx

SUBROUTINE calc_dy(nx,ny,lon_v,lat_v,dy)
    implicit none
    integer, intent(in) :: nx,ny
    real(kind=8), intent(in) :: lon_v(nx,ny-1),lat_v(nx,ny-1)
    real(kind=8), intent(out) :: dy(nx,ny)

    real(kind=8), dimension(nx,0:ny) :: mylon
    real(kind=8), dimension(nx,0:ny) :: mylat

    integer :: i,j

    mylon(1:nx,1:ny-1) = lon_v
    mylat(1:nx,1:ny-1) = lat_v

    mylon(:,0) = 2*lon_v(:,1)-lon_v(:,2)
    mylon(:,ny) = 2*lon_v(:,ny-1)-lon_v(:,ny-2)

    mylat(:,0) = 2*lat_v(:,1)-lat_v(:,2)
    mylat(:,ny) = 2*lat_v(:,ny-1)-lat_v(:,ny-2)

    do i = 1,nx
      do j = 1,ny
        CALL dislonlat(mylon(i,j-1),mylat(i,j-1),mylon(i,j),mylat(i,j),dy(i,j))
      end do
    end do
END SUBROUTINE calc_dy

SUBROUTINE calc_dp(nx,ny,nz,zeta,h,hc,s_r,Cs_r,mask,dp)
    implicit none
    integer, intent(in) :: nx,ny,nz
    real(kind=8), intent(in) :: zeta(nx,ny),h(nx,ny),hc,s_r(nz),Cs_r(nz)
    integer, intent(in) :: mask(nx,ny)
    real(kind=8), intent(out) :: dp(nx,ny,nz)

    integer :: i,j,k

    do i = 1,nx
    do j = 1,ny
    if (mask(i,j).eq.1) then
      do k = 1,nz
        dp(i,j,k) = -1*(zeta(i,j)+(zeta(i,j)+h(i,j))*(hc*s_r(k)+h(i,j)*Cs_r(k))/&
            (hc+h(i,j)))
      end do
    end if
    end do
    end do

END SUBROUTINE calc_dp

SUBROUTINE calc_v(nx,ny,nz,zeta,h,hc,s_w,Cs_w,mask,f,v)
    implicit none
    integer, intent(in) :: nx,ny,nz
    real(kind=8), intent(in) :: zeta(nx,ny),h(nx,ny),hc,s_w(nz+1),Cs_w(nz+1)
    real(kind=8), intent(in) :: f(nx,ny)
    integer, intent(in) :: mask(nx,ny)
    real(kind=8), intent(out) :: v(nx,ny,nz)

    real(kind=8) :: depth(nx,ny,nz+1)
    integer :: i,j,k

    v = 0
    CALL calc_dp(nx,ny,nz+1,zeta,h,hc,s_w,Cs_w,mask,depth)

    do i = 1,nx
      do j = 1,ny
        if (mask(i,j).eq.1) then
          do k = 1,nz
            v(i,j,k) = f(i,j)*abs(depth(i,j,k+1)-depth(i,j,k))
          end do
        end if
      end do
    end do

     
    
END SUBROUTINE calc_v
