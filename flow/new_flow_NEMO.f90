PROGRAM new_flow_NEMO
    USE NETCDF
    implicit none
    integer :: date0(3),date(3),date1(3)
    integer :: ncid,dimid,varid,nx,ny,nz,vid(9)
    integer :: xd,yd,zd,td,cnt,ncid2,varid2,julian
    integer :: oncid,i,j,k,vi(2),vj,ui,uj(2)
    real(kind=8), allocatable :: depth(:)
    real(kind=8), allocatable :: lon(:),lat(:),f(:,:),v(:,:,:),mass(:,:,:),dz(:)
    real(kind=8), allocatable :: dx(:,:),dy(:,:),uo(:,:,:),vo(:,:,:),zeta(:,:)
    real(kind=8), allocatable :: thetao(:,:,:),so(:,:,:),tmp4d(:,:,:,:)
    real(kind=8) :: time,tref,Cp,uflow,vflow,uflow_e,uflow_w,vflow_n,vflow_s
    real(kind=8) :: vinsalt,voutsalt,den,den1,mytime
    character(len=200) :: pref,path,fname
    integer, allocatable :: mask(:,:,:),tmp2d(:,:),mask2d(:,:)

    100 format(a,'/',i4.4,'/',i2.2,'/')
    101 format('BAL-NEMO_PHY-',i4.4,i2.2,i2.2,'00.nc')
    102 format('NEMO_rho_',i4.4,i2.2,i2.2,'.nc')
    103 format('BAL-NEMO_PHY-',i4.4,i2.2,i2.2,'12.nc')

    pref = '/users/work/mmuzyka/programs/ROMS_bc/new_files'
    date0 = (/2023,1,1/)
    date1 = (/2023,1,1/)
!    date1 = (/2023,9,10/)
    date = date0
    tref = -2.0
    Cp = 3985.0
    vi = (/45,130/)
    vj = 230
    ui = 156
    uj = (/100,145/)
    

    
    write(path,100) trim(pref),date(1),date(2)
    write(fname,101) date(1),date(2),date(3)

    CALL check(nf90_open(trim(path)//trim(fname),NF90_NOWRITE,ncid),10)
    CALL check(nf90_inq_dimid(ncid, "lon", dimid),11)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nx),12)
    CALL check(nf90_inq_dimid(ncid, "lat", dimid),13)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=ny),14)
    CALL check(nf90_inq_dimid(ncid, "depth", dimid),15)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nz),16)

    allocate( lon(nx), lat(ny), depth(nz), dz(nz), f(nx,ny) )
    allocate( mask(nx,ny,nz),tmp2d(nx,ny),v(nx,ny,nz),mass(nx,ny,nz) )
    allocate( mask2d(nx,ny),thetao(nx,ny,nz),so(nx,ny,nz) )
    allocate( dx(nx,ny),dy(nx,ny),uo(nx,ny,nz),vo(nx,ny,nz) )
    allocate( tmp4d(nx,ny,nz,24) ) 

    CALL check(nf90_inq_varid(ncid,"depth",varid),27)
    CALL check(nf90_get_var(ncid,varid,depth),28)
    CALL check(nf90_inq_varid(ncid,"lon",varid),27)
    CALL check(nf90_get_var(ncid,varid,lon),28)
    CALL check(nf90_inq_varid(ncid,"lat",varid),27)
    CALL check(nf90_get_var(ncid,varid,lat),28)
    CALL check(nf90_inq_varid(ncid,"uo",varid),27)
    CALL check(nf90_get_var(ncid,varid,uo,start=(/1,1,1,1/),&
        count=(/nx,ny,nz,1/)),28)

    CALL check(nf90_close(ncid),35)

!    write(*,*) nx,ny,nz

!    CALL dens(3.0,35.5,-5000,den)
!    write(*,*) 'den=',den

!    CALL calc_den(real(3.0,8),real(35.5,8),real(-5000,8),den1,den)
!    write(*,*) 'den=',den
!    write(*,*) 'den1=',den1

    mask2d = 1
    where(uo(:,:,1).lt.-900) mask2d = 0
    mask2d(1:150,285) = 0

    mask = 0
    mask(360,200,:) = 1
    do 
      cnt = 0 
      do i = 2,nx-1
        do j = 2,ny-1
          if (mask(i,j,1).eq.1) then
            if ( mask(i-1,j,1).eq.0 .and. mask2d(i-1,j).eq.1 ) then
              mask(i-1,j,:) = 1
              cnt = cnt+1
            end if
            if ( mask(i+1,j,1).eq.0 .and. mask2d(i+1,j).eq.1 ) then
              mask(i+1,j,:) = 1
              cnt = cnt+1
            end if
            if ( mask(i,j-1,1).eq.0 .and. mask2d(i,j-1).eq.1 ) then
              mask(i,j-1,:) = 1
              cnt = cnt+1
            end if
            if ( mask(i,j+1,1).eq.0 .and. mask2d(i,j+1).eq.1 ) then
              mask(i,j+1,:) = 1
              cnt = cnt+1
            end if
          end if
        end do
      end do
      if (cnt.eq.0) EXIT
    end do
    where(mask.eq.1 .and. uo.lt.-900) mask = 0

    dz(1) = depth(1)*2
!    write(*,*) 'depth1=',depth(1)
    do i = 2,nz
      dz(i) = (depth(i)-sum(dz(1:i-1)))*2 
    end do
    CALL field(nx,ny,lon,lat,f,dx,dy)
    do i = 1,nx
      do j = 1,ny
        do k = 1,nz
          v(i,j,k) = f(i,j)*dz(k)
        end do
      end do
    end do

!    do i = 1,nz
!        write(*,*) i,dz(i)
!    end do

    cnt = 0
    CALL check(nf90_create( "NEMO_flow_2023testowy.nc", NF90_NETCDF4, oncid ), 200)
    CALL check(nf90_def_dim( oncid, 'time', NF90_UNLIMITED, td ), 201)

    CALL check(nf90_def_var( oncid, 'time', NF90_DOUBLE, (/td/), vid(1)), 217)
    CALL check(nf90_put_att( oncid, vid(1), 'standard_name', 'time'),220)
    CALL check(nf90_put_att( oncid, vid(1), 'long_name', 'Validity time'),220)
!    CALL check(nf90_put_att( oncid, vid(1), 'units', 'days since 1900-01-01 00:00:00'),220)
    CALL check(nf90_put_att( oncid, vid(1), 'units', 'seconds since 1968-05-23 00:00:00'),220)
    CALL check(nf90_put_att( oncid, vid(1), 'calendar', 'proleptic_gregorian'),220)
    CALL check(nf90_put_att( oncid, vid(1), 'axis', 'T'),220)

!    CALL check(nf90_def_var( oncid, 'uflow', NF90_DOUBLE, (/td/), vid(2)), 217)
!    CALL check(nf90_put_att( oncid, vid(2), 'standard_name', 'u-section water mass'),220)
!    CALL check(nf90_put_att( oncid, vid(2), 'units', 'm3/s'),220)

!    CALL check(nf90_def_var( oncid, 'vflow', NF90_DOUBLE, (/td/), vid(3)), 221)
!    CALL check(nf90_put_att( oncid, vid(3), 'standard_name', 'v-section water mass'),222)
!    CALL check(nf90_put_att( oncid, vid(3), 'units', 'm3/s'),223)

    CALL check(nf90_def_var( oncid, 'uflow_e', NF90_DOUBLE, (/td/), vid(4)), 217)
    CALL check(nf90_put_att( oncid, vid(4), 'standard_name', 'u-section water mass in east dir'),220)
    CALL check(nf90_put_att( oncid, vid(4), 'units', 'm3/s'),220)

    CALL check(nf90_def_var( oncid, 'uflow_w', NF90_DOUBLE, (/td/), vid(5)), 217)
    CALL check(nf90_put_att( oncid, vid(5), 'standard_name', 'u-section water mass in west dir'),220)
    CALL check(nf90_put_att( oncid, vid(5), 'units', 'm3/s'),220)

    CALL check(nf90_def_var( oncid, 'vflow_n', NF90_DOUBLE, (/td/), vid(6)), 221)
    CALL check(nf90_put_att( oncid, vid(6), 'standard_name', 'v-section water mass in north dir'),222)
    CALL check(nf90_put_att( oncid, vid(6), 'units', 'm3/s'),223)

    CALL check(nf90_def_var( oncid, 'vflow_s', NF90_DOUBLE, (/td/), vid(7)), 221)
    CALL check(nf90_put_att( oncid, vid(7), 'standard_name', 'v-section water mass in south dir'),222)
    CALL check(nf90_put_att( oncid, vid(7), 'units', 'm3/s'),223)

    CALL check(nf90_def_var( oncid, 'vinsalt', NF90_DOUBLE, (/td/), vid(8)), 221)
    CALL check(nf90_put_att( oncid, vid(8), 'standard_name', 'v-section salt inflow'),222)
    CALL check(nf90_put_att( oncid, vid(8), 'units', 't/s'),223)

    CALL check(nf90_def_var( oncid, 'voutsalt', NF90_DOUBLE, (/td/), vid(9)), 221)
    CALL check(nf90_put_att( oncid, vid(9), 'standard_name', 'v-section salt outflow'),222)
    CALL check(nf90_put_att( oncid, vid(9), 'units', 't/s'),223)

    CALL check( nf90_enddef( oncid ), 208 )
    do
        write(fname,101) date(1),date(2),date(3)
        write(*,*) trim(path)//trim(fname)
        CALL check(nf90_open(trim(path)//trim(fname),NF90_NOWRITE,ncid),10)

        write(fname,103) date(1),date(2),date(3)
        write(*,*) trim(path)//trim(fname)
        CALL check(nf90_open(trim(path)//trim(fname),NF90_NOWRITE,ncid2),10)

        CALL check(nf90_inq_varid(ncid,"uo",varid),27)
        CALL check(nf90_get_var(ncid,varid,tmp4d(:,:,:,1:12)),28)
        CALL check(nf90_inq_varid(ncid2,"uo",varid2),27)
        CALL check(nf90_get_var(ncid2,varid2,tmp4d(:,:,:,13:24)),28)
        uo = sum(tmp4d,dim=4)/24.0
        write(*,*) uo(360,200,1)

        CALL check(nf90_inq_varid(ncid,"vo",varid),27)
        CALL check(nf90_get_var(ncid,varid,tmp4d(:,:,:,1:12)),28)
        CALL check(nf90_inq_varid(ncid2,"vo",varid2),27)
        CALL check(nf90_get_var(ncid2,varid2,tmp4d(:,:,:,13:24)),28)
        vo = sum(tmp4d,dim=4)/24.0
        write(*,*) vo(360,200,1)

        CALL check(nf90_inq_varid(ncid,"thetao",varid),27)
        CALL check(nf90_get_var(ncid,varid,tmp4d(:,:,:,1:12)),28)
        CALL check(nf90_inq_varid(ncid2,"thetao",varid2),27)
        CALL check(nf90_get_var(ncid2,varid2,tmp4d(:,:,:,13:24)),28)
        thetao = sum(tmp4d,dim=4)/24.0
        write(*,*) thetao(360,200,1)

        CALL check(nf90_inq_varid(ncid,"so",varid),27)
        CALL check(nf90_get_var(ncid,varid,tmp4d(:,:,:,1:12)),28)
        CALL check(nf90_inq_varid(ncid2,"so",varid2),27)
        CALL check(nf90_get_var(ncid2,varid2,tmp4d(:,:,:,13:24)),28)
        so = sum(tmp4d,dim=4)/24.0
        write(*,*) so(360,200,1)

        CALL check(nf90_inq_varid(ncid,"time",varid),27)
        CALL check(nf90_get_var(ncid,varid,time,start=(/12/)),28)

        CALL check(nf90_close(ncid),35)
        CALL check(nf90_close(ncid2),35)

!        where(mask.eq.0) temp=0
!        where(mask.eq.0) salt=0
!        CALL dens(nx,ny,nz,temp,salt,depth,den)
!        mass = salt*0.001*den*v*mask
!        heat = (temp-tref)*Cp*den*v*mask  
!        uflow = 0
        uflow_e = 0
        uflow_w = 0
!        vflow = 0
        vflow_n = 0
        vflow_s = 0
        vinsalt = 0
        voutsalt = 0

        do j = uj(1),uj(2)
          do k = 1,nz
            if (mask(ui,j,k).eq.1) then
!              uflow = uflow+uo(ui,j,k)*dy(ui,j)*dz(k)
              if (uo(ui,j,k).gt.0) then
                uflow_e = uflow_e+uo(ui,j,k)*dy(ui,j)*dz(k)
              else
                uflow_w = uflow_w+uo(ui,j,k)*dy(ui,j)*dz(k)
              end if
            end if
          end do 
        end do 
        do i = vi(1),vi(2)
          do k = 1,nz
            if (mask(i,vj,k).eq.1) then
              CALL calc_den(thetao(i,vj,k),so(i,vj,k),-1*depth(k),den1,den)
!              write(*,*) den*0.001*so(i,vj,k)
              if (vo(i,vj,k).gt.0) then
                vflow_n = vflow_n+vo(i,vj,k)*dx(i,vj)*dz(k)
                voutsalt = voutsalt+vo(i,vj,k)*dx(i,vj)*dz(k)*&
                    den*1e-6*so(i,vj,k)
              else
                vflow_s = vflow_s+vo(i,vj,k)*dx(i,vj)*dz(k)
                vinsalt = vinsalt+vo(i,vj,k)*dx(i,vj)*dz(k)*&
                    den*1e-6*so(i,vj,k)
              end if
            end if
          end do 
        end do 

        CALL elapsed(date,(/1968,5,23/),julian)
        mytime = real(julian,8)*86400+12*3600.0 
        cnt = cnt+1
        CALL check( nf90_put_var( oncid, vid(1), mytime, start = (/cnt/) ), 209 )
!        CALL check( nf90_put_var( oncid, vid(2), uflow, start = (/cnt/) ), 209 )
!        CALL check( nf90_put_var( oncid, vid(3), vflow, start = (/cnt/) ), 209 )
        CALL check( nf90_put_var( oncid, vid(4), uflow_e, start = (/cnt/) ), 209 )
        CALL check( nf90_put_var( oncid, vid(5), uflow_w, start = (/cnt/) ), 209 )
        CALL check( nf90_put_var( oncid, vid(6), vflow_n, start = (/cnt/) ), 209 )
        CALL check( nf90_put_var( oncid, vid(7), vflow_s, start = (/cnt/) ), 209 )
        CALL check( nf90_put_var( oncid, vid(8), vinsalt, start = (/cnt/) ), 209 )
        CALL check( nf90_put_var( oncid, vid(9), voutsalt, start = (/cnt/) ), 209 )
 
        if (all(date==date1)) EXIT
        CALL add_day(date,.true.)
        write(path,100) trim(pref),date(1),date(2)
        write(fname,101) date(1),date(2),date(3)
    end do
    CALL check(nf90_close( oncid ), 213 )

!    write(*,*) salt(200,200,1),den(200,200,1),v(200,200,1)
    
    
    CALL check(nf90_create( "field.nc", NF90_NETCDF4, ncid ), 200)
    CALL check(nf90_def_dim( ncid, 'x', nx, xd ), 201)
    CALL check(nf90_def_dim( ncid, 'y', ny, yd ), 202)
    CALL check(nf90_def_dim( ncid, 'z', nz, zd ), 202)
!    CALL check(nf90_def_dim( ncid, 'z', nz, zd ), 203)
    CALL check(nf90_def_var( ncid, 'f', NF90_DOUBLE, (/xd,yd/), vid(1)), 204) 
    CALL check(nf90_def_var( ncid, 'mask', NF90_INT, (/xd,yd,zd/), vid(2)), 204) 
    CALL check(nf90_def_var( ncid, 'dx', NF90_DOUBLE, (/xd,yd/), vid(3)), 204) 
    CALL check(nf90_def_var( ncid, 'dy', NF90_DOUBLE, (/xd,yd/), vid(4)), 204) 
!    CALL check(nf90_def_var( ncid, 'pden1', NF90_DOUBLE, (/xd,yd,zd/), vid(1)), 204) 
!    CALL check(nf90_def_var( ncid, 'den1', NF90_DOUBLE, (/xd,yd,zd/), vid(2)), 205) 
!    CALL check(nf90_def_var( ncid, 'pden2', NF90_DOUBLE, (/xd,yd,zd/), vid(3)), 206) 
!    CALL check(nf90_def_var( ncid, 'den2', NF90_DOUBLE, (/xd,yd,zd/), vid(4)), 207) 
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
!    CALL check( nf90_put_var( ncid, vid(1), pden ), 209 )
!    CALL check( nf90_put_var( ncid, vid(2), den ), 210 )
!    CALL dens_POP(43,43,nz,t,s,d,pden,den)
    CALL check( nf90_put_var( ncid, vid(1), mask(:,:,1)*f*1e-6 ), 211 )
    CALL check( nf90_put_var( ncid, vid(2), mask ), 211 )
    CALL check( nf90_put_var( ncid, vid(3), mask(:,:,1)*dx*1e-3 ), 211 )
    CALL check( nf90_put_var( ncid, vid(4), mask(:,:,1)*dy*1e-3 ), 211 )
!    CALL check( nf90_put_var( ncid, vid(3), pden ), 211 )
!    CALL check( nf90_put_var( ncid, vid(4), den ), 212 )
!    CALL check( nf90_put_var( ncid, vid(5), t(:,1,1) ), 212 )
!    CALL check( nf90_put_var( ncid, vid(6), s(1,:,1) ), 212 )
!    CALL check( nf90_put_var( ncid, vid(7), d(1,1,:) ), 212 )
    CALL check( nf90_close( ncid ), 213 )
     

    deallocate( lon, lat, depth, dz, f, dx, dy )
    deallocate( mask, tmp2d, v, mass, uo, vo ) 
    deallocate( mask2d, thetao, so, tmp4d )

END PROGRAM new_flow_NEMO

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

SUBROUTINE calc_den(temp,salt,Tp,den1,den)
    implicit none
    real(kind=8), intent(in) :: temp,salt,Tp
    real(kind=8), intent(out) :: den1,den
    real(kind=8) :: Tt,Ts,sqrtTs,Tpr10
    real(kind=8) :: C(0:2),cff
    real(kind=8) :: bulk0,bulk1,bulk2,bulk

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

    Tt=MAX(-2.0,temp)
    Ts=MAX(0.0,salt)
    sqrtTs=SQRT(Ts)
    Tpr10=0.1*Tp

    C(0)=Q00+Tt*(Q01+Tt*(Q02+Tt*(Q03+Tt*(Q04+Tt*Q05))))
    C(1)=U00+Tt*(U01+Tt*(U02+Tt*(U03+Tt*U04)))
    C(2)=V00+Tt*(V01+Tt*V02)

    den1=C(0)+Ts*(C(1)+sqrtTs*C(2)+Ts*W00)

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
    bulk=bulk0-Tp*(bulk1-Tp*bulk2)

    cff=1.0/(bulk+Tpr10)
    den=den1*bulk*cff
END SUBROUTINE calc_den

SUBROUTINE dens(temp,salt,depth,den)
    implicit none
    real(kind=8), intent(in) :: temp,salt,depth
    real(kind=8), intent(out) :: den

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

          Tt=MAX(-2.0,temp)
          Ts=MAX(0.0,salt)
          sqrtTs=SQRT(Ts)
          Tp=depth
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

SUBROUTINE field(nx,ny,lon,lat,f,dx,dy)
    implicit none
    integer, intent(in) :: nx,ny
    real(kind=8), intent(in) :: lon(nx),lat(ny)
    real(kind=8), intent(out) :: f(nx,ny),dx(nx,ny),dy(nx,ny)

    real(kind=8) :: tlon(2),tlat(2),a(2),b(2),c,s(2),p
    real(kind=8), dimension(0:nx+1) :: mylon
    real(kind=8), dimension(0:ny+1) :: mylat

    integer :: i,j

    mylon(1:nx) = lon
    mylat(1:ny) = lat

    mylon(0) = 2*lon(1)-lon(2)
    mylon(nx+1) = 2*lon(nx)-lon(nx-1)
    mylat(0) = 2*lat(1)-lat(2)
    mylat(ny+1) = 2*lat(ny)-lat(ny-1)

    do i = 1,nx
      tlon(1) = (mylon(i)+mylon(i-1))*0.5
      tlon(2) = (mylon(i+1)+mylon(i))*0.5
      do j = 1,ny
        tlat(1) = (mylat(j)+mylat(j-1))*0.5
        tlat(2) = (mylat(j+1)+mylat(j))*0.5
        CALL dislonlat(tlon(1),tlat(2),tlon(2),tlat(2),a(1)) 
        CALL dislonlat(tlon(1),tlat(1),tlon(2),tlat(1),a(2)) 
        CALL dislonlat(tlon(1),tlat(1),tlon(1),tlat(2),b(1)) 
        CALL dislonlat(tlon(2),tlat(1),tlon(2),tlat(2),b(2)) 
        CALL dislonlat(tlon(1),tlat(1),tlon(2),tlat(2),c)
        p = 0.5*(a(1)+b(1)+c)
        s(1) = sqrt(p*(p-a(1))*(p-b(1))*(p-c))
        p = 0.5*(a(2)+b(2)+c)
        s(2) = sqrt(p*(p-a(2))*(p-b(2))*(p-c))
        f(i,j) = sum(s)
        dx(i,j) = 0.5*sum(a)
        dy(i,j) = 0.5*sum(b)
      end do
    end do
END SUBROUTINE field

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
