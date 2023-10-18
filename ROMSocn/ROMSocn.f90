PROGRAM ROMSocn
    USE NETCDF
    implicit none
    integer :: i,j,k,t,f
    integer :: ncid,dimid,vid,nx,ny,nz,nt,nf
    integer :: oncid,ovid(6),cnt
    integer :: datestr(3),datestp(3),datetmp(3)
    integer :: xdimid,ydimid,zdimid,tdimid
    integer, allocatable :: date(:,:),mask_rho(:,:)
    real(kind=8), allocatable :: h(:,:),Cs_r(:),s_rho(:)
    real(kind=8), allocatable :: zeta(:,:),z_rho(:,:,:),temp(:,:,:),salt(:,:,:)
    real(kind=8), allocatable :: rho(:,:,:),Cs_w(:),s_w(:),Hz(:,:,:),pot(:,:)
    real(kind=8) :: hc,dt,dp(2),time,ncouple
    character(len=200) :: filename,path
    character(len=50) :: fout

    100 format(a,'/ocean_avg_',i4,'-',i2.2,'-',i2.2,'.nc')
    101 format('OCNtoICE_',i4,'-',i2.2,'-',i2.2,'.nc')
    path = '/users/work/mmuzyka/CSDIR/metro_560x600_50lvlsv5/run/baltic'

    ncouple = 15.0
    datestr = (/1999,1,10/)
    datestp = (/1999,1,11/)
!    datestp = (/1993,10,1/)

    datetmp = datestr
    nf = 0
    do 
      nf = nf+1

      if ( all(datetmp.eq.datestp) ) EXIT
      CALL add_day(datetmp,.true.)
    end do
    allocate( date(nf,3) )

    datetmp = datestr
    nf = 0
    do
      nf = nf+1
      date(nf,:) = datetmp

      if ( all(datetmp.eq.datestp) ) EXIT
      CALL add_day(datetmp,.true.)
    end do

    write(filename,100) trim(path),datestr(1),datestr(2),datestr(3)

    CALL check(nf90_open(trim(filename),NF90_NOWRITE,ncid),310)
    CALL check(nf90_inq_dimid(ncid, "xi_rho", dimid),311)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nx),312)

    CALL check(nf90_inq_dimid(ncid, "eta_rho", dimid),313)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=ny),314)

    CALL check(nf90_inq_dimid(ncid, "s_rho", dimid),315)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nz),316)

    CALL check(nf90_inq_dimid(ncid, "ocean_time", dimid),315)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nt),316)

    allocate(h(nx,ny),mask_rho(nx,ny),Cs_r(nz),s_rho(nz),zeta(nx,ny))
    allocate(z_rho(nx,ny,nz),temp(nx,ny,nz),salt(nx,ny,nz),pot(nx,ny))
    allocate(rho(nx,ny,nz),Cs_w(nz+1),s_w(nz+1),Hz(nx,ny,nz))

    CALL check(nf90_inq_varid(ncid,"h",vid),17)
    CALL check(nf90_get_var(ncid,vid,h),18)
    CALL check(nf90_inq_varid(ncid,"hc",vid),19)
    CALL check(nf90_get_var(ncid,vid,hc),20)
    CALL check(nf90_inq_varid(ncid,"Cs_w",vid),21)
    CALL check(nf90_get_var(ncid,vid,Cs_w),22)
    CALL check(nf90_inq_varid(ncid,"s_w",vid),23)
    CALL check(nf90_get_var(ncid,vid,s_w),24)
    CALL check(nf90_inq_varid(ncid,"Cs_r",vid),31)
    CALL check(nf90_get_var(ncid,vid,Cs_r),32)
    CALL check(nf90_inq_varid(ncid,"s_rho",vid),33)
    CALL check(nf90_get_var(ncid,vid,s_rho),34)
    CALL check(nf90_inq_varid(ncid,"mask_rho",vid),33)
    CALL check(nf90_get_var(ncid,vid,mask_rho),34)
    CALL check(nf90_inq_varid(ncid,"dt",vid),35)
    CALL check(nf90_get_var(ncid,vid,dt),36)

    CALL check(nf90_close(ncid),360)



    do f = 1,nf
    write(filename,100) trim(path),date(f,1),date(f,2),date(f,3)
    write(fout,101) date(f,1),date(f,2),date(f,3)
    CALL check( nf90_create( trim(fout),NF90_NETCDF4,oncid ), 500 )
    CALL check( nf90_def_dim( oncid, 'xi_rho', nx, xdimid ), 501 )
    CALL check( nf90_def_dim( oncid, 'eta_rho', ny, ydimid ), 502 )
    CALL check( nf90_def_dim( oncid, 's_rho', nz, zdimid ), 503 )
    CALL check( nf90_def_dim( oncid, 'time', NF90_UNLIMITED, tdimid ), 504 )

    CALL check(nf90_def_var( oncid, 'z_rho', NF90_DOUBLE,&
        (/xdimid,ydimid,zdimid,tdimid/), ovid(1) ), 508)
    CALL check(nf90_def_var( oncid, 'rho', NF90_DOUBLE,&
        (/xdimid,ydimid,zdimid,tdimid/), ovid(2) ), 508)
    CALL check(nf90_def_var( oncid, 'Hz', NF90_DOUBLE,&
        (/xdimid,ydimid,zdimid,tdimid/), ovid(3) ), 508)
    CALL check( nf90_def_var( oncid, 'time', NF90_DOUBLE, (/tdimid/), ovid(4) ), 504)
    CALL check(nf90_def_var( oncid, 'pot1', NF90_DOUBLE,&
        (/xdimid,ydimid,tdimid/), ovid(5) ), 508)
    CALL check(nf90_def_var( oncid, 'pot2', NF90_DOUBLE,&
        (/xdimid,ydimid,tdimid/), ovid(6) ), 508)

    CALL check(nf90_open(trim(filename),NF90_NOWRITE,ncid),310)
    CALL check(nf90_inq_varid(ncid,'ocean_time',vid),120)
    CALL check( nf90_copy_att(ncid, vid, 'long_name', oncid, ovid(4) ), 504)
    CALL check( nf90_copy_att(ncid, vid, 'units', oncid, ovid(4) ), 504)
    CALL check( nf90_copy_att(ncid, vid, 'calendar', oncid, ovid(4) ), 504)
    CALL check( nf90_copy_att(ncid, vid, 'field', oncid, ovid(4) ), 504)
 
    CALL check( nf90_enddef(oncid), 516 )
!OPEN INPUT

    do t = 1,nt
    CALL check(nf90_inq_varid(ncid,"zeta",vid),17)
    CALL check(nf90_get_var(ncid,vid,zeta,start=(/1,1,t/)),18)
    CALL check(nf90_inq_varid(ncid,"temp",vid),17)
    CALL check(nf90_get_var(ncid,vid,temp,start=(/1,1,1,t/)),18)
    CALL check(nf90_inq_varid(ncid,"salt",vid),17)
    CALL check(nf90_get_var(ncid,vid,salt,start=(/1,1,1,t/)),18)
    CALL check(nf90_inq_varid(ncid,"ocean_time",vid),17)
    CALL check(nf90_get_var(ncid,vid,time,start=(/t/)),18)

    do i = 1,nx
    do j = 1,ny
    dp(2) = zeta(i,j)+(zeta(i,j)+h(i,j))*&
            (hc*s_w(1)+h(i,j)*Cs_w(1))/(hc+h(i,j))
    do k = 1,nz
      z_rho(i,j,k) = zeta(i,j)+(zeta(i,j)+h(i,j))*&
            (hc*s_rho(k)+h(i,j)*Cs_r(k))/(hc+h(i,j))
      dp(1) = dp(2)
      dp(2) = zeta(i,j)+(zeta(i,j)+h(i,j))*&
            (hc*s_w(k+1)+h(i,j)*Cs_w(k+1))/(hc+h(i,j))
      Hz(i,j,k) = dp(2)-dp(1)
    end do
    end do
    end do

    CALL dens(nx,ny,nz,temp,salt,z_rho,rho)

    CALL check(nf90_put_var(oncid, ovid(1), z_rho, start=(/1,1,1,t/)), 517 )
    CALL check(nf90_put_var(oncid, ovid(2), rho, start=(/1,1,1,t/)), 517 )
    CALL check(nf90_put_var(oncid, ovid(3), Hz, start=(/1,1,1,t/)), 517 )
    CALL check(nf90_put_var(oncid, ovid(4), time, start=(/t/)), 517 )
    CALL frzmlt(nx,ny,nz,mask_rho,temp,salt,z_rho,Hz,rho,ncouple*dt,0,pot) 
    CALL check(nf90_put_var(oncid, ovid(5), pot, start=(/1,1,t/)), 517 )
    CALL frzmlt(nx,ny,nz,mask_rho,temp,salt,z_rho,Hz,rho,ncouple*dt,1,pot) 
    CALL check(nf90_put_var(oncid, ovid(6), pot, start=(/1,1,t/)), 517 )
    write(*,*) t
    end do

    CALL check(nf90_close(ncid),360)
!CLOSE INPUT

    CALL check( nf90_close(oncid), 600 )

    end do


    deallocate(date,h,mask_rho,Cs_r,s_rho,zeta,z_rho)
    deallocate(temp,salt,rho,Cs_w,s_w,Hz,pot)
   
END PROGRAM ROMSocn

SUBROUTINE check(status,label)
    USE NETCDF
    implicit none
    INTEGER, INTENT(in) :: status , label

    IF (status .ne. nf90_noerr) then
        print*, trim(nf90_strerror(status)),' label=',label
        stop "Stopped"
    END IF
END SUBROUTINE check

SUBROUTINE dens(nx,ny,nz,temp,salt,depth,rho)
    implicit none
    integer, intent(in) :: nx,ny,nz
    real(kind=8), intent(in) :: temp(nx,ny,nz),salt(nx,ny,nz),depth(nx,ny,nz)
    real(kind=8), intent(out) :: rho(nx,ny,nz)

    integer :: i,j,k
    real(kind=8) :: Tp, Tpr10, Ts, Tt, sqrtTs, cff
    real(kind=8), dimension(0:9) :: C
    real(kind=8) :: den,den1,bulk,bulk0,bulk1,bulk2

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
          den=den1*bulk*cff
!          den=den-1000.0

          rho(i,j,k)=den
        END DO
      END DO
    END DO

END SUBROUTINE dens

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

SUBROUTINE frzmlt(nx,ny,nz,mask,temp,salt,depth,Hz,rho,dt,ver,pot)
    implicit none
    integer, intent(in) :: nx,ny,nz,ver
    integer, intent(in) :: mask(nx,ny)
    real(kind=8), intent(in) :: temp(nx,ny,nz),salt(nx,ny,nz),depth(nx,ny,nz)
    real(kind=8), intent(in) :: Hz(nx,ny,nz),rho(nx,ny,nz),dt
    real(kind=8), intent(out) :: pot(nx,ny)

    integer :: i,j,k
    real(kind=8) :: Cp = 3985.0,  z_r_max= -5.0, rho0= 1025.0         ! Joules/kg/degC
    real(kind=8) :: Sold,t_fr,qfraz_prod,meltpot,meltheat,t_freeze
    real(kind=8) :: aver,avet,aves,thick

    pot = 0.0
    do i = 1, nx
    do j = 1, ny
      if (mask(i,j).eq.1) then
        pot(i,j) = 0.0
        do k = 1, nz
          if (depth(i,j,k) < z_r_max) cycle
          Sold = max(0.0,salt(i,j,k))
          t_fr = t_freeze(Sold)
          qfraz_prod= rho(i,j,k)*Cp*Hz(i,j,k) *                             &
     &             max(t_fr - temp(i,j,k),0.0)

          meltpot= rho(i,j,k)*Cp*Hz(i,j,k) *                                &
     &             min(t_fr - temp(i,j,k),0.0)

          pot(i,j) = pot(i,j) + qfraz_prod
          meltheat= max (meltpot, -pot(i,j) )
          pot(i,j) = pot(i,j) + meltheat
        end do
        pot(i,j) = pot(i,j)/dt

        if ( pot(i,j)<=0.0 ) then
        SELECT CASE (ver)
          CASE(0)
            Sold = max(0.0,salt(i,j,nz))
            t_fr = t_freeze(Sold)
            pot(i,j) = rho0*Cp*min(t_fr-temp(i,j,nz),0.0)*5.0/dt
          CASE(1)
            aver = 0.0
            avet = 0.0
            aves = 0.0
            thick = 0.0
            do k = 1, nz
              if (depth(i,j,k) < z_r_max) cycle
              thick = thick + Hz(i,j,k) 
              aver = aver+Hz(i,j,k)*rho(i,j,k)
              avet = avet+Hz(i,j,k)*temp(i,j,k)
              Sold = max(0.0,salt(i,j,k))
              aves = aves+Hz(i,j,k)*Sold
            end do
            if (thick.gt.0.0) then
              aver = aver/thick
              avet = avet/thick
              aves = aves/thick
            end if
            t_fr = t_freeze(aves) 
            pot(i,j) = aver*Cp*min(t_fr-avet,0.0)*thick/dt
          CASE DEFAULT
            write(*,*) 'melting potential method unsupported'
        END SELECT
        end if
      end if
    end do
    end do

END SUBROUTINE frzmlt

    real(kind=8) function t_freeze(s1)
      real(kind=8), intent(in) :: s1
      t_freeze = -0.054*s1
      return
    end function t_freeze
