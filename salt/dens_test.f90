PROGRAM dens_test
    USE NETCDF
    implicit none
    integer :: date0(3),date(3),date1(3)
    integer :: ncid,dimid,varid,nx,ny,nz,vid(8)
    integer :: xd,yd,zd
    integer :: oncid,i
    real, allocatable :: depth(:)
    real(kind=8), allocatable :: t(:,:,:),s(:,:,:),d(:,:,:)
    real(kind=8), allocatable :: pden(:,:,:),den(:,:,:)
    character(len=200) :: pref,path,fname,oname

    100 format(a,'/',i4.4,'/',i2.2,'/')
    101 format('CMEMS_BAL_PHY_reanalysis_dailymeans_',i4.4,i2.2,i2.2,'.nc')
    102 format('NEMO_rho_',i4.4,i2.2,i2.2,'.nc')

    pref = '/users/work/mmuzyka/programs/ROMS_bc/files'
    date0 = (/1999,1,1/)
    date1 = (/1999,1,5/)
    date = date0

    
    write(path,100) trim(pref),date(1),date(2)
    write(fname,101) date(1),date(2),date(3)

    CALL check(nf90_open(trim(path)//trim(fname),NF90_NOWRITE,ncid),10)
    CALL check(nf90_inq_dimid(ncid, "longitude", dimid),11)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nx),12)
    CALL check(nf90_inq_dimid(ncid, "latitude", dimid),13)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=ny),14)
    CALL check(nf90_inq_dimid(ncid, "depth", dimid),15)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nz),16)

    allocate( depth(nz) )
    allocate( t(43,43,nz),s(43,43,nz),d(43,43,nz) )
    allocate( pden(43,43,nz),den(43,43,nz) ) 

    CALL check(nf90_inq_varid(ncid,"depth",varid),27)
    CALL check(nf90_get_var(ncid,varid,depth),28)

    CALL check(nf90_close(ncid),35)

    write(*,*) nx,ny,nz

    do i = 1,43
      t(i,:,:) = i-3
    end do
    do i = 1,43
      s(:,i,:) = i-1
    end do
    do i = 1,nz
      d(:,:,i) = depth(i)
    end do

    
    CALL check(nf90_create( "test_dens.nc", NF90_NETCDF4, ncid ), 200)
    CALL check(nf90_def_dim( ncid, 'x', 43, xd ), 201)
    CALL check(nf90_def_dim( ncid, 'y', 43, yd ), 202)
    CALL check(nf90_def_dim( ncid, 'z', nz, zd ), 203)
    CALL check(nf90_def_var( ncid, 'pden1', NF90_DOUBLE, (/xd,yd,zd/), vid(1)), 204) 
    CALL check(nf90_def_var( ncid, 'den1', NF90_DOUBLE, (/xd,yd,zd/), vid(2)), 205) 
    CALL check(nf90_def_var( ncid, 'pden2', NF90_DOUBLE, (/xd,yd,zd/), vid(3)), 206) 
    CALL check(nf90_def_var( ncid, 'den2', NF90_DOUBLE, (/xd,yd,zd/), vid(4)), 207) 
    CALL check(nf90_def_var( ncid, 'x', NF90_DOUBLE, (/xd/), vid(5)), 217)
    CALL check(nf90_put_att( ncid, vid(5), 'standard_name', 'temperature' ),220)
    CALL check(nf90_put_att( ncid, vid(5), 'units', 'degrees_C' ),220)
    CALL check(nf90_put_att( ncid, vid(5), 'axis', 'X' ),220)
    CALL check(nf90_def_var( ncid, 'y', NF90_DOUBLE, (/yd/), vid(6)), 227)
    CALL check(nf90_put_att( ncid, vid(6), 'standard_name', 'salinity' ),220)
    CALL check(nf90_put_att( ncid, vid(6), 'units', '0.001' ),220)
    CALL check(nf90_put_att( ncid, vid(6), 'axis', 'Y' ),220)
    CALL check(nf90_def_var( ncid, 'z', NF90_DOUBLE, (/zd/), vid(7)), 237)
    CALL check(nf90_put_att( ncid, vid(7), 'standard_name', 'depth' ),220)
    CALL check(nf90_put_att( ncid, vid(7), 'units', 'm' ),220)
    CALL check(nf90_put_att( ncid, vid(7), 'axis', 'Z' ),220)
    CALL check(nf90_put_att( ncid, vid(7), 'positive', 'down' ),220)
    CALL check( nf90_enddef( ncid ), 208 )
    CALL dens(43,43,nz,t,s,d,pden,den)
    CALL check( nf90_put_var( ncid, vid(1), pden ), 209 )
    CALL check( nf90_put_var( ncid, vid(2), den ), 210 )
    CALL dens_POP(43,43,nz,t,s,d,pden,den)
    CALL check( nf90_put_var( ncid, vid(3), pden ), 211 )
    CALL check( nf90_put_var( ncid, vid(4), den ), 212 )
    CALL check( nf90_put_var( ncid, vid(5), t(:,1,1) ), 212 )
    CALL check( nf90_put_var( ncid, vid(6), s(1,:,1) ), 212 )
    CALL check( nf90_put_var( ncid, vid(7), d(1,1,:) ), 212 )
    CALL check(nf90_close( ncid ), 213 )
     

    deallocate( depth )
    deallocate( t,s,d ) 
    deallocate( pden,den )

    write(*,*)  0.059808*(exp(-0.025*100) - 1.0)     &
            + 0.100766*100 + 2.28405e-7*100**2   
 
END PROGRAM dens_test

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


SUBROUTINE dens(nx,ny,nz,temp,salt,depth,pden,den)
    implicit none
    integer, intent(in) :: nx,ny,nz
    real(kind=8), intent(in) :: temp(nx,ny,nz),salt(nx,ny,nz),depth(nx,ny,nz)
    real(kind=8), intent(out) :: pden(nx,ny,nz),den(nx,ny,nz)

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
          pden(i,j,k)=den1
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
