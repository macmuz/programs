PROGRAM ROMStoPOP
    USE NETCDF
    USE omp_lib
    implicit none
    integer :: nxR,nyR,nzR,nxP,nyP,nzP,id(2),nt,tmpdim(2)
    integer :: ncid,dimid,varid,varid2,reclen,i,t
    integer :: ncidR,ncidP,tdimid,varids(7),natt
    integer :: xdimid,ydimid,zdimid,ndimid,dimids3(3)
    integer, allocatable :: idx(:,:,:),tmp2(:,:)
    real(kind=8), allocatable :: lonR(:,:),latR(:,:),w(:,:,:)
    real(kind=8), allocatable :: lonP(:,:),latP(:,:),time(:),output3d(:,:,:,:)
    real(kind=8), allocatable :: input2d(:,:),input3d(:,:,:),z_t(:)
    real(kind=8), allocatable :: h(:,:),Cs_r(:),s_r(:),dpR(:,:,:),tmp1(:,:,:)
    real(kind=8) :: hc
    character(len=250) :: ROMSgrid,POPgrid,ROMSfile
    character(len=100) :: num(2,2),idxfile,wfile,output,timeatt,calendar
    character(len=50) :: att
    logical :: ex 
    logical, allocatable :: mask1(:,:),mask2(:,:)

    !SET STACKSIZE EQUAL 1024 MB per thread
    CALL KMP_SET_STACKSIZE_S(3221225472)
    !END SET

    101 format(a,'x',a,"to",a,'x',a,"_idx.bin")
    102 format(a,'x',a,"to",a,'x',a,"_w.bin")

    ROMSgrid = '/users/work/mmuzyka/CSDIR/input_wider/ROMS_grid_05NM_wider_filter.nc'
    POPgrid = '/users/work/mmuzyka/programs/ROMS_ini/PLGNG001.pop.h.2016-07-01-03600.nc'

    call get_command_argument(1,ROMSfile)
    write(*,*) trim(ROMSfile)
    id(1) = INDEX(ROMSfile, '/', BACK=.true.)
    id(2) = INDEX(ROMSfile, '.nc')

    !READ ROMSgrid
    CALL check(nf90_open(trim(ROMSfile),NF90_NOWRITE,ncid),10)
    CALL check(nf90_inq_dimid(ncid, "xi_rho", dimid),11)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nxR),12)
    CALL check(nf90_inq_dimid(ncid, "eta_rho", dimid),13)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nyR),14)
    CALL check(nf90_inq_dimid(ncid, "s_rho", dimid),15)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nzR),16)
    CALL check(nf90_inq_dimid(ncid, "ocean_time", dimid),15)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nt),16)
    allocate(lonR(nxR,nyR),latR(nxR,nyR),time(nt),mask1(nxR,nyR))
    allocate(h(nxR,nyR),CS_r(nzR),s_r(nzR),dpR(nxR,nyR,nzR))
    allocate(input2d(nxR,nyR),input3d(nxR,nyR,nzR))
    CALL check(nf90_inq_varid(ncid,"lon_rho",varid),17)
    CALL check(nf90_get_var(ncid,varid,lonR),18)
    CALL check(nf90_inq_varid(ncid,"lat_rho",varid),19)
    CALL check(nf90_get_var(ncid,varid,latR),20)
    CALL check(nf90_inq_varid(ncid,"h",varid),19)
    CALL check(nf90_get_var(ncid,varid,h),20)
    CALL check(nf90_inq_varid(ncid,"hc",varid),19)
    CALL check(nf90_get_var(ncid,varid,hc),20)
    CALL check(nf90_inq_varid(ncid,"Cs_r",varid),19)
    CALL check(nf90_get_var(ncid,varid,Cs_r),20)
    CALL check(nf90_inq_varid(ncid,"s_rho",varid),19)
    CALL check(nf90_get_var(ncid,varid,s_r),20)
    CALL check(nf90_inq_varid(ncid,"ocean_time",varid),19)
    CALL check(nf90_get_var(ncid,varid,time),20)
    CALL check(nf90_get_att(ncid,varid,"units",timeatt),20)
    CALL check(nf90_get_att(ncid,varid,"calendar",calendar),20)
    CALL check(nf90_inq_varid(ncid,"mask_rho",varid),19)
    CALL check(nf90_get_var(ncid,varid,input2d),20)
    CALL check(nf90_close(ncid),35)

    !READ POPgrid
    CALL check(nf90_open(trim(POPgrid),NF90_NOWRITE,ncid),110)
    CALL check(nf90_inq_dimid(ncid, "nlon", dimid),111)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nxP),112)
    CALL check(nf90_inq_dimid(ncid, "nlat", dimid),113)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nyP),114)
    CALL check(nf90_inq_dimid(ncid, "z_t", dimid),115)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nzP),116)
    allocate(lonP(nxP,nyP),latP(nxP,nyP),idx(nxP,nyP,2),output3d(nxP,nyP,nzP,nt))
    allocate(mask2(nxP,nyP),tmp2(nxP,nyP),w(nxP,nyP,4),z_t(nzP))
    allocate(tmp1(nxR,nyR,nzP))
    CALL check(nf90_inq_varid(ncid,"TLONG",varid),117)
    CALL check(nf90_get_var(ncid,varid,lonP),118)
    CALL check(nf90_inq_varid(ncid,"TLAT",varid),119)
    CALL check(nf90_get_var(ncid,varid,latP),120)
    CALL check(nf90_inq_varid(ncid,"KMT",varid),119)
    CALL check(nf90_get_var(ncid,varid,tmp2),120)
    CALL check(nf90_inq_varid(ncid,"z_t",varid),119)
    CALL check(nf90_get_var(ncid,varid,z_t),120)
    CALL check(nf90_close(ncid),135)

    z_t = z_t*0.01

    mask1 = .false.
    where(input2d.gt.0) mask1 = .true.
    mask2 = .false.
    where(tmp2.gt.0) mask2 = .true.

    write(num(1,1),"(i4)") nxR 
    write(num(1,2),"(i4)") nyR 
    write(num(2,1),"(i4)") nxP 
    write(num(2,2),"(i4)") nyP 
    write(idxfile,101) trim(ADJUSTL(num(1,1))),trim(ADJUSTL(num(1,2))),&
        trim(ADJUSTL(num(2,1))),trim(ADJUSTL(num(2,2)))
    write(wfile,102) trim(ADJUSTL(num(1,1))),trim(ADJUSTL(num(1,2))),&
        trim(ADJUSTL(num(2,1))),trim(ADJUSTL(num(2,2)))

    write(*,*) trim(idxfile)
    write(*,*) trim(wfile)
    write(output,*) ROMSfile(id(1)+1:id(2)-1)//'_'//trim(ADJUSTL(num(2,1)))//&
            'x'//trim(ADJUSTL(num(2,2)))//'_zeta.nc'
    output = adjustl(output)
    write(*,*) trim(output)

    INQUIRE(FILE=trim(idxfile),EXIST=ex)
    inquire(iolength = reclen) idx
    if (ex) then
      OPEN(10,FILE=trim(idxfile),ACCESS='DIRECT', recl=reclen, FORM='UNFORMATTED')
      read(10, rec=1) idx
      CLOSE(10)
    else
      CALL findCorner(nxR,nyR,nxP,nyP,lonR,latR,lonP,latP,mask2,idx)
      OPEN(10,FILE=trim(idxfile),ACCESS='DIRECT', recl=reclen, FORM='UNFORMATTED')
      write(10, rec=1) idx
      CLOSE(10) 
    end if

    CALL findcoef(nxR,nyR,nxP,nyP,lonR,latR,lonP,latP,idx,w)
   
    CALL check(nf90_create( "idx.nc", NF90_CLOBBER, ncid ), 200)
    CALL check(nf90_def_dim( ncid, 'x', nxP, xdimid ), 201)
    CALL check(nf90_def_dim( ncid, 'y', nyP, ydimid ), 202)
    CALL check(nf90_def_dim( ncid, 'z', 2, zdimid ), 202)
    CALL check(nf90_def_dim( ncid, 'n', 4, ndimid ), 202)
    dimids3 = (/ xdimid, ydimid, zdimid /)
    CALL check(nf90_def_var( ncid, 'idx', NF90_INT, dimids3, varid), 205)
    CALL check(nf90_def_var( ncid, 'w', NF90_DOUBLE, &
        (/xdimid,ydimid,ndimid/), varid2), 205)
    CALL check( nf90_enddef( ncid ), 207 )
    CALL check( nf90_put_var( ncid, varid, idx ), 209 )
    CALL check( nf90_put_var( ncid, varid2, w ), 209 )
    CALL check(nf90_close( ncid ), 210 )

    CALL check(nf90_open(trim(POPgrid),NF90_NOWRITE,ncidP),10)
    CALL check(nf90_create( trim(output), NF90_NETCDF4, ncid ), 200)
    CALL check(nf90_def_dim( ncid, 'nlon', nxP, xdimid ), 201)
    CALL check(nf90_def_dim( ncid, 'nlat', nyP, ydimid ), 202)
    CALL check(nf90_def_dim( ncid, 'z_t', nzP, zdimid ), 202)
    CALL check(nf90_def_dim( ncid, 'time', NF90_UNLIMITED, tdimid ), 202)

!    CALL check(nf90_def_dim( ncid, 'dim1', nxR, tmpdim(1) ), 202)
!    CALL check(nf90_def_dim( ncid, 'dim2', nyR, tmpdim(2) ), 202)

    CALL check(nf90_def_var( ncid, 'time', NF90_DOUBLE, (/tdimid/), varids(1)), 205)
    CALL check(nf90_put_att( ncid, varids(1), "long_name", "time" ), 205)
    CALL check(nf90_put_att( ncid, varids(1), "units", trim(timeatt) ), 205)
    CALL check(nf90_put_att( ncid, varids(1), "calendar", trim(calendar) ), 205)

    CALL check(nf90_def_var( ncid, 'z_t', NF90_FLOAT, (/zdimid/), varids(2)), 205)
    CALL check(nf90_inq_varid(ncidP,"z_t",varid),119)
    CALL check(NF90_INQUIRE_VARIABLE(ncidP,varid,natts = natt),119)
    do i = 1, natt
      CALL check(nf90_inq_attname(ncidP, varid, i, att),119)
      CALL check(nf90_copy_att(ncidP, varid, trim(att), ncid, varids(2)), 205)
    end do
    
    CALL check(nf90_def_var( ncid, 'TLONG', NF90_DOUBLE,&
        (/xdimid,ydimid/), varids(3)), 205)
    CALL check(nf90_inq_varid(ncidP,"TLONG",varid),119)
    CALL check(NF90_INQUIRE_VARIABLE(ncidP,varid,natts = natt),119)
    do i = 1, natt
      CALL check(nf90_inq_attname(ncidP, varid, i, att),119)
      CALL check(nf90_copy_att(ncidP, varid, trim(att), ncid, varids(3)), 205)
    end do

    CALL check(nf90_def_var( ncid, 'TLAT', NF90_DOUBLE,&
        (/xdimid,ydimid/), varids(4)), 205)
    CALL check(nf90_inq_varid(ncidP,"TLAT",varid),119)
    CALL check(NF90_INQUIRE_VARIABLE(ncidP,varid,natts = natt),119)
    do i = 1, natt
      CALL check(nf90_inq_attname(ncidP, varid, i, att),119)
      CALL check(nf90_copy_att(ncidP, varid, trim(att), ncid, varids(4)), 205)
    end do

    CALL check(nf90_close( ncidP ), 210 )

    CALL check(nf90_def_var( ncid, 'TEMP', NF90_FLOAT,&
        (/xdimid,ydimid,zdimid,tdimid/), varids(5)), 205)
    CALL check(nf90_put_att( ncid, varids(5), "long_name",&
        "Potential Temperature" ), 205)
    CALL check(nf90_put_att( ncid, varids(5), "units", "degC" ), 205)
    CALL check(nf90_put_att( ncid, varids(5), "coordinates",&
        "TLONG TLAT z_t time" ), 205)

    CALL check(nf90_def_var( ncid, 'SALT', NF90_FLOAT,&
        (/xdimid,ydimid,zdimid,tdimid/), varids(6)), 205)
    CALL check(nf90_put_att( ncid, varids(6), "long_name",&
        "Salinity" ), 205)
    CALL check(nf90_put_att( ncid, varids(6), "units",&
        "gram/kilogram" ), 205)
    CALL check(nf90_put_att( ncid, varids(6), "coordinates",&
        "TLONG TLAT z_t time" ), 205)

!    CALL check(nf90_def_var( ncid, 'tmp', NF90_FLOAT,&
!        (/tmpdim(1),tmpdim(2),zdimid,tdimid/), varids(7)), 205)

    CALL check( nf90_enddef( ncid ), 207 )

    CALL check( nf90_put_var( ncid, varids(1), time ), 209 )
    CALL check( nf90_put_var( ncid, varids(2), z_t ), 209 )
    CALL check( nf90_put_var( ncid, varids(3), lonP ), 209 )
    CALL check( nf90_put_var( ncid, varids(4), latP ), 209 )


    CALL check(nf90_open(trim(ROMSfile),NF90_NOWRITE,ncidR),10)
    do t = 1,nt
      CALL check(nf90_inq_varid(ncidR,"zeta",varid),17)
      CALL check(nf90_get_var(ncidR,varid,input2d,start=(/1,1,t/)),18)

      dpR = 0
      CALL calc_dp(nxR,nyR,nzR,input2d,h,hc,s_r,Cs_r,mask1,dpR)
 
      CALL check(nf90_inq_varid(ncidR,"temp",varid),17)
      CALL check(nf90_get_var(ncidR,varid,input3d,start=(/1,1,1,t/)),18)
      CALL vertical(nxR,nyR,nzR,nzP,z_t,dpR,mask1,input3d,tmp1)
!      CALL check( nf90_put_var( ncid, varids(7), tmp1,&
!            start=(/1,1,1,t/), count=(/nxR,nyR,nzP,1/) ), 209 )

      output3d = 0
      CALL interp(nxR,nyR,nxP,nyP,nzP,tmp2,idx,w,tmp1,output3d) 
      CALL check( nf90_put_var( ncid, varids(5), output3d,&
            start=(/1,1,1,t/), count=(/nxP,nyP,nzP,1/) ), 209 )

      CALL check(nf90_inq_varid(ncidR,"salt",varid),17)
      CALL check(nf90_get_var(ncidR,varid,input3d,start=(/1,1,1,t/)),18) 
      CALL vertical(nxR,nyR,nzR,nzP,z_t,dpR,mask1,input3d,tmp1)

      output3d = 0 
      CALL interp(nxR,nyR,nxP,nyP,nzP,tmp2,idx,w,tmp1,output3d) 
      CALL check( nf90_put_var( ncid, varids(6), output3d,&
            start=(/1,1,1,t/), count=(/nxP,nyP,nzP,1/) ), 209 )

    end do
    CALL check(nf90_close(ncidR),35)

    CALL check(nf90_close( ncid ), 210 )

!    do i = 1, nzP
!      write(*,*) z_t(i)*0.01
!    end do
!    do i = 1, nzR
!      write(*,*) dpR(700,200,i)
!    end do


    deallocate(lonR,latR,lonP,latP,idx,output3d)
    deallocate(mask2,tmp2,w,h,Cs_r,s_r,dpR,tmp1)
    deallocate(input2d,input3d,time,z_t,mask1)
    
END PROGRAM ROMStoPOP

SUBROUTINE check(status,label)
    USE NETCDF
    implicit none
    INTEGER, INTENT(in) :: status , label

    IF (status .ne. nf90_noerr) then
        print*, trim(nf90_strerror(status)),' label=',label
        stop "Stopped"
    END IF
END SUBROUTINE check

SUBROUTINE interp(nxin,nyin,nxout,nyout,nz,kmt,idx,w,input,output)
    implicit none
    integer, intent(in) :: nxin,nyin,nxout,nyout,nz
    integer, intent(in) :: kmt(nxout,nyout),idx(nxout,nyout,2)
    real(kind=8), intent(in) :: w(nxout,nyout,4),input(nxin,nyin,nz)
    real(kind=8), intent(out) :: output(nxout,nyout,nz)
    
    integer :: i,j,k,x,y

    do i = 1, nxout
    do j = 1, nyout
    do k = 1, kmt(i,j)
      x = idx(i,j,1)
      y = idx(i,j,2)
      output(i,j,k) = input(x,y,k)*w(i,j,1)+&
                      input(x+1,y,k)*w(i,j,2)+&
                      input(x+1,y+1,k)*w(i,j,3)+&
                      input(x,y+1,k)*w(i,j,4)
    end do
    end do
    end do

END SUBROUTINE interp

SUBROUTINE calc_dp(nx,ny,nz,zeta,h,hc,s_r,Cs_r,mask,dp)
    implicit none
    integer, intent(in) :: nx,ny,nz
    real(kind=8), intent(in) :: zeta(nx,ny),h(nx,ny),hc,s_r(nz),Cs_r(nz)
    logical, intent(in) :: mask(nx,ny)
    real(kind=8), intent(out) :: dp(nx,ny,nz)

    integer :: i,j,k

    do i = 1,nx 
    do j = 1,ny 
    if (mask(i,j)) then
      do k = 1,nz
        dp(i,j,k) = -1*(zeta(i,j)+(zeta(i,j)+h(i,j))*(hc*s_r(k)+h(i,j)*Cs_r(k))/&
            (hc+h(i,j)))
      end do 
    end if
    end do 
    end do 
 
END SUBROUTINE calc_dp

SUBROUTINE vertical(nx,ny,nzin,nzout,z_t,dp,mask,input,output)
    implicit none
    integer, intent(in) :: nx,ny,nzin,nzout
    real(kind=8), intent(in) :: z_t(nzout),dp(nx,ny,nzin),input(nx,ny,nzin)
    logical, intent(in) :: mask(nx,ny) 
    real(kind=8), intent(out) :: output(nx,ny,nzout)
    
    integer :: i,j,k
    real(kind=8) :: splint,y2(nzin)
   
!$OMP PARALLEL DO DEFAULT(PRIVATE),SHARED(dp,input,output),&
!$OMP& FIRSTPRIVATE(nx,ny,nzin,nzout,mask,z_t) SCHEDULE(DYNAMIC)
    do i = 1, nx
    do j = 1, ny
    if (mask(i,j)) then
      output(i,j,:) = 0
      CALL spline(dp(i,j,:),input(i,j,:),nzin,y2)
      do k = 1, nzout
        if (z_t(k).gt.dp(i,j,nzin) .and. z_t(k).lt.dp(i,j,1)) then
          output(i,j,k) = splint(dp(i,j,:),input(i,j,:),y2,nzin,z_t(k))
!          write(*,*) z_t(k),dp(i,j,nzin),dp(i,j,1)
        elseif(z_t(k).ge.dp(i,j,1)) then
          output(i,j,k) = input(i,j,1)
        else
          output(i,j,k) = input(i,j,nzin)
        end if
      end do
    end if
    end do
    end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO DEFAULT(PRIVATE),SHARED(output),&
!$OMP& FIRSTPRIVATE(nx,ny,nzout,mask) SCHEDULE(DYNAMIC)
    do k = 1, nzout
      CALL extrap(output(:,:,k),mask,nx,ny,100,.true.)
    end do
!$OMP END PARALLEL DO

END SUBROUTINE vertical

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

SUBROUTINE findCorner(nx1,ny1,nx2,ny2,lon1,lat1,lon2,lat2,mask2,idx)
    implicit none
    integer, intent(in) :: nx1,ny1,nx2,ny2
    real(kind=8), intent(in) :: lon1(nx1,ny1),lat1(nx1,ny1)
    real(kind=8), intent(in) :: lon2(nx2,ny2),lat2(nx2,ny2)
    integer, intent(out) :: idx(nx2,ny2,2)
    logical, intent(in) :: mask2(nx2,ny2)

    integer :: i,j,k,l,locTmp(2),x,y
    real(kind=8) :: dis_array(nx1,ny1),pt(2),cor(4,2)
    logical :: score

    idx = 0    

!$OMP PARALLEL DO DEFAULT(PRIVATE),SHARED(idx),&
!$OMP& FIRSTPRIVATE(nx1,ny1,nx2,ny2,lon1,lat1,lon2,lat2,mask2) SCHEDULE(DYNAMIC)
    do i = 1, nx2
    do j = 1, ny2
    if (mask2(i,j)) then
      pt = (/lon2(i,j),lat2(i,j)/)
      dis_array(:,:) = sqrt((pt(1)-lon1(:,:))**2+(pt(2)-lat1(:,:))**2)
      locTmp = minloc(dis_array)
      x = locTmp(1)
      y = locTmp(2)

      if (x.gt.1 .and. y.gt.1) then
        k = x-1
        l = y-1
        cor(1,:) = (/lon1(k,l),lat1(k,l)/)
        cor(2,:) = (/lon1(k,l+1),lat1(k,l+1)/)
        cor(3,:) = (/lon1(k+1,l+1),lat1(k+1,l+1)/)
        cor(4,:) = (/lon1(k+1,l),lat1(k+1,l)/)
        CALL in_convex_polygon(cor,pt,score)
        if (score) then
          idx(i,j,:) = (/k,l/)
          cycle
        end if 
      end if

      if (x.gt.1 .and. y.gt.0 .and. y.lt.ny1) then
        k = x-1
        l = y
        cor(1,:) = (/lon1(k,l),lat1(k,l)/)
        cor(2,:) = (/lon1(k,l+1),lat1(k,l+1)/)
        cor(3,:) = (/lon1(k+1,l+1),lat1(k+1,l+1)/)
        cor(4,:) = (/lon1(k+1,l),lat1(k+1,l)/)
        CALL in_convex_polygon(cor,pt,score)
        if (score) then
          idx(i,j,:) = (/k,l/)
          cycle
        end if
      end if

      if (x.gt.0 .and. x.lt.nx1 .and. y.gt.1) then
        k = x
        l = y-1
        cor(1,:) = (/lon1(k,l),lat1(k,l)/)
        cor(2,:) = (/lon1(k,l+1),lat1(k,l+1)/)
        cor(3,:) = (/lon1(k+1,l+1),lat1(k+1,l+1)/)
        cor(4,:) = (/lon1(k+1,l),lat1(k+1,l)/)
        CALL in_convex_polygon(cor,pt,score)
        if (score) then
          idx(i,j,:) = (/k,l/)
          cycle
        end if
      end if

      if (x.gt.0 .and. x.lt.nx1 .and. y.gt.0 .and. y.lt.ny1) then
        k = x
        l = y
        cor(1,:) = (/lon1(k,l),lat1(k,l)/)
        cor(2,:) = (/lon1(k,l+1),lat1(k,l+1)/)
        cor(3,:) = (/lon1(k+1,l+1),lat1(k+1,l+1)/)
        cor(4,:) = (/lon1(k+1,l),lat1(k+1,l)/)
        CALL in_convex_polygon(cor,pt,score)
        if (score) then
          idx(i,j,:) = (/k,l/)
          cycle
        end if
      end if
    
    end if
    end do
    write(*,*) 'i=',i
    end do
!$OMP END PARALLEL DO

END SUBROUTINE findCorner

SUBROUTINE findcoef(nx1,ny1,nx2,ny2,lon1,lat1,lon2,lat2,idx,w)
    implicit none
    integer, intent(in) :: nx1,ny1,nx2,ny2
    real(kind=8), intent(in) :: lon1(nx1,ny1),lat1(nx1,ny1)
    real(kind=8), intent(in) :: lon2(nx2,ny2),lat2(nx2,ny2)
    integer, intent(in) :: idx(nx2,ny2,2)
    real(kind=8), intent(out) :: w(nx2,ny2,4)

    integer :: i,j,k,l,iter
    real(kind=8) :: p(2),v(4,2),a(2),b(2),c(2),d(2)
    real(kind=8) :: x,y,dx(2),tol,f(2),Df(2,2),q(2)
    w = 0
!$OMP PARALLEL DO DEFAULT(PRIVATE),SHARED(w),&
!$OMP& FIRSTPRIVATE(nx1,ny1,nx2,ny2,lon1,lat1,lon2,lat2,idx) SCHEDULE(DYNAMIC)
    do i = 1, nx2
    do j = 1, ny2
    if (idx(i,j,1).gt.0) then
      k = idx(i,j,1)
      l = idx(i,j,2)
      p = (/lon2(i,j),lat2(i,j)/)
      v(1,:) = (/lon1(k,l),lat1(k,l)/)
      v(2,:) = (/lon1(k+1,l),lat1(k+1,l)/)
      v(3,:) = (/lon1(k+1,l+1),lat1(k+1,l+1)/)
      v(4,:) = (/lon1(k,l+1),lat1(k,l+1)/)
      a = v(1,:)-p
      b = v(2,:)-v(1,:)
      c = v(4,:)-v(1,:)
      d = v(1,:)-v(2,:)-v(4,:)+v(3,:)
      
      x = 0.5 
      y = 0.5
      dx = (/0.1,0.1/)
      iter = 0
      tol = 1.0e-12 
      do
        f = a+b*x+c*y+d*x*y
        Df(:,1)=(b + d*y)
        Df(:,2)=(c + d*x)
        CALL cramer(Df,dx,-1*f)
        x=x+dx(1)
        y=y+dx(2)
        if ( sqrt(dx(1)**2+dx(2)**2).gt.10.0 ) iter=20
        if ( sqrt(dx(1)**2+dx(2)**2).le.tol .or. iter.ge.20) EXIT
      end do
      if (iter.lt.20) then
        w(i,j,1) = (1-x)*(1-y) 
        w(i,j,2) = x*(1-y)
        w(i,j,3) = x*y
        w(i,j,4) = (1-x)*y

        q=w(i,j,1)*v(1,:)+w(i,j,2)*v(2,:)+w(i,j,3)*v(3,:)+w(i,j,4)*v(4,:)
!        w(i,j,1:2) = q
!        write(*,*) p-q
      end if
    end if
    end do
    end do 
!$OMP END PARALLEL DO
END SUBROUTINE findcoef

SUBROUTINE cramer(A,x,b)
    implicit none
    real(kind=8), intent(in) :: A(2,2),b(2)
    real(kind=8), intent(out) :: x(2)

    real(kind=8) :: W,Wx,Wy

    W = A(1,1)*A(2,2)-A(2,1)*A(2,2)
    Wx = b(1)*A(2,2)-b(2)*A(1,2)
    Wy = A(1,1)*b(2)-A(2,1)*b(1)

    if (W.ne.0.0) then
      x(1) = Wx/W
      x(2) = Wy/W
    else
      write(*,*) "nieoznaczony"
    end if
    
END SUBROUTINE cramer

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

SUBROUTINE extrap(a,mask,lon,lat,maxscn,fill)
    implicit none
    integer, intent(in) :: lon,lat,maxscn
    real(kind=8), intent(inout) :: a(lon,lat)
    logical, intent(in) :: mask(lon,lat)
    logical, intent(in) :: fill

    integer :: i,j,n,cnt,overall
    real(kind=8) :: relc,ave
    real(kind=8), dimension(lon,lat) :: sor,res
    logical :: mask_tmp(lon,lat),mask_tmp2(lon,lat)

    relc = 0.6
    sor = 0.0
    where(.not.mask) sor=relc

    if (fill) then
      mask_tmp2 = mask
      mask_tmp = mask_tmp2
      do

        overall = 0
        do i = 1, lon
        do j = 1, lat

          if (.not.mask_tmp(i,j)) then

            cnt = 0
            ave = 0.0

            if ( i.gt.1 ) then
              if ( mask_tmp(i-1,j) ) then
                ave = ave+a(i-1,j)
                cnt = cnt+1
              end if
            end if

            if ( j.gt.1 ) then
              if ( mask_tmp(i,j-1) ) then
                ave = ave+a(i,j-1)
                cnt = cnt+1
              end if
            end if

            if ( i.lt.lon ) then
              if ( mask_tmp(i+1,j) ) then
                ave = ave+a(i+1,j)
                cnt = cnt+1
              end if
            end if

            if ( j.lt.lat ) then
              if ( mask_tmp(i,j+1) ) then
                ave = ave+a(i,j+1)
                cnt = cnt+1
              end if
            end if

            if ( cnt.gt.0 ) then
              a(i,j) = ave/(real(cnt,8))
              overall = overall+cnt
              mask_tmp2(i,j) = .true.
            end if

          end if

        end do
        end do

        mask_tmp = mask_tmp2
        if ( overall.eq.0 ) EXIT

      end do
      
    end if

    do n=1,maxscn
      do i=2,lon-1
        do j=2,lat-1
          res(i,j)=0.25*(a(i-1,j)+a(i+1,j)+a(i,j-1)+a(i,j+1))-a(i,j)
        end do
      end do

      do i=2,lon-1
        res(i,1)=0.3333*(a(i-1,1)+a(i+1,1)+a(i,2))-a(i,1)
        res(i,lat)=0.3333*(a(i-1,lat)+a(i+1,lat)+a(i,lat-1))-a(i,lat)
      end do
      do j=2,lat-1
        res(1,j)=0.3333*(a(1,j-1)+a(1,j+1)+a(2,j))-a(1,j)
        res(lon,j)=0.3333*(a(lon,j-1)+a(lon,j+1)+a(lon-1,j))-a(lon,j)
      end do

      res(1,1)=0.5*(a(1,2)+a(2,1))-a(1,1)
      res(lon,1)=0.5*(a(lon,2)+a(lon-1,1))-a(lon,1)
      res(1,lat)=0.5*(a(1,lat-1)+a(2,lat))-a(1,lat)
      res(lon,lat)=0.5*(a(lon,lat-1)+a(lon-1,lat))-a(lon,lat)

      res=res*sor
      a=a+res
    end do
END SUBROUTINE extrap

    
