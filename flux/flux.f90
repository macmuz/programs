PROGRAM flux
    USE netcdf
    implicit none
    integer :: ncid, varid, dimid, nx, ny, nz
    integer :: i,j,k
    integer :: xdimid,ydimid,zdimid,dimids3(3) 
    real(kind=8) :: hc
    real(kind=8), allocatable :: h(:,:),cs_w(:),zeta(:,:),depth(:,:,:) 
    real(kind=8), allocatable :: mask(:,:),s_w(:)
    character(len=200) :: filename

    filename = '/users/work/mmuzyka/CSDIR/new_metro_560x600_02/run/baltic/ocean_avg_2018-08-10.nc'

!READ GRID
    CALL check(nf90_open(trim(filename),NF90_NOWRITE,ncid),310)

    CALL check(nf90_inq_dimid(ncid, "xi_rho", dimid),311)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nx),312)

    CALL check(nf90_inq_dimid(ncid, "eta_rho", dimid),313)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=ny),314)

    CALL check(nf90_inq_dimid(ncid, "s_rho", dimid),315)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nz),316)

    allocate( h(nx,ny), cs_w(nz+1), zeta(nx,ny), depth(nx,ny,nz+1) )
    allocate( mask(nx,ny), s_w(nz+1) )

    CALL check(nf90_inq_varid(ncid,"hc",varid),317)
    CALL check(nf90_get_var(ncid,varid,hc),318)

    CALL check(nf90_inq_varid(ncid,"Cs_w",varid),319)
    CALL check(nf90_get_var(ncid,varid,cs_w),320)

    CALL check(nf90_inq_varid(ncid,"s_w",varid),319)
    CALL check(nf90_get_var(ncid,varid,s_w),320)

    CALL check(nf90_inq_varid(ncid,"h",varid),319)
    CALL check(nf90_get_var(ncid,varid,h),320)

    CALL check(nf90_inq_varid(ncid,"zeta",varid),319)
    CALL check(nf90_get_var(ncid,varid,zeta,start=(/1,1,2/)),320)

    CALL check(nf90_inq_varid(ncid,"mask_rho",varid),319)
    CALL check(nf90_get_var(ncid,varid,mask),320)

    CALL check(nf90_close(ncid),360)
!END READ GRID

    
    where(mask.lt.0.5) zeta = 0.0
    do i = 1, nx   
      do j = 1, ny   
        do k = 1, nz+1   
          depth(i,j,k) = zeta(i,j)+(zeta(i,j)+h(i,j))*(hc*s_w(k)+h(i,j)*cs_w(k))/(hc+h(i,j))
        end do
      end do
    end do

    CALL check( nf90_create( 'depth.nc',NF90_NETCDF4,ncid ), 500 )
    CALL check( nf90_def_dim( ncid, 'nx', nx, xdimid ), 501 )
    CALL check( nf90_def_dim( ncid, 'ny', ny, ydimid ), 502 )
    CALL check( nf90_def_dim( ncid, 'nz', nz+1, zdimid ), 502 )
    dimids3 = (/ xdimid, ydimid, zdimid /)
    CALL check(nf90_def_var( ncid, 'depth', NF90_DOUBLE,&
        dimids3, varid ), 508)
    CALL check( nf90_enddef(ncid), 516 )
    CALL check( nf90_put_var( ncid, varid, depth ), 517 )
    CALL check( nf90_close(ncid), 600 ) 

    deallocate( h, cs_w, zeta, depth, mask, s_w )
END PROGRAM flux

SUBROUTINE check(status,label)
    USE NETCDF
    implicit none
    INTEGER, INTENT(in) :: status , label

    IF (status .ne. nf90_noerr) then
        print*, trim(nf90_strerror(status)),' label= ',label
        stop "Stopped"
    END IF
END SUBROUTINE check
