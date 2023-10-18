PROGRAM create1992
    USE NETCDF
    implicit none
    integer, parameter :: nx=600,ny=640,nz=66
    integer :: i,ncid1,ncid2,vid
    integer :: int2d(nx,ny)
    real(kind=4) :: tmp2d(nx,ny),tmp3d(nx,ny,nz),fv
    character(len=50) :: template,output,input
    character(len=200) :: command
    logical :: tmask(nx,ny,nz),umask(nx,ny,nz)

    101 format('cp ',a,' ',a)
    
    tmask = .true.
    umask = .true.

    template = 'PLGNG003.pop.h.2016-11-20-03600.nc'
    output = 'PLGNG003.pop.h.1992-07-01-03600.nc'
    input = 'RS001.pop.hs.RS001.1992-07-01_00.nc'

    write(command,101) trim(template),trim(output)
    write(*,*) trim(command)
    CALL SYSTEM(trim(command))

    CALL check(nf90_open(trim(input),NF90_NOWRITE,ncid1),10)
    CALL check(nf90_open(trim(output),NF90_WRITE,ncid2),11)

    CALL check(nf90_inq_varid(ncid2,"KMT",vid),17)
    CALL check(nf90_get_var(ncid2,vid,int2d),18)
    do i = 1,nz
      where(int2d.lt.i) tmask(:,:,i) = .false.
    end do
   
    CALL check(nf90_inq_varid(ncid2,"KMU",vid),17)
    CALL check(nf90_get_var(ncid2,vid,int2d),18)
    do i = 1,nz
      where(int2d.lt.i) umask(:,:,i) = .false.
    end do
   
    CALL check(nf90_inq_varid(ncid1,"SHGT",vid),17)
    CALL check(nf90_get_var(ncid1,vid,tmp2d),18)

    CALL check(nf90_inq_varid(ncid2,"SSH",vid),17) 
    CALL check(nf90_get_att(ncid2,vid,"_FillValue",fv),17)
    where(tmask(:,:,1).eq..false.) tmp2d = fv 
    do i = 1,24
      CALL check(nf90_put_var(ncid2,vid,tmp2d,start = (/1,1,i/),&
        count=(/nx,ny,1/)),244)
    end do
 
    CALL check(nf90_inq_varid(ncid1,"UBTROP",vid),17)
    CALL check(nf90_get_var(ncid1,vid,tmp2d),18)

    CALL check(nf90_inq_varid(ncid2,"SU",vid),17) 
    CALL check(nf90_get_att(ncid2,vid,"_FillValue",fv),17) 
    where(umask(:,:,1).eq..false.) tmp2d = fv 
    do i = 1,24
      CALL check(nf90_put_var(ncid2,vid,tmp2d,start = (/1,1,i/),&
        count=(/nx,ny,1/)),244)
    end do
 
    CALL check(nf90_inq_varid(ncid1,"VBTROP",vid),17)
    CALL check(nf90_get_var(ncid1,vid,tmp2d),18)

    CALL check(nf90_inq_varid(ncid2,"SV",vid),17) 
    CALL check(nf90_get_att(ncid2,vid,"_FillValue",fv),17) 
    where(umask(:,:,1).eq..false.) tmp2d = fv 
    do i = 1,24
      CALL check(nf90_put_var(ncid2,vid,tmp2d,start = (/1,1,i/),&
        count=(/nx,ny,1/)),244)
    end do
 
    CALL check(nf90_inq_varid(ncid1,"UVEL",vid),17)
    CALL check(nf90_get_var(ncid1,vid,tmp3d),18)

    CALL check(nf90_inq_varid(ncid2,"UVEL",vid),17) 
    CALL check(nf90_get_att(ncid2,vid,"_FillValue",fv),17) 
    where(umask.eq..false.) tmp3d = fv 
    do i = 1,24
      CALL check(nf90_put_var(ncid2,vid,tmp3d,start = (/1,1,1,i/),&
        count=(/nx,ny,nz,1/)),244)
    end do
 
    CALL check(nf90_inq_varid(ncid1,"VVEL",vid),17)
    CALL check(nf90_get_var(ncid1,vid,tmp3d),18)

    CALL check(nf90_inq_varid(ncid2,"VVEL",vid),17) 
    CALL check(nf90_get_att(ncid2,vid,"_FillValue",fv),17) 
    where(umask.eq..false.) tmp3d = fv 
    do i = 1,24
      CALL check(nf90_put_var(ncid2,vid,tmp3d,start = (/1,1,1,i/),&
        count=(/nx,ny,nz,1/)),244)
    end do
 
    CALL check(nf90_inq_varid(ncid1,"TEMP",vid),17)
    CALL check(nf90_get_var(ncid1,vid,tmp3d),18)

    CALL check(nf90_inq_varid(ncid2,"TEMP",vid),17) 
    CALL check(nf90_get_att(ncid2,vid,"_FillValue",fv),17) 
    where(tmask.eq..false.) tmp3d = fv 
    do i = 1,24
      CALL check(nf90_put_var(ncid2,vid,tmp3d,start = (/1,1,1,i/),&
        count=(/nx,ny,nz,1/)),244)
    end do
 
    CALL check(nf90_inq_varid(ncid1,"SALT",vid),17)
    CALL check(nf90_get_var(ncid1,vid,tmp3d),18)

    CALL check(nf90_inq_varid(ncid2,"SALT",vid),17) 
    CALL check(nf90_get_att(ncid2,vid,"_FillValue",fv),17) 
    where(tmask.eq..false.) tmp3d = fv 
    do i = 1,24
      CALL check(nf90_put_var(ncid2,vid,tmp3d,start = (/1,1,1,i/),&
        count=(/nx,ny,nz,1/)),244)
    end do
 
    CALL check(nf90_close(ncid2),35)
    CALL check(nf90_close(ncid1),36)
    

END PROGRAM create1992

SUBROUTINE check(status,label)
    USE NETCDF
    implicit none
    INTEGER, INTENT(in) :: status , label

    IF (status .ne. nf90_noerr) then
        print*, trim(nf90_strerror(status)),' label= ',label
        stop "Stopped"
    END IF
END SUBROUTINE check
