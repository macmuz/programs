PROGRAM move_forcing
    USE netcdf
    USE omp_lib
    implicit none
    integer :: ncid, vid, reclen, varid(3), ncid2
    integer :: x, y, i, j, k, l, narg, ntime, varid2(3)
    character(len=250) :: igrid, ogrid
    integer :: idx(560,600,4,2)
    integer :: x_dimid, y_dimid, z_dimid, t_dimid
    integer :: dimids3(3), dimids4(4)
    real(kind=8) :: point(2), corners(4,2), time
    real(kind=8) :: ilon(600,640), ilat(600,640), input(600,640)
    real(kind=8) :: olon(560,600), olat(560,600), output(560,600)
    real(kind=8) :: dis_array(600,640)
    real(kind=8) :: W(560,600,4)
    character(len=20) :: file1, file2, varname, timename, units
    character(len=20) :: varname2
    character(len=250) :: inputfile
    character(len=50) :: long_name(3), outputfile
    logical :: ex1, ex2, fexit, score, mask(560,600)

    file1 = "idx_forcing.bin"
    file2 = "W_forcing.bin"

    igrid = '/users/work/mmuzyka/CSDIR/metro_files/baltic_gridv4.nc'
    ogrid = 'ROMS_grid_2_3km_560x600_NetCDF4.nc'

!############# command line arguments #################################
    narg = command_argument_count()
    if (narg.lt.3 .or. narg.gt.4) then
        write(*,*) "Program must have 3 (4) arguments: &
                    &inputfile, timename, varname, (varname2 optional)"
        stop
    end if
    call get_command_argument(1,inputfile)
    call get_command_argument(2,timename)
    call get_command_argument(3,varname)
    if (narg.eq.4) call get_command_argument(4,varname2)

    i = INDEX(inputfile, '/', BACK=.true.)
    outputfile = inputfile(i+1:)
    i = INDEX(outputfile, '.nc', BACK=.true.)
    outputfile(i:) = '_560x600.nc'
!######################################################################

    CALL check(nf90_open(trim(igrid),NF90_NOWRITE,ncid),1)
    CALL check(nf90_inq_varid(ncid,"lon_rho",vid),7)
    CALL check(nf90_get_var(ncid,vid,ilon),8)
    CALL check(nf90_inq_varid(ncid,"lat_rho",vid),7)
    CALL check(nf90_get_var(ncid,vid,ilat),8)
    CALL check(nf90_close(ncid),9)

    CALL check(nf90_open(trim(ogrid),NF90_NOWRITE,ncid),1)
    CALL check(nf90_inq_varid(ncid,"lon_rho",vid),7)
    CALL check(nf90_get_var(ncid,vid,olon),8)
    CALL check(nf90_inq_varid(ncid,"lat_rho",vid),7)
    CALL check(nf90_get_var(ncid,vid,olat),8)
    CALL check(nf90_close(ncid),9)
    
!INTERP COEF
    INQUIRE(FILE=trim(file1),EXIST=ex1)
    INQUIRE(FILE=trim(file2),EXIST=ex2)
    if (ex1 .and. ex2) then

      inquire(iolength = reclen) idx

      OPEN(10,FILE=trim(file1),ACCESS='DIRECT', recl=reclen, FORM='UNFORMATTED')
      read(10, rec=1) idx
      CLOSE(10)

      inquire(iolength = reclen) W

      OPEN(10,FILE=trim(file2),ACCESS='DIRECT', recl=reclen, FORM='UNFORMATTED')
      read(10, rec=1) W
      CLOSE(10)

    else

    idx = 0
    W = 0
!$OMP PARALLEL DO DEFAULT(PRIVATE),SHARED(ilon,ilat,olon,olat),&
!$OMP& SHARED(idx,W), SCHEDULE(DYNAMIC)
    do j = 1, 600
    do i = 1, 560

        dis_array(:,:) = sqrt( (olon(i,j)-ilon(:,:))**2+(olat(i,j)-ilat(:,:))**2 )
        idx(i,j,1,:) = minloc( dis_array )
        x = idx(i,j,1,1)
        y = idx(i,j,1,2)
        if (x.gt.1 .and. x.lt.600 .and. y.gt.1 .and. y.lt.640) then
          point(1) = olon(i,j)
          point(2) = olat(i,j)
          fexit = .false.
          do k = x-1,x
            if (fexit) EXIT
          do l = y-1,y
            if (fexit) EXIT
            corners(1,:) = (/ilon(k,l),ilat(k,l)/)
            corners(2,:) = (/ilon(k,l+1),ilat(k,l+1)/)
            corners(3,:) = (/ilon(k+1,l+1),ilat(k+1,l+1)/)
            corners(4,:) = (/ilon(k+1,l),ilat(k+1,l)/)
            CALL in_convex_polygon(corners,point,score)
            if (score) then
              fexit = .true.
              idx(i,j,1,:) = (/k,l/)
              idx(i,j,2,:) = (/k,l+1/)
              idx(i,j,3,:) = (/k+1,l+1/)
              idx(i,j,4,:) = (/k+1,l/)
              call calc_w2( ilon(k:k+1,l:l+1), ilat(k:k+1,l:l+1), point, W(i,j,:) )
            end if
          end do
          end do
        else
          idx(i,j,1,:) = 0
        end if

    end do
    end do
!$OMP END PARALLEL DO
      inquire(iolength = reclen) idx

      OPEN(10,FILE=trim(file1),ACCESS='DIRECT', recl=reclen, FORM='UNFORMATTED')
      write(10, rec=1) idx
      CLOSE(10)

      inquire(iolength = reclen) W

      OPEN(10,FILE=trim(file2),ACCESS='DIRECT', recl=reclen, FORM='UNFORMATTED')
      write(10, rec=1) W
      CLOSE(10)

      CALL check(nf90_create( "forcing_interp_coef.nc", NF90_NETCDF4, ncid ), 200)
      CALL check(nf90_def_dim( ncid, 'x', 560, x_dimid ), 201)
      CALL check(nf90_def_dim( ncid, 'y', 600, y_dimid ), 202)
      CALL check(nf90_def_dim( ncid, 'corner', 4, z_dimid ), 202)
      CALL check(nf90_def_dim( ncid, 'index', 2, t_dimid ), 202)
      dimids3 = (/ x_dimid, y_dimid, z_dimid /)
      dimids4 = (/ x_dimid, y_dimid, z_dimid, t_dimid /)
      CALL check(nf90_def_var( ncid, 'idx', NF90_INT, dimids4, varid(1)), 205)
      CALL check(nf90_def_var( ncid, 'w', NF90_DOUBLE, dimids3, varid(2)), 205)
      CALL check( nf90_enddef( ncid ), 207 )
      CALL check( nf90_put_var( ncid, varid(1), idx ), 209 )
      CALL check( nf90_put_var( ncid, varid(2), W ), 209 )
      CALL check(nf90_close( ncid ), 210 )

    end if !files exist
!END INTERP COEF

    mask = .true.
    where(idx(:,:,1,1).eq.0) mask = .false.

    CALL check(nf90_open(trim(inputfile),NF90_NOWRITE,ncid),1)
    CALL check(nf90_inq_dimid(ncid, trim(timename), t_dimid),2)
    CALL check(nf90_inquire_dimension(ncid, t_dimid, len = ntime),3)
    CALL check(nf90_inq_varid(ncid, trim(timename), varid(1)),4)
    CALL check(nf90_get_att(ncid, varid(1), 'long_name', long_name(1)),5)
    CALL check(nf90_inq_varid(ncid, trim(varname), varid(2)),6)
    CALL check(nf90_get_att(ncid, varid(2), 'long_name', long_name(2)),7)
    CALL check(nf90_get_att(ncid, varid(2), 'units', units),8)
    if (narg.eq.4) then
    CALL check(nf90_inq_varid(ncid, trim(varname2), varid(3)),6)
    CALL check(nf90_get_att(ncid, varid(3), 'long_name', long_name(3)),7)
    end if

    CALL check(nf90_create( trim(outputfile), NF90_NETCDF4, ncid2 ), 200)
    CALL check(nf90_def_dim( ncid2, 'xi_rho', 560, x_dimid ), 201)
    CALL check(nf90_def_dim( ncid2, 'eta_rho', 600, y_dimid ), 202)
    CALL check(nf90_def_dim( ncid2, trim(timename), NF90_UNLIMITED, t_dimid ), 202)
    dimids3 = (/ x_dimid, y_dimid, t_dimid /)
    CALL check(nf90_def_var( ncid2, trim(timename), NF90_DOUBLE, (/t_dimid/), varid2(1)), 205)
    CALL check(nf90_put_att(ncid2, varid2(1), 'long_name', long_name(1)), 205)
    CALL check(nf90_put_att(ncid2, varid2(1), 'units', &
        'days since 1968-05-23 00:00:00 GMT'), 205)
    CALL check(nf90_put_att(ncid2, varid2(1), 'calendar', 'gregorian'), 205)
    CALL check(nf90_def_var( ncid2, trim(varname), NF90_REAL, dimids3, varid2(2)), 205)
    CALL check(nf90_put_att(ncid2, varid2(2), 'long_name', long_name(2)), 205)
    CALL check(nf90_put_att(ncid2, varid2(2), 'units', units), 205)
    CALL check(nf90_put_att(ncid2, varid2(2), 'time', trim(timename)), 205)
    if (narg.eq.4) then
    CALL check(nf90_def_var( ncid2, trim(varname2), NF90_REAL, dimids3, varid2(3)), 205)
    CALL check(nf90_put_att(ncid2, varid2(3), 'long_name', long_name(3)), 205)
    CALL check(nf90_put_att(ncid2, varid2(3), 'units', units), 205)
    CALL check(nf90_put_att(ncid2, varid2(3), 'time', trim(timename)), 205)
    end if
    CALL check( nf90_enddef( ncid2 ), 207 )

    do i = 1,ntime
      CALL check( nf90_get_var(ncid, varid(1), time, start = (/ i /)) ,300)
      CALL check( nf90_get_var(ncid, varid(2), input, start = (/ 1,1,i /)) ,300)

      CALL calcme(600,640,560,600,input,output,W,idx)
      CALL extrap(output,mask,560,600,100,2)

      CALL check( nf90_put_var( ncid2, varid2(1), time, start = (/ i /)), 209 )
      CALL check( nf90_put_var( ncid2, varid2(2), real(output,4), &
        start = (/ 1,1,i /) ), 209 )
      if (narg.eq.4) then
      CALL check( nf90_get_var(ncid, varid(3), input, start = (/ 1,1,i /)) ,300)
      CALL calcme(600,640,560,600,input,output,W,idx)
      CALL extrap(output,mask,560,600,100,2)
      CALL check( nf90_put_var( ncid2, varid2(3), real(output,4), &
        start = (/ 1,1,i /) ), 209 )
      end if
      write(*,*) 'i=',i,trim(varname)
    end do

    CALL check(nf90_close( ncid2 ), 210 )

    CALL check(nf90_close( ncid ), 210 )

END PROGRAM move_forcing

SUBROUTINE check(status,label)
    USE NETCDF
    implicit none
    INTEGER, INTENT(in) :: status , label

    IF (status .ne. nf90_noerr) then
        print*, trim(nf90_strerror(status)),' label=',label
        stop "Stopped"
    END IF
END SUBROUTINE check

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

SUBROUTINE calc_w2(lon,lat,point_xy,w)
    implicit none
    real(kind=8), intent(in) :: lon(2,2),lat(2,2)
    real(kind=8), intent(in) :: point_xy(2)
    real(kind=8), intent(out) :: w(4)

    integer :: idx(2)
    real(kind=8) :: x, y, dis_array(2,2)

    x = point_xy(1)
    y = point_xy(2)
    dis_array(:,:) = sqrt((x-lon(:,:))**2+(y-lat(:,:))**2)
    idx = minloc( dis_array )

    if ( dis_array(idx(1),idx(2)).eq.0 ) then
        w = 0
        if (x.eq.1 .and. y.eq.1) then
            w(1) = 1
        elseif (x.eq.1 .and. y.eq.2) then
            w(2) = 1
        elseif (x.eq.2 .and. y.eq.2) then
            w(3) = 1
        elseif (x.eq.2 .and. y.eq.1) then
            w(4) = 1
        end if
    else
        dis_array(:,:) = 1/dis_array(:,:)
        w(1) = dis_array(1,1)/sum(dis_array)
        w(2) = dis_array(1,2)/sum(dis_array)
        w(3) = dis_array(2,2)/sum(dis_array)
        w(4) = dis_array(2,1)/sum(dis_array)
    end if


END SUBROUTINE calc_w2

SUBROUTINE calcme(nxin,nyin,nxout,nyout,inarray,outarray,W,idx)
    implicit none
    integer, intent(in) :: nxin,nyin,nxout,nyout
    real(kind=8), intent(in) :: inarray(nxin,nxin),W(nxout,nyout,4)
    real(kind=8), intent(out) :: outarray(nxout,nyout)
    integer, intent(in) :: idx(nxout,nyout,4,2)

    real(kind=8) :: tmp
    integer :: i,j,k,x,y

    do i = 1,nxout
    do j = 1,nyout
      if (idx(i,j,1,1).ne.0) then
        tmp = 0
        do k = 1,4
          x = idx(i,j,k,1)
          y = idx(i,j,k,2)
          if (x.gt.0 .and. y.gt.0) then
            tmp = tmp+inarray(x,y)*W(i,j,k)
          end if
        end do
        outarray(i,j) = tmp
      end if
    end do
    end do
END SUBROUTINE calcme

SUBROUTINE extrap(a,mask,lon,lat,maxscn,met)
    implicit none
    integer, intent(in) :: lon,lat,maxscn,met
    real(kind=8), intent(inout) :: a(lon,lat)
    logical, intent(in) :: mask(lon,lat)

    integer :: i,j,n,cnt,overall
    real(kind=8) :: relc,ave
    real(kind=8), dimension(lon,lat) :: sor,res
    logical :: mask_tmp(lon,lat),mask_tmp2(lon,lat)

    relc=1.0
    sor = 0.0
    where(.not.mask) sor=relc

    select case(met)
    case(0)
      where(.not.mask) a=0.0
    case(1)
      cnt = 0
      ave = 0.0
      do i=1,lon
      do j=1,lat
        if (mask(i,j)) then
            ave=ave+a(i,j)
            cnt=cnt+1
        end if
      end do
      end do
      if ( cnt.GT.0 ) ave = ave/real(cnt,8)
      where(.not.mask) a=ave
    case(2)
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
    end select

    do n=1,maxscn

!$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j),SCHEDULE(DYNAMIC)
        do i=2,lon-1
          do j=2,lat-1
            res(i,j)=0.25*(a(i-1,j)+a(i+1,j)+a(i,j-1)+a(i,j+1))-a(i,j)
          end do
        end do
!$OMP END PARALLEL DO

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
