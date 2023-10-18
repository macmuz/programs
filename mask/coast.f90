PROGRAM coast
    USE NETCDF
    implicit none
    character(len=150) :: filename,outfile
    character(len=15) :: varname,dim1,dim2
    integer :: ncid,varid,dimids2(2)
    integer :: nlons,nlats,x_dimid,y_dimid
    integer :: i,j,cnt
    real(kind=8), allocatable :: bathy(:,:)
    real(kind=8) :: fv

!##################################
    filename = 'bathy2_out.nc'
    outfile = 'coast.nc'
    varname = 'h'
    fv = -999.0
!##################################

    CALL check(nf90_open(trim(filename), NF90_NOWRITE, ncid),100)

    CALL check(nf90_inq_varid(ncid, trim(varname), varid),101)
    CALL check(nf90_inquire_variable(ncid, varid, dimids=dimids2),102)
    CALL check(nf90_inquire_dimension(ncid, dimids2(1), len=nlons, &
                                                     name=dim1),103)
    CALL check(nf90_inquire_dimension(ncid, dimids2(2), len=nlats, &
                                                     name=dim2),104)
    allocate(bathy(nlons,nlats))
    CALL check(nf90_get_var(ncid, varid, bathy),105)


    
    do i = 3, nlons-2
    do j = 3, nlats-2
      if ( bathy(i,j).ne.fv .and. bathy(i-1,j).eq.fv .and. &
           bathy(i,j+1).eq.fv ) then
      if ( bathy(i,j-1).ne.fv .and. bathy(i-1,j-1).eq.fv .and. &
           bathy(i,j-2).eq.fv ) then
      if ( bathy(i+1,j).ne.fv .and. bathy(i+1,j+1).eq.fv .and. &
           bathy(i+2,j).eq.fv ) then
        bathy(i,j) = fv
        bathy(i,j-1) = fv
        bathy(i+1,j) = fv
      end if
      end if
      end if

      if ( bathy(i,j).ne.fv .and. bathy(i+1,j).eq.fv .and. &
           bathy(i,j+1).eq.fv ) then
      if ( bathy(i,j-1).ne.fv .and. bathy(i+1,j-1).eq.fv .and. &
           bathy(i,j-2).eq.fv ) then
      if ( bathy(i-1,j).ne.fv .and. bathy(i-2,j).eq.fv .and. &
           bathy(i-1,j+1).eq.fv ) then
        bathy(i,j) = fv
        bathy(i,j-1) = fv
        bathy(i-1,j) = fv
      end if
      end if
      end if

      if ( bathy(i,j).ne.fv .and. bathy(i+1,j).eq.fv .and. &
           bathy(i,j-1).eq.fv ) then
      if ( bathy(i-1,j).ne.fv .and. bathy(i-2,j).eq.fv .and. &
           bathy(i-1,j-1).eq.fv ) then
      if ( bathy(i,j+1).ne.fv .and. bathy(i,j+2).eq.fv .and. &
           bathy(i+1,j+1).eq.fv ) then
        bathy(i,j) = fv
        bathy(i-1,j) = fv
        bathy(i,j+1) = fv
      end if
      end if
      end if

      if ( bathy(i,j).ne.fv .and. bathy(i-1,j).eq.fv .and. &
           bathy(i,j-1).eq.fv ) then
      if ( bathy(i+1,j).ne.fv .and. bathy(i+2,j).eq.fv .and. &
           bathy(i+1,j-1).eq.fv ) then
      if ( bathy(i,j+1).ne.fv .and. bathy(i,j+2).eq.fv .and. &
           bathy(i-1,j+1).eq.fv ) then
        bathy(i,j) = fv
        bathy(i+1,j) = fv
        bathy(i,j+1) = fv
      end if
      end if
      end if
    end do
    end do


  
    do  
    cnt = 1 
    do i = 2, nlons-1
    do j = 2, nlats-1
      if ( bathy(i,j).ne.fv .and. bathy(i-1,j).eq.fv .and. &
           bathy(i+1,j).eq.fv ) then
        cnt = cnt+1
        bathy(i,j) = fv
      end if
      if ( bathy(i,j).ne.fv .and. bathy(i,j-1).eq.fv .and. &
           bathy(i,j+1).eq.fv ) then
        cnt = cnt+1
        bathy(i,j) = fv
      end if
    end do
    end do
    write(*,*) cnt
      if (cnt.le.1) EXIT
    end do
    
    bathy(865,230) = fv
    bathy(877,221) = fv
    bathy(876,222) = fv
    bathy(880,219) = fv


    CALL check(nf90_create( trim(outfile), NF90_CLOBBER, ncid ), 200)
    CALL check(nf90_def_dim( ncid, trim(dim1), nlons, x_dimid ), 201)
    CALL check(nf90_def_dim( ncid, trim(dim2), nlats, y_dimid ), 202)
    dimids2 = (/ x_dimid, y_dimid /)
    CALL check(nf90_def_var( ncid, trim(varname), NF90_DOUBLE, &
                            dimids2, varid), 203)
    CALL check( nf90_put_att( ncid, varid, '_FillValue', fv), 204)
    CALL check( nf90_enddef( ncid ), 205 )
    CALL check( nf90_put_var( ncid, varid, bathy ), 206)
    CALL check(nf90_close( ncid ), 207)


    deallocate(bathy)

END PROGRAM coast

SUBROUTINE check(status,label)
    USE NETCDF
    implicit none
    INTEGER, INTENT(in) :: status , label

    IF (status .ne. nf90_noerr) then
        print*, trim(nf90_strerror(status)),'label=',label
        stop "Stopped"
    END IF
END SUBROUTINE check
