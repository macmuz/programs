PROGRAM mask
    USE netcdf
    implicit none
    character(len=150) :: filename,outfile
    character(len=15) :: varname,dim1,dim2
    integer :: ncid,varid,dimid,nlons,nlats
    integer :: x_dimid,y_dimid,dimids2(2)
    integer :: i,j,cnt,itmp,jtmp,cnt2,k
    integer, allocatable :: list(:,:)
    real(kind=8) :: fv
    real(kind=8), allocatable :: bathy(:,:)
    logical, allocatable :: processed(:,:),mymask(:,:)

!##################################
    filename = 'coast.nc'
    varname = 'h'
    fv = -999.0
    outfile = 'out.nc'
!##################################

    CALL check(nf90_open(trim(filename), NF90_NOWRITE, ncid),100)
    CALL check(nf90_inq_varid(ncid, trim(varname), varid),101)
    CALL check(nf90_inquire_variable(ncid, varid, dimids=dimids2),102)
    CALL check(nf90_inquire_dimension(ncid, dimids2(1), len=nlons, &
                                                     name=dim1),103)
    CALL check(nf90_inquire_dimension(ncid, dimids2(2), len=nlats, &
                                                     name=dim2),104)
    allocate(bathy(nlons,nlats),processed(nlons,nlats))
    allocate(list(nlons*nlats,2),mymask(nlons,nlats))
    CALL check(nf90_get_var(ncid, varid, bathy),105)
    CALL check(nf90_close(ncid),106)

    mymask=.FALSE.
    where(bathy.GT.0) mymask=.TRUE.
    processed = .FALSE.
    do i = 1,nlons
    write(*,*) 'I=',i
    do j = 1,nlats
      list=0
      if (bathy(i,j).GT.0 .and. .not.processed(i,j)) then
        write(*,*) i,j
        cnt=1
        list(cnt,1)=i
        list(cnt,2)=j
        processed(i,j)=.TRUE.
        cnt2=cnt
        do
          itmp=list(cnt,1)
          jtmp=list(cnt,2)

          if(itmp.gt.1) then
            if (bathy(itmp-1,jtmp).GT.0 .and. &
            .not.processed(itmp-1,jtmp)) then
              cnt2=cnt2+1
              list(cnt2,1)=itmp-1
              list(cnt2,2)=jtmp
              processed(itmp-1,jtmp)=.TRUE.
              write(*,*) list(cnt2,:)
            end if
          end if

          if(itmp.lt.nlons) then
            if (bathy(itmp+1,jtmp).GT.0 .and. &
            .not.processed(itmp+1,jtmp)) then
              cnt2=cnt2+1
              list(cnt2,1)=itmp+1
              list(cnt2,2)=jtmp
              processed(itmp+1,jtmp)=.TRUE.
              write(*,*) list(cnt2,:)
            end if
          end if

          if(jtmp.gt.1) then
            if (bathy(itmp,jtmp-1).GT.0 .and. &
            .not.processed(itmp,jtmp-1)) then
              cnt2=cnt2+1
              list(cnt2,1)=itmp
              list(cnt2,2)=jtmp-1
              processed(itmp,jtmp-1)=.TRUE.
              write(*,*) list(cnt2,:)
            end if
          end if

          if(jtmp.lt.nlats) then
            if (bathy(itmp,jtmp+1).GT.0 .and. &
            .not.processed(itmp,jtmp+1)) then
              cnt2=cnt2+1
              list(cnt2,1)=itmp
              list(cnt2,2)=jtmp+1
              processed(itmp,jtmp+1)=.TRUE.
              write(*,*) list(cnt2,:)
            end if
          end if

          cnt=cnt+1
          if (cnt.GT.cnt2) then
            if (cnt2.LT.100000) then
              do k = 1,cnt2
                itmp=list(k,1)
                jtmp=list(k,2)
                write(*,*) 'removing',itmp,jtmp
                mymask(itmp,jtmp)=.FALSE.
              end do
            end if
            EXIT          
          end if
        end do
      end if    
    end do
    end do

    where(.not.mymask) bathy=fv
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

    deallocate(bathy,processed,list,mymask)
END PROGRAM mask

SUBROUTINE check(status,label)
    USE NETCDF
    implicit none
    INTEGER, INTENT(in) :: status , label

    IF (status .ne. nf90_noerr) then
        print*, trim(nf90_strerror(status)),'label=',label
        stop "Stopped"
    END IF
END SUBROUTINE check
