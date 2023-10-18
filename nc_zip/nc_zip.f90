PROGRAM nc_zip
    USE netcdf
    implicit none
    integer, parameter :: nb=1, deflate=9
    logical, parameter :: shuffle=.true.
    integer :: i, j, k, narg, sidx, ncid_in, ncid_out
    integer :: nDim, nVar, nGAtt, ulimitID, nAtt
    integer :: ndims, xtype, n, xtypezip, shp(4), tmp
    integer, allocatable :: dimids(:), dimlens(:), idx(:)
    integer, allocatable :: dimid(:), varid(:), dimlen(:), zip_varid(:)
    character(len=300) :: infile, outfile
    character(len=50), allocatable :: var(:)
    character(len=50), allocatable :: dimname(:), varname(:)
    character(len=50) :: attname
    logical :: fv, copy

    real(kind=4) :: float_fv, float_sf, float_af, float_min, float_max
    real(kind=8) :: double_fv, double_sf, double_af, double_min, double_max

    integer(kind=nb), allocatable :: zipped(:,:,:,:),zipped_1d(:)
    logical, allocatable :: mask(:,:,:,:),mask_1d(:)
    character(kind=1), allocatable :: data_char(:,:,:,:),data_char_1d(:)
    integer(kind=1), allocatable :: data_byte(:,:,:,:),data_byte_1d(:)
    integer(kind=2), allocatable :: data_short(:,:,:,:),data_short_1d(:)
    integer(kind=4), allocatable :: data_int(:,:,:,:),data_int_1d(:)
    real(kind=4), allocatable :: data_float(:,:,:,:),data_float_1d(:)
    real(kind=8), allocatable :: data_double(:,:,:,:),data_double_1d(:)

    character(kind=1) :: scalar_char
    integer(kind=1) :: scalar_byte
    integer(kind=2) :: scalar_short
    integer(kind=4) :: scalar_int
    real(kind=4) :: scalar_float
    real(kind=8) :: scalar_double
   
    select case(nb)
      CASE(1)
        xtypezip = NF90_BYTE
      CASE(2)
        xtypezip = NF90_SHORT
      CASE(4)
        xtypezip = NF90_INT
    end select
 
    narg = command_argument_count()
    allocate( var(narg-1), zip_varid(narg-1) )

    call get_command_argument(1,infile)
    do i = 1, narg-1
      call get_command_argument(i+1,var(i))
    end do    

    sidx = index(infile,'.nc')
    write(outfile,'(A,A)') trim(infile(1:sidx-1)),'_zipped.nc'


    call check( nf90_open(trim(infile),NF90_NOWRITE,ncid_in), 1 )
    call check( nf90_inquire(ncid_in, nDim, nVar, nGAtt, ulimitID), 2 )

    allocate( dimid(nDim), dimlen(nDim), varid(nVar) )
    allocate( dimname(nDim), varname(nVar) )

    do i = 1, narg-1
      call check( nf90_inq_varid( ncid_in, trim(var(i)), zip_varid(i)), 3 )
    end do

    do i = 1, nDim
      call check( nf90_inquire_dimension(ncid_in, i, name=dimname(i), len=dimlen(i)), 3 )
    end do

    do i = 1, nVar
      call check( nf90_inquire_variable(ncid_in, i, name=varname(i) ), 3 )
    end do    


    call check( nf90_create(trim(outfile), NF90_NETCDF4, ncid_out), 300 )

    do i = 1, nDim
      n = dimlen(i)
      if ( i.eq.ulimitID ) n=NF90_UNLIMITED
      call check( nf90_def_dim(ncid_out, trim(dimname(i)), n, dimid(i) ), 301 )
    end do

    do i = 1, nGAtt
      call check( nf90_inq_attname(ncid_in, NF90_GLOBAL, i, attname), 3 )
      call check( nf90_copy_att(ncid_in, NF90_GLOBAL, trim(attname), ncid_out, NF90_GLOBAL), 302 )
    end do
    
    call check( nf90_enddef(ncid_out), 101 )

    do i = 1, nVar
    
!      write(*,*) trim(varname(i))

      copy = .false.

      call check( nf90_inquire_variable(ncid_in, i, ndims=ndims), 103 )
      allocate( dimids(ndims), dimlens(ndims), idx(ndims) )
      call check( nf90_inquire_variable(ncid_in, i, xtype=xtype, dimids=dimids, nAtts=nAtt), 103 )
      do j = 1, ndims
        dimids(j) = dimid(dimids(j))
        dimlens(j) = dimlen(dimids(j))
      end do
      do j = 1, nAtt
        call check( nf90_inq_attname(ncid_in, i, j, attname), 3 )
        if ( trim(attname).eq.'add_offset' .or. trim(attname).eq.'scale_factor' ) then
          write(*,*) 'Variable: ',trim(varname(i)), ' is already packed. Omitting.'
          copy = .true. 
          EXIT
        end if
      end do
      if ( ANY(zip_varid.eq.i) ) then
      if (xtype.eq.NF90_BYTE) then
        write(*,*) 'Cannot pack "BYTE" variable. Omitting'
        copy = .true. 
      end if
      if (xtype.eq.NF90_CHAR) then
        write(*,*) 'Cannot pack "CHARACTER" variable. Omitting'
        copy = .true.
      end if
      if (xtype.eq.NF90_SHORT) then
        write(*,*) 'Cannot pack "SHORT" variable. Omitting'
        copy = .true.
      end if
      if (xtype.eq.NF90_INT) then
        write(*,*) 'Cannot pack "INT" variable. Omitting'
        copy = .true.
      end if
      if (ndims.eq.0) then
        write(*,*) trim(varname(i)), ' is a scalar. Omitting'
        copy = .true.
      end if
      end if

      shp(:) = 1
      if ( ndims.gt.4 ) then
        shp(:) = dimlens(1:4)
      else
        shp(1:ndims) = dimlens(:)
      end if

      if ( ALL(zip_varid.ne.i) .or. copy ) then

        call check( nf90_redef(ncid_out), 101 )
        if (ndims.eq.0) then
          call check( nf90_def_var(ncid_out, trim(varname(i)), xtype, &
            varid(i) ), 4 )
        else
          call check( nf90_def_var(ncid_out, trim(varname(i)), xtype, &
            dimids, varid(i), shuffle=shuffle, deflate_level = deflate), 4 )
        end if
        do j = 1, nAtt
          call check( nf90_inq_attname(ncid_in, i, j, attname), 3 )
          call check( nf90_copy_att(ncid_in, i, trim(attname), ncid_out, varid(i)), 302 )
        end do
        call check( nf90_enddef(ncid_out), 101 )

        select case (xtype)
          CASE (NF90_BYTE)
            if (ndims.eq.0) then
              call check( nf90_get_var(ncid_in, i, scalar_byte), 99 )
              call check( nf90_put_var(ncid_out, varid(i), scalar_byte), 100 )
            else
            allocate( data_byte(shp(1),shp(2),shp(3),shp(4)) )

            if ( ndims.gt.4 ) then
              idx(:) = 1
              do j = 1, product(dimlens(5:ndims))
                do k = 4,ndims-1
                  if (idx(k).gt.dimlens(k)) then
                    idx(k+1) = idx(k+1)+1
                    idx(k) = 1
                  end if
                end do

                call check( nf90_get_var(ncid_in, i, data_byte, start = idx, &
                count = (/dimlens(1), dimlens(2), dimlens(3), dimlens(4),&
                (k**0, k=1,ndims-4)/) ), 105 )

                call check( nf90_put_var(ncid_out, varid(i), data_byte, start = idx, &
                count = (/dimlens(1), dimlens(2), dimlens(3), dimlens(4),&
                (k**0, k=1,ndims-4)/) ), 205 )

                idx(4) = idx(4)+dimlens(4)
              end do
            else
              call check( nf90_get_var(ncid_in, i, &
                data_byte(1:shp(1),1:shp(2),1:shp(3),1:shp(4))), 105 )
              call check( nf90_put_var(ncid_out, varid(i), &
                data_byte(1:shp(1),1:shp(2),1:shp(3),1:shp(4))), 205 )
            end if

            deallocate( data_byte )
            end if
          CASE (NF90_CHAR)
            if (ndims.eq.0) then
              call check( nf90_get_var(ncid_in, i, scalar_char), 99 )
              call check( nf90_put_var(ncid_out, varid(i), scalar_char), 100 )
            else
            allocate( data_char(shp(1),shp(2),shp(3),shp(4)) )

            if ( ndims.gt.4 ) then
              idx(:) = 1
              do j = 1, product(dimlens(5:ndims))
                do k = 4,ndims-1
                  if (idx(k).gt.dimlens(k)) then
                    idx(k+1) = idx(k+1)+1
                    idx(k) = 1
                  end if
                end do

                call check( nf90_get_var(ncid_in, i, data_char, start = idx, &
                count = (/dimlens(1), dimlens(2), dimlens(3), dimlens(4),&
                (k**0, k=1,ndims-4)/) ), 105 )

                call check( nf90_put_var(ncid_out, varid(i), data_char, start = idx, &
                count = (/dimlens(1), dimlens(2), dimlens(3), dimlens(4),&
                (k**0, k=1,ndims-4)/) ), 205 )

                idx(4) = idx(4)+dimlens(4)
              end do
            else
              call check( nf90_get_var(ncid_in, i, &
                data_char(1:shp(1),1:shp(2),1:shp(3),1:shp(4))), 105 )
              call check( nf90_put_var(ncid_out, varid(i), &
                data_char(1:shp(1),1:shp(2),1:shp(3),1:shp(4))), 205 )
            end if

            deallocate( data_char )
            end if
          CASE (NF90_SHORT)
            if (ndims.eq.0) then
              call check( nf90_get_var(ncid_in, i, scalar_short), 99 )
              call check( nf90_put_var(ncid_out, varid(i), scalar_short), 100 )
            else
            allocate( data_short(shp(1),shp(2),shp(3),shp(4)) )

            if ( ndims.gt.4 ) then
              idx(:) = 1
              do j = 1, product(dimlens(5:ndims))
                do k = 4,ndims-1
                  if (idx(k).gt.dimlens(k)) then
                    idx(k+1) = idx(k+1)+1
                    idx(k) = 1
                  end if
                end do

                call check( nf90_get_var(ncid_in, i, data_short, start = idx, &
                count = (/dimlens(1), dimlens(2), dimlens(3), dimlens(4),&
                (k**0, k=1,ndims-4)/) ), 105 )

                call check( nf90_put_var(ncid_out, varid(i), data_short, start = idx, &
                count = (/dimlens(1), dimlens(2), dimlens(3), dimlens(4),&
                (k**0, k=1,ndims-4)/) ), 205 )

                idx(4) = idx(4)+dimlens(4)
              end do
            else
              call check( nf90_get_var(ncid_in, i, &
                data_short(1:shp(1),1:shp(2),1:shp(3),1:shp(4))), 105 )
              call check( nf90_put_var(ncid_out, varid(i), &
                data_short(1:shp(1),1:shp(2),1:shp(3),1:shp(4))), 205 )
            end if

            deallocate( data_short )
            end if
          CASE (NF90_INT)
            if (ndims.eq.0) then
              call check( nf90_get_var(ncid_in, i, scalar_int), 99 )
              call check( nf90_put_var(ncid_out, varid(i), scalar_int), 100 )
            else
            allocate( data_int(shp(1),shp(2),shp(3),shp(4)) )

            if ( ndims.gt.4 ) then
              idx(:) = 1
              do j = 1, product(dimlens(5:ndims))
                do k = 4,ndims-1
                  if (idx(k).gt.dimlens(k)) then
                    idx(k+1) = idx(k+1)+1
                    idx(k) = 1
                  end if
                end do

                call check( nf90_get_var(ncid_in, i, data_int, start = idx, &
                count = (/dimlens(1), dimlens(2), dimlens(3), dimlens(4),&
                (k**0, k=1,ndims-4)/) ), 105 )

                call check( nf90_put_var(ncid_out, varid(i), data_int, start = idx, &
                count = (/dimlens(1), dimlens(2), dimlens(3), dimlens(4),&
                (k**0, k=1,ndims-4)/) ), 205 )

                idx(4) = idx(4)+dimlens(4)
              end do
            else
              call check( nf90_get_var(ncid_in, i, &
                data_int(1:shp(1),1:shp(2),1:shp(3),1:shp(4))), 105 )
              call check( nf90_put_var(ncid_out, varid(i), &
                data_int(1:shp(1),1:shp(2),1:shp(3),1:shp(4))), 205 )
            end if

            deallocate( data_int )
            end if
          CASE (NF90_FLOAT)
            if (ndims.eq.0) then
              call check( nf90_get_var(ncid_in, i, scalar_float), 99 )
              call check( nf90_put_var(ncid_out, varid(i), scalar_float), 100 )
            else
            allocate( data_float(shp(1),shp(2),shp(3),shp(4)) )

            if ( ndims.gt.4 ) then
              idx(:) = 1
              do j = 1, product(dimlens(5:ndims))
                do k = 4,ndims-1
                  if (idx(k).gt.dimlens(k)) then
                    idx(k+1) = idx(k+1)+1
                    idx(k) = 1
                  end if
                end do

                call check( nf90_get_var(ncid_in, i, data_float, start = idx, &
                count = (/dimlens(1), dimlens(2), dimlens(3), dimlens(4),&
                (k**0, k=1,ndims-4)/) ), 105 )

                call check( nf90_put_var(ncid_out, varid(i), data_float, start = idx, &
                count = (/dimlens(1), dimlens(2), dimlens(3), dimlens(4),&
                (k**0, k=1,ndims-4)/) ), 205 )

                idx(4) = idx(4)+dimlens(4)
              end do
            else
              call check( nf90_get_var(ncid_in, i, &
                data_float(1:shp(1),1:shp(2),1:shp(3),1:shp(4))), 105 )
              call check( nf90_put_var(ncid_out, varid(i), &
                data_float(1:shp(1),1:shp(2),1:shp(3),1:shp(4))), 205 )
            end if

            deallocate( data_float )
            end if
          CASE (NF90_DOUBLE)
            if (ndims.eq.0) then
              call check( nf90_get_var(ncid_in, i, scalar_double), 99 )
              call check( nf90_put_var(ncid_out, varid(i), scalar_double), 100 )
            else
            allocate( data_double(shp(1),shp(2),shp(3),shp(4)) )

            if ( ndims.gt.4 ) then
              idx(:) = 1
              do j = 1, product(dimlens(5:ndims))
                do k = 4,ndims-1
                  if (idx(k).gt.dimlens(k)) then
                    idx(k+1) = idx(k+1)+1
                    idx(k) = 1
                  end if
                end do

                call check( nf90_get_var(ncid_in, i, data_double, start = idx, &
                count = (/dimlens(1), dimlens(2), dimlens(3), dimlens(4),&
                (k**0, k=1,ndims-4)/) ), 105 )

                call check( nf90_put_var(ncid_out, varid(i), data_double, start = idx, &
                count = (/dimlens(1), dimlens(2), dimlens(3), dimlens(4),&
                (k**0, k=1,ndims-4)/) ), 205 )

                idx(4) = idx(4)+dimlens(4)
              end do
            else
              call check( nf90_get_var(ncid_in, i, &
                data_double(1:shp(1),1:shp(2),1:shp(3),1:shp(4))), 105 )
              call check( nf90_put_var(ncid_out, varid(i), &
                data_double(1:shp(1),1:shp(2),1:shp(3),1:shp(4))), 205 )
            end if

            deallocate( data_double )
            end if
        end select

      else

        fv = .false.
 
        select case (xtype)
          CASE (NF90_FLOAT)
  !##########################################################          
            if ( ndims.gt.4 ) then
                allocate( data_float_1d(product(dimlens)), zipped_1d(product(dimlens)) )
                allocate( mask_1d(product(dimlens)) )
                mask_1d(:) = .true.
            else
                allocate( mask(shp(1),shp(2),shp(3),shp(4)) )
                mask(:,:,:,:) = .true.
            end if
            allocate( data_float(shp(1),shp(2),shp(3),shp(4)) )
            allocate( zipped(shp(1),shp(2),shp(3),shp(4)) )

            do j = 1, nAtt
              call check( nf90_inq_attname(ncid_in, i, j, attname), 3 )
              if ( trim(attname).eq.'missing_value' ) then
                fv = .true.
                call check( nf90_get_att(ncid_in, i, trim(attname), float_fv), 4 )
              end if
              if ( trim(attname).eq.'_FillValue' ) then
                fv = .true.
                call check( nf90_get_att(ncid_in, i, trim(attname), float_fv), 4 )
              end if
            end do

            call check( nf90_redef(ncid_out), 101 )
            call check( nf90_def_var(ncid_out, trim(varname(i)), &
                xtypezip, dimids, varid(i), shuffle=shuffle, &
                deflate_level = deflate), 4 )
            do j = 1, nAtt
              call check( nf90_inq_attname(ncid_in, i, j, attname), 3 )
              if ( trim(attname).ne.'missing_value' .and. &
                  trim(attname).ne.'_FillValue') then
                call check( nf90_copy_att(ncid_in, i, trim(attname), &
                    ncid_out, varid(i)), 302 )
              end if
            end do
            if (fv) then
              call check( nf90_put_att(ncid_out, varid(i), '_FillValue', &
                int(-2**(8*nb-1),nb) ), 37 )
            end if

            if ( ndims.gt.4 ) then
              tmp = product(dimlens(1:4))
              idx(:) = 1
              do j = 1, product(dimlens(5:ndims))
                do k = 4,ndims-1
                  if (idx(k).gt.dimlens(k)) then
                    idx(k+1) = idx(k+1)+1
                    idx(k) = 1
                  end if
                end do

                call check( nf90_get_var(ncid_in, i, data_float, start = idx, &
                count = (/dimlens(1), dimlens(2), dimlens(3), dimlens(4), &
                (k**0, k=1,ndims-4)/) ), 105 )

                data_float_1d(1+(j-1)*tmp:j*tmp) = reshape(data_float,(/tmp/))                    

                idx(4) = idx(4)+dimlens(4)
              end do
            else
              call check( nf90_get_var(ncid_in, i, &
                data_float(1:shp(1),1:shp(2),1:shp(3),1:shp(4))), 105 )
            end if
  
120         continue         
            if (fv) then
              if ( ndims.gt.4 ) then 
                if (float_fv.eq.float_fv) then
                  where(data_float_1d.eq.float_fv) mask_1d = .false.
                else
                  where(data_float_1d.ne.data_float_1d) mask_1d = .false.
                end if
                zipped_1d = int(-2**(8*nb-1),nb)
              else 
                if (float_fv.eq.float_fv) then
                  where(data_float.eq.float_fv) mask = .false.
                else
                  where(data_float.ne.data_float) mask = .false.
                end if 
                zipped = int(-2**(8*nb-1),nb)
              end if
            end if

            if ( ndims.gt.4 ) then
              float_min = minval(data_float_1d,mask=mask_1d)
              float_max = maxval(data_float_1d,mask=mask_1d)
            else
              float_min = minval(data_float,mask=mask)
              float_max = maxval(data_float,mask=mask)
            end if

            if ( float_min.lt.-1e30 ) then
                fv = .true.
                float_fv = float_min
                call check( nf90_put_att(ncid_out, varid(i), '_FillValue', &
                int(-2**(8*nb-1),nb) ), 37 )
                goto 120
            end if
            if ( float_max.gt.1e30 ) then
                fv = .true.
                float_fv = float_max
                call check( nf90_put_att(ncid_out, varid(i), '_FillValue', &
                int(-2**(8*nb-1),nb) ), 37 )
                goto 120
            end if

            if (fv) then
               float_sf = (float_max-float_min)/(2**(8*nb)-2)
               float_af = float_min+(2**(8*nb-1)-1)*float_sf 
            else
               float_sf = (float_max-float_min)/(2**(8*nb)-1) 
               float_af = float_min+(2**(8*nb-1))*float_sf 
            end if

            write(*,*) trim(varname(i))
            write(*,*) 'scale_factor:',float_sf
            write(*,*) 'add_offset:',float_af
            write(*,*) 'min:',nint((float_min-float_af)/float_sf,nb)
            write(*,*) 'max:',nint((float_max-float_af)/float_sf,nb)

            call check( nf90_put_att(ncid_out, varid(i), 'add_offset', float_af ), 38 ) 
            call check( nf90_put_att(ncid_out, varid(i), 'scale_factor', float_sf ), 39 ) 
            call check( nf90_enddef(ncid_out), 101 )
            

            if ( ndims.gt.4 ) then
              where( mask_1d ) zipped_1d = nint((data_float_1d-float_af)/float_sf,nb) 
              tmp = product(dimlens(1:4))
              idx(:) = 1
              do j = 1, product(dimlens(5:ndims))
                do k = 4,ndims-1
                  if (idx(k).gt.dimlens(k)) then
                    idx(k+1) = idx(k+1)+1
                    idx(k) = 1
                  end if
                end do

                zipped = reshape(zipped_1d(1+(j-1)*tmp:j*tmp),&
                    (/dimlens(1),dimlens(2),dimlens(3),dimlens(4)/))

                call check( nf90_put_var(ncid_out, varid(i), zipped, start = idx, &
                count = (/dimlens(1), dimlens(2), dimlens(3), dimlens(4), &
                  (k**0, k=1,ndims-4)/) ), 105 )

                idx(4) = idx(4)+dimlens(4)
              end do
            else
              where( mask ) zipped = nint((data_float-float_af)/float_sf,nb) 
              call check( nf90_put_var(ncid_out, varid(i), &
                zipped(1:shp(1),1:shp(2),1:shp(3),1:shp(4))), 105 )
            end if

            deallocate( data_float, zipped )
            if ( ndims.gt.4 ) then
                deallocate( data_float_1d, zipped_1d, mask_1d )
            else
                deallocate( mask )
            end if
  !##########################################################          
          CASE (NF90_DOUBLE)
!##########################################################
            if ( ndims.gt.4 ) then
                allocate( data_double_1d(product(dimlens)), zipped_1d(product(dimlens)) )
                allocate( mask_1d(product(dimlens)) )
                mask_1d(:) = .true.
            else
                allocate( mask(shp(1),shp(2),shp(3),shp(4)) )
                mask(:,:,:,:) = .true.
            end if
            allocate( data_double(shp(1),shp(2),shp(3),shp(4)) )
            allocate( zipped(shp(1),shp(2),shp(3),shp(4)) )

            do j = 1, nAtt
              call check( nf90_inq_attname(ncid_in, i, j, attname), 3 )
              if ( trim(attname).eq.'missing_value' ) then
                fv = .true.
                call check( nf90_get_att(ncid_in, i, trim(attname), double_fv), 4 )
              end if
              if ( trim(attname).eq.'_FillValue' ) then
                fv = .true.
                call check( nf90_get_att(ncid_in, i, trim(attname), double_fv), 4 )
              end if
            end do

            call check( nf90_redef(ncid_out), 101 )
            call check( nf90_def_var(ncid_out, trim(varname(i)), &
                xtypezip, dimids, varid(i), shuffle=shuffle, &
                deflate_level = deflate), 4 )
            do j = 1, nAtt
              call check( nf90_inq_attname(ncid_in, i, j, attname), 3 )
              if ( trim(attname).ne.'missing_value' .and. &
                  trim(attname).ne.'_FillValue') then
                call check( nf90_copy_att(ncid_in, i, trim(attname), &
                    ncid_out, varid(i)), 302 )
              end if
            end do
            if (fv) then
              call check( nf90_put_att(ncid_out, varid(i), '_FillValue', &
                int(-2**(8*nb-1),nb) ), 37 )
            end if

            if ( ndims.gt.4 ) then
              tmp = product(dimlens(1:4))
              idx(:) = 1
              do j = 1, product(dimlens(5:ndims))
                do k = 4,ndims-1
                  if (idx(k).gt.dimlens(k)) then
                    idx(k+1) = idx(k+1)+1
                    idx(k) = 1
                  end if
                end do

                call check( nf90_get_var(ncid_in, i, data_double, start = idx, &
                count = (/dimlens(1), dimlens(2), dimlens(3), dimlens(4), &
                (k**0, k=1,ndims-4)/) ), 105 )

                data_double_1d(1+(j-1)*tmp:j*tmp) = reshape(data_double,(/tmp/))

                idx(4) = idx(4)+dimlens(4)
              end do
            else
              call check( nf90_get_var(ncid_in, i, &
                data_double(1:shp(1),1:shp(2),1:shp(3),1:shp(4))), 105 )
            end if

220         continue
            if (fv) then
              if ( ndims.gt.4 ) then
                if (double_fv.eq.double_fv) then
                  where(data_double_1d.eq.double_fv) mask_1d = .false.
                else
                  where(data_double_1d.ne.data_double_1d) mask_1d = .false.
                end if
                zipped_1d = int(-2**(8*nb-1),nb)
              else
                if (double_fv.eq.double_fv) then
                  where(data_double.eq.double_fv) mask = .false.
                else
                  where(data_double.ne.data_double) mask = .false.
                end if
                zipped = int(-2**(8*nb-1),nb)
              end if
            end if

            if ( ndims.gt.4 ) then
              double_min = minval(data_double_1d,mask=mask_1d)
              double_max = maxval(data_double_1d,mask=mask_1d)
            else
              double_min = minval(data_double,mask=mask)
              double_max = maxval(data_double,mask=mask)
            end if

            if ( double_min.lt.-1e30 ) then
                fv = .true.
                double_fv = double_min
                call check( nf90_put_att(ncid_out, varid(i), '_FillValue', &
                int(-2**(8*nb-1),nb) ), 37 )
                goto 220
            end if
            if ( double_max.gt.1e30 ) then
                fv = .true.
                double_fv = double_max
                call check( nf90_put_att(ncid_out, varid(i), '_FillValue', &
                int(-2**(8*nb-1),nb) ), 37 )
                goto 220
            end if

            if (fv) then
               double_sf = (double_max-double_min)/(2**(8*nb)-2)
               double_af = double_min+(2**(8*nb-1)-1)*double_sf
            else
               double_sf = (double_max-double_min)/(2**(8*nb)-1)
               double_af = double_min+(2**(8*nb-1))*double_sf
            end if

            write(*,*) trim(varname(i))
            write(*,*) 'scale_factor:',double_sf
            write(*,*) 'add_offset:',double_af
            write(*,*) 'min:',nint((double_min-double_af)/double_sf,nb)
            write(*,*) 'max:',nint((double_max-double_af)/double_sf,nb)

            call check( nf90_put_att(ncid_out, varid(i), 'add_offset', double_af ), 38 )
            call check( nf90_put_att(ncid_out, varid(i), 'scale_factor', double_sf ), 39 )
            call check( nf90_enddef(ncid_out), 101 )


            if ( ndims.gt.4 ) then
              where( mask_1d ) zipped_1d = nint((data_double_1d-double_af)/double_sf,nb)
              tmp = product(dimlens(1:4))
              idx(:) = 1
              do j = 1, product(dimlens(5:ndims))
                do k = 4,ndims-1
                  if (idx(k).gt.dimlens(k)) then
                    idx(k+1) = idx(k+1)+1
                    idx(k) = 1
                  end if
                end do

                zipped = reshape(zipped_1d(1+(j-1)*tmp:j*tmp),&
                    (/dimlens(1),dimlens(2),dimlens(3),dimlens(4)/))

                call check( nf90_put_var(ncid_out, varid(i), zipped, start = idx, &
                count = (/dimlens(1), dimlens(2), dimlens(3), dimlens(4), &
                  (k**0, k=1,ndims-4)/) ), 105 )

                idx(4) = idx(4)+dimlens(4)
              end do
            else
              where( mask ) zipped = nint((data_double-double_af)/double_sf,nb)
              call check( nf90_put_var(ncid_out, varid(i), &
                zipped(1:shp(1),1:shp(2),1:shp(3),1:shp(4))), 105 )
            end if

            deallocate( data_double, zipped )
            if ( ndims.gt.4 ) then
                deallocate( data_double_1d, zipped_1d, mask_1d )
            else
                deallocate( mask )
            end if
  !##########################################################          
        end select

      end if

      deallocate( dimids, dimlens, idx )
    end do

    call check( nf90_close(ncid_out), 100 )
    call check( nf90_close(ncid_in), 100 )

    deallocate( var, zip_varid, dimid, dimlen, varid, dimname, varname )
END PROGRAM nc_zip

SUBROUTINE check(status,label)
    USE NETCDF
    implicit none
    INTEGER, INTENT(in) :: status , label

    IF (status .ne. nf90_noerr) then
        print*, trim(nf90_strerror(status)),' label=',label
        stop "Stopped"
    END IF
END SUBROUTINE check
