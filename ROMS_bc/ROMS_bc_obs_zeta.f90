PROGRAM ROMS_BC
    USE NETCDF
    implicit none
    integer, parameter :: nx_zeta=36,ny_zeta=22
    integer, parameter :: xin=50,yin=70,zin=56
    integer :: narg, year, id, ncid, dimid, vid, varid(8)
    integer :: nx, ny, nz, i, j, k, l, x, y
    integer :: xi_d,s_d,temp_t_d,zeta_t_d,salt_t_d
    integer :: dimids4(4), eta_d
    integer :: dimids3temp(3),dimids2zeta(2),dimids3salt(3)
    integer :: minpos(2), date(3), date_0(3),date_tmp(3)
    integer :: jbc, ibcstr, ibcend, cnt, julian
    integer :: timeid, nt_zeta, nt_tr,minz
    integer :: istr,iend,span
    integer, allocatable :: mask(:,:), idx(:,:,:)
    real(kind=8) :: hc, gote_loc(2), zeta_ave, splint
    real(kind=8) :: point(2), corners(4,2), dzdx
    real(kind=8) :: lon_tr(xin), lat_tr(yin), depth_tr(zin)
    real(kind=8) :: tmpzin(zin), tmp_tr(xin,yin,zin)  
    real(kind=8), allocatable :: Cs_r(:), s_rho(:), h(:,:)
    real(kind=8), allocatable :: lon(:,:), lat(:,:), depth(:,:)
    real(kind=8), allocatable :: dis_array(:,:), skag_zeta(:)
    real(kind=8), allocatable :: gote_zeta(:), W(:,:), zeta(:,:)
    real(kind=8), allocatable :: temp(:,:,:), salt(:,:,:), time(:)
    real(kind=8), allocatable :: tempin(:,:,:,:), saltin(:,:,:,:)
    real(kind=8), allocatable :: tr_tmp(:,:),temp_day(:,:,:),salt_day(:,:,:)
    real(kind=8), allocatable :: tmp_temp(:,:), tmp_salt(:,:), y2a(:,:)
    real(kind=8), allocatable :: y2_1d(:), lon_val(:)
    real(kind=8), allocatable :: y2a_2(:,:), y2_1d_2(:), dp_val(:)
    character(len=150) :: buffer, gridfile, fileout, zetafile, trfile
    logical :: leapy, fexit, score, mask_zeta(nx_zeta,ny_zeta), first
    logical :: mask_tr(xin,yin,zin),ex
   
    zetafile = 'zeta_bc_v5.nc'
    gote_loc(1) = 11.791
    gote_loc(2) = 57.685
    date_0=(/1968,5,23/)

!READ command line arg: YEAR
    narg = command_argument_count()
    if (narg.ne.2) then
      write(*,*) "Program must have YEAR and gridfile as arguments"
      stop
    end if
    call get_command_argument(1,buffer)
    read(buffer, *) year
    call get_command_argument(2,buffer)
    write(gridfile, *) buffer
    gridfile = adjustl(gridfile)
    id = index(buffer,'.nc')
    write(fileout,'(A,A,i4.4,A)') trim(buffer(1:id-1)),'_bc_obs_zetav5_',year,'.nc'
!READ command line

!CFRETE outfile name
    write(*,'(a,i4.4)') 'year: ',year
    write(*,'(a,a)') 'fileout: ',trim(fileout)
!END CREATE

!FIND timeid for zetafile
    timeid = 0
    do i = 1992,year-1
        CALL leap(i,leapy)
        if (leapy) then
            timeid = timeid+24*366
        else
            timeid = timeid+24*365
        end if
    end do
    timeid=timeid+1
    write(*,'(A,i)') 'timeid: ', timeid

    CALL leap(year,leapy)
    nt_tr = 365
    if (leapy) nt_tr = 366
!    nt_tr = 1 !!!!!!TESTOWE!!!!!!!!!!!!!!!
    nt_zeta = 24*nt_tr
    allocate( skag_zeta(nt_zeta), dis_array(nx_zeta,ny_zeta) )
    allocate( gote_zeta(nt_zeta), y2a(xin,yin) ) 
    allocate( tempin(xin,yin,zin,nt_tr+2), saltin(xin,yin,zin,nt_tr+2) ) 
!END FIND

!READ grid
    CALL check(nf90_open(trim(gridfile),NF90_NOWRITE,ncid),10)
    CALL check(nf90_inq_dimid(ncid, "xi_rho", dimid),11)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nx),12)
    CALL check(nf90_inq_dimid(ncid, "eta_rho", dimid),13)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=ny),14)
    CALL check(nf90_inq_dimid(ncid, "s_rho", dimid),15)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nz),16)
    allocate( Cs_r(nz), s_rho(nz), h(nx,ny), depth(nx,nz) )
    allocate( lon(nx,ny), lat(nx,ny), mask(nx,ny) )
    CALL check(nf90_inq_varid(ncid,"h",vid),17)
    CALL check(nf90_get_var(ncid,vid,h),18)
    CALL check(nf90_inq_varid(ncid,"hc",vid),19)
    CALL check(nf90_get_var(ncid,vid,hc),20)
    CALL check(nf90_inq_varid(ncid,"lon_rho",vid),21)
    CALL check(nf90_get_var(ncid,vid,lon),22)
    CALL check(nf90_inq_varid(ncid,"lat_rho",vid),27)
    CALL check(nf90_get_var(ncid,vid,lat),28)
    CALL check(nf90_inq_varid(ncid,"Cs_r",vid),31)
    CALL check(nf90_get_var(ncid,vid,Cs_r),32)
    CALL check(nf90_inq_varid(ncid,"s_rho",vid),33)
    CALL check(nf90_get_var(ncid,vid,s_rho),34)
    CALL check(nf90_inq_varid(ncid,"mask_rho",vid),33)
    CALL check(nf90_get_var(ncid,vid,mask),34)
    CALL check(nf90_close(ncid),35)
!END READ grid

!FIND NORTH BC
    jbc = 0
    cnt = 0
    do i = 1, nx
      do j = ny, 1, -1
        if ( mask(i,j).eq.1 ) then
          if ( j.eq.jbc ) then
            cnt = cnt+1
          else
            cnt = 1
            jbc = j
            ibcstr = i
          end if
          EXIT
        end if 
      end do
      if ( cnt.ge.20 ) EXIT
    end do
    do i = ibcstr, nx
      if ( mask(i,jbc).eq.0 ) then
        ibcend = i-1
        EXIT
      end if
    end do
    write(*,*) 'jbc,ibcstr,ibcend: ',jbc, ibcstr, ibcend
    ibcstr = ibcstr
    ibcend = ibcend
    allocate( idx(1+ibcend-ibcstr,4,2), W(1+ibcend-ibcstr,4) )
    allocate( tr_tmp(11+ibcend-ibcstr,zin+1) )
    allocate( y2a_2(11+ibcend-ibcstr,zin+1) )
    allocate( y2_1d_2(11+ibcend-ibcstr), dp_val(11+ibcend-ibcstr) )
    allocate( y2_1d(xin), lon_val(xin) )
    allocate( zeta(nx,nt_zeta), temp(nx,nz,24*nt_tr) )
    allocate( salt(nx,nz,24*nt_tr), time(nt_zeta) )

    do i = 1,nx
      do k = 1,nz
        depth(i,k) = h(i,jbc)*(hc*s_rho(k)+h(i,jbc)*&
                                Cs_r(k))/(hc+h(i,jbc))
      end do
    end do

    !%%%%%%%%%%
    istr = ibcstr-4
    iend = ibcend+4
    span = 1+iend-istr
    !%%%%%%%%%%

    allocate( salt_day(span,nz,nt_tr+2), temp_day(span,nz,nt_tr+2))
    allocate( tmp_temp(span,nz), tmp_salt(span,nz) )
!END FIND

!READ ZETA FROM FILE
    CALL check(nf90_open(trim(zetafile),NF90_NOWRITE,ncid),1015)
    CALL check(nf90_inq_varid(ncid,"skag",vid),21)
    CALL check(nf90_get_var(ncid,vid,skag_zeta,&
         start = (/timeid/),count = (/nt_zeta/)),22)
    CALL check(nf90_inq_varid(ncid,"gote",vid),31)
    CALL check(nf90_get_var(ncid,vid,gote_zeta,&
         start = (/timeid/),count = (/nt_zeta/)),32)
    CALL check(nf90_close(ncid),23)
!END READ

!ZETA interpolation
    do k = 1, nt_zeta
      zeta(ibcstr,k) = skag_zeta(k)
      zeta(ibcend,k) = gote_zeta(k)
      dzdx = (gote_zeta(k)-skag_zeta(k))/real(ibcend-ibcstr,8)
      do i = ibcstr+1,ibcend-1
        zeta(i,k) = zeta(i-1,k)+dzdx
      end do
      zeta(1:ibcstr-1,k)=zeta(ibcstr,k)
      zeta(ibcend+1:nx,k)=zeta(ibcend,k) 
    end do
!END ZETA

!go to 2222 !TESTOWE!


!READ TRACERS FROM FILES
    100 format('files/'i4.4,'/',i2.2,'/CMEMS_BAL_PHY_reanalysis_dailymeans_',i4,i2.2,i2.2,'.nc')
    date = (/year-1,12,31/)
    first = .true.
    do i = 1,nt_tr+2
      write(trfile,100) date(1),date(2),date(1),date(2),date(3)
      if (first) then
        INQUIRE(FILE=trim(trfile),EXIST=ex)
        if (.not.ex) then
          date_tmp = date
          CALL add_day(date_tmp,.true.)
          write(trfile,100) date_tmp(1),date_tmp(2),date_tmp(1),&
            date_tmp(2),date_tmp(3) 
        end if
      end if    
 
      CALL check(nf90_open(trim(trfile),NF90_NOWRITE,ncid),1016)
      
      if (first) then  
        CALL check(nf90_inq_varid(ncid,"depth",vid),17)
        CALL check(nf90_get_var(ncid,vid,tmpzin),18)

        do k = 1,zin
          depth_tr(k) = (-1)*tmpzin(1+zin-k) !-1 BECAUSE DEPTH IN ROMS IS NEGATIVE
        end do
        depth_tr(zin) = 0
        do k = 1,zin
          if ( depth_tr(k).gt.minval(depth) ) exit
        end do
        minz = k-1
        

        CALL check(nf90_inq_varid(ncid,"longitude",vid),19)
        CALL check(nf90_get_var(ncid,vid,lon_tr,start=(/20/),&
            count=(/xin/)),18)

        do k = 1,xin
          if ( lon(ibcstr,jbc).lt.lon_tr(k) ) exit
        end do
        write(*,*) k
        write(*,*) lon(ibcstr,jbc),lon_tr(k)

        CALL check(nf90_inq_varid(ncid,"latitude",vid),20)
        CALL check(nf90_get_var(ncid,vid,lat_tr,start=(/250/),&
            count=(/yin/)),18)

        do k = xin,1,-1
          if ( lon(ibcend,jbc).gt.lon_tr(k) ) exit
        end do
        write(*,*) k
        write(*,*) lon(ibcend,jbc),lon_tr(k)

      end if

      CALL check(nf90_inq_varid(ncid,"thetao",vid),17)
      CALL check(nf90_get_var(ncid,vid,tmp_tr,start=(/20,250,1,1/),&
          count=(/xin,yin,zin,1/)),18)

      do k = 1,56
        tempin(:,:,k,i) = tmp_tr(:,:,57-k)
      end do

      if (first) then
        mask_tr = .true.    
        where( tempin(:,:,:,1).ne.tempin(:,:,:,1) ) mask_tr = .false.   

!        mask_tr(:12,:,:) = .false.
!        mask_tr(27:,:,:) = .false. 
        mask_tr(:10,:,:) = .false.
        mask_tr(36:,:,:) = .false. 

        first = .false.
      end if

      CALL check(nf90_inq_varid(ncid,"so",vid),17)
      CALL check(nf90_get_var(ncid,vid,tmp_tr,start=(/20,250,1,1/),&
          count=(/xin,yin,zin,1/)),18)

      do k = 1,56
        saltin(:,:,k,i) = tmp_tr(:,:,57-k)
      end do

      CALL check(nf90_close(ncid),23)
      write(*,*) 'read day:',i

      CALL extrap3d(tempin(:,:,:,i),mask_tr,xin,yin,zin,500,1)
      CALL extrap3d(saltin(:,:,:,i),mask_tr,xin,yin,zin,500,2)

      CALL add_day(date,.true.)
    end do 
!END READ TRACERS


!    CALL check( nf90_create( 'diagnostic.nc',NF90_NETCDF4,ncid ), 300 )
!    CALL check( nf90_def_dim( ncid, 'xi_rho', xin, xi_d ), 301 )
!    CALL check( nf90_def_dim( ncid, 'eta_rho', yin, eta_d ), 301 )
!    CALL check( nf90_def_dim( ncid, 's_rho', zin, s_d ), 302 )
!    CALL check( nf90_def_dim( ncid, 'temp_time', nt_tr+2, temp_t_d ), 303 )
!    dimids4 = (/ xi_d,eta_d,s_d,temp_t_d /)
!    CALL check(nf90_def_var( ncid, 'temp_north', NF90_FLOAT, dimids4,&
!               varid(2) ), 316)
!    CALL check(nf90_def_var( ncid, 'salt_north', NF90_FLOAT, dimids4,&
!               varid(3) ), 316)
!    CALL check( nf90_enddef(ncid), 340 )
!    CALL check( nf90_put_var( ncid, varid(2), real(tempin,4) ), 342 )
!    CALL check( nf90_put_var( ncid, varid(3), real(saltin,4) ), 342 )
!    CALL check(nf90_close(ncid),350)
    

    
!NEW SPLINE TRACER INTERPOLATION

    do l = 1, nt_tr+2
      call spline3d(xin,yin,zin,lon_tr,lat_tr,depth_tr,&
        span,nz,lon(istr:iend,jbc),lat(istr:iend,jbc),&
        depth(istr:iend,:),tempin(:,:,:,l),temp_day(:,:,l))
 
      call spline3d(xin,yin,zin,lon_tr,lat_tr,depth_tr,&
        span,nz,lon(istr:iend,jbc),lat(istr:iend,jbc),&
        depth(istr:iend,:),saltin(:,:,:,l),salt_day(:,:,l)) 
    end do

!END NEW SPLINE TRACER INTERPOLATION

goto 1000
!NEW TRACER INTERPOLATION

    CALL check( nf90_create( 'after_1st_spline.nc',NF90_NETCDF4,ncid ), 300 )
    CALL check( nf90_def_dim( ncid, 'xi_rho', 11+ibcend-ibcstr, xi_d ), 301 )
    CALL check( nf90_def_dim( ncid, 's_rho', zin+1, s_d ), 302 )
    CALL check( nf90_def_dim( ncid, 'temp_time', nt_tr+2, temp_t_d ), 303 )
    dimids3temp = (/ xi_d,s_d,temp_t_d /)
    CALL check(nf90_def_var( ncid, 'temp_north', NF90_FLOAT, dimids3temp,&
               varid(1) ), 316)
    CALL check(nf90_def_var( ncid, 'salt_north', NF90_FLOAT, dimids3temp,&
               varid(2) ), 316)
    CALL check( nf90_enddef(ncid), 340 )

    do l = 1,nt_tr+2

      do k = 1,zin
        CALL splie2(lat_tr,xin,yin,tempin,y2a)
        do x = ibcstr-5, ibcend+5
          do i = 1, xin
            lon_val(i) = splint(lat_tr,tempin(i,:,k,l),y2a(i,:),yin,lat(x,jbc))
          end do
          call spline(lon_tr,lon_val,xin,y2_1d)
          tr_tmp(x+6-ibcstr,k) = splint(lon_tr,lon_val,y2_1d,xin,lon(x,jbc))
        end do
      end do
      tr_tmp(:,zin+1) = tr_tmp(:,zin)

      CALL check( nf90_put_var( ncid, varid(1), real(tr_tmp,4), &
        start=(/1,1,l/) ),342 )

      !temp_day
!      CALL splie2(depth_tr,11+ibcend-ibcstr,zin+1,tr_tmp,y2a_2)
!      do x = ibcstr-4, ibcend+4
!        do y = 1, nz
!          do i = 1, 11+ibcend-ibcstr
!            dp_val(i) = splint(depth_tr,tr_tmp(i,:),y2a_2(i,:),&
!                zin+1,depth(x,y)) 
!          end do
!          call spline(lon(ibcstr-5:ibcend+5,jbc),dp_val,11+ibcend-ibcstr,y2_1d_2)
!          temp_day(x,y,l) = splint(lon(ibcstr-5:ibcend+5,jbc),dp_val,y2_1d_2,&
!            11+ibcend-ibcstr,lon(x,jbc)) 
!        end do
!      end do
      

      do k = 1,zin
        CALL splie2(lat_tr,xin,yin,saltin,y2a)
        do x = ibcstr-5, ibcend+5
          do i = 1, xin
            lon_val(i) = splint(lat_tr,saltin(i,:,k,l),y2a(i,:),yin,lat(x,jbc))
          end do
          call spline(lon_tr,lon_val,xin,y2_1d)
          tr_tmp(x+6-ibcstr,k) = splint(lon_tr,lon_val,y2_1d,xin,lon(x,jbc))
        end do
      end do
      tr_tmp(:,zin+1) = tr_tmp(:,zin)

      CALL check( nf90_put_var( ncid, varid(2), real(tr_tmp,4), &
        start=(/1,1,l/) ),343 )

      !salt_day
      CALL splie2(depth_tr,11+ibcend-ibcstr,zin+1,tr_tmp,y2a_2)
      do x = ibcstr-4, ibcend+4
        do y = 1, nz
          do i = 1, 11+ibcend-ibcstr
            dp_val(i) = splint(depth_tr,tr_tmp(i,:),y2a_2(i,:),&
                zin+1,depth(x,y))
          end do
          call spline(lon(ibcstr-5:ibcend+5,jbc),dp_val,11+ibcend-ibcstr,y2_1d_2)
          salt_day(x,y,l) = splint(lon(ibcstr-5:ibcend+5,jbc),dp_val,y2_1d_2,&
            11+ibcend-ibcstr,lon(x,jbc))
        end do
      end do

    end do

    CALL check(nf90_close(ncid),350) 
!END OF NEW TRACER INTERPOLATION


!TRACER INTERPOLATION
    deallocate( dis_array )
    allocate( dis_array(xin,yin) )
    idx = 0
    W = 0
    do i = ibcstr, ibcend
      do k = 1,xin
        do l = 1,yin
          dis_array(k,l) = sqrt((lon_tr(k)-lon(i,jbc))**2+&
              (lat_tr(l)-lat(i,jbc))**2)
        end do
      end do
      minpos = minloc(dis_array)
      x = minpos(1)
      y = minpos(2)
      if (x.gt.1 .and. x.lt.xin .and. y.gt.1 .and. y.lt.yin) then
          point(1) = lon(i,jbc)
          point(2) = lat(i,jbc)
          fexit = .false.
          do k = x-1,x
            if (fexit) EXIT
          do l = y-1,y
            if (fexit) EXIT
            corners(1,:) = (/lon_tr(k),lat_tr(l)/)
            corners(2,:) = (/lon_tr(k),lat_tr(l+1)/)
            corners(3,:) = (/lon_tr(k+1),lat_tr(l+1)/)
            corners(4,:) = (/lon_tr(k+1),lat_tr(l)/)
            CALL in_convex_polygon(corners,point,score)
            if (score) then
              if ( .not.mask_tr(k,l,zin) .or. &
                   .not.mask_tr(k+1,l,zin) .or. &
                   .not.mask_tr(k,l+1,zin) .or. &
                   .not.mask_tr(k+1,l+1,zin) ) then
                write(*,*) 'interpolation from land, i=',i
                fexit = .true.
                idx(1+i-ibcstr,1,:) = 0
              else
                fexit = .true.
                idx(1+i-ibcstr,1,:) = (/k,l/)
                idx(1+i-ibcstr,2,:) = (/k,l+1/)
                idx(1+i-ibcstr,3,:) = (/k+1,l+1/)
                idx(1+i-ibcstr,4,:) = (/k+1,l/)
                call calc_w( lon_tr(k:k+1), lat_tr(l:l+1),&
                  point, W(1+i-ibcstr,:) )
              end if
            end if
          end do
          end do
        else
          idx(1+i-ibcstr,1,:) = 0
        end if
    end do

    do i = 1,nt_tr+2
      tr_tmp = -999
      do k = 1,zin
        CALL calcme(xin,yin,1+ibcend-ibcstr,tempin(:,:,k,i),&
          tr_tmp(:,k),W,idx)
      end do
      tr_tmp(:,zin+1) = tr_tmp(:,zin)
      do j = 1,ibcend-ibcstr
        if ( tr_tmp(j+1,zin).eq.-999 .and. tr_tmp(j,zin).ne.-999 ) then
          tr_tmp(j+1,:) = tr_tmp(j,:)
        end if
      end do
      do j = 1+ibcend-ibcstr,2,-1
        if ( tr_tmp(j-1,zin).eq.-999 .and. tr_tmp(j,zin).ne.-999 ) then
          tr_tmp(j-1,:) = tr_tmp(j,:)
        end if
      end do
      
      CALL vertical_interp(1+ibcend-ibcstr,nz,zin+1,depth_tr,&
          depth(ibcstr:ibcend,:),tr_tmp,temp_day(:,:,i))
      
      tr_tmp = -999
      do k = 1,zin
        CALL calcme(xin,yin,1+ibcend-ibcstr,saltin(:,:,k,i),&
          tr_tmp(:,k),W,idx)
      end do
      tr_tmp(:,zin+1) = tr_tmp(:,zin)
      do j = 1,ibcend-ibcstr
        if ( tr_tmp(j+1,zin).eq.-999 .and. tr_tmp(j,zin).ne.-999 ) then
          tr_tmp(j+1,:) = tr_tmp(j,:)
        end if
      end do
      do j = 1+ibcend-ibcstr,2,-1
        if ( tr_tmp(j-1,zin).eq.-999 .and. tr_tmp(j,zin).ne.-999 ) then
          tr_tmp(j-1,:) = tr_tmp(j,:)
        end if
      end do

      CALL vertical_interp(1+ibcend-ibcstr,nz,zin+1,depth_tr,&
          depth(ibcstr:ibcend,:),tr_tmp,salt_day(:,:,i)) 
    end do 
!END TRACER INTERPOLATION


1000 continue
!****************LINEAR_INTERPOLATION_OF_TRACERS***********
    cnt = 0
    do i = 1,nt_tr
      tmp_temp(:,:) = (temp_day(:,:,i+1)-temp_day(:,:,i))/24
      tmp_salt(:,:) = (salt_day(:,:,i+1)-salt_day(:,:,i))/24
      do j = 1,12
        temp(ibcstr-4:ibcend+4,:,24*cnt+j) = temp_day(:,:,i)+(j+11)*tmp_temp(:,:)
        salt(ibcstr-4:ibcend+4,:,24*cnt+j) = salt_day(:,:,i)+(j+11)*tmp_salt(:,:)
        do k = 1,nz
          temp(1:ibcstr-5,k,24*cnt+j)=temp(ibcstr-4,k,24*cnt+j)
          salt(1:ibcstr-5,k,24*cnt+j)=salt(ibcstr-4,k,24*cnt+j)
          temp(ibcend+5:nx,k,24*cnt+j)=temp(ibcend+4,k,24*cnt+j)
          salt(ibcend+5:nx,k,24*cnt+j)=salt(ibcend+4,k,24*cnt+j)
        end do
      end do
      tmp_temp(:,:) = (temp_day(:,:,i+2)-temp_day(:,:,i+1))/24
      tmp_salt(:,:) = (salt_day(:,:,i+2)-salt_day(:,:,i+1))/24
      do j = 13,24
        temp(ibcstr-4:ibcend+4,:,24*cnt+j) = temp_day(:,:,i+1)+(j-13)*tmp_temp(:,:)
        salt(ibcstr-4:ibcend+4,:,24*cnt+j) = salt_day(:,:,i+1)+(j-13)*tmp_salt(:,:)
        do k = 1,nz
          temp(1:ibcstr-5,k,24*cnt+j)=temp(ibcstr-4,k,24*cnt+j)
          salt(1:ibcstr-5,k,24*cnt+j)=salt(ibcstr-4,k,24*cnt+j)
          temp(ibcend+5:nx,k,24*cnt+j)=temp(ibcend+4,k,24*cnt+j)
          salt(ibcend+5:nx,k,24*cnt+j)=salt(ibcend+4,k,24*cnt+j)
        end do
      end do

      cnt = cnt+1
    end do
!**********************************************************

!2222 continue !TESTOWE!

!FILL TIME VALUE
    cnt = 0
    date = (/year,1,1/)
    do i = 1,nt_tr
        cnt = cnt+1
        CALL elapsed(date,date_0,julian)
        do j = 1,24
            time((cnt-1)*24+j) = real(julian,8)+real(j-1,8)/24
        end do
        CALL add_day(date,.true.)
    end do
!END FILL

    
    CALL check( nf90_create( trim(fileout),NF90_NETCDF4,ncid ), 300 )
    CALL check( nf90_def_dim( ncid, 'xi_rho', nx, xi_d ), 301 )
    CALL check( nf90_def_dim( ncid, 's_rho', nz, s_d ), 302 )
    CALL check( nf90_def_dim( ncid, 'temp_time', 24*nt_tr, temp_t_d ), 303 )
    CALL check( nf90_def_dim( ncid, 'zeta_time', nt_zeta, zeta_t_d ), 403 )
    CALL check( nf90_def_dim( ncid, 'salt_time', 24*nt_tr, salt_t_d ), 404 )
    dimids3temp = (/ xi_d,s_d,temp_t_d /)
    dimids2zeta = (/ xi_d,zeta_t_d /)
    dimids3salt = (/ xi_d,s_d,salt_t_d /)

    CALL check(nf90_def_var( ncid, 'temp_time', NF90_DOUBLE,&
        (/temp_t_d/),varid(1) ), 304)
    CALL check( nf90_put_att( ncid, varid(1), 'long_name',&
        'potential temperature time' ), 305 )
    CALL check( nf90_put_att( ncid, varid(1), 'units',&
        'days since 1968-05-23 00:00:00 GMT' ), 306 )
    CALL check( nf90_put_att( ncid, varid(1), 'calendar',&
        'gregorian' ), 307 )

    CALL check(nf90_def_var( ncid, 'zeta_time', NF90_DOUBLE,&
        (/zeta_t_d/),varid(3) ), 308)
    CALL check( nf90_put_att( ncid, varid(3), 'long_name',&
        'free-surface time' ), 309 )
    CALL check( nf90_put_att( ncid, varid(3), 'units',&
        'days since 1968-05-23 00:00:00 GMT' ), 310 )
    CALL check( nf90_put_att( ncid, varid(3), 'calendar',&
        'gregorian' ), 311 )

    CALL check(nf90_def_var( ncid, 'salt_time', NF90_DOUBLE,&
        (/salt_t_d/),varid(5) ), 312)
    CALL check( nf90_put_att( ncid, varid(5), 'long_name',&
        'surface net heat flux time' ), 313 )
    CALL check( nf90_put_att( ncid, varid(5), 'units',&
        'days since 1968-05-23 00:00:00 GMT' ), 314 )
    CALL check( nf90_put_att( ncid, varid(5), 'calendar',&
        'gregorian' ), 315 )

    CALL check(nf90_def_var( ncid, 'temp_north', NF90_FLOAT, dimids3temp,&
               varid(2) ), 316)
    CALL check( nf90_put_att( ncid, varid(2), 'long_name',&
        'potential temperature northern boundary condition' ), 317 )
    CALL check( nf90_put_att( ncid, varid(2), 'units',&
        'Celsius' ), 318 )
    CALL check( nf90_put_att( ncid, varid(2), 'time',&
        'temp_time' ), 319 )

    CALL check(nf90_def_var( ncid, 'zeta_north', NF90_FLOAT, dimids2zeta,&
               varid(4) ), 320)
    CALL check( nf90_put_att( ncid, varid(4), 'long_name',&
        'free-surface northern boundary condition' ), 321 )
    CALL check( nf90_put_att( ncid, varid(4), 'units',&
        'meter' ), 322 )
    CALL check( nf90_put_att( ncid, varid(4), 'time',&
        'zeta_time' ), 323 )

    CALL check(nf90_def_var( ncid, 'salt_north', NF90_FLOAT, dimids3salt,&
               varid(6) ), 324)
    CALL check( nf90_put_att( ncid, varid(6), 'long_name',&
        'salinity northern boundary condition' ), 325 )
    CALL check( nf90_put_att( ncid, varid(6), 'time',&
        'salt_time' ), 326 )

    CALL check(nf90_def_var( ncid, 'lon', NF90_DOUBLE, (/xi_d/),&
               varid(7) ), 327)
    CALL check( nf90_put_att( ncid, varid(7), 'long_name',&
        'longitude of RHO-points' ), 328 )
    CALL check( nf90_put_att( ncid, varid(7), 'units',&
        'degree_east' ), 329 )

    CALL check(nf90_def_var( ncid, 'dp', NF90_DOUBLE, (/xi_d,s_d/),&
               varid(8) ), 330)
    CALL check( nf90_put_att( ncid, varid(8), 'long_name',&
        'depth of RHO-points' ), 331 )
    CALL check( nf90_put_att( ncid, varid(8), 'units',&
        'meter' ), 332 )

    CALL check( nf90_enddef(ncid), 340 )
    CALL check( nf90_put_var( ncid, varid(1), time ), 341 )
    CALL check( nf90_put_var( ncid, varid(2), real(temp,4) ), 342 )
    CALL check( nf90_put_var( ncid, varid(3), time ), 343 )
    CALL check( nf90_put_var( ncid, varid(4), real(zeta,4) ), 344 )
    CALL check( nf90_put_var( ncid, varid(5), time ), 345 )
    CALL check( nf90_put_var( ncid, varid(6), real(salt,4) ), 346 )
    CALL check( nf90_put_var( ncid, varid(7), lon(:,jbc) ), 347 )
    CALL check( nf90_put_var( ncid, varid(8), depth ), 348 )
    CALL check(nf90_close(ncid),350)

    deallocate( skag_zeta, dis_array, gote_zeta, tr_tmp )
    deallocate( Cs_r, s_rho, h, lon, lat, mask)
    deallocate( idx, W, zeta, temp, salt, y2a )
    deallocate( y2_1d, lon_val )
    deallocate( y2a_2, y2_1d_2, dp_val )
    deallocate( time, depth, tempin, saltin)
    deallocate( salt_day, temp_day, tmp_temp, tmp_salt )
END PROGRAM ROMS_BC

SUBROUTINE check(status,label)
    USE NETCDF
    implicit none
    INTEGER, INTENT(in) :: status , label

    IF (status .ne. nf90_noerr) then
        print*, trim(nf90_strerror(status)),' label= ',label
        stop "Stopped"
    END IF
END SUBROUTINE check

SUBROUTINE leap(year,answer)
    implicit none
    integer, intent(in) :: year
    logical, intent(out) :: answer

    answer = .false.
    if( mod(year,4).eq.0 ) answer = .true.
    if( mod(year,4).eq.100 ) answer = .false.
    if( mod(year,4).eq.400 ) answer = .true.
END SUBROUTINE leap

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

!!$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j),SCHEDULE(DYNAMIC)
        do i=2,lon-1
          do j=2,lat-1
            res(i,j)=0.25*(a(i-1,j)+a(i+1,j)+a(i,j-1)+a(i,j+1))-a(i,j)
          end do
        end do
!!$OMP END PARALLEL DO

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

SUBROUTINE calc_w(lonin,latin,point_xy,w)
    implicit none
    real(kind=8), intent(in) :: lonin(2),latin(2)
    real(kind=8), intent(in) :: point_xy(2)
    real(kind=8), intent(out) :: w(4)

    real(kind=8) :: x, y, lon(2,2), lat(2,2)
    real(kind=8) :: A, B, C, s, t
    real(kind=8) :: y1,y2,x1,x2

    x = point_xy(1)
    y = point_xy(2)
    lon(1,:) = lonin
    lon(2,:) = lonin
    lat(:,1) = latin
    lat(:,2) = latin

    A = ( lon(1,1)-lon(1,2) )*( lat(2,1)-lat(2,2) )-&
        ( lat(1,1)-lat(1,2) )*( lon(2,1)-lon(2,2) )

    B = y*( ( lon(2,1)-lon(2,2) )-( lon(1,1)-lon(1,2) ) )-&
        x*( ( lat(2,1)-lat(2,2) )-( lat(1,1)-lat(1,2) ) )+&
        ( lon(1,1)-lon(1,2) )*lat(2,2)-&
        ( lat(1,1)-lat(1,2) )*lon(2,2)+&
        ( lat(2,1)-lat(2,2) )*lon(1,2)-&
        ( lon(2,1)-lon(2,2) )*lat(1,2)

    C = y*( lon(2,2)-lon(1,2) )-&
        x*( lat(2,2)-lat(1,2) )+&
        lon(1,2)*lat(2,2)-&
        lat(1,2)*lon(2,2)

    if ( A.ne.0 ) then 
      t = ( -B+sqrt(B**2-4*A*C) )/(2*A)
      s = (  y-lat(1,2)-( lat(1,1)-lat(1,2) )*t  )/&
        (  lat(2,2)+( lat(2,1)-lat(2,2) )*t-&
        lat(1,2)-( lat(1,1)-lat(1,2) )*t  )

      if ( (t<0).OR.(t>1) ) then
        t = ( -B-sqrt(B**2-4*A*C) )/(2*A)
        s = (  y-lat(1,2)-( lat(1,1)-lat(1,2) )*t  )/&
          (  lat(2,2)+( lat(2,1)-lat(2,2) )*t-&
          lat(1,2)-( lat(1,1)-lat(1,2) )*t  )
      end if

      w(1) = (1-s)*t
      w(2) = (1-s)*(1-t)
      w(3) = s*(1-t)
      w(4) = s*t
    else
      x1 = lonin(1)
      x2 = lonin(2)
      y1 = latin(1)
      y2 = latin(2)
      w(1) = ((x2-x)/(x2-x1))*((y2-y)/(y2-y1)) 
      w(2) = ((x2-x)/(x2-x1))*((y-y1)/(y2-y1)) 
      w(3) = ((x-x1)/(x2-x1))*((y-y1)/(y2-y1)) 
      w(4) = ((x-x1)/(x2-x1))*((y2-y)/(y2-y1)) 
    end if

END SUBROUTINE calc_w

SUBROUTINE calcme(nxin,nyin,nxout,inarray,outarray,W,idx)
    implicit none
    integer, intent(in) :: nxin,nyin,nxout
    real(kind=8), intent(in) :: inarray(nxin,nxin),W(nxout,4)
    real(kind=8), intent(out) :: outarray(nxout)
    integer, intent(in) :: idx(nxout,4,2)

    real(kind=8) :: tmp
    integer :: i,k,x,y

    do i = 1,nxout
      if (idx(i,1,1).ne.0) then
        tmp = 0
        do k = 1,4
          x = idx(i,k,1)
          y = idx(i,k,2)
          if (x.gt.0 .and. y.gt.0) then
            tmp = tmp+inarray(x,y)*W(i,k)
          end if
        end do
        outarray(i) = tmp
      end if
    end do
END SUBROUTINE calcme

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

SUBROUTINE extrap3d(a,mask,lon,lat,n,maxscn,met)
    USE omp_lib
    implicit none
    integer, intent(in) :: lon,lat,n,maxscn,met
    real(kind=8), intent(inout) :: a(lon,lat,n)
    logical, intent(in) :: mask(lon,lat,n)

    integer :: i,j,k,l,cnt,overall
    real(kind=8) :: relc,ave
    real(kind=8), dimension(lon,lat,n) :: sor,res
    logical :: mask_tmp(lon,lat,n),mask_tmp2(lon,lat,n)

    relc = 0.6
    sor = 0.0
    where(.not.mask) sor=relc

    !FILLING LAND WITH LAYER AVERAGE
    select case(met)
    case(0)
      where(.not.mask) a=0.0
    case(1)
      do k = 1,n

      cnt = 0
      ave = 0.0

      do i=1,lon
      do j=1,lat
        if (mask(i,j,k)) then
            ave=ave+a(i,j,k)
            cnt=cnt+1
        end if
      end do
      end do
      if ( cnt.GT.0 ) ave = ave/real(cnt,8)
      where(.not.mask(:,:,k)) a(:,:,k)=ave

      end do
    case(2)
      mask_tmp2 = mask
      mask_tmp = mask_tmp2
      do

      overall = 0
      do i = 1, lon
      do j = 1, lat
      do k = 1, n

        if (.not.mask_tmp(i,j,k)) then

          cnt = 0
          ave = 0.0

          if ( i.gt.1 ) then
            if ( mask_tmp(i-1,j,k) ) then
              ave = ave+a(i-1,j,k)
              cnt = cnt+1
            end if
          end if

          if ( j.gt.1 ) then
            if ( mask_tmp(i,j-1,k) ) then
              ave = ave+a(i,j-1,k)
              cnt = cnt+1
            end if
          end if

          if ( k.gt.1 ) then
            if ( mask_tmp(i,j,k-1) ) then
              ave = ave+a(i,j,k-1)
              cnt = cnt+1
            end if
          end if

          if ( i.lt.lon ) then
            if ( mask_tmp(i+1,j,k) ) then
              ave = ave+a(i+1,j,k)
              cnt = cnt+1
            end if
          end if

          if ( j.lt.lat ) then
            if ( mask_tmp(i,j+1,k) ) then
              ave = ave+a(i,j+1,k)
              cnt = cnt+1
            end if
          end if

          if ( k.lt.n ) then
            if ( mask_tmp(i,j,k+1) ) then
              ave = ave+a(i,j,k+1)
              cnt = cnt+1
            end if
          end if

          if ( cnt.gt.0 ) then
            a(i,j,k) = ave/(real(cnt,8))
            overall = overall+cnt
            mask_tmp2(i,j,k) = .true.
          end if

        end if

      end do
      end do
      end do
      mask_tmp = mask_tmp2
      if ( overall.eq.0 ) EXIT

      end do
    end select
    !END FILLING

    do l = 1,maxscn

!$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k),SCHEDULE(DYNAMIC)
      do i=2,lon-1
        do j=2,lat-1
          res(i,j,n)=0.2*(a(i-1,j,n)+a(i+1,j,n)+a(i,j-1,n)+a(i,j+1,n)+&
            a(i,j,n-1))-a(i,j,n)
          res(i,j,1)=0.2*(a(i-1,j,1)+a(i+1,j,1)+a(i,j-1,1)+a(i,j+1,1)+&
            a(i,j,2))-a(i,j,1)
        do k = 2,n-1
          res(i,j,k)=0.1667*(a(i-1,j,k)+a(i+1,j,k)+a(i,j-1,k)+a(i,j+1,k)+&
            a(i,j,k+1)+a(i,j,k-1))-a(i,j,k)
        end do
        end do
      end do
!$OMP END PARALLEL DO

      do k = 2,n-1
      do i=2,lon-1
        res(i,1,k)=0.2*(a(i-1,1,k)+a(i+1,1,k)+a(i,1,k+1)+a(i,1,k-1)+a(i,2,k))-a(i,1,k)
        res(i,lat,k)=0.2*(a(i-1,lat,k)+a(i+1,lat,k)+a(i,lat,k+1)+a(i,lat,k-1)+a(i,lat-1,k))-a(i,lat,k)
      end do
      do j=2,lat-1
        res(1,j,k)=0.2*(a(1,j-1,k)+a(1,j+1,k)+a(1,j,k+1)+a(1,j,k-1)+a(2,j,k))-a(1,j,k)
        res(lon,j,k)=0.2*(a(lon,j-1,k)+a(lon,j+1,k)+a(lon,j,k+1)+a(lon,j,k-1)+a(lon-1,j,k))-a(lon,j,k)
      end do
      end do

      do i=2,lon-1
        res(i,1,1)=0.25*(a(i-1,1,1)+a(i+1,1,1)+a(i,1,2)+a(i,2,1))-a(i,1,1)
        res(i,lat,1)=0.25*(a(i-1,lat,1)+a(i+1,lat,1)+a(i,lat,2)+a(i,lat-1,1))-a(i,lat,1)
        res(i,1,n)=0.25*(a(i-1,1,n)+a(i+1,1,n)+a(i,1,n-1)+a(i,2,n))-a(i,1,n)
        res(i,lat,n)=0.25*(a(i-1,lat,n)+a(i+1,lat,n)+a(i,lat-1,n)+a(i,lat,n-1))-a(i,lat,n)
      end do

      do j=2,lat-1
        res(1,j,1)=0.25*(a(1,j-1,1)+a(1,j+1,1)+a(1,j,2)+a(2,j,1))-a(1,j,1)
        res(lon,j,1)=0.25*(a(lon,j-1,1)+a(lon,j+1,1)+a(lon,j,2)+a(lon-1,j,1))-a(lon,j,1)
        res(1,j,n)=0.25*(a(1,j-1,n)+a(1,j+1,n)+a(1,j,n-1)+a(2,j,n))-a(1,j,n)
        res(lon,j,n)=0.25*(a(lon,j-1,n)+a(lon,j+1,n)+a(lon,j,n-1)+a(lon-1,j,n))-a(lon,j,n)
      end do

      do k = 2,n-1
       res(1,1,k)=0.25*(a(1,1,k-1)+a(1,1,k+1)+a(2,1,k)+a(1,2,k))-a(1,1,k)
       res(lon,1,k)=0.25*(a(lon,1,k-1)+a(lon,1,k+1)+a(lon-1,1,k)+a(lon,2,k))-a(lon,1,k)
       res(1,lat,k)=0.25*(a(1,lat,k-1)+a(1,lat,k+1)+a(2,lat,k)+a(1,lat-1,k))-a(1,lat,k)
       res(lon,lat,k)=0.25*(a(lon,lat,k-1)+a(lon,lat,k+1)+a(lon-1,lat,k)+a(lon,lat-1,k))-a(lon,lat,k)
      end do

      res(1,1,1)=0.3333*(a(2,1,1)+a(1,2,1)+a(1,1,2))-a(1,1,1)
      res(lon,1,1)=0.3333*(a(lon-1,1,1)+a(lon,2,1)+a(lon,1,2))-a(lon,1,1)
      res(1,lat,1)=0.3333*(a(2,lat,1)+a(1,lat-1,1)+a(1,lat,2))-a(1,lat,1)
      res(lon,lat,1)=0.3333*(a(lon-1,lat,1)+a(lon,lat-1,1)+a(lon,lat,2))-a(lon,lat,1)
      res(1,1,n)=0.3333*(a(2,1,n)+a(1,2,n)+a(1,1,n-1))-a(1,1,n)
      res(lon,1,n)=0.3333*(a(lon-1,1,n)+a(lon,2,n)+a(lon,1,n-1))-a(lon,1,n)
      res(1,lat,n)=0.3333*(a(2,lat,n)+a(1,lat-1,n)+a(1,lat,n-1))-a(1,lat,n)
      res(lon,lat,n)=0.3333*(a(lon-1,lat,n)+a(lon,lat-1,n)+a(lon,lat,n-1))-a(lon,lat,n)


      res=res*sor
      a=a+res
    end do

END SUBROUTINE extrap3d

SUBROUTINE vertical_interp(nx,Nlvl,zin,depth,z_rho,input,output)
    implicit none
    integer, intent(in) :: nx,Nlvl,zin
    real(kind=8), intent(in) :: depth(zin),z_rho(nx,Nlvl)
    real(kind=8), intent(in) :: input(nx,zin)
    real(kind=8), intent(out) :: output(nx,Nlvl)

    integer :: i,k,m
    real(kind=8) :: y0,y1,x0,x1,x

    do i = 1, nx
    do k = Nlvl,1,-1
      do m = zin,2,-1
        if (z_rho(i,k).le.depth(m) .and. z_rho(i,k).gt.depth(m-1)) then
          x = z_rho(i,k)
          x0 = depth(m-1)
          x1 = depth(m)
          y0 = input(i,m-1)
          y1 = input(i,m)
          output(i,k) = y0+(y1-y0)*(x-x0)/(x1-x0)
          EXIT
        else
           continue
        end if
      end do
    end do
    end do
END SUBROUTINE vertical_interp

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

SUBROUTINE splie2(x2a,m,n,ya,y2a)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: m,n
    REAL(kind=8), DIMENSION(n), INTENT(IN) :: x2a
    REAL(kind=8), DIMENSION(m,n), INTENT(IN) :: ya
    REAL(kind=8), DIMENSION(m,n), INTENT(OUT) :: y2a

    INTEGER :: j

    DO j=1,m
        call spline(x2a,ya(j,:),n,y2a(j,:))
    END DO
END SUBROUTINE splie2

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


!spline3d(xin,yin,zin+1,lon_tr,lat_tr,depth_tr,&
!   span,nz,lon(istr:iend,jbc),lat(istr:iend,jbc),
!   depth(istr:iend,:),tempin(:,:,:,l),temp_day(:,:,l))
SUBROUTINE spline3d(nxin,nyin,nzin,inlon,inlat,indp,&
            nxout,nzout,outlon,outlat,outdp,&
            input, output)
    implicit none
    integer, intent(in) :: nxin, nyin, nzin, nxout, nzout
    real(kind=8), intent(in) :: inlon(nxin), inlat(nyin),indp(nzin)
    real(kind=8), intent(in) :: outlon(nxout), outlat(nxout),outdp(nxout,nzout)
    real(kind=8), intent(in) :: input(nxin,nyin,nzin)
    real(kind=8), intent(out) :: output(nxout,nzout)

    integer :: i,j,k,x,z
    real(kind=8) :: splint
    real(kind=8) :: y2a(nxin,nyin,nzin)
    real(kind=8) :: slice(nxin,nzin), y2_slice(nxin,nzin)
    real(kind=8) :: line(nzin), y2_line(nzin)

    do i = 1, nxin
      do k = 1, nzin
        call spline(inlat,input(i,:,k),nyin,y2a(i,:,k)) 
      end do
    end do 

    do x = 1, nxout

      do i = 1, nxin
        do k = 1, nzin
          slice(i,k) = splint(inlat,input(i,:,k),y2a(i,:,k),nyin,outlat(x)) 
        end do
      end do
      
      do k = 1, nzin
        call spline(inlon,slice(:,k),nxin,y2_slice(:,k))
      end do

      do k = 1, nzin
        line(k) = splint(inlon,slice(:,k),y2_slice(:,k),nxin,outlon(x))
      end do

      call spline(indp,line,nzin,y2_line)
      do z = 1, nzout
        output(x,z) = splint(indp,line,y2_line,nzin,outdp(x,z))
      end do

    end do

END SUBROUTINE spline3d
