PROGRAM robert
    USE NETCDF
    implicit none
    integer :: ncid,dimid,vid,nt 
    character(len=20) :: input,output
    real(kind=8), allocatable :: skag(:),gote(:)
    real(kind=8), allocatable :: skag2(:),gote2(:)

    input = 'zeta_bc_v3.nc'
    output = 'zeta_bc_v5.nc'

    CALL check(nf90_open(trim(input),NF90_NOWRITE,ncid),310)
    CALL check(nf90_inq_dimid(ncid, "time", dimid),311)
    CALL check(nf90_inquire_dimension(ncid, dimid, len=nt),312)

    allocate( skag(nt), gote(nt) )
    allocate( skag2(nt), gote2(nt) )

    CALL check(nf90_inq_varid(ncid,"skag",vid),17)
    CALL check(nf90_get_var(ncid,vid,skag),18)
    CALL check(nf90_inq_varid(ncid,"gote",vid),19)
    CALL check(nf90_get_var(ncid,vid,gote),20)
    CALL check(nf90_close(ncid),21)

    CALL filter(nt,skag,skag2,0.50)
    CALL filter(nt,gote,gote2,0.50)

    CALL check(nf90_open(trim(output),NF90_WRITE,ncid),210)
    CALL check(nf90_inq_varid(ncid,"skag",vid),117)
    CALL check(nf90_put_var(ncid,vid,skag2),118)
    CALL check(nf90_inq_varid(ncid,"gote",vid),119)
    CALL check(nf90_put_var(ncid,vid,gote2),120)
    CALL check(nf90_close(ncid),121)
    
    deallocate( skag, gote )
    deallocate( skag2, gote2 )
    
END PROGRAM robert

SUBROUTINE check(status,label)
    USE NETCDF
    implicit none
    INTEGER, INTENT(in) :: status , label

    IF (status .ne. nf90_noerr) then
        print*, trim(nf90_strerror(status)),' label=',label
        stop "Stopped"
    END IF
END SUBROUTINE check

SUBROUTINE filter(nt,input,output,alpha)
    implicit none
    integer,intent(in) :: nt
    real(kind=4), intent(in) :: alpha
    real(kind=8), intent(in) :: input(nt)
    real(kind=8), intent(out) :: output(nt)

    integer :: i

    output(1) = input(1)
    write(*,*) alpha
    
    do i = 2,nt-1
       output(i) = input(i)+alpha*(input(i+1)-2*input(i)+output(i-1)) 
    end do
    
    output(nt) = input(nt)+alpha*(output(nt-1)-input(nt))
    

END SUBROUTINE filter
