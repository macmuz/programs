PROGRAM TEST
    implicit none
    integer :: narg,i
    character(len=10) :: input

    100 format('task number: ', a, ' step: ',i2.2)

    narg = command_argument_count()
    if (narg.ne.1) then
      write(*,*) "Program must have 1 argument"
      stop
    end if

    call get_command_argument(1,input)

    do i = 1,5
        call sleep(5)
        write(*,100) trim(input),i
    end do
END PROGRAM TEST
