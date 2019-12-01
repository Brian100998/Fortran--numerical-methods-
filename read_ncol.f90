    subroutine readncol(ncol,x, n)

    use nrtype

    implicit none
    real(DP), dimension(:,:), allocatable, intent(out)   :: x        ! array to fill data in
    integer(I4B), intent(out)                                :: n      ! number of readed points
    integer(I4B), intent(in)                                 :: ncol

    real(DP), dimension(:,:), allocatable                :: oldx     ! auxiliady arrays
    character (len=50)          :: datafile
    integer(I4B)                :: state

    ! ask user the name of the file
    write(*,'(a)',advance="no") "Enter name of the n-column data file: "
    read(*,*) datafile

    ! open file to read
    open(75,file=datafile,action="read",status="old")

    ! initialization
    if (allocated(x)) deallocate( x )
    allocate(x(0,0))
    n = 0

    do
        allocate(oldx(n+1, ncol))
        oldx(1:n,:) = x(1:n,:)

        read(75,*,iostat = state) oldx(n+1,:)
        if (state /= 0 ) exit
        n = n + 1

        deallocate(x)
        allocate(x(n,ncol))
        x(1:n,:) = oldx(1:n,:)  ! x = oldx doesn't work anymore (it used to)
        deallocate(oldx)
    end do

    close(75)

    end subroutine