program IAM

  implicit none
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
  INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)

  real(kind=dp) :: a(4,151), b(4,151), c(151), minq,maxq
  real(kind=dp), allocatable :: q(:), aff(:,:), iat(:), imol(:), itot(:), sincq(:), geom(:,:)
  integer(kind=ikind) :: nq, natoms, i,j
  character(len=10) :: atoms(151)
  character(len=50) :: xyzfile
  character(len=10), allocatable :: list_of_atoms(:)
  character(len=8) :: filename="affl.txt"


  IF (COMMAND_ARGUMENT_COUNT() .LT. 1) THEN
    WRITE (*, *) 'ONE INPUTS REQUIRED, ABORTING...'
    STOP
  else 
    CALL GET_COMMAND_ARGUMENT(1, xyzfile) ! file to be checked
  end if

  call table_of_ff(filename, atoms, a, b, c) 

  minq = 0.000000000001
  maxq = 5
  nq = 200
  call read_xyz(geom, xyzfile, list_of_atoms, natoms)
  allocate(q(nq))
  allocate(sincq(nq))
  allocate(iat(nq))
  allocate(imol(nq))
  allocate(itot(nq))
  allocate(aff(nq,natoms))

  call linspace(minq,maxq,nq,q)

  do i = 1, natoms, 1
    call obtain_form_factors(aff(:,i), list_of_atoms(i), q, atoms, a, b, c)
  end do


  imol(:) = 0

  iat = SUM(aff**2, DIM=2)
  do i = 1, natoms
  do j = i+1, natoms
  sincq = sinc(q * dist(geom(i,:),geom(j,:),3))
  imol = imol+(2*aff(:,i)*aff(:,j)*sincq)
  end do
  end do 

  itot=imol+iat

  call write_signal(itot,'alal') 

  contains

    function dist(geom1,geom2,n)
            INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
            INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
            real(kind=dp), external :: dnrm2
            real(kind=dp), intent(in):: geom1(:), geom2(:)
            real(kind=dp) :: dist
            integer(kind=ikind), intent(in) :: n

            dist = dnrm2(n, geom1-geom2, 1)
    end function
            



    elemental function sinc (a)
             INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
            INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
            real(kind=dp), intent(in) :: a
             real(kind=dp):: sinc
             if (abs(a) < 1.0d-10) then
                sinc = 1
             else
                sinc = sin(a) / (a)
             end if
    end function
    subroutine linspace(from, to, n,array)
     INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
     INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
     integer,intent(in) :: n
    real(kind=dp), intent(in) :: from, to
    real(kind=dp), intent(out),dimension(n) :: array
    real(kind=dp) :: range

    integer :: i

    range = to - from

    do i=1, n
        array(i) = from + range * (i - 1) / (n - 1)
    end do
end subroutine
    subroutine obtain_form_factors(aff1,atom,q1,atoms,a,b,c)
        implicit none
        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        real(kind=dp), parameter :: pi=dacos(-1.00_dp)
        character(len=10), intent(in) :: atom
        real(kind=dp), intent(in), dimension(:) :: q1
        real(kind=dp),intent(out),dimension(:):: aff1
        character(len=10),intent(in),dimension(151):: atoms
        real(kind=dp), intent(in),dimension(4,151) :: a,b
        real(kind=dp), intent(in),dimension(151):: c
        integer(kind=ikind) :: index,i





        ! allocate(aff1(size(q1)))

        do i=1,151
           if (atom==atoms(i)) then
               index=i
            end if
        end do

        do i=1,size(q1)
            aff1(i)=c(index)+sum(a(:,index)*exp(-b(:,index)*(q1(i)/(4*pi))**2.0))
        enddo

    end subroutine obtain_form_factors
  
    subroutine table_of_ff(filename,atoms,a,b,c)
        implicit none
        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        integer :: i
        character(len=*), intent(in) :: filename

        character(len=10),intent(out), dimension(151):: atoms
        real(kind=dp), intent(out),dimension(4,151) :: a,b
        real(kind=dp), intent(out),dimension(151):: c

        open(unit=15,file=TRIM(filename))

        do i=1,151
            read(15,*)atoms(i),a(1,i) ,b(1,i),a(2,i),b(2,i),a(3,i),b(3,i),a(4,i),b(4,i),c(i)

        end do
        close(unit=15)
    end subroutine table_of_ff

    subroutine write_signal(signal,filename)
      real(kind=dp), intent(in) :: signal(:)
      character(len=*), intent(in) :: filename

      ! open(17,file=trim(filename))
      ! do i = 1,size(signal)
      ! write(17,'(ES12.5)') signal(i)
      ! close(17)
      do i = 1,size(signal)
      write(*,'(ES12.5)') signal(i)
      end do
      end subroutine

  subroutine read_xyz(GEOM, filename, atoms, natoms)
    implicit none

  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
  INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
    real(kind=dp),  ALLOCATABLE, intent(OUT) :: GEOM(:,:)
    CHARACTER(len=50),intent(IN):: filename
    CHARACTER(LEN=10), DIMENSION(:), ALLOCATABLE , intent(OUT):: atoms
    integer(kind=ikind), intent(out) :: natoms
    integer(kind=ikind) :: i, fid, io

    fid=89
    OPEN (fid, file=TRIM(filename))
    READ (fid, *,IOSTAT=io) natoms 
    !This catches the end of a file to stop the calculation
    if (io .LT. 0) then
      return
    end if 
    read(fid,*)

    if (.not. allocated(GEOM)) then
    allocate (GEOM(natoms, 3))
    allocate (atoms(natoms))
    end if

    DO i = 1, natoms
    READ (fid, *,IOSTAT=io) atoms(i), GEOM(i, :)
    end do

  end subroutine
end program
