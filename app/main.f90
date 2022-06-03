!
! 2021 Advent of Code - Solution
!
program main
  use day19
  implicit none

  type(scanner_t), allocatable :: scanners(:)
  integer :: i, j, manh, manh_max
  integer, allocatable :: rms(:,:,:), buoys(:,:)
  real :: time0, time1

  !character(len=*), parameter :: fileinp = 'test.inp'
  character(len=*), parameter :: fileinp = 'input.inp'

  ! Input from file
  call buoys_read(fileinp, scanners)
  do i = 1, size(scanners)
    call scanners(i) % print(.false.)
  enddo

  ! Calibrate scanners
  call cpu_time(time0)
  rms = rotation_matrices()
  call calibrate_all(scanners, rms)
  do i=1,size(scanners)
    call scanners(i) % print(.false.)
  enddo

  ! Extract buoys
  buoys = extract_all(scanners, rms)

  ! Part One - Number of buoys
  !print '(3(i6))', buoys
  print '("Number of buoys = ",i0)', size(buoys,2)

  ! Part Two - Manhattan distance
  manh_max = -huge(manh)
  do i=1,size(scanners)-1
  do j=i+1,size(scanners)
    manh = scanners(i) % manhattan(scanners(j))
    if (manh > manh_max) manh_max = manh
  end do
  end do
  print '("Maximum Manhattan distance = ",i0)', manh_max

  call cpu_time(time1)
  print '("Time taken = ",f7.2," seconds")', time1-time0
end program main
