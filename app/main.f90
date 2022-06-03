program main
  use day19
  implicit none

  type(scanner_t), allocatable :: scanners(:)
  integer :: i, j, npoints, manh, manh_max
  integer, allocatable :: rm(:,:,:), buoys(:,:), buoysr(:,:)
  real :: time0, time1

  !character(len=*), parameter :: fileinp = 'test.inp'
  character(len=*), parameter :: fileinp = 'input.inp'

  ! Input from file
  call buoys_read(fileinp, scanners)
  npoints = 0
  do i = 1, size(scanners)
    call scanners(i) % print(.false.)
    npoints = npoints + scanners(i)%nb()
  enddo

  ! Calibrate scanners
  call cpu_time(time0)
  rm = rotation_matrices()
  call calibrate_all(scanners, rm)
  do i=1,size(scanners)
    call scanners(i) % print(.false.)
  enddo

  ! Extract buoys
  allocate(buoys(3, npoints))
  j = 0
  do i=1, size(scanners)
    call scanners(i) % extract_buoys(rm, buoys(:, j+1:j+scanners(i)%nb()))
    j = j+scanners(i)%nb()
  end do
  call sort(buoys, 3)
  call sort(buoys, 2)
  call sort(buoys, 1)
  buoysr = remove_duplicates(buoys)

  ! Results of part One
  !print '(3(i6))', buoysr
  print '("Number of buoys = ",i0)', size(buoysr,2)

  ! Part Two
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
