module day19
  implicit none
  private
  public scanner_t, buoys_read, calibrate_all, extract_all
  public rotation_matrices

  integer, parameter :: EYE_INDEX = 1, ROT_NULL = -99

  type scanner_t
    integer, allocatable :: buoys(:,:)
    integer :: id = -1
    integer :: rot = ROT_NULL
    integer :: pos(3) = 0
  contains
    procedure :: nb => scanner_nbuoys
    procedure :: manhattan => scanner_manhattan
    procedure :: extract_buoys => scanner_extract_buoys
    procedure :: print => scanner_print
  end type

contains

  pure integer function scanner_nbuoys(this) result(nb)
    class(scanner_t), intent(in) :: this
    if (allocated(this%buoys)) then
      nb = size(this%buoys,dim=2)
    else
      nb = -1
    endif
  end function scanner_nbuoys



  pure integer function scanner_manhattan(this, other) result(d)
    class(scanner_t), intent(in) :: this, other
    integer :: diss(size(this%pos,1))
    diss = abs(this%pos-other%pos) 
    d = sum(diss)
  end function scanner_manhattan



  pure subroutine scanner_extract_buoys(this, rms, buoys)
    class(scanner_t), intent(in) :: this
    integer, intent(in) :: rms(:,:,:)
    integer, intent(out) :: buoys(:,:)
!
! Get buoys positions in global co-ordinates
!
    integer :: i

    if (this%rot == ROT_NULL) &
    &   error stop 'extract_buoys - scanner not calibrated' 
    if (size(buoys,2) /= this%nb() .or. size(buoys,1) /= size(this%buoys,1)) &
    &   error stop 'extract_buoys - output array shape inconsistent'
    do i = 1, this%nb()
      buoys(:,i) = matmul(rms(:,:,this%rot), this%buoys(:,i)) + this%pos
    end do
  end subroutine scanner_extract_buoys



  subroutine scanner_print(this, buoys_printed)
    class(scanner_t), intent(in) :: this
    logical, intent(in), optional :: buoys_printed

    logical buoys_printed0
    buoys_printed0 = .true.
    if (present(buoys_printed)) buoys_printed0 = buoys_printed

    if (this % rot /= ROT_NULL) then
      print '("Scanner ",i2,"  nb=",i2,"  rot=",i2,"  pos= [",3(i6),"]")', &
      &   this%id, this%nb(), this%rot, this%pos
    else
      print '("Scanner ",i2,"  nb=",i2,"  uncalibrated")', &
      &   this%id, this%nb()
    endif
    if (buoys_printed0) then
      print '(3(i6))', this%buoys
      print *
    end if
  end subroutine scanner_print



  pure function extract_all(scs,rms) result(res)
    type(scanner_t), intent(in) :: scs(:)
    integer, intent(in) :: rms(:,:,:)
    integer, allocatable :: res(:,:)
!
! Get buoys list from all scanners. Remove duplicates.
!
    integer, allocatable :: buoys(:,:)
    integer :: i, j, nb

    nb = 0
    do i=1,size(scs)
      nb = nb + scs(i) % nb()
    end do
    allocate(buoys(3, nb))
    j = 0
    do i=1, size(scs)
      call scs(i) % extract_buoys(rms, buoys(:, j+1:j+scs(i)%nb()))
      j = j+scs(i)%nb()
    end do
    call sort(buoys, 3)
    call sort(buoys, 2)
    call sort(buoys, 1)
    res = remove_duplicates(buoys)
    deallocate(buoys)
  end function extract_all



  subroutine calibrate_all(scs, rms)
    type(scanner_t), intent(inout) :: scs(:)
    integer, intent(in) :: rms(:,:,:)

    integer :: list(2,size(scs)*(size(scs)-1)/2)
    integer :: newlist(2,size(scs)*(size(scs)-1)/2)
    integer :: i, j, k, n_list, n_newlist
    logical :: skipped

    ! prepare working list
    n_list = 0
    do i=1, size(scs)-1
    do j=i+1, size(scs)
      n_list = n_list + 1
      list(:,n_list) = [i, j]
    end do
    end do

    ! un-calibrate all scanners, mark first scanner as a reference [0,0,0]
    scs % rot = ROT_NULL
    scs(1) % rot = EYE_INDEX
    scs(1) % pos = 0

    do
      ! test all scanner pairs combinations
      print '("Scanners left to calibrate = ",i0,"  pairs to test = ",i0)', &
      &   count(scs%rot==ROT_NULL), n_list
      if (n_list == 0) exit
      n_newlist = 0
      print *
      do k=1,n_list
        call calibrate_pair(scs(list(1,k)), scs(list(2,k)), rms, skipped)
        if (.not. skipped) cycle
        n_newlist = n_newlist + 1
        newlist(:,n_newlist) = list(:,k)
      end do
      n_list = n_newlist
      list(:,1:n_list) = newlist(:,1:n_list)
    end do
  end subroutine calibrate_all



  subroutine calibrate_pair(s1, s2, rms, skipped)
    type(scanner_t), intent(inout) :: s1, s2
    integer, intent(in) :: rms(:,:,:)
    logical, intent(out) :: skipped
!
! Orient the calibrated scanner from the pair accordingly and test for all 24
! orientations of the uncalibrated scanner. If both scanners are uncalibrated,
! then return for now.
!
    integer, allocatable :: dis(:,:), map(:,:), mat(:,:)
    integer :: k

    skipped = .false.
    do k=1, size(rms,dim=3)
      if (s1%rot /= ROT_NULL .and. s2%rot == ROT_NULL) then
        dis = scanner_pair(s1, rms(:,:,s1%rot), s2, rms(:,:,k))
      elseif (s1%rot == ROT_NULL .and. s2%rot /= ROT_NULL) then
        dis = scanner_pair(s1, rms(:,:,k), s2, rms(:,:,s2%rot))
      elseif (s1%rot /= ROT_NULL .and. s2%rot /= ROT_NULL) then
        return ! both scanners calibrated, nothing to do
      else
        skipped = .true. ! both positions uknown, we will call back later
        return
      endif
      map = scanner_pair_map(s1,s2)
      call sort(dis,3,map)
      call sort(dis,2,map)
      call sort(dis,1,map)
      mat = extract_match(dis,map)
      if (size(mat,dim=2) > 0) exit
    end do
    if (size(mat,dim=2) <= 0) return

    ! match found - update scanner position and orientation
    if (s1%rot /= ROT_NULL .and. s2%rot == ROT_NULL) then
      s2 % rot = k
      s2 % pos = s1 % pos + mat(1:3,1)
      !call s2 % print(.false.)
    elseif (s1%rot == ROT_NULL .and. s2%rot /= ROT_NULL) then
      s1 % rot = k
      s1 % pos = s2 % pos - mat(1:3,1)
      !call s1 % print(.false.)
    else
      error stop 'calibrate_pair - impossible branch'
    endif
  end subroutine calibrate_pair



  pure function scanner_pair(s1,rm1,s2,rm2) result(dis)
    type(scanner_t), intent(in) :: s1, s2
    integer, intent(in) :: rm1(:,:), rm2(:,:)
    integer dis(3, s1%nb()*s2%nb())
!
! For every pair {[a];[b]} where [a] is buoy from S1 and [b] is buoy from S2:
! 1. orient [a] and [b] to a common reference coordinate system
!     [a'] = [[rot A]]*[a]  and  [b'] = [[rot B]]*[b]
! 2. create vector 
!     [v]_ab = [a']-[b']
! 3. store this vector to the list of vectors "dis"
!
! If both [a] and [b] reference the same buoy with position [w'], and
! [S1] and [S2] are positions of scanners 1 and 2, then
!   [a'] = [w']-[S1]  
!   [b'] = [w']-[S2]
!   [v]_ab = [a']-[b'] = [w']-[S1]-[w']+[S2] = [S2]-[S1]
! and therefore, the vector is the relative distance from S1 to S2.
!
! When the same vector repeats in the list "dis", this indicates the overlap of
! several buoys and a correct relative orientation between scanners.
!
! This function should be called together with "scanner_pair_map" if the
! identification of [a],[b] entries from list "dis" will be later needed.
!
    integer :: i, j, ij
    integer :: b(3), apos(3,s1%nb())

    do i=1, s1%nb()
      apos(:,i) = matmul(rm1, s1%buoys(:,i))
    end do

    ij = 0
    do j=1, s2%nb()
      b = matmul(rm2, s2%buoys(:,j))
      do i=1, s1%nb()
        ij = ij+1
        dis(:, ij) = apos(:,i) - b
      end do
    end do
  end function



  pure function scanner_pair_map(s1,s2) result(map)
    type(scanner_t), intent(in) :: s1, s2
    integer map(2, s1%nb()*s2%nb())
!
! Allows to track corresponding entries for vectors [a],[b] in the list "dis"
! This function should be called together with "scanner_pair".
!
    integer :: i, j, ij
    ij = 0
    do j=1, s2%nb()   
      do i=1, s1%nb() 
        ij = ij+1
        map(:, ij) = [i, j]
      end do
    end do
  end function scanner_pair_map



  pure subroutine sort(arr,ind,map)
    integer, intent(inout) :: arr(:,:)
    integer, intent(in) :: ind
    integer, intent(inout), optional :: map(:,:)
!
! Insertion sort according to "ind" column in the vector list "arr"
!
    integer :: n, i, j
    integer :: sav(size(arr,dim=1))
    integer, allocatable :: sav_map(:)
    
    if (ind<1 .or. ind>size(arr,dim=1)) error stop 'sort index out of bounds'
    n = size(arr,dim=2)
    if (present(map)) then
      if (size(map,dim=2) /= n) error stop 'sort - arr and map not same size'
      allocate(sav_map(size(map,dim=1)))
    end if

    do i=2,n
      sav = arr(:,i)
      if (present(map)) sav_map = map(:,i)
      do j=i-1,1,-1
        if (arr(ind,j) <= sav(ind)) exit
        arr(:,j+1) = arr(:,j) 
        if (present(map)) map(:,j+1) = map(:,j)
      end do
      arr(:,j+1) = sav
      if (present(map)) map(:,j+1) = sav_map
    end do
    if (present(map)) deallocate(sav_map)
  end subroutine sort



  function extract_match(dis,map) result(match)
    integer, intent(in) :: dis(:,:), map(:,:)
    integer, allocatable :: match(:,:)

    integer, parameter :: INROW_REQ = 12
    integer :: i, inrow, blocks_found, nf
    integer :: tmp(size(dis,1)+size(map,1), size(dis,2))

    if (size(dis,2) /= size(map,2)) error stop 'extract - arr and map not same size'
    nf = 0
    blocks_found = 0
    inrow = 1
    do i=2, size(dis,2)
      if (all(dis(:,i-1)==dis(:,i))) then
        inrow = inrow + 1
        cycle
      endif
      ! inrow broken
      if (inrow >= INROW_REQ) then
        tmp(:size(dis,1),  nf+1:nf+inrow) = dis(:,i-inrow:i-1)
        tmp(size(dis,1)+1:,nf+1:nf+inrow) = map(:,i-inrow:i-1)
        nf = nf + inrow
        blocks_found = blocks_found + 1
      end if
      inrow = 1
    end do

    ! if the last item is part of block (this will be rarely run)
    if (inrow >= INROW_REQ) then
print *, 'rare - rare'
        tmp(:size(dis,1),  nf+1:nf+inrow) = dis(:,i-inrow:i-1)
        tmp(size(dis,1)+1:,nf+1:nf+inrow) = map(:,i-inrow:i-1)
        nf = nf + inrow
        blocks_found = blocks_found + 1
    end if

    allocate(match(size(dis,1)+size(map,1),nf))
    match = tmp(:,1:nf)
    if (nf>0) print '("Matching entries = ",i0," in ",i0," blocks")', &
    &  nf, blocks_found
  end function extract_match



  pure function remove_duplicates(pos) result(cpos)
    integer, intent(in) :: pos(:,:)
    integer, allocatable :: cpos(:,:)

    integer :: i, uniq
    integer :: tmp(size(pos,1),size(pos,2))

    uniq = 0
    do i=2, size(pos,2)
      if (all(pos(:,i-1)==pos(:,i))) then
        cycle
      endif
      ! inrow broken
      uniq = uniq + 1
      tmp(:,uniq) = pos(:,i-1)
    end do
    uniq = uniq + 1
    tmp(:,uniq) = pos(:,i-1)
    allocate(cpos(size(pos,1),uniq))
    cpos = tmp(:,1:uniq)
  end function remove_duplicates



  pure function rotation_matrices() result(rms)
    integer :: rms(3,3,24)
!
! Vector [v] is transformed to alternative co-ordinate system [valt] using
! rotation matrix [[A]] 
!   [valt] = [[A]] * [v]
!
! Rows of rotation matrix are base vectors defining the coordinate system
! There are six directions where the first base vector [a] can point:
!   [a] = { [1,0,0], [-1,0,0], [0,1,0], [0,-1,0], [0,0,1], [0,0,-1] }
! Then there are four directions where the second base vector [b] can point.
! For example, if [a] = [1,0,0] these posibilities are left:
!   [b] = { [0,1,0], [0,-1,0], [0,0,1], [0,0,-1] }
! The third base vector [c] is then fixed as the cross-product:
!   [c] = [a] x [b]
!
! There are therefore 6*4 = 24 possible rotation matrices.
!
    integer :: imat, i, j, ii, jj

    rms = 0
    imat = 0
    do i=1, 3
    do j=1, 3
      if (j==i) cycle
      do ii=1,-1,-2
      do jj=1,-1,-2
        imat = imat + 1
        rms(1,i,imat) = ii
        rms(2,j,imat) = jj
        rms(3,:,imat) = cross_product(rms(1,:,imat), rms(2,:,imat))
      enddo
      enddo
    enddo
    enddo

    ! Assert the identity matrix is at EYE_INDEX position
    associate (E=>EYE_INDEX)
      if (rms(1,1,E)/=1 .or. rms(2,2,E)/=1 .or. rms(3,3,E)/=1) &
      &  error stop 'rotation_matrices - matrix at EYE_INDEX is not identity matrix'
    end associate
  end function rotation_matrices



  pure function cross_product(a,b) result(c)
    integer, intent(in) :: a(:), b(:)
    integer :: c(3)
    c(1) = a(2)*b(3)-a(3)*b(2)
    c(2) = a(3)*b(1)-a(1)*b(3)
    c(3) = a(1)*b(2)-a(2)*b(1)
  end function



  subroutine buoys_read(file, scanners)
    character(len=*), intent(in) :: file
    type(scanner_t), allocatable, intent(out) :: scanners(:)
!
! Parse input
!
    integer, parameter :: NSMAX = 100, BMAX = 100, DIM=3
    type(scanner_t) :: tmp_scanners(NSMAX)
    integer :: tmp_buoys(DIM,BMAX)

    integer :: fid, ns, nb, iostat
    character(len=1000) :: line

    open(newunit=fid, file=file, status='old')
    ns = 0
    do
      read(fid,'(a)',iostat=iostat) line
      if (iostat /= 0) exit
      if (line(1:11) == '--- scanner') call scanner_read()
      if (line=='') cycle
    enddo

    allocate(scanners(ns))
    scanners = tmp_scanners(:ns)
  contains
    subroutine scanner_read()
      ns = ns + 1
      nb = 0
      read(line(12:),*) tmp_scanners(ns) % id
      do
        read(fid,*,iostat=iostat) tmp_buoys(:,nb+1)
        if (iostat /= 0) exit
        nb = nb + 1
      enddo
      backspace(fid)
      tmp_scanners(ns) % buoys = tmp_buoys(:,:nb)
    end subroutine
  end subroutine buoys_read

end module day19
