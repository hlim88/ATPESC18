function MLIFE_nextstate(matrix, lni, lnj, row, col)
  use mlife2dtypes, only : BORN, DIES
  integer MLIFE_nextstate
  integer, allocatable, intent(in) :: matrix(:,:)
  integer :: row, col, lni, lnj
  integer :: sum

  ! add values of all eight neighbors
  sum = matrix(row-1,col-1) + matrix(row-1,col) + &
        matrix(row-1,col+1) + matrix(row,col-1) + &
        matrix(row,col+1)   + matrix(row+1,col-1) + &
        matrix(row+1,col)   + matrix(row+1,col+1)

  if (sum < 2 .or. sum > 3) then
     MLIFE_nextstate = DIES
  else if (sum == 3) then
     MLIFE_nextstate = BORN
  else
     MLIFE_nextstate = matrix(row,col)
  endif
end function MLIFE_nextstate

subroutine MLIFE_TimeIterations(patch, nIter, doCheckpoint, m1, m2, &
     exchangeInit, exchange, exchangeEnd, timedata)
  use mpi
  use mlife2dtypes
  use mlife2dfuncs
  use mlife2dio
  use mlifetesting
  implicit none
  Type(MLIFEPatchDesc), intent(in) :: patch
  integer, intent(in) :: nIter
  logical, intent(in) :: doCheckpoint
  integer, allocatable, intent(inout) :: m1(:,:), m2(:,:)
  interface
     function MLIFE_nextstate(matrix, lni, lnj, row, col)
       integer, allocatable, intent(in) :: matrix(:,:)
       integer row, col, lni, lnj
       integer MLIFE_nextstate
     end function MLIFE_nextstate
  end interface
  interface
     function exchangeInit(patch, m1, m2)
       use mlife2dtypes
       Type(MLIFEPatchDesc), intent(in) :: patch
       integer, allocatable, intent(inout) :: m1(:,:), m2(:,:)
       integer exchangeInit
     end function exchangeInit
  end interface
  interface
     function exchangeEnd()
       use mlife2dtypes
       integer exchangeEnd
     end function exchangeEnd
  end interface
  interface
     function exchange(patch, matrix, timedata, phase)
       use mlife2dtypes
       use mpi
       Type(MLIFEPatchDesc) :: patch
       integer, allocatable, intent(inout) :: matrix(:,:)
       Type(MLIFETiming)    :: timedata
       integer, intent(in)  :: phase
       integer              :: exchange
     end function exchange
  end interface
  Type(MLIFETiming), intent(out) :: timedata

  double precision t1, t2, t3
  integer i, j, k, ierr
  integer LCols, LRows

  LCols = patch%lnj
  LRows = patch%lni

  ! Initialize the timedata
  timedata%packtime   = 0.0
  timedata%unpacktime = 0.0
  timedata%exchtime   = 0.0
  timedata%itertime   = 0.0

  ! Initialize mesh
  !print *, 'init: ', lbound(m1,1), ubound(m1,1), lbound(m1,2), ubound(m1,2)
  if (.not.allocated(m1)) then
     print *, 'Panic: m1 not allocated'
  endif
  call MLIFE_InitLocalMesh(patch, m1, m2)
  if (doCheckpoint) then
     call MLIFEIO_Checkpoint(patch, m1, 0, MPI_INFO_NULL)
  endif

  ! Initialize exchange
  ierr = exchangeInit(patch, m1, m2)

  ! Time the iterations
  t1 = 0

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  t2 = MPI_Wtime()
  do k=0, nIter, 2
     t3 = MPI_Wtime()
     ierr = exchange(patch, m1, timedata, 1)
     t1 = t1 + MPI_Wtime() - t3

     ! calculate new state for all non-boundary elements
     ! Change this loop to compute a different solution (e.g., an
     ! explicit method for the wave or heat equation or the matrix-vector
     ! product for an implicit method.  Also change the second loop below.
     if (testing >=0) then
     do i=1, LRows
        do j=1, LCols
           m2(i,j) = m1(i-1,j)
        end do
     end do
     else
     do i=1, LRows
        do j=1, LCols
           m2(i,j) = MLIFE_nextstate(m1, LRows, LCols, i, j)
        end do
     end do
     endif

     if (doCheckpoint) then
        call MLIFEIO_Checkpoint(patch, m1, k+1, MPI_INFO_NULL)
     end if

     if (k+1 > nIter) then
        exit   ! Don't do the second half of the 2-step advance
     endif

     ! Do a second iteration using m2 instead of m1, so that we don't
     ! have to swap pointers or copy arrays

     t3 = MPI_Wtime()
     ierr = exchange(patch, m2, timedata, 2)
     t1 = t1 + MPI_Wtime() - t3

     ! calculate new state for all non-boundary elements
     if (testing >=0) then
     do i=1, LRows
        do j=1, LCols
           m1(i,j) = m2(i-1,j)
        end do
     end do
     else
     do i=1, LRows
        do j=1, LCols
           m1(i,j) = MLIFE_nextstate(m2, LRows, LCols, i, j)
        end do
     end do
     endif

     if (doCheckpoint) then
        call MLIFEIO_Checkpoint(patch, m2, k+2, MPI_INFO_NULL)
     end if
  enddo
  t2 = MPI_Wtime() - t2

  ierr = exchangeEnd()

  ! Get the maximum time over all involved processes
  timedata%packtime   = timedata%packtime / nIter;
  timedata%unpacktime = timedata%unpacktime / nIter;
  timedata%exchtime   = t1 / nIter;
  timedata%itertime   = t2 / nIter;
  call MPI_Allreduce(MPI_IN_PLACE, timedata, 4, &
       MPI_DOUBLE, MPI_MAX, patch%comm, ierr)

end subroutine MLIFE_TimeIterations
