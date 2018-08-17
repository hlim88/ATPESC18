!
!  (C) 2004 by University of Chicago.
!      See COPYRIGHT in top-level directory.
!
function MLIFE_ExchangeInitPt2pt9(patch, m1, m2)
  use mlife2dtypes
  Type(MLIFEPatchDesc) :: patch
  integer MLIFE_ExchangeInitPt2pt9
  integer, allocatable, intent(inout) :: m1(:,:), m2(:,:)
  MLIFE_ExchangeInitPt2pt9 = 0
end function MLIFE_ExchangeInitPt2pt9

function MLIFE_ExchangeEndPt2pt9()
  use mlife2dtypes
  integer MLIFE_ExchangeEndPt2pt9
  MLIFE_ExchangeEndPt2pt9 = 0
end function MLIFE_ExchangeEndPt2pt9

function MLIFE_ExchangePt2pt9(patch, matrix, timedata, phase)
  use mlife2dtypes
  use mpi
  implicit none
  Type(MLIFEPatchDesc) :: patch
  integer, allocatable, intent(inout) :: matrix(:,:)
  Type(MLIFETiming)    :: timedata
  integer, intent(in)  :: phase
  integer MLIFE_ExchangePt2pt9
  !
  integer  :: reqs(16)
  integer  :: comm
  integer  :: dtype = MPI_DATATYPE_NULL  ! Sets the initial value, not on each entry
  integer  :: LRows, LCols, ierr

  comm  = patch%comm
  LRows = patch%lni
  LCols = patch%lnj

  ! Send and receive boundary information

  if (dtype == MPI_DATATYPE_NULL) then
     call MPI_Type_vector(LCols, 1, LRows+2, MPI_INTEGER, dtype, ierr )
     call MPI_Type_commit(dtype, ierr)
  endif

  ! first, move the left, right edges */
  call MPI_Isend(matrix(1,1), LRows, MPI_INTEGER, patch%left, &
       0, comm, reqs(1), ierr)
  call MPI_Irecv(matrix(1,0), LRows, MPI_INTEGER, patch%left, &
       0, comm, reqs(2), ierr)
  call MPI_Isend(matrix(1,LCols), LRows, MPI_INTEGER, patch%right, &
       0, comm, reqs(3), ierr)
  call MPI_Irecv(matrix(0,LCols+1), LRows, MPI_INTEGER, patch%right, &
       0, comm, reqs(4), ierr)

  ! move the top, bottom edges
  ! move the top, bottom edges (including diagonals)
  call MPI_Isend(matrix(1,1), 1, dtype, patch%up, 0, comm, reqs(5), ierr)
  call MPI_Irecv(matrix(0,1), 1, dtype, patch%up, 0, comm, reqs(6), ierr)
  call MPI_Isend(matrix(Lrows,1), 1, dtype, patch%down, 0, comm, reqs(7), ierr)
  call MPI_Irecv(matrix(Lrows+1,1), 1, dtype, patch%down, &
       0, comm, reqs(8), ierr)

  ! Move the diagonals.  If we know that we have an eager send limit
  ! of at least one word, these could send/recv instead of isend/irecv.
  call MPI_Isend( matrix(1,1), 1, MPI_INTEGER, patch%ul, 0, comm, reqs(9), &
       ierr )
  call MPI_Isend( matrix(LRows,1), 1, MPI_INTEGER, patch%ll, 0, comm, &
       reqs(10), ierr)
  call MPI_Isend( matrix(1,LCols), 1, MPI_INTEGER, patch%ur, 0, comm, &
       reqs(11), ierr)
  call MPI_Isend( matrix(LRows,LCols), 1, MPI_INTEGER, patch%lr, 0, comm, &
       reqs(12), ierr)
  call MPI_Irecv( matrix(0,0), 1, MPI_INTEGER, patch%ul, 0, comm, reqs(13), &
       ierr)
  call MPI_Irecv( matrix(0,LCols+1), 1, MPI_INTEGER, patch%ur, 0, comm, &
       reqs(14), ierr)
  call MPI_Irecv( matrix(LRows+1,0), 1, MPI_INTEGER, patch%ll, 0, comm, &
       reqs(15), ierr)
  call MPI_Irecv( matrix(LRows+1,LCols+1), 1, MPI_INTEGER, patch%lr, 0, comm,&
       reqs(16), ierr)

  call MPI_Waitall(16, reqs, MPI_STATUSES_IGNORE, ierr)

  MLIFE_ExchangePt2pt9 = MPI_SUCCESS
  return
end function MLIFE_ExchangePt2pt9
