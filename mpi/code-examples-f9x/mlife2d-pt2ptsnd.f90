!
!  (C) 2004 by University of Chicago.
!      See COPYRIGHT in top-level directory.
!
function MLIFE_ExchangeInitPt2ptSnd(patch, m1, m2)
  use mlife2dtypes
  Type(MLIFEPatchDesc) :: patch
  integer MLIFE_ExchangeInitPt2ptSnd
  integer, allocatable, intent(inout) :: m1(:,:), m2(:,:)
  MLIFE_ExchangeInitPt2ptSnd = 0
end function MLIFE_ExchangeInitPt2ptSnd

function MLIFE_ExchangeEndPt2ptSnd()
  use mlife2dtypes
  integer MLIFE_ExchangeEndPt2ptSnd
  MLIFE_ExchangeEndPt2ptSnd = 0
end function MLIFE_ExchangeEndPt2ptSnd

function MLIFE_ExchangePt2ptSnd(patch, matrix, timedata, phase)
  use mlife2dtypes
  use mpi
  implicit none
  Type(MLIFEPatchDesc) :: patch
  integer, allocatable, intent(inout) :: matrix(:,:)
  Type(MLIFETiming)    :: timedata
  integer, intent(in)  :: phase
  integer MLIFE_ExchangePt2ptSnd
  !
  integer  :: reqs(2)
  integer  :: comm
  integer  :: dtype = MPI_DATATYPE_NULL  ! Sets the initial value, not on each entry
  integer  :: LRows, LCols, ierr

  comm = patch%comm
  LRows = patch%lni
  LCols = patch%lnj

  ! Send and receive boundary information

  if (dtype == MPI_DATATYPE_NULL) then
     call MPI_Type_vector(LCols+2, 1, LRows+2, MPI_INTEGER, dtype, ierr)
     call MPI_Type_commit(dtype, ierr)
  endif

  ! first, move the left, right edges */
  call MPI_Irecv(matrix(1,0), LRows, MPI_INTEGER, patch%left, &
       0, comm, reqs(1), ierr)
  call MPI_Irecv(matrix(0,LCols+1), LRows, MPI_INTEGER, patch%right, &
       0, comm, reqs(2), ierr)
  call MPI_Send(matrix(1,1), LRows, MPI_INTEGER, patch%left, &
       0, comm, ierr)
  call MPI_Send(matrix(1,LCols), LRows, MPI_INTEGER, patch%right, &
       0, comm, ierr)

  ! We need to wait on these for the trick that we use to move
  ! the diagonal terms to work
  call MPI_Waitall(2, reqs, MPI_STATUSES_IGNORE, ierr)

  ! move the top, bottom edges (including diagonals)
  ! move the top, bottom edges (including diagonals)
  call MPI_Irecv(matrix(0,0), 1, dtype, patch%up, 0, comm, reqs(1), ierr)
  call MPI_Irecv(matrix(Lrows+1,0), 1, dtype, patch%down, &
       0, comm, reqs(2), ierr)
  call MPI_Send(matrix(1,0), 1, dtype, patch%up, 0, comm, ierr)
  call MPI_Send(matrix(Lrows,0), 1, dtype, patch%down, 0, comm, ierr)

  call MPI_Waitall(2, reqs, MPI_STATUSES_IGNORE, ierr)

  MLIFE_ExchangePt2ptSnd = MPI_SUCCESS
  return
end function MLIFE_ExchangePt2ptSnd
