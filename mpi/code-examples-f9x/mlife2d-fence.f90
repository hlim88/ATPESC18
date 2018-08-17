!
!  (C) 2004 by University of Chicago.
!      See COPYRIGHT in top-level directory.
!

module mlife2dfence
    integer wins(2)
    integer above_LRows, below_LRows, left_LCols
end module mlife2dfence

function MLIFE_ExchangeInitFence(patch, m1, m2)
  use mpi
  use mlife2dtypes
  use mlife2dfence
  implicit none
  Type(MLIFEPatchDesc) :: patch
  integer MLIFE_ExchangeInitFence
  integer, allocatable, intent(inout) :: m1(:,:), m2(:,:)
  integer nprocs, sizeofint
  integer LRows, LCols, ierr
  integer (kind=MPI_ADDRESS_KIND) winsize

  LRows = patch%lni;
  LCols = patch%lnj;

! Create windows
! We could also set the window info to include no_locks=true
  call MPI_Type_size(MPI_INTEGER, sizeofint, ierr)
  ! We can't depend on the MPI module correctly casting the window size
  winsize = (LRows+2)*(LCols+2)*sizeofint
  call MPI_Win_create(m1, winsize, &
       sizeofint, MPI_INFO_NULL, patch%comm, wins(1), ierr)

  call MPI_Win_create(m2, winsize, &
       sizeofint, MPI_INFO_NULL, patch%comm, wins(2), ierr)

!   for one-sided communication, we need to know the number of
!   local rows in rank above and the number of local
!   columns in rank left and right in order to do
!   the puts into the right locations in memory.

  call MPI_Comm_size(patch%comm, nprocs, ierr)

  if (patch%up == MPI_PROC_NULL) then
     above_LRows = 0
  else
     call MPI_Recv(above_LRows, 1, MPI_INTEGER, patch%up, 0, patch%comm, &
          MPI_STATUS_IGNORE, ierr)
  endif
  if (patch%down /= MPI_PROC_NULL) then
     ! Could just call MPI_Send with a dest of patch%down, since a
     ! send to MPI_PROC_NULL is a no-op
     call MPI_Send(patch%lni, 1, MPI_INTEGER, patch%down, 0, patch%comm, &
          ierr)
  endif

  if (patch%left == MPI_PROC_NULL) then
     left_LCols = 0
  else
     call MPI_Recv(left_LCols, 1, MPI_INTEGER, patch%left, 0, patch%comm, &
          MPI_STATUS_IGNORE, ierr)
  endif
  if (patch%right /= MPI_PROC_NULL) then
     call MPI_Send(patch%lnj, 1, MPI_INT, patch%right, 0, patch%comm, ierr)
  endif

  if (patch%down == MPI_PROC_NULL) then
     below_LRows = 0
  else
     call MPI_Recv(below_LRows, 1, MPI_INTEGER, patch%down, 0, patch%comm, &
          MPI_STATUS_IGNORE, ierr)
  endif
  if (patch%up /= MPI_PROC_NULL) then
     call MPI_Send(patch%lni, 1, MPI_INTEGER, patch%up, 0, patch%comm, ierr)
  endif

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  MLIFE_ExchangeInitFence = 0
end function MLIFE_ExchangeInitFence

function MLIFE_ExchangeEndFence()
  use mlife2dtypes
  use mlife2dfence
  integer MLIFE_ExchangeEndFence

  call MPI_Win_free(wins(1), ierr)
  call MPI_Win_free(wins(2), ierr)
  MLIFE_ExchangeEndFence = 0
end function MLIFE_ExchangeEndFence

function MLIFE_ExchangeFence(patch, matrix, timedata, phase)
  use mlife2dtypes
  use mlife2dfence
  use mpi
  implicit none
  Type(MLIFEPatchDesc) :: patch
  integer, allocatable, intent(inout) :: matrix(:,:)
  Type(MLIFETiming)    :: timedata
  integer, intent(in)  :: phase
  integer MLIFE_ExchangeFence
  !
  integer  :: comm
  integer  :: mytype = MPI_DATATYPE_NULL  ! Sets the initial value, not on each entry
  integer  :: up_type = MPI_DATATYPE_NULL
  integer  :: down_type = MPI_DATATYPE_NULL
  integer  :: LRows, LCols, ierr
  integer(kind=MPI_ADDRESS_KIND) :: disp

  comm  = patch%comm
  LRows = patch%lni
  LCols = patch%lnj

  ! create datatype if not already created
  if (mytype == MPI_DATATYPE_NULL) then
     call MPI_Type_vector(LCols+2, 1, LRows+2, MPI_INTEGER, mytype, ierr)
     call MPI_Type_commit(mytype, ierr)
  endif
  if (up_type == MPI_DATATYPE_NULL) then
     call MPI_Type_vector(LCols+2, 1, above_Lrows+2, MPI_INTEGER, up_type, &
          ierr)
     call MPI_Type_commit(up_type, ierr)
  endif
  if (down_type == MPI_DATATYPE_NULL) then
     call MPI_Type_vector(LCols+2, 1, below_LRows+2, MPI_INTEGER, &
          down_type, ierr)
     call MPI_Type_commit(down_type, ierr)
  endif

  call MPI_Win_fence(MPI_MODE_NOPRECEDE, wins(phase), ierr)

  ! first put the left, right edges
  disp = (left_LCols + 1) * (LRows + 2) + 1
  call MPI_Put(matrix(1,1), LRows, MPI_INTEGER, patch%left, disp, LRows, &
       MPI_INTEGER, wins(phase), ierr)

  disp = 1
  call MPI_Put(matrix(1,LCols), LRows, MPI_INTEGER, patch%right, disp, LRows, &
       MPI_INTEGER, wins(phase), ierr)

  ! Complete the right/left transfers for the diagonal trick
  call MPI_Win_fence( 0, wins(phase), ierr)

  ! now put the top, bottom edges (including the diagonal points)
  disp = (above_LRows+1)
  call MPI_Put(matrix(1,0), 1, mytype, patch%up, &
       disp, 1, up_type, wins(phase), ierr)

  disp = 0
  call MPI_Put(matrix(LRows,0), 1, mytype, patch%down, disp, &
       1, down_type, wins(phase), ierr)

  call MPI_Win_fence(MPI_MODE_NOSTORE + MPI_MODE_NOPUT + &
       MPI_MODE_NOSUCCEED, wins(phase), ierr)

  MLIFE_ExchangeFence = MPI_SUCCESS
  return
end function MLIFE_ExchangeFence
