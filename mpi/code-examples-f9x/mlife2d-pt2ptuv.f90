!
!  (C) 2004 by University of Chicago.
!      See COPYRIGHT in top-level directory.
!

module mlife2dpt2ptuv
    integer, dimension(:), allocatable :: usbuf, dsbuf, urbuf, drbuf
    integer                            :: buflen
end module mlife2dpt2ptuv

function MLIFE_ExchangeInitPt2ptUV(patch, m1, m2)
  use mlife2dtypes
  use mlife2dpt2ptuv
  implicit none
  Type(MLIFEPatchDesc) :: patch
  integer MLIFE_ExchangeInitPt2ptUV
  integer, allocatable, intent(inout) :: m1(:,:), m2(:,:)

  allocate(usbuf(0:patch%lnj+1))
  allocate(dsbuf(0:patch%lnj+1))
  allocate(urbuf(0:patch%lnj+1))
  allocate(drbuf(0:patch%lnj+1))
  buflen = patch%lnj+2

  MLIFE_ExchangeInitPt2ptUV = 0
end function MLIFE_ExchangeInitPt2ptUV

function MLIFE_ExchangeEndPt2ptUV()
  use mlife2dtypes
  use mlife2dpt2ptuv
  integer MLIFE_ExchangeEndPt2ptUV

  deallocate(usbuf,dsbuf,urbuf,drbuf)
  MLIFE_ExchangeEndPt2ptUV = 0
end function MLIFE_ExchangeEndPt2ptUV

function MLIFE_ExchangePt2ptUV(patch, matrix, timedata, phase)
  use mlife2dtypes
  use mlife2dpt2ptuv
  use mpi
 implicit none
  Type(MLIFEPatchDesc) :: patch
  integer, allocatable, intent(inout) :: matrix(:,:)
  Type(MLIFETiming)    :: timedata
  integer, intent(in)  :: phase
  integer MLIFE_ExchangePt2ptUV
  !
  integer  :: reqs(4)
  integer  :: comm, ierr
  integer  :: LRows, LCols, j
  double precision :: t1, t2

  comm = patch%comm
  LRows = patch%lni
  LCols = patch%lnj

  ! Send and receive boundary information

  ! first, move the left, right edges */
  call MPI_Isend(matrix(1,1), LRows, MPI_INTEGER, patch%left, &
       0, comm, reqs(1), ierr)
  call MPI_Irecv(matrix(1,0), LRows, MPI_INTEGER, patch%left, &
       0, comm, reqs(2), ierr)
  call MPI_Isend(matrix(1,LCols), LRows, MPI_INTEGER, patch%right, &
       0, comm, reqs(3), ierr)
  call MPI_Irecv(matrix(0,LCols+1), LRows, MPI_INTEGER, patch%right, &
       0, comm, reqs(4), ierr)
  ! We need to wait on these for the trick that we use to move
  ! the diagonal terms to work
  call MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE, ierr)

  ! move the top, bottom edges (including diagonals)

  call MPI_Irecv(urbuf, buflen, MPI_INTEGER, patch%up, 0, comm, reqs(2), ierr)
  call MPI_Irecv(drbuf, buflen, MPI_INTEGER, &
       patch%down, 0, comm, reqs(4), ierr)
  ! User pack of the buffers to send
  t1 = MPI_Wtime();
  if (patch%up .ne. MPI_PROC_NULL) then
     do j=0,LCols+1
        usbuf(j) = matrix(1,j)
     enddo
  endif
  if (patch%down .ne. MPI_PROC_NULL) then
     do j=0, LCols+1
        dsbuf(j) = matrix(LRows,j)
     enddo
  endif
  t1 = MPI_Wtime() - t1
  call MPI_Isend(usbuf, buflen, MPI_INTEGER, patch%up, 0, comm, reqs(1), ierr)
  call MPI_Isend(dsbuf, buflen, MPI_INTEGER, &
       patch%down, 0, comm, reqs(3), ierr)
  ! We need to wait on these for the trick that we use to move
  ! the diagonal terms to work
  call MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE, ierr)
  ! Need to unpack
  t2 = MPI_Wtime();
  if (patch%up .ne. MPI_PROC_NULL) then
     do j=0, LCols+1
        matrix(0,j)       = urbuf(j)
     enddo
  endif
  if (patch%down .ne. MPI_PROC_NULL) then
     do j=0, LCols+1
        matrix(LRows+1,j) = drbuf(j)
     enddo
  endif
  t2 = MPI_Wtime() - t2

  timedata%packtime = timedata%packtime + t1
  timedata%unpacktime = timedata%unpacktime + t2
  MLIFE_ExchangePt2ptUV = MPI_SUCCESS

  return
end function MLIFE_ExchangePt2ptUV
