module mlife2dio
  integer, private :: mlifeio_comm
  contains

!
subroutine MLIFEIO_Init(comm)
  use mpi
  implicit none
  integer comm, ierr

  call MPI_Comm_dup(comm, mlifeio_comm, ierr)
end subroutine MLIFEIO_Init

subroutine MLIFEIO_Finalize()
  use mpi
  implicit none
  integer ierr
  if (mlifeio_comm .ne. MPI_COMM_NULL) then
     call MPI_Comm_free(mlifeio_comm, ierr)
  endif
end subroutine MLIFEIO_Finalize

subroutine MLIFEIO_Restart(patch, matrix, info)
  use mlife2dtypes
end subroutine MLIFEIO_Restart

function MLIFEIO_Can_restart()
  logical MLIFEIO_Can_restart
  MLIFEIO_Can_restart = .false.
end function MLIFEIO_Can_restart

subroutine MLIFEIO_Checkpoint(patch, matrix, iter, info)
  use mlife2dtypes
  use msleep
  use mpi
  implicit none
  Type(MLIFEPatchDesc)             :: patch
  integer, allocatable, intent(in) :: matrix(:,:)
  integer                          :: iter
  integer                          :: info
  integer                          :: mrank, nprocs, ierr
  integer MAX_SIZE
  parameter (MAX_SIZE=256)
  character(MAX_SIZE) cbuf
  integer  buf(0:MAX_SIZE)
  integer  np, row, i, j, npcols, maxcols

  call MPI_Comm_size(mlifeio_comm, nprocs, ierr)
  call MPI_Comm_rank(mlifeio_comm, mrank, ierr)

  ! Check first that the data is not too big for this output
  call MPI_Allreduce(patch%lnj, maxcols, 1, MPI_INTEGER, MPI_MAX, &
       mlifeio_comm, ierr)
  if (maxcols + 2 > MAX_SIZE .or. patch%gNJ > MAX_SIZE) then
     if (mrank == 0) then
        print *, "Maximum width(y) exceeded for stdout output.  Max is ", &
             MAX_SIZE, ", this run has a max of ", patch%gNJ
     endif
     call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  endif

  if (mrank == 0) then
     npcols = patch%pNJ

     print *, "[H[2J# Iteration ", iter
     do row=1, patch%gNI
        ! Clear the cbuf
        do i=1, patch%gNJ
           cbuf(i:i) = ' '
        enddo
        ! We know how many processes there are in each row
        np = 0
        ! Does this process contribute to this row?
        if (row >= patch%gI .and. row < patch%gI + patch%lni) then
           do i=1, patch%lnj
              if (matrix(row-patch%gI+1,i) == 1) then
                 cbuf(i+patch%gJ-1:i+patch%gJ-1) = '*'
              else
                 cbuf(i+patch%gJ-1:i+patch%gJ-1) = ' '
              endif
           enddo
           np = np + 1
        endif
        ! Receive from all of the processes for this row
        do while (np < npcols)
           call MPI_Recv(buf, MAX_SIZE, MPI_INTEGER, MPI_ANY_SOURCE, row, &
                mlifeio_comm, MPI_STATUS_IGNORE, ierr)
           ! For each entry in buf that is set, set the
	   ! corresponding element in cbuf.  buf(0) is the
	   ! first col index, buf(1) is the number of elements
           do i=1, buf(1)
              if (buf(i+1) == 1) then
                 cbuf(buf(0)-1+i:buf(0)-1+i) = '*'
              else
                 cbuf(buf(0)-1+i:buf(0)-1+i) = ' '
              endif
           enddo
           np = np + 1
        enddo
        ! no need to null cbuf !        cbuf[patch->gNJ] = 0;
        ! The odd characters are commands to an xterm/vt100 window
        !printf("[%03d;%03dH%3d: %s", row+1, 1, row, cbuf );
        print *, cbuf(1:patch%gNJ)
     enddo
     ! The odd characters are commands to an xterm/vt100 window
     !printf("[%03d;%03dH", patch%gNJ+2, 1);
     !fflush(stdout);
  else
     buf(0) = patch%gJ
     buf(1) = patch%lnj
     do i=1, patch%lni
        row = (i-1) + patch%gI
        do j=1,patch%lnj
           buf(1+j) = matrix(i,j)
        enddo
        call MPI_Send(buf, patch%lnj + 2, MPI_INTEGER, 0, row, mlifeio_comm, &
             ierr)
     enddo
  endif

  call MLIFEIO_msleep(250)

end subroutine MLIFEIO_Checkpoint

subroutine mprint(iter, lni, lnj, matrix)
  use msleep
  integer matrix(0:lni+1,0:lnj+1)
  integer lni, lnj, iter, i, j
  character(100) cbuf
  character(150) outbuf
  integer, parameter :: BORN = 1

  outbuf(1:19) = "[H[2J# Iteration "
  write(outbuf(20:24), "(I5)") iter
  print *, outbuf(1:29)

  outbuf(1:) = " "

  do i=1, lni
     cbuf(1:) = " "
     do j=1, lnj
        if (matrix(i,j) == BORN) then
           cbuf(j:j) = "*"
        else
           cbuf(j:j) = " "
        endif
     end do
     outbuf(1:2)   = "["
     outbuf(6:6)   = ";"
     outbuf(10:10) = "H"
     write(outbuf(3:5),"(I3.3)") i+1
     write(outbuf(7:9),"(I3.3)") j
     outbuf(11:11+lnj) = cbuf(1:lnj)
     print *, outbuf
  enddo
  call MLIFEIO_msleep(1)
end subroutine mprint

end module mlife2dio
