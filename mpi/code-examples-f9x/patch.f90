! SLIDE: 2D Life Code Walkthrough */
!
!  (C) 2013 by University of Chicago.
!      See COPYRIGHT in top-level directory.
!

subroutine MLIFE_PatchCreateProcessMesh(options, patch)
  use mlife2dtypes
  use mpi
  implicit none
  Type(MLIFEOptions) :: options
  Type(MLIFEPatchDesc) :: patch
  integer dims(2), up, down, left, right, ul, ur, ll, lr
  integer prow, pcol, nprocs, mrank, ierr

  dims(1) = options%pNI
  dims(2) = options%pNJ;

  patch%comm = MPI_COMM_WORLD
  call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, mrank, ierr)

  call MPI_Dims_create(nprocs, 2, dims, ierr)

  patch%pNI = dims(1)
  patch%pNJ = dims(2)

  patch%gNI = options%gNI
  patch%gNJ = options%gNJ

  ! compute the cartesian coords of this process; number across
  ! rows changing column by 1 changes rank by 1)
  prow = mrank / dims(2)
  pcol = mod(mrank,dims(2))
  patch%patchI = prow
  patch%patchJ = pcol

  ! Compute the neighbors
  left  = MPI_PROC_NULL
  right = MPI_PROC_NULL
  up    = MPI_PROC_NULL
  down  = MPI_PROC_NULL
  ul    = MPI_PROC_NULL
  ur    = MPI_PROC_NULL
  ll    = MPI_PROC_NULL
  lr    = MPI_PROC_NULL
  if (prow > 0) then
     up   = mrank - dims(2)
     if (pcol > 0) ul = up - 1
     if (pcol < dims(2) - 1) ur = up + 1
  endif
  if (pcol > 0) then
     left = mrank - 1
  endif
  if (prow < dims(1)-1) then
     down = mrank + dims(2);
     if (pcol > 0) ll = down - 1;
     if (pcol < dims(2) - 1) lr = down + 1;
  endif
  if (pcol < dims(2)-1) then
     right = mrank + 1;
  endif
  patch%left  = left;
  patch%right = right;
  patch%up    = up;
  patch%down  = down;
  patch%ul    = ul;
  patch%ur    = ur;
  patch%ll    = ll;
  patch%lr    = lr;

end subroutine MLIFE_PatchCreateProcessMesh

subroutine MLIFE_PatchCreateProcessMeshWithCart(options, patch)
  use mlife2dtypes
  use mpi
  implicit none
  Type(MLIFEOptions) :: options
  Type(MLIFEPatchDesc) :: patch
  integer dims(2), coords(2)
  logical periods(2)
  integer up, down, left, right, ul, ur, ll, lr, ierr
  integer nprocs, mrank

  dims(1) = options%pNI;
  dims(2) = options%pNJ;

  call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

  call MPI_Dims_create(nprocs, 2, dims, ierr)

  patch%pNI = dims(1);
  patch%pNJ = dims(2);

  patch%gNI = options%gNI;
  patch%gNJ = options%gNJ;

  ! Create a Cartesian communicator using the recommended (by dims_create)
  ! sizes
  periods(1) = .false.
  periods(2) = .false.
  call MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, .true., &
       patch%comm, ierr)

  ! The ordering of processes, relative to the rank in the new
  ! communicator, is defined by the MPI standard.  However, there are
  ! useful convenience functions

  call MPI_Comm_rank(patch%comm, mrank, ierr)
  call MPI_Cart_coords(patch%comm, mrank, 2, coords, ierr)
  ! compute the cartesian coords of this process; number across
  ! rows changing column by 1 changes rank by 1)
  patch%patchI = coords(1)  ! prow
  patch%patchJ = coords(2)  ! pcol

  ! Get the neighbors.  Shift can be used for the coordinate directions
  call MPI_Cart_shift(patch%comm, 0, 1, up, down, ierr)
  call MPI_Cart_shift(patch%comm, 1, 1, left, right, ierr)

  ! For the diagonal processes, we can either:
  ! 1. Use the defined layout to compute them
  ! 2. Communicate with the neighbors in the coordinate directions, who
  !    know those ranks (e.g., as the up neighbor for the ranks of the up
  !    neighbors left and right neighbors).
  ul    = MPI_PROC_NULL
  ur    = MPI_PROC_NULL
  ll    = MPI_PROC_NULL
  lr    = MPI_PROC_NULL
  if (up /= MPI_PROC_NULL) then
     if (left /= MPI_PROC_NULL) ul = up - 1
     if (right /= MPI_PROC_NULL) ur = up + 1
  endif
  if (down /= MPI_PROC_NULL) then
     if (left /= MPI_PROC_NULL) ll = down - 1
     if (right /= MPI_PROC_NULL) lr = down + 1
  endif

  patch%left  = left;
  patch%right = right;
  patch%up    = up;
  patch%down  = down;
  patch%ul    = ul;
  patch%ur    = ur;
  patch%ll    = ll;
  patch%lr    = lr;

end subroutine MLIFE_PatchCreateProcessMeshWithCart


subroutine MLIFE_PatchCreateDataMeshDesc(options, patch)
  use mlife2dtypes
  use mpi
  implicit none
  Type(MLIFEOptions) :: options
  Type(MLIFEPatchDesc) :: patch
  integer firstcol, firstrow, lastcol, lastrow

  ! compute the decomposition of the global mesh.
  ! these are for the "active" part of the mesh, and range from
  ! 1 to GRows by 1 to GCols
  firstcol = 1 + patch%patchJ * (patch%gNJ / patch%pNJ);
  firstrow = 1 + patch%patchI * (patch%gNI / patch%pNI);
  if (patch%patchJ == patch%pNJ - 1) then
     lastcol = patch%gNJ
  else
     lastcol  = 1 + (patch%patchJ + 1) * (patch%gNJ / patch%pNJ) - 1
  endif

  if (patch%patchI == patch%pNI - 1) then
     lastrow = patch%gNI
  else
     lastrow = 1 + (patch%patchI + 1) * (patch%gNI / patch%pNI) - 1
  endif

  patch%gI  = firstrow;
  patch%gJ  = firstcol;
  patch%lni = lastrow - firstrow + 1;
  patch%lnj = lastcol - firstcol + 1;

end subroutine MLIFE_PatchCreateDataMeshDesc

! Allocate a 2-D array; allocate
! local mesh with ghost cells and as a contiguous block so that
! strided access may be used for the left/right edges
!
! For simplicity, all patches have halo cells on all sides, even
! if the process shares a physical boundary.
subroutine MLIFE_AllocateLocalMesh(patch, m1, m2)
  use mlife2dtypes
  use mpi
  implicit none
  Type(MLIFEPatchDesc) :: patch
  integer, allocatable, intent(out) :: m1(:,:), m2(:,:)

  allocate(m1(0:patch%lni+1,0:patch%lnj+1))
  allocate(m2(0:patch%lni+1,0:patch%lnj+1))
end subroutine MLIFE_AllocateLocalMesh

subroutine MLIFE_FreeLocalMesh(patch, m1, m2)
  use mlife2dtypes
  use mpi
  implicit none
  Type(MLIFEPatchDesc) :: patch
  integer, allocatable, intent(inout) :: m1(:,:), m2(:,:)

  deallocate(m1)
  deallocate(m2)
end subroutine MLIFE_FreeLocalMesh

subroutine MLIFE_InitLocalMesh(patch, m1, m2)
  use mlife2dtypes
  use mlifetesting
  implicit none
  Type(MLIFEPatchDesc)  :: patch
  integer, allocatable, intent(inout) :: m1(:,:), m2(:,:)
  integer i, j, lni, lnj, sn
  integer, dimension(:), allocatable :: seed
  real dummy, r

  lni = patch%lni
  lnj = patch%lnj

  ! initialize the boundaries of the life matrix
  do j=0, lnj+1
     m1(0,j)     = DIES
     m1(lni+1,j) = DIES
     m2(0,j)     = DIES
     m2(lni+1,j) = DIES
  enddo
  do i=0, lni+1
     m1(i,0)     = DIES
     m1(i,lnj+1) = DIES
     m2(i,0)     = DIES
     m2(i,lnj+1) = DIES
  enddo

  ! Initialize the life matrix
  if (testing >= 0) then
     do j=0, lnj+1
        do i=0,lni+1
           m1(i,j) = DIES
           m2(i,j) = DIES
        enddo
     enddo
     if (patch%gI == 1) then
        do j=1,lnj
           m1(0,j) = BORN
            m2(0,j) = BORN
        enddo
     endif
  else
     ! Use the Fortran standard random function.  This is likely to give
     ! different results than the C random functions used in the C versions
     ! of this code.
     call random_seed(size = sn)
     allocate(seed(sn))
     do i=1, lni
        seed = IEOR(1000,(i + patch%gI-1)) * (/ (j-1, j=1,sn) /)
        call random_seed(put = seed)
        ! advance to the random number generator to the
        ! first *owned* cell in this row
        do j=1, patch%gJ
           call random_number(dummy)
        enddo

        do j=1, lnj
           call random_number(r)
           if (r .gt. 0.5) then
              m1(i,j) = BORN;
           else
              m1(i,j) = DIES;
           endif
        end do
     end do
     deallocate(seed)
  endif
end subroutine MLIFE_InitLocalMesh
