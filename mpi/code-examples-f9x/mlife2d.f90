module mlifetesting
  integer testing   ! set this to enable testing rather than the life game
end module mlifetesting

program main
  use mpi
  use mlife2dtypes
  use mlife2dfuncs
  use mlife2dio
  use mlifetesting
  implicit none
  Type(MLIFEOptions) :: options
  Type(MLIFEPatchDesc) :: patch
  Type(MLIFETiming) :: timedataPt2pt, timedataPt2pt9, timedataPt2ptSnd,     &
       timedataPt2ptUv, timedataFence
!  Type(MLIFETiming) :: timedataExample
  integer rank, provided, ierr, commtype, kk
  logical verbose, doCheckpoint
  integer maxcomtype, maxtestrun
  character(20) commstr
  integer, allocatable :: m1(:,:), m2(:,:)

  interface
     subroutine MLIFE_TimeIterations(patch, nIter, doCheck, m1, m2, &
          exchangeInit, exchange, exchangeEnd, timedata)
       use mlife2dtypes
       Type(MLIFEPatchDesc) :: patch
       integer, intent(in)  :: nIter
       logical, intent(in)  :: doCheck
       integer, allocatable, intent(inout) :: m1(:,:), m2(:,:)
       interface
         function exchangeInit(patch, m1, m2)
           use mlife2dtypes
           Type(MLIFEPatchDesc), intent(in) :: patch
           integer, allocatable, intent(inout) :: m1(:,:), m2(:,:)
           integer exchangeInit
         end function exchangeInit
         function exchangeEnd()
           integer exchangeEnd
         end function exchangeEnd
         function exchange(patch, matrix, timedata, phase)
           use mlife2dtypes
           Type(MLIFEPatchDesc) :: patch
           integer, allocatable, intent(inout) :: matrix(:,:)
           Type(MLIFETiming)    :: timedata
           integer, intent(in)  :: phase
           integer exchange
         end function exchange
       end interface
       Type(MLIFETiming)    :: timedata
     end subroutine MLIFE_TimeIterations
  end interface

  call MPI_Init_thread(MPI_THREAD_SINGLE, provided, ierr)
  call MLIFE_ParseArgs(options)
  verbose      = options%verbose
  doCheckpoint = options%doIO
  testing      = options%testing

  call MLIFEIO_Init(MPI_COMM_WORLD)

  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

  if (testing >= 0) then
     maxcomtype = 0
     maxtestrun = 0
  else
     maxcomtype = 1
     maxtestrun = 1
  endif

  do commtype=0,maxcomtype
     select case (commtype)
     case(0)
        call MLIFE_PatchCreateProcessMesh(options, patch)
        commstr = "COMM_WORLD"
     case(1)
        call MLIFE_PatchCreateProcessMeshWithCart(options, patch)
        commstr = "CART_CREATE"
     case default
        print *, "Internal error - commtype = ", commtype
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
     end select
     call MLIFE_PatchCreateDataMeshDesc(options, patch)

     if (verbose) then
        print *, "[", rank, "] Proc with ", commstr, " is [", &
             patch%patchI, ",", patch%patchJ, "], in [", &
             patch%pNI, ",", patch%pNJ, ", mesh [", &
             patch%gI, ",", patch%gI+patch%lni-1, "]x[", &
             patch%gJ, ",", patch%gJ+patch%lnj-1, "] within [", &
             patch%gNI, ",", patch%gnJ, "]"
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
     endif

     call MLIFE_AllocateLocalMesh(patch, m1, m2)

     ! Run these tests multiple times.  In many cases, the first
     ! run will take significantly longer, as library setup actions
     ! take place.  To illustrate this, write out the timing data for
     ! each run, so that the first and second times can be compared
     do kk=0, maxtestrun
        ! For each communication approach, perform these steps
	!       1: Initialize the mesh (initial data)
	!       2: Initialize the exchange
	!       3: Time the iterations, including any communication steps.
        if (verbose) print *, 'About to do pt2pt'
        call MLIFE_TimeIterations(patch, options%nIter, doCheckpoint, &
                                  m1, m2, &
                                  MLIFE_ExchangeInitPt2pt, &
                                  MLIFE_ExchangePt2pt, &
                                  MLIFE_ExchangeEndPt2pt, &
                                  timedataPt2pt)
        ! Uncomment the next line and the matching endif to use the -c option
        ! if (.false.) then
        if (verbose) print *, 'About to do pt2ptsnd'
        call MLIFE_TimeIterations(patch, options%nIter, doCheckpoint, &
                                  m1, m2, &
                                  MLIFE_ExchangeInitPt2ptSnd, &
                                  MLIFE_ExchangePt2ptSnd, &
                                  MLIFE_ExchangeEndPt2ptSnd, &
                                  timedataPt2ptSnd)
        if (verbose) print *, 'About to do pt2ptuv'
        call MLIFE_TimeIterations(patch, options%nIter, doCheckpoint, &
                                  m1, m2, &
                                  MLIFE_ExchangeInitPt2ptUV, &
                                  MLIFE_ExchangePt2ptUV, &
                                  MLIFE_ExchangeEndPt2ptUV, &
                                  timedataPt2ptUv)
        if (verbose) print *, 'About to do fence'
        call MLIFE_TimeIterations(patch, options%nIter, doCheckpoint, &
                                  m1, m2, &
                                  MLIFE_ExchangeInitFence, &
                                  MLIFE_ExchangeFence, &
                                  MLIFE_ExchangeEndFence, &
                                  timedataFence)
        if (verbose) print *, 'About to do pt2pt9'
        call MLIFE_TimeIterations(patch, options%nIter, doCheckpoint, &
                                  m1, m2, &
                                  MLIFE_ExchangeInitPt2pt9, &
                                  MLIFE_ExchangePt2pt9, &
                                  MLIFE_ExchangeEndPt2pt9, &
                                  timedataPt2pt9);

!        if (verbose) print *, 'About to do example'
!        call MLIFE_TimeIterations(patch, options%nIter, doCheckpoint, &
!                                  m1, m2, &
!                                  MLIFE_ExchangeInitExample, &
!                                  MLIFE_ExchangeExample, &
!                                  MLIFE_ExchangeEndExample, &
!                                  timedataExample);

        ! Print the total time taken
        if (rank == 0) then
           ! This allows meshes upto 99999 elements on each size, on
           ! process arrays that are upto 9999 processes in each dimension
           print '("Mesh size [",I5,",",I5,"] on process array [",I4,",",I4,"] for ",A)', &
                patch%gNI, patch%gNJ, patch%pNI, patch%pNJ, commstr
           print *, 'Exchange type	Per iter	Exchange	Packtime	Unpacktime'
           print '("[",I4,"] Pt2pt: 	", 1PE12.3, "	", 1PE12.3)', rank, &
                        timedataPt2pt%itertime, &
                        timedataPt2pt%exchtime
           print '("[",I4,"] Pt2ptSnd: 	", 1PE12.3, "	", 1PE12.3)', rank, &
                        timedataPt2ptSnd%itertime, &
                        timedataPt2ptSnd%exchtime
           print '("[",I4,"] Pt2ptUv: 	", 1PE12.3, "	", 1PE12.3, "	", 1PE12.3, "	", 1PE12.3)', rank, &
                        timedataPt2ptUv%itertime, &
                        timedataPt2ptUv%exchtime, &
                        timedataPt2ptUv%packtime, &
                        timedataPt2ptUv%unpacktime
           print '("[",I4,"] Fence: 	", 1PE12.3, "	", 1PE12.3)', rank, &
                        timedataFence%itertime, &
                        timedataFence%exchtime
           print '("[",I4,"] Pt2pt9: 	", 1PE12.3, "	", 1PE12.3)', rank, &
                        timedataPt2pt9%itertime, &
                        timedataPt2pt9%exchtime
!           print '("[",I4,"] Example: 	", 1PE12.3, "	", 1PE12.3)', rank, &
!                        timedataExample%itertime, &
!                        timedataExample%exchtime
        endif
     enddo
     if (verbose) print *, 'About to free mesh'
     call MLIFE_FreeLocalMesh(patch, m1, m2)
  enddo
  call MLIFEIO_Finalize()
  call MPI_Finalize(ierr)
end program main
