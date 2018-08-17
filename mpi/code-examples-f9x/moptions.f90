! SLIDE: 2D Life Code Walkthrough */
!
!  (C) 2013 by University of Chicago.
!      See COPYRIGHT in top-level directory.
!

! MLIFEParseArgs
!
! Note: Command line arguments are not guaranteed in the MPI
!       environment to be passed to all processes.  To be
!       portable, we must process on rank 0 and distribute
!       results.
!
subroutine MLIFE_ParseArgs(options)
  use mpi
  use mlife2dtypes
  use mlife2dfuncs, only : MLIFE_Abort
  Type(MLIFEOptions) options
  integer wrank, ierr, blklens(2), dtypes(2), argstype
  integer (kind=MPI_ADDRESS_KIND) displs(2)
  integer i, argc, cchar
  character(32) arg

  call MPI_Comm_rank(MPI_COMM_WORLD, wrank, ierr)
  blklens(1) = 9
  blklens(2) = len(options%prefix)
  displs(1)  = 0
  call MPI_Get_address(options, displs(1), ierr)
  call MPI_Get_address(options%prefix, displs(2), ierr)
  displs(2) = displs(2) - displs(1)
  displs(1) = 0
  dtypes(1) = MPI_INTEGER
  dtypes(2) = MPI_CHARACTER
  call MPI_Type_create_struct(2, blklens, displs, dtypes, argstype, ierr)
  call MPI_Type_commit(argstype, ierr)

  if (wrank == 0) then
     ! Initialize defaults
     options%gNI     = -20
     options%gNJ     = -20
     options%pNI     = 0
     options%pNJ     = 0
     options%nIter   = 1
     options%verbose = .false.
     options%doIO    = .false.
     options%restartIter = -1
     options%prefix  = "mlife"
     options%testing = -1
     argc = command_argument_count()
     i    = 1
     do while (i <= argc)
        call get_command_argument(i, arg)
        arglen = len_trim(arg)
        if (arg(1:1) == "-" .and. len_trim(arg) == 2) then
           cchar = ichar(arg(2:2))
           select case (cchar)
           case (ichar("a"))
              i = i + 1
              call get_command_argument(i, arg)
              read(arg,*) options%pNJ
           case (ichar("b"))
              i = i + 1
              call get_command_argument(i, arg)
              read(arg,*) options%pNI
           case (ichar("c"))
              options%doIO = .true.
           case (ichar("x"))
              i = i + 1
              call get_command_argument(i, arg)
              read(arg,*) options%gNJ
           case (ichar("y"))
              i = i + 1
              call get_command_argument(i, arg)
              read(arg,*) options%gNI
           case (ichar("i"))
              i = i + 1
              call get_command_argument(i, arg)
              read(arg,*) options%nIter
           case (ichar("p"))
              i = i + 1
              call get_command_argument(i, arg)
              options%prefix = arg(1:len_trim(arg))
           case (ichar("r"))
              i = i + 1
              call get_command_argument(i, arg)
              read(arg,*) options%restartIter
           case (ichar("v"))
              options%verbose = .true.
           case (ichar("t"))
              options%testing = 0
           case default
              print *, "-a <pj> - Number of processes in x (j) direction"
              print *, "-b <pi> - Number of processes in y (i) direction"
              print *, "-c      - Enable I/O (checkpoint)"
              print *, "-x <nj> - Size of mesh in x (j) direction"
              print *, "-y <ni> - Size of mesh in y (i) direction"
              print *, "-i <n>  - Number of iterations"
              print *, "-r <i>  - Iteration where restart begins (and read restart file)"
              print *, "-p <pre>- Filename prefix for I/O"
              print *, "-v      - Turn on verbose output"
              print *, "-t      - Testing mode instead of game of life"
              call MLIFE_Abort("")
           end select
        end if
        i = i + 1
     end do
     if (.not. options%doIO) then
        if (options%gNI .lt. 0 .and. options%gNJ .lt. 0) then
           options%gNI = 2048
           options%gNJ = 2048
        endif
     endif
     if (options%gNI .lt. 0) then
        options%gNI = - options%gNI
     endif
     if (options%gNJ .lt. 0) then
        options%gNJ = - options%gNJ
     endif
  end if

  call MPI_Bcast(options, 1, argstype, 0, MPI_COMM_WORLD, ierr)
  call MPI_Type_free(argstype, ierr)
end subroutine MLIFE_ParseArgs

subroutine MLIFE_Abort(str)
  use mpi
  character(*), intent(in) :: str
  print *, "MLIFE Aborting: ", str
  call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
end subroutine MLIFE_Abort
