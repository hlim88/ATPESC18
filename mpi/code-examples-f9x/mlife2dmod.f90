!
!  (C) 2004 by University of Chicago.
!      See COPYRIGHT in top-level directory.
!
module mlife2dtypes
  ! use mpi
type MLIFEPatchDesc
   integer comm               ! Communicator for processor mesh
   integer patchI, patchJ;    ! (I,J) location of patch in processor mesh
   integer pNI, pNJ;          ! Size of processor mesh
   integer left, right, up, down;  ! Neighbors to this patch: left = (I,J-1),
				   ! up = (I-1,J), etc.  (standard matrix
				   ! ordering)
   integer ul, ur, ll, lr;    ! Upper left, upper right, lower left,
			      ! lower right (needed for one option of the
			      ! 9 point stencil)
   integer gNI, gNJ;          ! Full size of data mesh
   integer gI, gJ;            ! Global I,J for the upper left corner of the
                              ! local patch
   integer lni, lnj;          ! Size of local patch
end type MLIFEPatchDesc

type MLIFETiming
    double precision packtime, unpacktime ! Time to pack/unpack data if separate
    double precision exchtime             ! Time for exchange
    double precision itertime             ! Time for the iteration loop
end type MLIFETiming

type MLIFEOptions
    integer gNI, gNJ           ! Size of global mesh
    integer pNI, pNJ           ! Size of processor mesh
    integer nIter              ! Number of iterations
    integer restartIter        ! Used to determine when to create restart files
    integer testing            ! < 0 if not testing.  >= 0 to select which test
    logical verbose            ! Whether verbose output used
    logical doIO               ! Whether I/O used
    character(64) prefix       ! Name for output file prefix
end type MLIFEOptions

integer, parameter :: BORN = 1, DIES=0

end module mlife2dtypes

module mlife2dfuncs
interface
   subroutine MLIFE_ParseArgs(options)
     use mlife2dtypes
     type(MLIFEOptions), intent(out) :: options
   end subroutine MLIFE_ParseArgs
end interface

interface
   subroutine MLIFE_PatchCreateProcessMesh(options, patch)
     use mlife2dtypes
     type(MLIFEOptions), intent(in)    ::options
     type(MLIFEPatchDesc), intent(out) ::patch
   end subroutine MLIFE_PatchCreateProcessMesh
end interface

interface
   subroutine MLIFE_PatchCreateProcessMeshWithCart(options, patch)
     use mlife2dtypes
     type(MLIFEOptions), intent(in)    ::options
     type(MLIFEPatchDesc), intent(out) ::patch
   end subroutine MLIFE_PatchCreateProcessMeshWithCart
end interface

interface
   subroutine MLIFE_PatchCreateDataMeshDesc(options, patch)
     use mlife2dtypes
     type(MLIFEOptions), intent(in)  ::options
     type(MLIFEPatchDesc)            ::patch
   end subroutine MLIFE_PatchCreateDataMeshDesc
end interface

interface
   subroutine MLIFE_AllocateLocalMesh(patch, m1, m2)
     use mlife2dtypes
     Type(MLIFEPatchDesc) :: patch
     integer, allocatable, intent(out) :: m1(:,:), m2(:,:)
   end subroutine MLIFE_AllocateLocalMesh
end interface

interface
   subroutine MLIFE_FreeLocalMesh(patch, m1, m2)
     use mlife2dtypes
     Type(MLIFEPatchDesc) :: patch
     integer, allocatable, intent(inout) :: m1(:,:), m2(:,:)
   end subroutine MLIFE_FreeLocalMesh
end interface

interface
   subroutine MLIFE_InitLocalMesh(patch, m1, m2)
     use mlife2dtypes
     Type(MLIFEPatchDesc) :: patch
     integer, allocatable, intent(inout) :: m1(:,:), m2(:,:)
   end subroutine MLIFE_InitLocalMesh
end interface

interface
   subroutine MLIFE_Abort(str)
     character(*), intent(in) :: str
   end subroutine MLIFE_Abort
end interface

! pt-2-pt with isend/irecv
interface
   function MLIFE_ExchangeInitPt2pt(patch, m1, m2)
     use mlife2dtypes
     Type(MLIFEPatchDesc), intent(in)    :: patch
     integer, allocatable, intent(inout) :: m1(:,:), m2(:,:)
     integer MLIFE_ExchangeInitPt2pt
   end function MLIFE_ExchangeInitPt2pt
end interface

interface
   function MLIFE_ExchangeEndPt2pt()
     integer MLIFE_ExchangeEndPt2pt
   end function MLIFE_ExchangeEndPt2pt
end interface

interface
   function MLIFE_ExchangePt2pt(patch, matrix, timedata, phase)
     use mlife2dtypes
     Type(MLIFEPatchDesc) :: patch
     integer, allocatable, intent(inout) :: matrix(:,:)
     Type(MLIFETiming)    :: timedata
     integer, intent(in)  :: phase
     integer MLIFE_ExchangePt2pt
   end function MLIFE_ExchangePt2pt
end interface

! pt-2-pt with isend/irecv, including diagonals
interface
   function MLIFE_ExchangeInitPt2pt9(patch, m1, m2)
     use mlife2dtypes
     Type(MLIFEPatchDesc), intent(in)    :: patch
     integer, allocatable, intent(inout) :: m1(:,:), m2(:,:)
     integer MLIFE_ExchangeInitPt2pt9
   end function MLIFE_ExchangeInitPt2pt9
end interface

interface
   function MLIFE_ExchangeEndPt2pt9()
     integer MLIFE_ExchangeEndPt2pt9
   end function MLIFE_ExchangeEndPt2pt9
end interface

interface
   function MLIFE_ExchangePt2pt9(patch, matrix, timedata, phase)
     use mlife2dtypes
     Type(MLIFEPatchDesc) :: patch
     integer, allocatable, intent(inout) :: matrix(:,:)
     Type(MLIFETiming)    :: timedata
     integer, intent(in)  :: phase
     integer MLIFE_ExchangePt2pt9
   end function MLIFE_ExchangePt2pt9
end interface

! pt-2-pt with irecv/Send
interface
   function MLIFE_ExchangeInitPt2ptSnd(patch, m1, m2)
     use mlife2dtypes
     Type(MLIFEPatchDesc), intent(in)    :: patch
     integer, allocatable, intent(inout) :: m1(:,:), m2(:,:)
     integer MLIFE_ExchangeInitPt2ptSnd
   end function MLIFE_ExchangeInitPt2ptSnd
end interface

interface
   function MLIFE_ExchangeEndPt2ptSnd()
     integer MLIFE_ExchangeEndPt2ptSnd
   end function MLIFE_ExchangeEndPt2ptSnd
end interface

interface
   function MLIFE_ExchangePt2ptSnd(patch, matrix, timedata, phase)
     use mlife2dtypes
     Type(MLIFEPatchDesc) :: patch
     integer, allocatable, intent(inout) :: matrix(:,:)
     Type(MLIFETiming)    :: timedata
     integer, intent(in)  :: phase
     integer MLIFE_ExchangePt2ptSnd
   end function MLIFE_ExchangePt2ptSnd
end interface

! pt-2-pt with irecv/isend but user pack/unpack of vector
interface
   function MLIFE_ExchangeInitPt2ptUV(patch, m1, m2)
     use mlife2dtypes
     Type(MLIFEPatchDesc), intent(in)    :: patch
     integer, allocatable, intent(inout) :: m1(:,:), m2(:,:)
     integer MLIFE_ExchangeInitPt2ptUV
   end function MLIFE_ExchangeInitPt2ptUV
end interface

interface
   function MLIFE_ExchangeEndPt2ptUV()
     integer MLIFE_ExchangeEndPt2ptUV
   end function MLIFE_ExchangeEndPt2ptUV
end interface

interface
   function MLIFE_ExchangePt2ptUV(patch, matrix, timedata, phase)
     use mlife2dtypes
     Type(MLIFEPatchDesc) :: patch
     integer, allocatable, intent(inout) :: matrix(:,:)
     Type(MLIFETiming)    :: timedata
     integer, intent(in)  :: phase
     integer MLIFE_ExchangePt2ptUV
   end function MLIFE_ExchangePt2ptUV
end interface


! RMA with Fence
interface
   function MLIFE_ExchangeInitFence(patch, m1, m2)
     use mlife2dtypes
     Type(MLIFEPatchDesc), intent(in)    :: patch
     integer, allocatable, intent(inout) :: m1(:,:), m2(:,:)
     integer MLIFE_ExchangeInitFence
   end function MLIFE_ExchangeInitFence
end interface

interface
   function MLIFE_ExchangeEndFence()
     integer MLIFE_ExchangeEndFence
   end function MLIFE_ExchangeEndFence
end interface

interface
   function MLIFE_ExchangeFence(patch, matrix, timedata, phase)
     use mlife2dtypes
     Type(MLIFEPatchDesc) :: patch
     integer, allocatable, intent(inout) :: matrix(:,:)
     Type(MLIFETiming)    :: timedata
     integer, intent(in)  :: phase
     integer MLIFE_ExchangeFence
   end function MLIFE_ExchangeFence
end interface

! Student example: uncomment for "mylife2d.f90" routines
! interface
!    function MLIFE_ExchangeInitExample(patch, m1, m2)
!      use mlife2dtypes
!      Type(MLIFEPatchDesc), intent(in)    :: patch
!      integer, allocatable, intent(inout) :: m1(:,:), m2(:,:)
!      integer MLIFE_ExchangeInitExample
!    end function MLIFE_ExchangeInitExample
! end interface

! interface
!    function MLIFE_ExchangeEndExample()
!      integer MLIFE_ExchangeEndExample
!    end function MLIFE_ExchangeEndExample
! end interface

! interface
!    function MLIFE_ExchangeExample(patch, matrix, timedata, phase)
!      use mlife2dtypes
!      Type(MLIFEPatchDesc) :: patch
!      integer, allocatable, intent(inout) :: matrix(:,:)
!      Type(MLIFETiming)    :: timedata
!      integer, intent(in)  :: phase
!      integer MLIFE_ExchangeExample
!    end function MLIFE_ExchangeExample
! end interface


end module mlife2dfuncs

