module msleep
  use, intrinsic :: iso_c_binding, only : c_int
  implicit none
  interface
     function MLIFEIO_sleep(sec) bind(C,name="sleep")
       import
       integer (c_int) :: MLIFEIO_sleep
       integer (c_int), intent(in), VALUE :: sec
     end function MLIFEIO_sleep
  end interface
  contains
    subroutine MLIFEIO_msleep(msec)
      integer msec
      integer rv
      if (msec < 1000) then
         rv = MLIFEIO_sleep(1)
      else
         rv = MLIFEIO_sleep(msec/1000)
      endif
    end subroutine MLIFEIO_msleep
end module msleep
