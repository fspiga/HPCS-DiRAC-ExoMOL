module accuracy
  implicit none
  !private

  public ik, hik, rk, trk, cl

  integer, parameter :: ik          = selected_int_kind(8)       ! "Normal" integers. This must map on
  integer, parameter :: hik         = selected_int_kind(8)       ! "Pointer" integers - sufficient to store
  integer, parameter :: rk          = selected_real_kind(12,25)  ! "Normal" reals and complex (complexi? :-)
  integer, parameter :: cl          = 80                         ! Max character string length
  integer, parameter :: trk        = selected_real_kind(12)

end module accuracy

