module hamiltonian
  implicit none
  private
  public :: HamiltonianSystem, kinetic_iface, potential_iface, dTdp_iface, dVdq_iface

  abstract interface
    real(8) function kinetic_iface(pp, ndof)
      real(8), intent(in) :: pp(:)
      integer, intent(in) :: ndof
    end function
    real(8) function potential_iface(qq, ndof)
      real(8), intent(in) :: qq(:)
      integer, intent(in) :: ndof
    end function
    subroutine dTdp_iface(pp, ndof, grad)
      real(8), intent(in) :: pp(:)
      integer, intent(in) :: ndof
      real(8), intent(out) :: grad(:)
    end subroutine
    subroutine dVdq_iface(qq, ndof, grad)
      real(8), intent(in) :: qq(:)
      integer, intent(in) :: ndof
      real(8), intent(out) :: grad(:)
    end subroutine
  end interface

  type :: HamiltonianSystem
    integer :: ndof = 0
    procedure(kinetic_iface), pointer, nopass :: T => null()
    procedure(potential_iface), pointer, nopass :: V => null()
    procedure(dTdp_iface), pointer, nopass :: dT_dp_ptr => null()
    procedure(dVdq_iface), pointer, nopass :: dV_dq_ptr => null()
  end type

end module
