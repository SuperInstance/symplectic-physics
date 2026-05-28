module conservation
  use hamiltonian
  implicit none
  private
  public :: energy_drift, volume_conserved, angular_momentum

contains

  subroutine energy_drift(sys, q_traj, p_traj, nsteps, ndof, drift)
    type(HamiltonianSystem), intent(in) :: sys
    real(8), intent(in) :: q_traj(:,:), p_traj(:,:)
    integer, intent(in) :: nsteps, ndof
    real(8), intent(out) :: drift(:)
    real(8) :: h0, hi
    real(8) :: qq(ndof), pp(ndof)
    integer :: ii

    qq = q_traj(:, 1)
    pp = p_traj(:, 1)
    h0 = sys%T(pp, ndof) + sys%V(qq, ndof)

    do ii = 1, nsteps + 1
      qq = q_traj(:, ii)
      pp = p_traj(:, ii)
      hi = sys%T(pp, ndof) + sys%V(qq, ndof)
      if (abs(h0) > 1.0d-15) then
        drift(ii) = abs(hi - h0) / abs(h0)
      else
        drift(ii) = abs(hi - h0)
      end if
    end do
  end subroutine

  logical function volume_conserved(mat_i, mat_f, nn, tol)
    real(8), intent(in) :: mat_i(:,:), mat_f(:,:)
    integer, intent(in) :: nn
    real(8), intent(in) :: tol
    real(8) :: det_i, det_f
    integer :: sz

    sz = 2 * nn
    call matrix_determinant(mat_i, sz, det_i)
    call matrix_determinant(mat_f, sz, det_f)

    volume_conserved = (abs(det_i - det_f) < tol)
  end function

  real(8) function angular_momentum(qq, pp, ndof)
    real(8), intent(in) :: qq(:), pp(:)
    integer, intent(in) :: ndof
    integer :: ii

    angular_momentum = 0.0d0
    do ii = 1, ndof/2
      angular_momentum = angular_momentum + qq(2*ii-1) * pp(2*ii) - qq(2*ii) * pp(2*ii-1)
    end do
  end function

  subroutine matrix_determinant(aa, nn, det)
    real(8), intent(in) :: aa(:,:)
    integer, intent(in) :: nn
    real(8), intent(out) :: det
    real(8) :: lu(nn, nn)
    integer :: ipiv(nn), info, ii

    lu = aa
    call dgetrf(nn, nn, lu, nn, ipiv, info)

    det = 1.0d0
    do ii = 1, nn
      det = det * lu(ii, ii)
      if (ipiv(ii) /= ii) det = -det
    end do
  end subroutine

end module
