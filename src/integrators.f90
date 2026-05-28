module integrators
  use hamiltonian
  implicit none
  private
  public :: symplectic_euler, stormer_verlet, yoshida4

contains

  subroutine symplectic_euler(sys, qq, pp, dt, nsteps, q_traj, p_traj)
    type(HamiltonianSystem), intent(in) :: sys
    real(8), intent(inout) :: qq(:), pp(:)
    real(8), intent(in) :: dt
    integer, intent(in) :: nsteps
    real(8), intent(out) :: q_traj(:,:), p_traj(:,:)
    real(8) :: grad(sys%ndof)
    integer :: kk

    q_traj(:, 1) = qq
    p_traj(:, 1) = pp

    do kk = 1, nsteps
      call sys%dV_dq_ptr(qq, sys%ndof, grad)
      pp = pp - dt * grad

      call sys%dT_dp_ptr(pp, sys%ndof, grad)
      qq = qq + dt * grad

      q_traj(:, kk+1) = qq
      p_traj(:, kk+1) = pp
    end do
  end subroutine

  subroutine stormer_verlet(sys, qq, pp, dt, nsteps, q_traj, p_traj)
    type(HamiltonianSystem), intent(in) :: sys
    real(8), intent(inout) :: qq(:), pp(:)
    real(8), intent(in) :: dt
    integer, intent(in) :: nsteps
    real(8), intent(out) :: q_traj(:,:), p_traj(:,:)
    real(8) :: grad(sys%ndof)
    integer :: kk

    q_traj(:, 1) = qq
    p_traj(:, 1) = pp

    do kk = 1, nsteps
      call sys%dV_dq_ptr(qq, sys%ndof, grad)
      pp = pp - 0.5d0 * dt * grad

      call sys%dT_dp_ptr(pp, sys%ndof, grad)
      qq = qq + dt * grad

      call sys%dV_dq_ptr(qq, sys%ndof, grad)
      pp = pp - 0.5d0 * dt * grad

      q_traj(:, kk+1) = qq
      p_traj(:, kk+1) = pp
    end do
  end subroutine

  subroutine yoshida4(sys, qq, pp, dt, nsteps, q_traj, p_traj)
    type(HamiltonianSystem), intent(in) :: sys
    real(8), intent(inout) :: qq(:), pp(:)
    real(8), intent(in) :: dt
    integer, intent(in) :: nsteps
    real(8), intent(out) :: q_traj(:,:), p_traj(:,:)
    real(8) :: w1, wm, grad(sys%ndof)
    integer :: kk

    ! Yoshida coefficients: w1 = 1/(2 - 2^(1/3)), wm = -2^(1/3)/(2-2^(1/3))
    w1 = 1.0d0 / (2.0d0 - 2.0d0**(1.0d0/3.0d0))
    wm = -2.0d0**(1.0d0/3.0d0) * w1

    q_traj(:, 1) = qq
    p_traj(:, 1) = pp

    do kk = 1, nsteps
      ! ========= First Verlet with step w1*dt =========
      call sys%dV_dq_ptr(qq, sys%ndof, grad)
      pp = pp - 0.5d0 * w1 * dt * grad
      call sys%dT_dp_ptr(pp, sys%ndof, grad)
      qq = qq + w1 * dt * grad
      ! ========= Second Verlet with step wm*dt =========
      ! (share the boundary half-kick with first and third)
      call sys%dV_dq_ptr(qq, sys%ndof, grad)
      pp = pp - 0.5d0 * (w1 + wm) * dt * grad
      call sys%dT_dp_ptr(pp, sys%ndof, grad)
      qq = qq + wm * dt * grad
      ! ========= Third Verlet with step w1*dt =========
      call sys%dV_dq_ptr(qq, sys%ndof, grad)
      pp = pp - 0.5d0 * (wm + w1) * dt * grad
      call sys%dT_dp_ptr(pp, sys%ndof, grad)
      qq = qq + w1 * dt * grad
      call sys%dV_dq_ptr(qq, sys%ndof, grad)
      pp = pp - 0.5d0 * w1 * dt * grad

      q_traj(:, kk+1) = qq
      p_traj(:, kk+1) = pp
    end do
  end subroutine

end module
