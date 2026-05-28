module nbody
  use hamiltonian
  implicit none
  private
  public :: nbody_gravity_init, nbody_step, nbody_run

  real(8), allocatable, save :: nb_masses(:)
  integer, save :: nb_n_bodies = 0

contains

  function nbody_T(pp, ndof) result(val)
    real(8), intent(in) :: pp(:)
    integer, intent(in) :: ndof
    real(8) :: val
    integer :: ii
    val = 0.0d0
    do ii = 1, ndof
      val = val + 0.5d0 * pp(ii) * pp(ii)
    end do
  end function

  function nbody_V(qq, ndof) result(val)
    real(8), intent(in) :: qq(:)
    integer, intent(in) :: ndof
    real(8) :: val
    val = 0.0d0
  end function

  subroutine nbody_dT_dp(pp, ndof, grad)
    real(8), intent(in) :: pp(:)
    integer, intent(in) :: ndof
    real(8), intent(out) :: grad(:)
    grad = pp
  end subroutine

  subroutine nbody_dV_dq(qq, ndof, grad)
    real(8), intent(in) :: qq(:)
    integer, intent(in) :: ndof
    real(8), intent(out) :: grad(:)
    integer :: ii, jj, xi, yi, xj, yj
    real(8) :: dx, dy, r, r3, force

    grad = 0.0d0
    do ii = 1, nb_n_bodies
      xi = 2*(ii-1) + 1
      yi = 2*(ii-1) + 2
      do jj = ii+1, nb_n_bodies
        xj = 2*(jj-1) + 1
        yj = 2*(jj-1) + 2

        dx = qq(xj) - qq(xi)
        dy = qq(yj) - qq(yi)
        r = sqrt(dx*dx + dy*dy)
        r = max(r, 1.0d-4)
        r3 = r * r * r

        force = nb_masses(ii) * nb_masses(jj) / r3

        grad(xi) = grad(xi) - force * dx
        grad(yi) = grad(yi) - force * dy
        grad(xj) = grad(xj) + force * dx
        grad(yj) = grad(yj) + force * dy
      end do
    end do
  end subroutine

  subroutine nbody_gravity_init(n_bodies, masses, positions, velocities, sys)
    integer, intent(in) :: n_bodies
    real(8), intent(in) :: masses(:)
    real(8), intent(in) :: positions(:,:), velocities(:,:)
    type(HamiltonianSystem), intent(out) :: sys

    nb_n_bodies = n_bodies
    if (allocated(nb_masses)) deallocate(nb_masses)
    allocate(nb_masses(n_bodies))
    nb_masses = masses

    sys%ndof = 2 * n_bodies
    sys%T => nbody_T
    sys%V => nbody_V
    sys%dT_dp_ptr => nbody_dT_dp
    sys%dV_dq_ptr => nbody_dV_dq
  end subroutine

  subroutine nbody_step(sys, qq, pp, dt, method)
    type(HamiltonianSystem), intent(in) :: sys
    real(8), intent(inout) :: qq(:), pp(:)
    real(8), intent(in) :: dt
    integer, intent(in) :: method
    real(8) :: grad(sys%ndof)

    select case(method)
    case(0)
      call sys%dV_dq_ptr(qq, sys%ndof, grad)
      pp = pp - dt * grad
      call sys%dT_dp_ptr(pp, sys%ndof, grad)
      qq = qq + dt * grad
    case(1)
      call sys%dV_dq_ptr(qq, sys%ndof, grad)
      pp = pp - 0.5d0 * dt * grad
      call sys%dT_dp_ptr(pp, sys%ndof, grad)
      qq = qq + dt * grad
      call sys%dV_dq_ptr(qq, sys%ndof, grad)
      pp = pp - 0.5d0 * dt * grad
    case(2)
      call yoshida_single(sys, qq, pp, dt)
    end select
  end subroutine

  subroutine yoshida_single(sys, qq, pp, dt)
    type(HamiltonianSystem), intent(in) :: sys
    real(8), intent(inout) :: qq(:), pp(:)
    real(8), intent(in) :: dt
    real(8) :: w1, wm, grad(sys%ndof)

    w1 = 1.0d0 / (2.0d0 - 2.0d0**(1.0d0/3.0d0))
    wm = -2.0d0**(1.0d0/3.0d0) * w1

    call sys%dV_dq_ptr(qq, sys%ndof, grad);  pp = pp - 0.5d0*w1*dt*grad
    call sys%dT_dp_ptr(pp, sys%ndof, grad);  qq = qq + w1*dt*grad
    call sys%dV_dq_ptr(qq, sys%ndof, grad);  pp = pp - 0.5d0*(w1+wm)*dt*grad
    call sys%dT_dp_ptr(pp, sys%ndof, grad);  qq = qq + wm*dt*grad
    call sys%dV_dq_ptr(qq, sys%ndof, grad);  pp = pp - 0.5d0*(wm+w1)*dt*grad
    call sys%dT_dp_ptr(pp, sys%ndof, grad);  qq = qq + w1*dt*grad
    call sys%dV_dq_ptr(qq, sys%ndof, grad);  pp = pp - 0.5d0*w1*dt*grad
  end subroutine

  subroutine nbody_run(sys, qq, pp, dt, nsteps, q_out, p_out)
    type(HamiltonianSystem), intent(in) :: sys
    real(8), intent(inout) :: qq(:), pp(:)
    real(8), intent(in) :: dt
    integer, intent(in) :: nsteps
    real(8), intent(out) :: q_out(:,:), p_out(:,:)
    integer :: kk

    q_out(:, 1) = qq
    p_out(:, 1) = pp

    do kk = 1, nsteps
      call nbody_step(sys, qq, pp, dt, 1)
      q_out(:, kk+1) = qq
      p_out(:, kk+1) = pp
    end do
  end subroutine

end module
