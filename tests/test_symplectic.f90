program test_symplectic
  use symplectic
  use hamiltonian
  use integrators
  use conservation
  use nbody
  implicit none

  integer :: npassed, nfailed, ntotal
  npassed = 0; nfailed = 0; ntotal = 0

  call test_symplectic_matrix_verification()
  call test_identity_is_symplectic()
  call test_compose_symplectic()
  call test_random_symplectic()
  call test_euler_energy_drift()
  call test_verlet_energy_oscillation()
  call test_yoshida_energy_drift()
  call test_pendulum_conservation()
  call test_phase_space_volume()
  call test_kepler_orbit()
  call test_three_body_stability()
  call test_angular_momentum_conservation()

  print '(a)', ''
  print '(a,i0,a,i0,a,i0,a)', 'Results: ', npassed, ' passed, ', nfailed, ' failed (', ntotal, ' total)'
  if (nfailed > 0) then
    print '(a)', 'SOME TESTS FAILED'
    stop 1
  else
    print '(a)', 'ALL TESTS PASSED'
  end if

contains

  subroutine check(name, condition)
    character(len=*), intent(in) :: name
    logical, intent(in) :: condition
    ntotal = ntotal + 1
    if (condition) then
      npassed = npassed + 1
      print '(a,a)', '  PASS: ', trim(name)
    else
      nfailed = nfailed + 1
      print '(a,a)', '  FAIL: ', trim(name)
    end if
  end subroutine

  ! Test 1: Symplectic matrix verification
  subroutine test_symplectic_matrix_verification()
    real(8) :: mat(4, 4)
    print '(a)', 'Test 1: Symplectic matrix verification'

    ! Shear matrix: [[I, A], [0, I]] in 2x2 block form
    mat = 0.0d0
    mat(1,1) = 1.0d0; mat(1,3) = 0.5d0
    mat(2,2) = 1.0d0; mat(2,4) = 0.5d0
    mat(3,3) = 1.0d0; mat(4,4) = 1.0d0

    call check('shear matrix is symplectic', is_symplectic(mat, 2, 1.0d-10))
  end subroutine

  ! Test 2: Identity is symplectic
  subroutine test_identity_is_symplectic()
    real(8) :: id4(4, 4)
    integer :: ii
    print '(a)', 'Test 2: Identity is symplectic'

    id4 = 0.0d0
    do ii = 1, 4
      id4(ii, ii) = 1.0d0
    end do

    call check('4x4 identity is symplectic', is_symplectic(id4, 2, 1.0d-14))
  end subroutine

  ! Test 3: Compose two symplectic matrices
  subroutine test_compose_symplectic()
    real(8) :: m1(4,4), m2(4,4), m3(4,4)
    print '(a)', 'Test 3: Compose two symplectic matrices'

    m1 = 0.0d0
    m1(1,1) = 1.0d0; m1(1,3) = 0.3d0
    m1(2,2) = 1.0d0; m1(2,4) = 0.3d0
    m1(3,3) = 1.0d0; m1(4,4) = 1.0d0

    m2 = 0.0d0
    m2(1,1) = 1.0d0; m2(2,2) = 1.0d0
    m2(3,1) = 0.7d0; m2(3,3) = 1.0d0
    m2(4,2) = 0.7d0; m2(4,4) = 1.0d0

    call symplectic_compose(m1, m2, m3, 2)
    call check('composition is symplectic', is_symplectic(m3, 2, 1.0d-10))
  end subroutine

  ! Test 4: Random symplectic matrix
  subroutine test_random_symplectic()
    real(8), allocatable :: mat(:,:)
    print '(a)', 'Test 4: Random symplectic matrix'

    call random_symplectic(mat, 2, 42)
    call check('random symplectic passes verification', is_symplectic(mat, 2, 1.0d-8))
  end subroutine

  ! Harmonic oscillator helpers
  function ho_T(pp, ndof) result(val)
    real(8), intent(in) :: pp(:)
    integer, intent(in) :: ndof
    real(8) :: val
    val = 0.5d0 * pp(1)**2
  end function

  function ho_V(qq, ndof) result(val)
    real(8), intent(in) :: qq(:)
    integer, intent(in) :: ndof
    real(8) :: val
    val = 0.5d0 * qq(1)**2
  end function

  subroutine ho_dTdp(pp, ndof, grad)
    real(8), intent(in) :: pp(:)
    integer, intent(in) :: ndof
    real(8), intent(out) :: grad(:)
    grad(1) = pp(1)
  end subroutine

  subroutine ho_dVdq(qq, ndof, grad)
    real(8), intent(in) :: qq(:)
    integer, intent(in) :: ndof
    real(8), intent(out) :: grad(:)
    grad(1) = qq(1)
  end subroutine

  ! Test 5: Euler energy drift
  subroutine test_euler_energy_drift()
    type(HamiltonianSystem) :: sys
    real(8) :: qq(1), pp(1), dt
    real(8) :: q_traj(1, 1001), p_traj(1, 1001), drift(1001)
    real(8) :: max_drift

    print '(a)', 'Test 5: Euler energy drift < 1e-2 over 1000 steps'

    sys%ndof = 1
    sys%T => ho_T; sys%V => ho_V
    sys%dT_dp_ptr => ho_dTdp; sys%dV_dq_ptr => ho_dVdq

    qq = 1.0d0; pp = 0.0d0; dt = 0.01d0

    call symplectic_euler(sys, qq, pp, dt, 1000, q_traj, p_traj)
    call energy_drift(sys, q_traj, p_traj, 1000, 1, drift)

    max_drift = maxval(drift)
    call check('Euler energy drift < 1e-2', max_drift < 1.0d-2)
  end subroutine

  ! Test 6: Verlet energy oscillation
  subroutine test_verlet_energy_oscillation()
    type(HamiltonianSystem) :: sys
    real(8) :: qq(1), pp(1), dt
    real(8) :: q_traj(1, 1001), p_traj(1, 1001), drift(1001)
    real(8) :: max_drift

    print '(a)', 'Test 6: Verlet energy oscillates but doesn''t drift'

    sys%ndof = 1
    sys%T => ho_T; sys%V => ho_V
    sys%dT_dp_ptr => ho_dTdp; sys%dV_dq_ptr => ho_dVdq

    qq = 1.0d0; pp = 0.0d0; dt = 0.01d0

    call stormer_verlet(sys, qq, pp, dt, 1000, q_traj, p_traj)
    call energy_drift(sys, q_traj, p_traj, 1000, 1, drift)

    max_drift = maxval(drift)
    call check('Verlet energy bounded < 1e-4', max_drift < 1.0d-4)
  end subroutine

  ! Test 7: Yoshida energy drift
  subroutine test_yoshida_energy_drift()
    type(HamiltonianSystem) :: sys
    real(8) :: qq(1), pp(1), dt
    real(8) :: q_traj(1, 1001), p_traj(1, 1001), drift(1001)
    real(8) :: max_drift

    print '(a)', 'Test 7: Yoshida energy drift < 1e-9 over 1000 steps'

    sys%ndof = 1
    sys%T => ho_T; sys%V => ho_V
    sys%dT_dp_ptr => ho_dTdp; sys%dV_dq_ptr => ho_dVdq

    qq = 1.0d0; pp = 0.0d0; dt = 0.01d0

    call yoshida4(sys, qq, pp, dt, 1000, q_traj, p_traj)
    call energy_drift(sys, q_traj, p_traj, 1000, 1, drift)

    max_drift = maxval(drift)
    call check('Yoshida energy drift < 1e-9', max_drift < 1.0d-9)
  end subroutine

  ! Pendulum helpers
  function pend_V(qq, ndof) result(val)
    real(8), intent(in) :: qq(:)
    integer, intent(in) :: ndof
    real(8) :: val
    val = 1.0d0 - cos(qq(1))
  end function

  subroutine pend_dVdq(qq, ndof, grad)
    real(8), intent(in) :: qq(:)
    integer, intent(in) :: ndof
    real(8), intent(out) :: grad(:)
    grad(1) = sin(qq(1))
  end subroutine

  ! Test 8: Pendulum conservation
  subroutine test_pendulum_conservation()
    type(HamiltonianSystem) :: sys
    real(8) :: qq(1), pp(1), dt
    real(8) :: q_traj(1, 10001), p_traj(1, 10001), drift(10001)
    real(8) :: max_drift

    print '(a)', 'Test 8: Pendulum energy conservation over 10000 steps'

    sys%ndof = 1
    sys%T => ho_T; sys%V => pend_V
    sys%dT_dp_ptr => ho_dTdp; sys%dV_dq_ptr => pend_dVdq

    qq = 1.0d0; pp = 0.0d0; dt = 0.001d0

    call stormer_verlet(sys, qq, pp, dt, 10000, q_traj, p_traj)
    call energy_drift(sys, q_traj, p_traj, 10000, 1, drift)

    max_drift = maxval(drift)
    call check('Pendulum energy drift < 1e-4 over 10000 steps', max_drift < 1.0d-4)
  end subroutine

  ! Test 9: Phase space volume
  subroutine test_phase_space_volume()
    real(8), allocatable :: m1(:,:), m2(:,:)
    real(8) :: id4(4,4)
    integer :: ii

    print '(a)', 'Test 9: Phase space volume preserved'

    id4 = 0.0d0
    do ii = 1, 4
      id4(ii,ii) = 1.0d0
    end do
    call check('identity preserves volume', volume_conserved(id4, id4, 2, 1.0d-10))

    call random_symplectic(m1, 2, 10)
    call random_symplectic(m2, 2, 20)
    call check('random symplectic preserves volume', volume_conserved(m1, m2, 2, 1.0d-6))
  end subroutine

  ! Test 10: Kepler orbit
  subroutine test_kepler_orbit()
    type(HamiltonianSystem) :: sys
    real(8) :: qq(4), pp(4), dt, r_init, r_final, r_min, r_max, r
    integer :: kk

    print '(a)', 'Test 10: 2-body Kepler orbit stability over 10000 steps'

    ! Use center-of-mass frame. Equal masses.
    ! Body 1 at (-0.5, 0), body 2 at (0.5, 0)
    ! For circular orbit: v = 0.5*sqrt(1/separation) per body in CoM frame
    ! F = m1*m2/r^2 = 1/1 = 1, a = v^2/(r/2) for each body, v = 0.5
    call nbody_gravity_init(2, [1.0d0, 1.0d0], &
      reshape([-0.5d0, 0.0d0, 0.5d0, 0.0d0], [2,2]), &
      reshape([0.0d0, -0.5d0, 0.0d0, 0.5d0], [2,2]), sys)

    qq = [-0.5d0, 0.0d0, 0.5d0, 0.0d0]
    pp = [0.0d0, -0.5d0, 0.0d0, 0.5d0]

    r_init = sqrt((qq(3)-qq(1))**2 + (qq(4)-qq(2))**2)
    r_min = r_init; r_max = r_init
    dt = 0.01d0

    do kk = 1, 10000
      call nbody_step(sys, qq, pp, dt, 1)
      r = sqrt((qq(3)-qq(1))**2 + (qq(4)-qq(2))**2)
      r_min = min(r_min, r)
      r_max = max(r_max, r)
    end do

    r_final = sqrt((qq(3)-qq(1))**2 + (qq(4)-qq(2))**2)
    call check('Kepler orbit radius bounded', abs(r_final - r_init) < 0.5d0)
    call check('Kepler orbit didn''t blow up', maxval(abs(qq)) < 10.0d0)
    call check('Kepler orbit stays elliptical', r_max / r_min < 5.0d0)
  end subroutine

  ! Test 11: 3-body stability
  subroutine test_three_body_stability()
    type(HamiltonianSystem) :: sys
    real(8) :: qq(6), pp(6), dt
    real(8) :: max_disp
    integer :: kk

    print '(a)', 'Test 11: 3-body stability test'

    call nbody_gravity_init(3, [1.0d0, 1.0d0, 1.0d0], &
      reshape([0.0d0, 0.0d0, 2.0d0, 0.0d0, 1.0d0, 1.73d0], [2,3]), &
      reshape([0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0], [2,3]), sys)

    qq = [0.0d0, 0.0d0, 2.0d0, 0.0d0, 1.0d0, 1.73d0]
    pp = [0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0]
    dt = 0.0005d0

    do kk = 1, 2000
      call nbody_step(sys, qq, pp, dt, 1)
    end do

    max_disp = maxval(abs(qq))
    call check('3-body didn''t diverge', max_disp < 20.0d0)
  end subroutine

  ! Test 12: Angular momentum conservation
  subroutine test_angular_momentum_conservation()
    type(HamiltonianSystem) :: sys
    real(8) :: qq(4), pp(4), dt, l0, lf
    integer :: kk

    print '(a)', 'Test 12: Angular momentum conservation in N-body'

    call nbody_gravity_init(2, [1.0d0, 1.0d0], &
      reshape([0.0d0, 0.0d0, 1.0d0, 0.0d0], [2,2]), &
      reshape([0.0d0, -1.0d0, 0.0d0, 1.0d0], [2,2]), sys)

    qq = [0.0d0, 0.0d0, 1.0d0, 0.0d0]
    pp = [0.0d0, -1.0d0, 0.0d0, 1.0d0]

    l0 = angular_momentum(qq, pp, 4)
    dt = 0.001d0

    do kk = 1, 5000
      call nbody_step(sys, qq, pp, dt, 1)
    end do

    lf = angular_momentum(qq, pp, 4)
    call check('angular momentum conserved < 1e-4', abs(lf - l0) < 1.0d-4)
  end subroutine

end program
