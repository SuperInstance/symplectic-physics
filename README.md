# symplectic-physics

Fortran 2008 library implementing symplectic integrators for Hamiltonian mechanics. Cache-friendly, BLAS-compatible, and close to the metal.

## Modules

| Module | Description |
|--------|-------------|
| `symplectic` | Symplectic matrix operations — verify, compose, random generation |
| `hamiltonian` | Hamiltonian system type for separable H = T + V |
| `integrators` | Symplectic Euler, Störmer-Verlet, 4th-order Yoshida |
| `conservation` | Energy drift, phase space volume, angular momentum tracking |
| `nbody` | N-body gravitational simulation using symplectic integrators |

## Building

Requires `gfortran` with Fortran 2008 support, LAPACK, and BLAS.

```bash
make        # Build library and run tests
make lib    # Build static library only
make test   # Build and run tests
make clean  # Clean build artifacts
```

Output: `build/libsymplectic.a` static library.

## Quick Start

```fortran
use hamiltonian
use integrators

type(HamiltonianSystem) :: sys
real(8) :: q(1), p(1)
real(8), allocatable :: q_traj(:,:), p_traj(:,:)

! Define your system (e.g., harmonic oscillator)
sys%ndof = 1
sys%T => my_kinetic
sys%V => my_potential
sys%dT_dp_ptr => my_dT_dp
sys%dV_dq_ptr => my_dV_dq

q = 1.0d0; p = 0.0d0
allocate(q_traj(1, 1001), p_traj(1, 1001))

! Integrate with Störmer-Verlet (2nd order, time-reversible)
call stormer_verlet(sys, q, p, 0.01d0, 1000, q_traj, p_traj)
```

## Integrators

### Symplectic Euler (1st order)
Fast, simple, symplectic. Good for quick prototyping.

### Störmer-Verlet (2nd order)
Time-reversible, excellent energy conservation. The workhorse integrator.

### Yoshida 4th order
Six-stage composition method. Near-machine-precision energy conservation for smooth systems.

## N-Body Simulation

```fortran
use nbody

type(HamiltonianSystem) :: sys
real(8) :: q(4), p(4), q_out(4, 1001), p_out(4, 1001)

call nbody_gravity_init(2, [1.0d0, 1.0d0], &
    reshape([0.0d0, 0.0d0, 1.0d0, 0.0d0], [2,2]), &
    reshape([0.0d0, -1.0d0, 0.0d0, 1.0d0], [2,2]), sys)

q = [0.0d0, 0.0d0, 1.0d0, 0.0d0]
p = [0.0d0, -1.0d0, 0.0d0, 1.0d0]

call nbody_run(sys, q, p, 0.001d0, 1000, q_out, p_out)
```

## Design Principles

- **Symplectic by construction** — integrators preserve the symplectic 2-form
- **Separable Hamiltonians** — exploits H = T(p) + V(q) for efficient splitting
- **BLAS/LAPACK compatible** — uses standard linear algebra routines
- **Cache-friendly** — contiguous array layout, no unnecessary allocations

## Testing

12 tests covering:
- Symplectic matrix verification and composition
- Energy conservation for harmonic oscillator, pendulum
- Phase space volume preservation
- N-body Kepler orbit stability
- 3-body stability
- Angular momentum conservation

## License

MIT

Part of the [SuperInstance OpenConstruct](https://github.com/SuperInstance/OpenConstruct) ecosystem.
