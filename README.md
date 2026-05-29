# symplectic-physics

**Symplectic integrators for Hamiltonian mechanics — Fortran 2008, BLAS-compatible, cache-friendly, close to the metal.**

Fortran implementation of the same integrators found in [symplectic-spin](https://github.com/SuperInstance/symplectic-spin) (Rust), optimized for high-performance physics simulations requiring LAPACK/BLAS acceleration.

## What This Gives You

- **Three integrators** — Symplectic Euler (1st order), Störmer-Verlet (2nd order), Yoshida 4th order
- **Symplectic matrix ops** — verify, compose, and generate random symplectic matrices
- **Conservation tracking** — energy drift, phase space volume, angular momentum
- **N-body simulation** — gravitational N-body with symplectic integration, Kepler orbit stability
- **BLAS/LAPACK compatible** — uses standard linear algebra routines for performance
- **Cache-friendly** — contiguous array layout, no unnecessary allocations

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

## Modules

| Module | Description |
|--------|-------------|
| `symplectic` | Symplectic matrix operations — verify, compose, random generation |
| `hamiltonian` | Hamiltonian system type for separable H = T + V |
| `integrators` | Symplectic Euler, Störmer-Verlet, 4th-order Yoshida |
| `conservation` | Energy drift, phase space volume, angular momentum tracking |
| `nbody` | N-body gravitational simulation using symplectic integrators |

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

## Building

Requires `gfortran` with Fortran 2008 support, LAPACK, and BLAS.

```bash
make        # Build library and run tests
make lib    # Build static library only → build/libsymplectic.a
make test   # Build and run tests
make clean  # Clean build artifacts
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

## How It Fits

Part of the SuperInstance symplectic ecosystem:

- **[symplectic-spin](https://github.com/SuperInstance/symplectic-spin)** — Same integrators in pure Rust (zero deps)
- **symplectic-physics** — Fortran 2008 with BLAS/LAPACK (this repo)
- **[spectral-mechanics](https://github.com/SuperInstance/spectral-mechanics)** — Graphs as spring-mass systems using Verlet

## License

MIT

Part of the [SuperInstance](https://github.com/SuperInstance) ecosystem.
