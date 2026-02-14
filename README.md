# DSMC Suite

A parallel DSMC (Direct Simulation Monte Carlo) framework for planetary exosphere simulation, converted from MATLAB research code to C + OpenMP.

Originally developed for modeling Ceres' transient H2O atmosphere during PhD research (advisor: Prof. Wing-Huen Ip), now generalized for any airless body.

## Simulation Pipeline

```
[1] Ray Tracer (MC shadow detection)
     -> Solar flux per surface face
          |
[2] Thermal Model (1D Crank-Nicolson)
     -> Surface temperature field
          |
[3] Ballistic DSMC (collisionless)
     -> Atmospheric density + velocity records
          |
[4] Advanced DSMC (with collisions)
     -> Final density + collision statistics
```

## Building

Requires: C11 compiler, CMake 3.12+, OpenMP

```bash
# macOS (Apple Clang + Homebrew libomp)
brew install libomp cmake
mkdir build && cd build
cmake .. && make -j$(nproc)

# Linux (GCC)
mkdir build && cd build
cmake .. && make -j$(nproc)
```

Run tests:
```bash
./test_verlet
```

## Usage

```bash
# Run individual modules
./dsmc_suite raytrace -m vertices.dat -f faces.dat \
    --npoint 642 --nface 1280 --za 27.5 --tt 48 --nmonte 80 \
    -o output_prefix -t 8

./dsmc_suite thermal -o output_prefix -t 8

./dsmc_suite dsmc -n 10000 -o output_prefix -t 8

./dsmc_suite dsmc-adv -n 10000 -o output_prefix -t 8

# Run full pipeline (raytrace -> thermal -> ballistic -> advanced)
./dsmc_suite pipeline -n 10000 -o output_prefix -t 8
```

### Vesta Ray Tracing Example

Generate a rotation GIF using the NASA Dawn mission Vesta mesh:

```bash
# Downloads the GLB model from NASA, extracts mesh, runs ray tracer
python3 scripts/vesta_real_gif.py
```

## Project Structure

```
dsmc_suite/
├── include/
│   ├── types.h          # Vec3, Particle, SphericalGrid, Mesh structs
│   ├── constants.h      # Physical constants (Ceres defaults)
│   ├── grid.h           # Spherical grid with uniform/exponential height
│   ├── velocity.h       # Maxwell-Boltzmann velocity table sampling
│   ├── integrator.h     # Stormer-Verlet orbital integrator
│   ├── collision.h      # Hard-sphere DSMC collisions
│   ├── polyfit.h        # Least-squares polynomial fitting
│   ├── thermal.h        # 1D thermal diffusion + sublimation
│   ├── raytracer.h      # MC shadow/flux on triangle meshes
│   ├── fileio.h         # .dat file I/O
│   └── particle.h       # Particle initialization
├── src/
│   ├── main.c           # CLI entry point
│   ├── common/          # Shared physics modules
│   ├── dsmc/            # Ballistic + advanced DSMC
│   ├── raytracer/       # MC ray tracing
│   └── thermal/         # Crank-Nicolson thermal solver
├── scripts/
│   ├── vesta_raytrace_gif.py   # Synthetic Vesta mesh + GIF
│   └── vesta_real_gif.py       # NASA Dawn mesh + GIF
├── tests/
│   └── test_verlet.c    # Orbital + polynomial fit tests
└── docs/
    └── simulation_report.md
```

## Physics Modules

### Ray Tracer (`shadow_maker.c`)
Monte Carlo shadow detection on triangulated asteroid surfaces. For each face at each rotational time step, samples random points and traces rays toward the sun to determine occlusion by neighboring faces.

### Thermal Model (`thermal_model.c`)
1D radial heat diffusion with nonlinear surface boundary condition: solar absorption - thermal radiation (sigma T^4) - sublimation cooling. Solved with Crank-Nicolson implicit scheme + Thomas algorithm for the tridiagonal system. Surface temperature found by bisection on the energy balance equation. H2O vapor pressure from Antoine equation.

### Ballistic DSMC (`ballistic.c`)
Collisionless particle transport. Launches particles from the surface with Maxwell-Boltzmann velocities weighted by local temperature. Integrates orbits with Stormer-Verlet under the body's gravitational field. Tracks photodissociation lifetime decay and day/night surface freeze-out.

### Advanced DSMC (`advanced.c`)
Two-phase simulation with molecular collisions:
1. **Phase 1**: Collisionless run to collect velocity statistics per grid cell
2. **Phase 2**: Polynomial CDF fit to velocity distributions, then full DSMC with accumulative collision probability and elastic hard-sphere scattering

## Parallelization

All compute-intensive loops use OpenMP:
- Particle trajectories: `#pragma omp parallel for schedule(dynamic)`
- Ray tracing: `#pragma omp parallel for collapse(2)`
- Thermal solver: `#pragma omp parallel for collapse(2)`
- Density accumulation: `#pragma omp atomic`
- Thread-safe RNG via per-thread seeds (`rand_r`)

## Applicable Bodies

Change physical constants in `constants.h` to simulate different targets:

| Body | Atmosphere | Data Source |
|------|-----------|-------------|
| Ceres | H2O (transient) | Dawn mission |
| Europa | H2O, O2 | Europa Clipper (2030) |
| Enceladus | H2O (plumes) | Cassini legacy |
| Moon | Na, K, H2O | Chang'e, Artemis |
| Mercury | Na, Ca | BepiColombo |

## Original MATLAB Sources

| C Module | MATLAB File |
|----------|-------------|
| `shadow_maker.c` | `1_Shadow_maker(NObug0720).m` |
| `thermal_model.c` | `m0-1 thermal model.m` |
| `ballistic.c` | `m2_ballistic_motion_correct0227.m` |
| `advanced.c` | `test_ballistic.m` |

## License

Research code. Contact author for usage.
