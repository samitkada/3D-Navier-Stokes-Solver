## Physics & Engineering Models

This project implements a three-dimensional computational fluid dynamics (CFD) solver for the incompressible Navier-Stokes equations using finite difference methods and pressure-velocity coupling. The simulation models the classical lid-driven cavity problem, a benchmark flow configuration widely used to validate numerical fluid dynamics algorithms and study vortex formation in confined flows.

The model simulates a cubic cavity filled with fluid where the top wall moves at a constant velocity while all other walls remain stationary. The moving lid transfers momentum into the fluid, generating complex three-dimensional vortex structures driven by the interaction of inertial and viscous forces.

---

### Incompressible Flow Modeling

The simulation assumes an incompressible fluid, meaning fluid density remains constant throughout the domain. Under this assumption, the velocity field must satisfy mass conservation at every point within the cavity.

The solver simultaneously computes:

- Three-dimensional velocity components
- Pressure distribution
- Vorticity field
- Velocity magnitude
- Recirculating vortex structures

This allows the complete flow field to be reconstructed and analyzed throughout the computational domain.

---

### Lid-Driven Cavity Benchmark

The lid-driven cavity is one of the most widely studied benchmark problems in computational fluid dynamics.

The physical setup consists of:

- A cubic fluid domain
- Stationary side walls
- Stationary bottom wall
- A moving top lid with prescribed velocity

As the lid moves, momentum is transferred into the fluid, producing primary and secondary vortex structures that evolve toward a steady-state flow solution.

Because no analytical solution exists for the general three-dimensional case, the problem is commonly used to verify numerical CFD methods and compare solver performance.

---

### Reynolds Number Effects

Flow behavior is controlled by the Reynolds number, which represents the ratio of inertial forces to viscous forces.

Lower Reynolds numbers produce:

- Smooth flow structures
- Strong viscous damping
- Stable vortex formation

Higher Reynolds numbers produce:

- Stronger recirculation regions
- Sharper velocity gradients
- More complex vortex interactions
- Increased numerical difficulty

The simulation allows investigation of how changing Reynolds number affects overall flow structure and vortex development.

---

### Numerical Discretization

The governing fluid equations are discretized using second-order finite difference methods on a structured Cartesian grid.

The computational domain is divided into thousands of control points where velocity and pressure are solved numerically.

Spatial discretization includes:

- Central difference gradients
- Finite difference Laplacian operators
- Three-dimensional structured grid generation
- Explicit time integration

This approach provides a balance between computational efficiency and solution accuracy.

---

### Pressure-Velocity Coupling

One of the primary challenges in incompressible CFD is enforcing mass conservation while solving for pressure and velocity simultaneously.

The solver uses a projection-based pressure correction approach similar to SIMPLE and Chorin fractional-step methods.

Each time step consists of:

1. Computing intermediate velocities from momentum equations
2. Solving a pressure Poisson equation
3. Correcting the velocity field using pressure gradients
4. Enforcing incompressibility throughout the domain

This process ensures the final velocity field satisfies both momentum conservation and continuity requirements.

---

### Vorticity Analysis

To characterize rotational flow structures, the solver computes the full three-dimensional vorticity field.

Vorticity is a measure of local fluid rotation and is commonly used to identify:

- Primary cavity vortices
- Secondary recirculation regions
- Shear-layer development
- Vortex core structures

The simulation calculates vorticity magnitude throughout the domain and visualizes coherent vortex structures using three-dimensional isosurfaces.

---

### Computational Methods

The solver combines several fundamental CFD techniques, including:

- Incompressible Navier-Stokes modeling
- Finite difference discretization
- Fractional-step pressure correction
- Pressure Poisson equation solution
- Explicit time marching
- Three-dimensional vortex identification
- Structured Cartesian mesh generation

These methods form the foundation of many modern CFD solvers used in aerospace, mechanical engineering, and fluid mechanics research.

---

### Engineering Applications

Although the lid-driven cavity is a benchmark problem, the numerical methods used in this project are directly applicable to real engineering systems involving fluid flow.

Relevant applications include:

- Aircraft aerodynamic analysis
- Turbomachinery flow modeling
- Combustion chamber flow studies
- HVAC and ventilation systems
- Internal flow passages
- Mixing and transport processes
- General CFD solver development

The project demonstrates core computational fluid dynamics concepts used in professional aerospace and mechanical engineering simulation environments.

---

### Visualization & Post-Processing

The solver generates multiple flow-field diagnostics for engineering analysis, including:

- Velocity magnitude contours
- Velocity vector fields
- Centerline velocity profiles
- Three-dimensional velocity isosurfaces
- Vorticity magnitude contours
- Vortex-core visualizations

These outputs provide insight into flow development, momentum transport, and vortex formation within the cavity and enable qualitative and quantitative evaluation of solver performance.
