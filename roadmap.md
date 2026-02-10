# SDF Geometry Kernel Roadmap
*picoGK-inspired, C++-based computational geometry kernel*

---

## Vision

Build a **field-driven implicit geometry kernel** capable of generating
manufacturable solids for advanced engineering workflows
(additive manufacturing, WAAM, motors, aerospace, robotics, etc.).

This is **not a CAD system**.
It is a **computational engineering kernel**.

---

## Core Principles

- Geometry = result of equations and fields
- Signed Distance Functions (SDF) as the core representation
- Deterministic, robust, and automatable
- Manufacturing-aware by design
- No sketches, no feature trees, no GUI dependency

---

## Technology Choices

- **Language:** C++ (core kernel)
- **Paradigm:** Implicit geometry (SDF)
- **Output:** Triangle mesh (STL / OBJ)
- **Acceleration:** Multithreading, SIMD, spatial pruning
- **Future:** Python bindings, GPU acceleration

---

## Phase 0 — Architectural Decisions

**Goal:** Lock design constraints early.

### Tasks
- Define SDF-based representation
- Define coordinate system & units
- Decide resolution & tolerance strategy
- Define kernel responsibility boundaries

### Deliverables
- Architecture document
- Core interfaces defined

---

## Phase 1 — Minimal Viable Kernel (MVP)

**Goal:** Generate a solid and export it.

### Tasks
- Define `SDF` interface
- Implement basic primitives:
  - Sphere
  - Box
  - Cylinder
  - Plane
- Implement Boolean operations:
  - Union
  - Intersection
  - Difference
- Uniform grid sampling
- Marching Cubes meshing
- STL export

### Done When
- STL opens correctly in a mesh viewer
- Booleans work reliably

---

## Phase 2 — Engineering-Grade SDF Behavior

**Goal:** Make geometry numerically robust.

### Tasks
- Smooth min/max blending
- Offset operations
- Shell and thickness control
- Distance normalization / re-distancing
- Gradient evaluation

### Done When
- Fillets and offsets behave predictably
- No self-intersections from offsets

---

## Phase 3 — Performance & Scaling

**Goal:** Evaluate millions of points efficiently.

### Tasks
- Bounding volumes for SDF nodes
- Early-out distance checks
- Multithreaded evaluation
- Cache repeated evaluations

### Done When
- Linear scaling with CPU cores
- Large scenes evaluate efficiently

---

## Phase 4 — Adaptive & Sparse Sampling

**Goal:** High resolution only where needed.

### Tasks
- Octree or adaptive grid
- Zero-crossing detection
- Dual Contouring (preferred)
- Crack-free meshing
- Feature-preserving refinement

### Done When
- Sharp features survive
- Memory usage scales well

---

## Phase 5 — Field-Driven Geometry

**Goal:** Geometry reacts to physics and constraints.

### Tasks
- Define `Field` interface
- Support spatially varying parameters
- Couple geometry to:
  - Temperature fields
  - Stress fields
  - Distance fields
- Smooth parameter transitions

### Done When
- Geometry changes continuously with fields
- No discontinuities from field coupling

---

## Phase 6 — Manufacturing Awareness (WAAM / AM)

**Goal:** Generate buildable geometry.

### Tasks
- Overhang detection
- Build-direction constraints
- Minimum curvature radius
- Bead-width / layer-height constraints
- Layer-wise validation

### Done When
- Geometry respects process limits
- No unbuildable regions are produced

---

## Phase 7 — ShapeKernel Layer

**Goal:** Abstract raw SDF math into engineering primitives.

### Tasks
- Implement reusable shapes:
  - Pipes
  - Shells
  - Smooth junctions
  - Structural ribs
  - Cooling channels
- Enforce safe defaults
- Hide raw Boolean logic from users

### Done When
- Users no longer write raw SDF equations

---

## Phase 8 — Optimization & Differentiation (Advanced)

**Goal:** Enable design by optimization.

### Tasks
- Sensitivity computation
- Finite differences / autodiff hooks
- Batch geometry generation
- Parallel design sweeps

### Done When
- Geometry is used inside optimization loops

---

## Stretch Goals

- GPU evaluation (CUDA / OpenCL)
- Python bindings
- Lattice & TPMS library
- Quasicrystal structures
- Differentiable geometry kernel

---

## Final Definition of Success

> A deterministic, field-driven geometry kernel that produces
> manufacturable, physics-informed solids automatically.

---

