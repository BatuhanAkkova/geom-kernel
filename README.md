# 🌀 GeomKernel: High-Performance SDF Geometry Kernel

[![C++20](https://img.shields.io/badge/C%2B%2B-20-blue.svg)](https://en.cppreference.com/w/cpp/20)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A modern, high-performance, header-only C++ library for **Signed Distance Function (SDF)** based geometry modeling. Designed for precision engineering, physics-informed structure design, and advanced manufacturing (WAAM/AM).

---

## Key Features

- **Engineering Precision**: Uses `double` (64-bit float) for all geometric calculations to ensure high-fidelity modeling.
- **Physics-Informed Design**: Seamlessly integrate analytical physics fields (heat transfer, pressure vessels) into your geometric primitives.
- **Advanced Modeling**:
  - **Primitives**: Spheres, Boxes, Cylinders, Toruses, and more.
  - **Smooth Booleans**: Create high-quality junctions with controllable blending.
  - **ShapeKernel Layer**: High-level engineering primitives (Pipes, Ribs, Cooling Channels, Shells).
- **Manufacturing Ready**: Built-in validation for Wire Arc Additive Manufacturing (WAAM), including overhang checks and curvature constraints.
- **Differentiable Core**: Supports Automatic Differentiation (Autodiff) for optimization-driven geometry and gradient-based solvers.
- **Header-Only & Lightweight**: Easily integrated as a git submodule with zero external dependencies (uses standard C++20).

---

## Architecture

`GeomKernel` represents geometry as a continuous scalar field where the value at any point `p` indicates the distance to the nearest surface:
- $f(p) < 0$: Inside the object
- $f(p) > 0$: Outside the object
- $f(p) = 0$: Surface interface

### Primary Components:
- **`SDF`**: Abstract base class for all signed distance functions.
- **`Primitives`**: Standard leaf nodes for the CSG tree.
- **`Booleans`**: Union, Intersection, and Difference operations.
- **`Meshers`**: High-performance Dual Contouring and Marching Cubes implementation for high-quality mesh generation from SDF fields.

---

## Quick Start

### 1. Integration (Git Submodule)

Add `GeomKernel` to your project's `extern` directory:

```bash
git submodule add https://github.com/your-org/geom-kernel extern/geom-kernel
```

In your `CMakeLists.txt`:

```cmake
add_subdirectory(extern/geom-kernel)
target_link_libraries(your_project PRIVATE GeomKernel)
```

### 2. Basic Example

Creating a simple box with a spherical cutout:

```cpp
#include "GeomKernel/Primitives.h"
#include "GeomKernel/Booleans.h"
#include <memory>

int main() {
    using namespace Geom;

    // 1. Define primitives
    auto box = std::make_shared<Box>(Point3(0,0,0), Vec3(10, 10, 10)); // 20x20x20 box
    auto sphere = std::make_shared<Sphere>(Point3(0,0,0), 12.0);      // r=12 sphere

    // 2. Perform Boolean operations
    auto result = std::make_shared<Difference>(box, sphere);

    // 3. Evaluate distance at a point
    Scalar dist = result->eval(Point3(5, 5, 5));
    
    return 0;
}
```

---

## Manufacturing Awareness (WAAM)

GeomKernel enforces real-world manufacturability:

```cpp
#include "GeomKernel/Manufacturing.h"

// Validate a geometry for WAAM printing
WAAMValidator validator;
validator.min_bead_width = 4.0; 
validator.max_overhang_angle = 45.0; // Degrees

auto report = validator.analyze(my_geometry_sdf);
if (!report.is_buildable) {
    std::cout << "Design contains unbuildable regions: " << report.reason << std::endl;
}
```

---

## Roadmap

- [ ] **GPU Acceleration**: OpenCL/CUDA backend for massively parallel evaluation.
- [ ] **Python Bindings**: Expose core kernel functionality to the Python ecosystem.
- [ ] **Lattice & TPMS**: Specialized library for Triply Periodic Minimal Surfaces.
- [ ] **Hessians**: Add second-order derivatives for advanced optimization.
- [ ] **Differentiable Kernel**: Make the kernel differentiable for optimization-driven geometry and gradient-based solvers.

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
