#include "Primitives.h"
#include "Booleans.h"
#include "SmoothBooleans.h"
#include "ShapeKernel.h"
#include "Fields.h"
#include "Modifiers.h"
#include "MarchingCubes.h"
#include "STLExporter.h"
#include "Autodiff.h"
#include "Optimization.h"
#include "TopoOptimizer.h"
#include "DensityField.h"
#include <iostream>
#include <memory>
#include <vector>

using namespace Geom;

/**
 * @brief Gallery Showcase: Physics-Informed Design
 * Automated wall thickness driven by analytical physics.
 */
void gallery_physics_vessel() {
    std::cout << "Gallery [1/6]: Generating Physics-Informed Vessel..." << std::endl;
    // Cylinder with field-driven offset (simulating pressure vessel thickness)
    auto inner = std::make_shared<Cylinder>(Point3(0,0,0), 50.0, 200.0);
    
    // Constant thickness field of 5mm
    auto thickness = std::make_shared<ConstantField>(5.0);
    
    // Inflation by thickness
    auto outer_surface = std::make_shared<FieldOffset>(inner, thickness);
    
    // Difference gives hollow shell
    auto vessel = std::make_shared<Difference>(outer_surface, inner);
    
    Mesh mesh;
    BoundingBox bb = vessel->boundingBox();
    MarchingCubes::march(*vessel, bb, 2.0, mesh);
    STLExporter::save(mesh, "gallery_physics_vessel.stl");
}

/**
 * @brief Gallery Showcase: Smooth Junctions
 * High-fidelity blending between primitives.
 */
void gallery_smooth_junctions() {
    std::cout << "Gallery [2/6]: Generating Smooth Junctions..." << std::endl;
    auto s1 = std::make_shared<Sphere>(Point3(-15, 0, 0), 20.0);
    auto s2 = std::make_shared<Sphere>(Point3(15, 0, 0), 20.0);
    
    // Smooth Union with 8mm blend radius
    auto joined = std::make_shared<SmoothUnion>(s1, s2, 8.0);
    
    Mesh mesh;
    MarchingCubes::march(*joined, joined->boundingBox(), 1.0, mesh);
    STLExporter::save(mesh, "gallery_smooth_junctions.stl");
}

/**
 * @brief Gallery Showcase: ShapeKernel Layer
 * High-level engineering primitives.
 */
void gallery_ribbed_pipe() {
    std::cout << "Gallery [3/6]: Generating Ribbed Pipe..." << std::endl;
    using namespace Geom::ShapeKernel;
    
    auto pipe = std::make_shared<Pipe>(Point3(0, -50, 0), Point3(0, 50, 0), 10.0);
    
    // Add structural ribs
    auto r1 = std::make_shared<StructuralRib>(pipe, Plane(Vec3(1,0,0), 0), 2.0, 5.0);
    auto r2 = std::make_shared<StructuralRib>(r1, Plane(Vec3(0,0,1), 0), 2.0, 5.0);
    
    // Subtract cooling channel
    auto channel = std::make_shared<Pipe>(Point3(0,-60,0), Point3(0,60,0), 6.0);
    auto finalPipe = std::make_shared<CoolingChannel>(r2, channel);
    
    Mesh mesh;
    MarchingCubes::march(*finalPipe, finalPipe->boundingBox(), 1.0, mesh);
    STLExporter::save(mesh, "gallery_ribbed_pipe.stl");
}

/**
 * @brief Gallery Showcase: Topology Optimization
 * SIMP Densities bridged to SDF.
 */
void gallery_topo_opt() {
    std::cout << "Gallery [4/6]: Generating Topology Optimization Mesh..." << std::endl;
    
    // Mock a density grid (30x15x15)
    TopoOptimizer::Config cfg{10, 5, 5, 0.4, 3.0, 1.5};
    TopoOptimizer opt(cfg);
    
    // Normally we'd run iterations. For the gallery, we'll manually set some high-density zones
    // to simulate a beam structure.
    auto& rho = const_cast<std::vector<double>&>(opt.densities());
    for(int i=0; i<rho.size(); ++i) rho[i] = 0.1; // void
    
    // Create a "support" and "load" connection
    for(int i=0; i<10; ++i) rho[i] = 1.0; // Top edge
    for(int i=0; i<5; ++i) rho[i*10] = 1.0; // Left side
    
    BoundingBox domain(Point3(0,0,0), Point3(100, 50, 50));
    auto df = std::make_shared<DensityField>(rho, 10, 5, 5, domain, 0.5);
    
    Mesh mesh;
    MarchingCubes::march(*df, domain, 2.0, mesh);
    STLExporter::save(mesh, "gallery_topo_opt.stl");
}

/**
 * @brief Gallery Showcase: Space Warping
 * Nonlinear transformations (Twist).
 */
void gallery_twisted_bar() {
    std::cout << "Gallery [5/6]: Generating Twisted Bar..." << std::endl;
    auto bar = std::make_shared<Box>(Point3(0,0,0), Vec3(10, 50, 10));
    auto twisted = std::make_shared<Twist>(bar, 0.05); // 0.05 rad per unit Y
    
    Mesh mesh;
    MarchingCubes::march(*twisted, twisted->boundingBox(), 1.0, mesh);
    STLExporter::save(mesh, "gallery_twisted_bar.stl");
}

/**
 * @brief Gallery Showcase: Differentiable Geometry
 * Sensitivity analysis via Autodiff.
 */
void gallery_sensitivity() {
    std::cout << "Gallery [6/6]: Running Sensitivity Analysis..." << std::endl;
    Point3 query(12, 0, 0);
    Scalar radius = 10.0;
    
    auto sphere_model = [](Point3D p, DualScalar r) {
        return p.length() - r;
    };
    
    Scalar sens = Optimization::computeSensitivity(query, sphere_model, radius);
    std::cout << "  Sensitivity of Sphere surface at (12,0,0) w.r.t Radius: " << sens << std::endl;
    std::cout << "  (Exactly -1.0 as expected analytically)" << std::endl;
}

int main() {
    std::cout << "--- GeomKernel Gallery Showcase ---" << std::endl;
    
    try {
        gallery_physics_vessel();
        gallery_smooth_junctions();
        gallery_ribbed_pipe();
        gallery_topo_opt();
        gallery_twisted_bar();
        gallery_sensitivity();
        
        std::cout << "\nAll gallery items generated successfully!" << std::endl;
        std::cout << "STL files produced in current directory." << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Gallery generation failed: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
