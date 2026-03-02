#include <iostream>
#include <iomanip>
#include <vector>
#include <cassert>
#include "Primitives.h"
#include "Booleans.h"
#include "Autodiff.h"
#include "Optimization.h"

using namespace Geom;

void test_autodiff_sphere() {
    std::cout << "Testing Autodiff: Sphere Sensitivity..." << std::endl;
    
    Point3 center(0, 0, 0);
    Scalar radius = 5.0;
    Sphere sphere(center, radius);

    Point3 query(3, 4, 0); // Distance should be sqrt(3^2 + 4^2) - 5 = 0
    
    // Sensitivity w.r.t radius
    // f(r) = length(p-c) - r
    // df/dr = -1
    
    auto differentiable_sphere = [&](Point3D p, DualScalar r) {
        DualScalar dx = p.x - center.x;
        DualScalar dy = p.y - center.y;
        DualScalar dz = p.z - center.z;
        return sqrt(dx*dx + dy*dy + dz*dz) - r;
    };

    Scalar sensitivity = Optimization::computeSensitivity(query, differentiable_sphere, radius);
    
    std::cout << "  Query Point: (3, 4, 0), Radius: 5.0" << std::endl;
    std::cout << "  Calculated Sensitivity (df/dr): " << sensitivity << " (Expected: -1.0)" << std::endl;
    
    assert(std::abs(sensitivity + 1.0) < EPSILON);
    std::cout << "  [PASS] Sphere sensitivity correct." << std::endl << std::endl;
}

void test_parallel_sweep() {
    std::cout << "Testing Parallel Design Sweep..." << std::endl;

    std::vector<Scalar> radii = {1.0, 2.0, 3.0, 4.0, 5.0};
    
    auto generator = [](Scalar r) {
        return std::make_shared<Sphere>(Point3(0,0,0), r);
    };

    auto objective = [](SDFPtr sdf) {
        // Use distance at a fixed point as objective
        return sdf->eval(Point3(10, 0, 0));
    };

    auto results = Optimization::sweep(radii, generator, objective);

    std::cout << std::setw(10) << "Radius" << std::setw(15) << "Objective" << std::endl;
    for (const auto& res : results) {
        std::cout << std::setw(10) << res.parameterValue << std::setw(15) << res.objectiveValue << std::endl;
        // Objective at (10,0,0) for sphere at (0,0,0) is 10 - r
        assert(std::abs(res.objectiveValue - (10.0 - res.parameterValue)) < EPSILON);
    }

    std::cout << "  [PASS] Parallel sweep produced correct results." << std::endl << std::endl;
}

int main() {
    try {
        test_autodiff_sphere();
        test_parallel_sweep();
        std::cout << "All Phase 8 tests passed!" << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Test failed: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
