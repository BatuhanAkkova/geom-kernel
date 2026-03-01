#include "ShapeKernel.h"
#include "STLExporter.h"
#include "MarchingCubes.h"
#include "Manufacturing.h"
#include <iostream>
#include <cassert>
#include <iomanip>

using namespace Geom;
using namespace Geom::ShapeKernel;

int main() {
    std::cout << "Starting Phase 7 — ShapeKernel Layer Tests..." << std::endl;

    // 1. Test Pipe
    auto pipe = std::make_shared<Pipe>(Point3(0, 0, 0), Point3(10, 10, 0), 2.0);
    assert(pipe->eval(Point3(0, 0, 0)) < 0);
    assert(pipe->eval(Point3(5, 5, 0)) < 0);
    assert(pipe->eval(Point3(11, 11, 0)) > 0);
    std::cout << "  - Pipe test passed." << std::endl;

    // 2. Test SmoothJunction
    auto sphere = std::make_shared<Sphere>(Point3(0, 0, 0), 5.0);
    auto cylinder = std::make_shared<Cylinder>(Point3(0, 0, 0), 2.0, 15.0);
    auto junction = std::make_shared<SmoothJunction>(sphere, cylinder, 2.0);
    // Smooth junction should be "larger" than raw union
    Scalar d_union = std::min(sphere->eval(Point3(4, 4, 0)), cylinder->eval(Point3(4, 4, 0)));
    Scalar d_smooth = junction->eval(Point3(4, 4, 0));
    assert(d_smooth <= d_union); 
    std::cout << "  - SmoothJunction test passed." << std::endl;

    // 3. Test StructuralRib
    auto baseBox = std::make_shared<Box>(Point3(0,0,0), Vec3(10, 10, 2));
    auto rib = std::make_shared<StructuralRib>(baseBox, Plane(Vec3(1,0,0), 0), 2.0, 5.0);
    // Point on the rib (x=0, y=0, z=5)
    assert(rib->eval(Point3(0, 0, 4)) < 0);
    // Point outside the rib height (z=10)
    assert(rib->eval(Point3(0, 0, 10)) > 0);
    std::cout << "  - StructuralRib test passed." << std::endl;

    // 4. Test CoolingChannel
    auto body = std::make_shared<Box>(Point3(0,0,0), Vec3(10,10,10));
    auto channelPath = std::make_shared<Pipe>(Point3(-15, 0, 0), Point3(15, 0, 0), 2.0);
    auto block = std::make_shared<CoolingChannel>(body, channelPath);
    // Center should now be empty (outside) due to subtraction
    assert(block->eval(Point3(0,0,0)) > 0);
    std::cout << "  - CoolingChannel test passed." << std::endl;

    // 5. Test EngineeringShell
    auto shell = std::make_shared<EngineeringShell>(sphere, 1.0, true);
    // Inside the sphere but deeper than 1.0 (e.g. at origin) should be OUTSIDE the shell wall
    assert(shell->eval(Point3(0,0,0)) > 0);
    // On the wall (dist -0.5 from surface)
    assert(shell->eval(Point3(4.5, 0, 0)) < 0);
    std::cout << "  - EngineeringShell test passed." << std::endl;

    // 6. Integration: Ribbed Pipe with Cooling Channel
    std::cout << "Generating Integration Model: Ribbed Pipe..." << std::endl;
    auto mainPipe = std::make_shared<Pipe>(Point3(0, -20, 0), Point3(0, 20, 0), 5.0);
    auto rib1 = std::make_shared<StructuralRib>(mainPipe, Plane(Vec3(1, 0, 0), 0), 1.0, 3.0);
    auto rib2 = std::make_shared<StructuralRib>(rib1, Plane(Vec3(0, 0, 1), 0), 1.0, 3.0);
    
    auto coolant = std::make_shared<Pipe>(Point3(0, -25, 0), Point3(0, 25, 0), 3.0);
    auto finalModel = std::make_shared<CoolingChannel>(rib2, coolant);

    // Mesh and Export
    Mesh mesh;
    MarchingCubes::march(*finalModel, finalModel->boundingBox(), 0.5, mesh);
    
    STLExporter::save(mesh, "phase7_ribbed_pipe.stl");
    std::cout << "  - Integration model exported to phase7_ribbed_pipe.stl" << std::endl;

    // 7. Manufacturing Validation of ShapeKernel Primitives
    std::cout << "Running Manufacturing Validation on ShapeKernel Primitives..." << std::endl;
    ManufacturingConfig config;
    config.maxOverhangAngle = 45.0;
    config.minCurvatureRadius = 5.0;
    config.minBeadWidth = 3.0;
    ManufacturingValidator validator(config);

    // Validate the ribbed pipe
    auto issues = validator.validate(finalModel, finalModel->boundingBox(), 1.0);
    std::cout << "  - Ribbed Pipe: Found " << issues.size() << " manufacturing issues." << std::endl;
    // We expect some overhang issues due to the horizontal ribs if they were oriented that way, 
    // but the point is it works with the validator.

    std::cout << "All Phase 7 tests passed!" << std::endl;
    return 0;
}
