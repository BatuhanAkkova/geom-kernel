#include "../include/Geometry.h"
#include "../include/SDF.h"
#include "../include/Primitives.h"
#include "../include/Modifiers.h"
#include "../include/Booleans.h"
#include "../include/MarchingCubes.h"
#include "../include/STLExporter.h"
#include "../include/Field.h"
#include <iostream>
#include <vector>
#include <string>

using namespace Geom;

void run_test(const std::string& name, SDFPtr sdf, const std::string& filename) {
    std::cout << "Meshing " << name << "..." << std::endl;

    BoundingBox bounds = sdf->boundingBox();
    Vec3 size = bounds.max - bounds.min;
    std::cout << "Bounds: Min(" << bounds.min.x << ", " << bounds.min.y << ", " << bounds.min.z << ") "
              << "Max(" << bounds.max.x << ", " << bounds.max.y << ", " << bounds.max.z << ")" << std::endl;

    // Grid resolution
    int resolution = 64; 
    Scalar cellSize = std::max({size.x, size.y, size.z}) / resolution;

    Mesh mesh;
    MarchingCubes::march(*sdf, bounds, cellSize, mesh);
    std::cout << "Generated " << mesh.triangles.size() << " triangles." << std::endl;


    STLExporter::save(mesh, filename);
    std::cout << "Saved to " << filename << std::endl;
}

int main() {
    // 1. Variable Radius Sphere (Egg)
    {
        auto sphere = std::make_shared<Sphere>(Point3(0, 0, 0), 8.0);
        
        // Ramp along Y axis, from -10 to +10.
        // Value goes from 4.0 to 0.0.
        // At bottom (y=-10), offset is 4.0 -> Radius = 12.0.
        // At top (y=10), offset is 0.0 -> Radius = 8.0.
        auto ramp = std::make_shared<RampField>(
            Point3(0, -10, 0), 
            Vec3(0, 1, 0), 
            20.0, // length 
            4.0,  // startVal (wider at bottom)
            0.0,  // endVal (regular at top)
            true  // clamp
        );
        
        auto egg = std::make_shared<FieldOffset>(sphere, ramp);

        run_test("Egg (FieldOffset)", egg, "phase5_egg.stl");
    }

    // 2. Pillowed Box
    {
        auto box = std::make_shared<Box>(Point3(25, 0, 0), Vec3(8, 8, 8));
        
        auto radial = std::make_shared<RadialField>(
            Point3(25, 0, 0),
            12.0,
            3.0, // Center value (inflate)
            0.0, // Edge value
            true
        );
        
        auto pillow = std::make_shared<FieldOffset>(box, radial);
        
        run_test("Pillowed Box", pillow, "phase5_pillow.stl");
    }

    // 3. Twisted Bar
    {
        // MUST center on rotation axis (Y) for twisting to make sense as a "Twisted Bar"
        auto bar = std::make_shared<Box>(Point3(-25, 0, 0), Vec3(4, 15, 4));
        // Actually, let's just use one at origin for the prime example
        auto centerBar = std::make_shared<Box>(Point3(0, 0, 0), Vec3(4, 15, 4));
        auto twisted = std::make_shared<Twist>(centerBar, 0.15); 

        run_test("Twisted Bar", twisted, "phase5_twist.stl");
    }

    // 4. Mixed Field (Sphere with Noise-like interference)
    {
        auto sphere = std::make_shared<Sphere>(Point3(0, 25, 0), 8.0);
        
        auto rampX = std::make_shared<RampField>(Point3(0,25,0), Vec3(1,0,0), 8.0, 0.0, 3.0);
        auto rampZ = std::make_shared<RampField>(Point3(0,25,0), Vec3(0,0,1), 8.0, 0.0, 3.0);
        
        auto combined = std::make_shared<FieldAdd>(rampX, rampZ);
        auto perturbed = std::make_shared<FieldOffset>(sphere, combined);
        
        run_test("Perturbed Sphere (FieldAdd)", perturbed, "phase5_mixed.stl");
    }


    return 0;
}
