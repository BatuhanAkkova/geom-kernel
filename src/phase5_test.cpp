#include "../include/Geometry.h"
#include "../include/SDF.h"
#include "../include/Primitives.h"
#include "../include/Modifiers.h"
#include "../include/Booleans.h"
#include "../include/MarchingCubes.h"
#include "../include/STLExporter.h"
#include "../include/Fields.h"
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
    // Sphere with a RampField applied to offset.
    {
        auto sphere = std::make_shared<Sphere>(Point3(0, 0, 0), 10.0);
        
        // Ramp along Y axis, from -15 to +15.
        // Value goes from 0.0 to 5.0.
        // So at bottom (y=-15), offset is 0. Radius = 10.
        // At top (y=15), offset is 5. Radius = 15.
        auto ramp = std::make_shared<RampField>(
            Point3(0, -15, 0), 
            Vec3(0, 1, 0), 
            30.0, // length 
            0.0,  // startVal
            5.0,  // endVal
            true  // clamp
        );
        
        auto egg = std::make_shared<FieldOffset>(sphere, ramp);

        run_test("Egg (FieldOffset)", egg, "phase5_egg.stl");
    }

    // 2. Pillowed Box
    // Box with RadialField offset.
    {
        auto box = std::make_shared<Box>(Point3(30, 0, 0), Vec3(10, 10, 10));
        
        // Radial field from center of box.
        // Radius 15.
        // Center: offset = 0.
        // Edge: offset = -2. (Shrink at edges? No, usually pillow means inflate center)
        // Let's try: Center offset = 2, Edge offset = 0.
        // Actually, FieldOffset is: dist - field.
        // New dist = dist - field.
        // If field > 0, dist becomes smaller -> surface moves outward (inflate).
        auto radial = std::make_shared<RadialField>(
            Point3(30, 0, 0),
            15.0,
            2.0, // Center value (inflate by 2)
            0.0, // Edge value
            true
        );
        
        auto pillow = std::make_shared<FieldOffset>(box, radial);
        
        run_test("Pillowed Box", pillow, "phase5_pillow.stl");
    }

    // 3. Twisted Bar
    {
        auto bar = std::make_shared<Box>(Point3(-30, 0, 0), Vec3(5, 20, 5));
        auto twisted = std::make_shared<Twist>(bar, 0.1); // Twist 0.1 rad per unit Y

        run_test("Twisted Bar", twisted, "phase5_twist.stl");
    }

    // 4. Mixed Field (Sphere with Noise-like interference)
    // Actually, let's use FieldAdd to combine two RampFields for a cross-gradient.
    {
        auto sphere = std::make_shared<Sphere>(Point3(0, 30, 0), 10.0);
        
        auto rampX = std::make_shared<RampField>(Point3(0,30,0), Vec3(1,0,0), 10.0, 0.0, 2.0);
        auto rampZ = std::make_shared<RampField>(Point3(0,30,0), Vec3(0,0,1), 10.0, 0.0, 2.0);
        
        auto combined = std::make_shared<FieldAdd>(rampX, rampZ);
        auto perturbed = std::make_shared<FieldOffset>(sphere, combined);
        
        run_test("Perturbed Sphere (FieldAdd)", perturbed, "phase5_mixed.stl");
    }


    return 0;
}
