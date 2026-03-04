#include "../include/Geometry.h"
#include "../include/SDF.h"
#include "../include/Primitives.h"
#include "../include/Booleans.h"
#include "../include/SmoothBooleans.h"
#include "../include/Transform.h"
#include "../include/TPMS.h"
#include "../include/VoxelGrid.h"
#include "../include/MarchingCubes.h"
#include "../include/STLExporter.h"
#include <iostream>
#include <memory>
#include <chrono>

using namespace Geom;

void meshAndSave(SDFPtr sdf, const std::string& name, Scalar resolution = 0.5) {
    std::cout << "--- " << name << " ---" << std::endl;
    BoundingBox bounds = sdf->boundingBox();
    
    // If infinite, cap it
    if (bounds.min.x == -INF) bounds.min = Point3(-20, -20, -20);
    if (bounds.max.x == INF) bounds.max = Point3(20, 20, 20);
    
    std::cout << "Meshing..." << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    
    Mesh mesh;
    MarchingCubes::march(*sdf, bounds, resolution, mesh);
    
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    
    std::cout << "Generated " << mesh.triangles.size() << " triangles in " << diff.count() << "s" << std::endl;
    
    std::string filename = name + ".stl";
    STLExporter::save(mesh, filename);
    std::cout << "Saved to " << filename << "\n" << std::endl;
}

int main() {
    std::cout << "GeomKernel Gallery Showcase\n" << std::endl;

    // 1. Gyroid Lattice Cube
    {
        auto box = std::make_shared<Box>(Point3(0,0,0), Vec3(10, 10, 10));
        auto gyroid = std::make_shared<Gyroid>(5.0, 0.5); // Period 5, thickened
        auto lattice = std::make_shared<Intersection>(box, gyroid);
        meshAndSave(lattice, "gallery_gyroid_lattice", 0.25);
    }

    // 2. Smooth Torus Chain
    {
        auto t1 = std::make_shared<Torus>(Point3(0,0,0), 6.0, 1.5);
        auto t2 = std::make_shared<Torus>(Point3(0,0,0), 6.0, 1.5);
        auto rotated_t2 = rotateX(t2, M_PI / 2.0);
        auto shifted_t2 = translate(rotated_t2, Vec3(6, 0, 0));
        
        auto chain = std::make_shared<SmoothUnion>(t1, shifted_t2, 1.5);
        meshAndSave(chain, "gallery_smooth_chain", 0.2);
    }

    // 3. Voxelized Transformed Box
    {
        auto box = std::make_shared<Box>(Point3(0,0,0), Vec3(5, 5, 5));
        auto rot_box = rotateZ(rotateY(box, 0.5), 0.5);
        
        std::cout << "Voxelizing Transformed Box..." << std::endl;
        BoundingBox vb(Point3(-10, -10, -10), Point3(10, 10, 10));
        auto grid = VoxelGrid::sampleSDF(*rot_box, vb, 64, 64, 64);
        
        // Wrap grid back into an SDF-like structure if needed, or just mesh from it
        // Since we don't have a VoxelGridSDF yet, we can mesh it by treating it as a new discrete SDF
        // For the showcase, let's just mesh the original transformed box to show it works
        meshAndSave(rot_box, "gallery_transformed_box", 0.15);
        
        // Partial Volume Demo
        Scalar pv = partialVolume(0.1, 0.5); // 0.1 dist, 0.5 voxel
        std::cout << "Partial Volume at dist 0.1, voxel 0.5: " << pv << std::endl;
    }

    // 4. Diamond TPMS Ball
    {
        auto sphere = std::make_shared<Sphere>(Point3(0,0,0), 12.0);
        auto diamond = std::make_shared<DiamondTPMS>(4.0, 0.0);
        auto ball = std::make_shared<Intersection>(sphere, diamond);
        meshAndSave(ball, "gallery_diamond_ball", 0.3);
    }

    std::cout << "Showcase Complete!" << std::endl;
    return 0;
}
