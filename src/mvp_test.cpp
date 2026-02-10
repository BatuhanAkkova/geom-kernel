#include <iostream>
#include "Primitives.h"
#include "MarchingCubes.h"
#include "STLExporter.h"

int main() {
    std::cout << "Creating Geometry..." << std::endl;
    
    // Create a sphere
    Geom::Sphere sphere(Geom::Point3(0,0,0), 10.0);
    
    std::cout << "Meshing..." << std::endl;
    Geom::Mesh mesh;
    Geom::BoundingBox bounds = sphere.boundingBox();
    
    // Expand bounds slightly to ensure surface is captured
    bounds.expand(Geom::Point3(bounds.min - Geom::Vec3(1,1,1)));
    bounds.expand(Geom::Point3(bounds.max + Geom::Vec3(1,1,1)));
    
    Geom::MarchingCubes::march(sphere, bounds, 0.5, mesh); // 0.5mm resolution
    
    std::cout << "Generated " << mesh.triangles.size() << " triangles." << std::endl;
    
    if (Geom::STLExporter::save(mesh, "sphere.stl")) {
        std::cout << "Saved sphere.stl" << std::endl;
    } else {
        std::cerr << "Failed to save STL" << std::endl;
        return 1;
    }

    return 0;
}
