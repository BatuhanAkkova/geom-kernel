#include "Primitives.h"
#include "MarchingCubes.h"
#include "STLExporter.h"
#include <iostream>
#include <vector>
#include <memory>

using namespace Geom;

void save_primitive(const std::string& name, SDF& sdf, Scalar resolution = 0.5) {
    std::cout << "Meshing " << name << "..." << std::endl;
    Mesh mesh;
    BoundingBox bounds = sdf.boundingBox();
    // Expand bounds slightly to ensure surface is captured
    bounds.expand(bounds.min - Vec3(1, 1, 1));
    bounds.expand(bounds.max + Vec3(1, 1, 1));
    
    MarchingCubes::march(sdf, bounds, resolution, mesh);
    std::string filename = "primitive_" + name + ".stl";
    if (STLExporter::save(mesh, filename)) {
        std::cout << "Saved " << filename << " (" << mesh.triangles.size() << " triangles)" << std::endl;
    } else {
        std::cerr << "Failed to save " << filename << std::endl;
    }
}

int main() {
    std::cout << "Running Primitive Generation Tests..." << std::endl;

    // 1. Sphere
    Sphere sphere(Point3(0, 0, 0), 10.0);
    save_primitive("sphere", sphere);

    // 2. Box
    Box box(Point3(0, 0, 0), Vec3(8, 5, 3));
    save_primitive("box", box);

    // 3. Cylinder
    Cylinder cylinder(Point3(0, 0, 0), 5.0, 15.0);
    save_primitive("cylinder", cylinder);

    // 4. Torus
    Torus torus(Point3(0, 0, 0), 10.0, 3.0);
    save_primitive("torus", torus);

    // 5. Plane (Note: Plane is infinite, so we use a bounding box for meshing)
    Plane plane(Vec3(0, 0, 1), 0.0); // Z=0 plane
    BoundingBox plane_bounds(Point3(-10, -10, -1), Point3(10, 10, 1));
    Mesh plane_mesh;
    MarchingCubes::march(plane, plane_bounds, 0.5, plane_mesh);
    STLExporter::save(plane_mesh, "primitive_plane.stl");
    std::cout << "Saved primitive_plane.stl" << std::endl;

    std::cout << "Primitive tests completed." << std::endl;
    return 0;
}
