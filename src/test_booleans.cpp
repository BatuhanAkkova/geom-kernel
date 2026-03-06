#include "Primitives.h"
#include "Booleans.h"
#include "MarchingCubes.h"
#include "STLExporter.h"
#include <iostream>
#include <memory>

using namespace Geom;

void save_boolean(const std::string& name, SDF& sdf, Scalar resolution = 0.5) {
    std::cout << "Meshing " << name << "..." << std::endl;
    Mesh mesh;
    BoundingBox bounds = sdf.boundingBox();
    bounds.expand(bounds.min - Vec3(1, 1, 1));
    bounds.expand(bounds.max + Vec3(1, 1, 1));
    
    MarchingCubes::march(sdf, bounds, resolution, mesh);
    std::string filename = "boolean_" + name + ".stl";
    STLExporter::save(mesh, filename);
    std::cout << "Saved " << filename << std::endl;
}

int main() {
    std::cout << "Running Boolean Calculation Tests..." << std::endl;

    auto s1 = std::make_shared<Sphere>(Point3(0, 0, 0), 10.0);
    auto s2 = std::make_shared<Sphere>(Point3(10, 0, 0), 10.0);

    // 1. Union
    Union u(s1, s2);
    save_boolean("union", u);

    // 2. Intersection
    Intersection i(s1, s2);
    save_boolean("intersection", i);

    // 3. Difference
    Difference d(s1, s2);
    save_boolean("difference", d);

    std::cout << "Boolean tests completed." << std::endl;
    return 0;
}
