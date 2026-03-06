#include "Primitives.h"
#include "Booleans.h"
#include "ShapeKernel.h"
#include "MarchingCubes.h"
#include "STLExporter.h"
#include <iostream>
#include <memory>

using namespace Geom;

void save_smooth(const std::string& name, SDF& sdf, Scalar resolution = 0.5) {
    std::cout << "Meshing Smooth " << name << "..." << std::endl;
    Mesh mesh;
    BoundingBox bounds = sdf.boundingBox();
    bounds.expand(bounds.min - Vec3(1, 1, 1));
    bounds.expand(bounds.max + Vec3(1, 1, 1));
    
    MarchingCubes::march(sdf, bounds, resolution, mesh);
    std::string filename = "smooth_" + name + ".stl";
    STLExporter::save(mesh, filename);
    std::cout << "Saved " << filename << std::endl;
}

int main() {
    std::cout << "Running Smooth Boolean Calculation Tests..." << std::endl;

    auto s1 = std::make_shared<Sphere>(Point3(0, 0, 0), 10.0);
    auto s2 = std::make_shared<Sphere>(Point3(12, 0, 0), 10.0);
    Scalar k = 4.0;

    // 1. Smooth Union
    SmoothUnion su(s1, s2, k);
    save_smooth("union", su);

    // 2. Smooth Intersection
    SmoothIntersection si(s1, s2, k);
    save_smooth("intersection", si);

    // 3. Smooth Difference
    SmoothDifference sd(s1, s2, k);
    save_smooth("difference", sd);

    // 4. Smooth Junction (ShapeKernel)
    auto cyl = std::make_shared<Cylinder>(Point3(0, 0, 0), 4.0, 30.0);
    Geom::ShapeKernel::SmoothJunction junction(s1, cyl, 3.0);
    save_smooth("junction", junction);

    std::cout << "Smooth Boolean tests completed." << std::endl;
    return 0;
}
