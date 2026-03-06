#include "TPMS.h"
#include "MarchingCubes.h"
#include "STLExporter.h"
#include <iostream>
#include <memory>

using namespace Geom;

void save_tpms(const std::string& name, SDF& sdf, const BoundingBox& bounds, Scalar resolution = 0.5) {
    std::cout << "Meshing " << name << "..." << std::endl;
    Mesh mesh;
    MarchingCubes::march(sdf, bounds, resolution, mesh);
    std::string filename = "tpms_" + name + ".stl";
    STLExporter::save(mesh, filename);
    std::cout << "Saved " << filename << std::endl;
}

int main() {
    std::cout << "Running TPMS Calculation Tests..." << std::endl;

    BoundingBox bounds(Point3(-10, -10, -10), Point3(10, 10, 10));
    Scalar period = 10.0;

    // 1. Gyroid
    Gyroid gyroid(period);
    save_tpms("gyroid", gyroid, bounds);

    // 2. Diamond
    DiamondTPMS diamond(period);
    save_tpms("diamond", diamond, bounds);

    // 3. Schwarz P
    SchwarzP schwarz(period);
    save_tpms("schwarzp", schwarz, bounds);

    std::cout << "TPMS tests completed." << std::endl;
    return 0;
}
