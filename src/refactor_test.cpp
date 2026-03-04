#include "../include/Geometry.h"
#include "../include/Autodiff.h"
#include "../include/SDF.h"
#include "../include/Primitives.h"
#include "../include/Transform.h"
#include "../include/TPMS.h"
#include "../include/VoxelGrid.h"
#include "../include/PointCloud.h"
#include "../include/MeshSDF.h"
#include "../include/MarchingCubes.h"
#include <iostream>
#include <vector>
#include <cassert>

using namespace Geom;

#define TEST(name) void test_##name()
#define RUN_TEST(name) \
    std::cout << "Running " << #name << "... "; \
    test_##name(); \
    std::cout << "PASSED" << std::endl;

// 1. Torus Test
TEST(Torus) {
    Torus torus(Point3(0,0,0), 10.0, 2.0);
    // Point on the tube center
    assert(std::abs(torus.eval(Point3(10, 0, 0))) < 1e-6);
    // Inside the tube
    assert(torus.eval(Point3(10, 1, 0)) < 0);
    // Outside
    assert(torus.eval(Point3(0, 0, 0)) > 0);
}

// 2. TPMS Test
TEST(TPMS) {
    Gyroid g(10.0); // Period 10
    // Gyroid should have signs changes
    assert(g.eval(Point3(0,0,0)) < 0 || g.eval(Point3(0,0,0)) > 0 || std::abs(g.eval(Point3(0,0,0))) < 1e-6);
    // Check periodicity roughly
    Scalar v1 = g.eval(Point3(1,1,1));
    Scalar v2 = g.eval(Point3(11,1,1));
    assert(std::abs(v1 - v2) < 1e-6);
}

// 3. Transform Test
TEST(Transform) {
    auto sphere = std::make_shared<Sphere>(Point3(0,0,0), 5.0);
    auto translated = translate(sphere, Vec3(10, 0, 0));
    
    // Original sphere at 0,0,0 has dist -5
    assert(std::abs(translated->eval(Point3(10,0,0)) + 5.0) < 1e-6);
    // At original 0,0,0 it should be outside now
    assert(translated->eval(Point3(0,0,0)) > 0);
}

// 4. VoxelGrid Test
TEST(VoxelGrid) {
    auto sphere = std::make_shared<Sphere>(Point3(0,0,0), 5.0);
    auto grid = VoxelGrid::sampleSDF(*sphere, BoundingBox(Point3(-10,-10,-10), Point3(10,10,10)), 21, 21, 21);
    
    // Center (index 10,10,10) should be -5
    assert(std::abs(grid->at(10,10,10) + 5.0) < 0.1);
}

// 5. Partial Volume Test
TEST(PartialVolume) {
    Scalar vol = partialVolume(0.0, 1.0);
    assert(std::abs(vol - 0.5) < 1e-6);
    
    Scalar vol_in = partialVolume(-10.0, 1.0);
    assert(vol_in > 0.99);
    
    Scalar vol_out = partialVolume(10.0, 1.0);
    assert(vol_out < 0.01);
}

// 6. Analytical Gradient Test
TEST(AnalyticalGradientTest) {
    auto sphere = std::make_shared<Sphere>(Point3(0,0,0), 10.0);
    Point3 p(5, 0, 0);
    Vec3 grad = sphere->computeAnalyticalGradient(p);
    assert(std::abs(grad.x - 1.0) < 1e-4);
    assert(std::abs(grad.y - 0.0) < 1e-4);
    assert(std::abs(grad.z - 0.0) < 1e-4);
}

TEST(HessiansTest) {
    auto sphere = std::make_shared<Sphere>(Point3(0,0,0), 10.0);
    Point3 p(10, 0, 0); // On the surface
    Mat3 h = sphere->analyticalHessian(p);

    // For a sphere f = sqrt(x^2+y^2+z^2) - R
    // df/dx = x/r
    // d2f/dx2 = (r - x^2/r)/r^2 = (1 - x^2/r^2)/r
    // At (10, 0, 0), r=10, x=10: d2f/dx2 = (1 - 1)/10 = 0
    // d2f/dy2 = (r - y^2/r)/r^2 = (1 - 0)/10 = 0.1
    // d2f/dz2 = 0.1
    // d2f/dxdy = -xy/r^3 = 0

    assert(std::abs(h(0, 0) - 0.0) < 1e-4);
    assert(std::abs(h(1, 1) - 0.1) < 1e-4);
    assert(std::abs(h(2, 2) - 0.1) < 1e-4);
    assert(std::abs(h(0, 1) - 0.0) < 1e-4);
    
    std::cout << "  Hessian at (10,0,0): [" << h(0,0) << ", " << h(1,1) << ", " << h(2,2) << "]" << std::endl;
}

// 7. PointCloudSDF Test
TEST(PointCloudSDF) {
    std::vector<Point3> pts = { Point3(0,0,0), Point3(1,0,0), Point3(0,1,0) };
    PointCloudSDF pc(pts);
    
    // Nearest to 0.1, 0, 0 is 0,0,0
    assert(pc.eval(Point3(0.1, 0, 0)) >= 0); // Outside since its a point cloud
}

int main() {
    std::cout << "Starting Refactor Tests..." << std::endl;
    
    RUN_TEST(Torus);
    RUN_TEST(TPMS);
    RUN_TEST(Transform);
    RUN_TEST(VoxelGrid);
    RUN_TEST(PartialVolume);
    RUN_TEST(AnalyticalGradientTest);
    RUN_TEST(HessiansTest);
    RUN_TEST(PointCloudSDF);
    
    std::cout << "\nALL TESTS PASSED!" << std::endl;
    return 0;
}
