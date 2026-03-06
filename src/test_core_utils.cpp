#include "Primitives.h"
#include "Transform.h"
#include "VoxelGrid.h"
#include "Octree.h"
#include "DualContouring.h"
#include "Autodiff.h"
#include "PointCloud.h"
#include "STLExporter.h"
#include <iostream>
#include <vector>
#include <cassert>

using namespace Geom;

void test_transforms() {
    std::cout << "Testing Transformations..." << std::endl;
    auto sphere = std::make_shared<Sphere>(Point3(0, 0, 0), 5.0);
    
    // Test Translation
    auto translated = translate(sphere, Vec3(10, 0, 0));
    assert(std::abs(translated->eval(Point3(10, 0, 0)) + 5.0) < 1e-6);
    
    // Test Scale (using uniformScale)
    auto scaled = uniformScale(sphere, 2.0);
    // Original radius 5, scaled by 2 -> radius 10.
    // At (10,0,0) distance should be 0.
    assert(std::abs(scaled->eval(Point3(10, 0, 0))) < 1e-6);

    std::cout << "Transform tests passed." << std::endl;
}

void test_voxelization() {
    std::cout << "Testing Voxelization..." << std::endl;
    Sphere sphere(Point3(0, 0, 0), 5.0);
    auto grid = VoxelGrid::sampleSDF(sphere, BoundingBox(Point3(-10, -10, -10), Point3(10, 10, 10)), 21, 21, 21);
    // Center point (index 10,10,10) is world (0,0,0). SDF at 0,0,0 is -5.
    assert(std::abs(grid->at(10, 10, 10) + 5.0f) < 0.1f);
    std::cout << "Voxelization tests passed." << std::endl;
}

void test_autodiff() {
    std::cout << "Testing Autodiff..." << std::endl;
    Sphere sphere(Point3(0, 0, 0), 10.0);
    Point3 p(10, 0, 0);
    
    Vec3 grad = sphere.gradient(p);
    assert(std::abs(grad.x - 1.0) < 1e-4);
    
    Field::enableHessian = true;
    Mat3 hess = sphere.hessian(p);
    assert(std::abs(hess(1,1) - 0.1) < 1e-4);
    std::cout << "Autodiff tests passed." << std::endl;
}

void test_dual_contouring() {
    std::cout << "Testing Dual Contouring..." << std::endl;
    Sphere sphere(Point3(0, 0, 0), 5.0);
    BoundingBox bounds(Point3(-6, -6, -6), Point3(6, 6, 6));
    Octree octree(sphere, bounds, 0.5);
    Mesh mesh;
    DualContouring::generateMesh(octree, sphere, mesh);
    STLExporter::save(mesh, "utils_dual_contouring.stl");
    std::cout << "Dual Contouring test passed and exported." << std::endl;
}

void test_point_cloud() {
    std::cout << "Testing Point Cloud SDF..." << std::endl;
    std::vector<Point3> pts;
    for(int i=0; i<360; i+=10) {
        Scalar r = i * M_PI / 180.0;
        pts.push_back(Point3(std::cos(r)*5.0, std::sin(r)*5.0, 0));
    }
    PointCloudSDF pc(pts);
    // Point on the circle should be approx 0
    assert(std::abs(pc.eval(Point3(5, 0, 0))) < 1.0);
    std::cout << "Point Cloud test passed." << std::endl;
}

int main() {
    try {
        test_transforms();
        test_voxelization();
        test_autodiff();
        test_dual_contouring();
        test_point_cloud();
        std::cout << "ALL CORE UTILS TESTS PASSED!" << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Test failed with exception: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
