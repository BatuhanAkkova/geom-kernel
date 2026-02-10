#include "Octree.h"
#include "DualContouring.h"
#include "Primitives.h"
#include "STLExporter.h"
#include <iostream>
#include <fstream>
#include <chrono>

using namespace Geom;

// Helper SDF for Union
class UnionSDF : public SDF {
    Sphere sa, sb;
    Vec3 offset;
public:
    UnionSDF(double r1, double r2, Vec3 off) : sa(Point3(0,0,0), r1), sb(Point3(0,0,0), r2), offset(off) {}
    double eval(const Point3& p) const override {
        double d1 = sa.eval(p);
        double d2 = sb.eval(p - offset);
        return std::min(d1, d2);
    }
    
    BoundingBox boundingBox() const override {
        BoundingBox b1 = sa.boundingBox();
        BoundingBox b2 = sb.boundingBox();
        // Shift b2
        b2.min = b2.min + offset;
        b2.max = b2.max + offset;
        
        return BoundingBox(
            Point3(std::min(b1.min.x, b2.min.x), std::min(b1.min.y, b2.min.y), std::min(b1.min.z, b2.min.z)),
            Point3(std::max(b1.max.x, b2.max.x), std::max(b1.max.y, b2.max.y), std::max(b1.max.z, b2.max.z))
        );
    }
};



int main() {
    std::cout << "Starting Phase 4: Adaptive & Sparse Sampling Test" << std::endl;

    // 1. Setup SDF (Sphere)
    Sphere sphere(Point3(0,0,0), 5.0);
    BoundingBox bounds(Point3(-6, -6, -6), Point3(6, 6, 6));

    // 2. Build Octree
    // Resolution approx 0.1
    // Bounds size 12. 
    // 12 / 0.1 = 120. Depth ~7 (128).
    // Timer start
    auto start = std::chrono::high_resolution_clock::now();
    Octree octree(sphere, bounds, 0.1); 
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;
    std::cout << "Total Octree Nodes: " << octree.nodes.size() << std::endl;

    // 3. Generate Mesh
    Mesh mesh;
    start = std::chrono::high_resolution_clock::now();
    DualContouring::generateMesh(octree, sphere, mesh);
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    std::cout << "Meshing Time: " << elapsed.count() << " ms" << std::endl;
              
    std::cout << "Vertices: " << mesh.triangles.size() * 3 << " (approx, unindexed)" << std::endl;
    std::cout << "Triangles: " << mesh.triangles.size() << std::endl;

    // 4. Export
    STLExporter::save(mesh, "octree_sphere.stl");
    std::cout << "Exported to octree_sphere.stl" << std::endl;
    
    // Test 2: Two Spheres (Union)
    std::cout << "\nTest 2: Union of two spheres" << std::endl;
    Sphere s1(Point3(0,0,0), 4.0);
    Sphere s2(Point3(0,0,0), 3.0);
    // Offset s2
    // We don't have Transform modifier locally here easily without full hierarchy,
    // but we can make a custom SDF lambda or class if needed.
    // Let's just use a simple approach: logic in a lambda wrapper.

    
    UnionSDF unionSdf(4.0, 3.0, Vec3(5, 0, 0));
    BoundingBox bounds2(Point3(-5, -5, -5), Point3(10, 5, 5));
    
    Octree octree2(unionSdf, bounds2, 0.15);
    Mesh mesh2;
    DualContouring::generateMesh(octree2, unionSdf, mesh2);
    
    STLExporter::save(mesh2, "octree_union.stl");
    std::cout << "Exported octree_union.stl (" << mesh2.triangles.size() << " triangles)" << std::endl;

    return 0;
}
