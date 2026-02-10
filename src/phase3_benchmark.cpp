#include "../include/SDF.h"
#include "../include/Primitives.h"
#include "../include/Booleans.h"
#include "../include/MarchingCubes.h"
#include <iostream>
#include <vector>
#include <chrono>
#include <memory>

using namespace Geom;

// specific complex scene for benchmarking
class ComplexScene : public SDF {
    std::vector<SDFPtr> spheres;
    BoundingBox bounds;
public:
    ComplexScene(int gridSize, double spacing, double radius) {
        bounds = BoundingBox();
        for(int z=0; z<gridSize; ++z) {
            for(int y=0; y<gridSize; ++y) {
                for(int x=0; x<gridSize; ++x) {
                    Point3 center(x * spacing, y * spacing, z * spacing);
                    auto s = std::make_shared<Sphere>(center, radius);
                    spheres.push_back(s);
                    
                    // Simple bounding box expansion
                    bounds.expand(s->boundingBox());
                }
            }
        }
    }

    Scalar eval(const Point3& p) const override {
        // Naive evaluation: Union of all spheres.
        // This is intentionally slow O(N) to simulate complex SDF graphs
        // or we can use a more optimized structure if we want to test just the meshing overhead.
        // But for "Scaling", we want to see how it handles many primitives.
        // However, standard Union is min().
        
        Scalar min_dist = 1e9;
        for(const auto& s : spheres) {
            min_dist = std::min(min_dist, s->eval(p));
        }
        return min_dist;
    }

    BoundingBox boundingBox() const override {
        return bounds;
    }
    
    // Helper to get raw sphere count
    size_t count() const { return spheres.size(); }
};

int main() {
    // 1. Setup Scene
    // 5x5x5 grid = 125 spheres.
    // If we do naive eval, valid point checks are 125 ops.
    // Grid resolution will determine total eval calls.
    int gridSize = 5; 
    double spacing = 2.5;
    double radius = 1.0;
    
    std::cout << "Building Scene with " << (gridSize*gridSize*gridSize) << " spheres..." << std::endl;
    ComplexScene scene(gridSize, spacing, radius);
    BoundingBox bounds = scene.boundingBox();
    
    // Add some padding to bounds
    bounds.min = bounds.min - Vec3(1,1,1);
    bounds.max = bounds.max + Vec3(1,1,1);
    
    std::cout << "Bounds: " << bounds.min.x << "," << bounds.min.y << "," << bounds.min.z << " to " 
              << bounds.max.x << "," << bounds.max.y << "," << bounds.max.z << std::endl;

    double resolutions[] = { 0.2, 0.1 }; // Coarse and Fine
    
    for(double res : resolutions) {
        std::cout << "\n-----------------------------------" << std::endl;
        std::cout << "Benchmarking Resolution: " << res << std::endl;
        
        Mesh mesh;
        auto start = std::chrono::high_resolution_clock::now();
        
        MarchingCubes::march(scene, bounds, res, mesh);
        
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        
        std::cout << "Time: " << elapsed.count() << " s" << std::endl;
        std::cout << "Triangles: " << mesh.triangles.size() << std::endl;
    }

    return 0;
}
