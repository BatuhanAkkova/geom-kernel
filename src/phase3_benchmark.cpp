#include "../include/SDF.h"
#include "../include/Primitives.h"
#include "../include/Booleans.h"
#include "../include/MarchingCubes.h"
#include <iostream>
#include <vector>
#include <chrono>
#include <memory>
#include <iomanip>

using namespace Geom;

// specific complex scene for benchmarking
class ComplexScene : public SDFNode<ComplexScene> {
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
                    bounds.expand(s->boundingBox());
                }
            }
        }
    }

    template <typename T>
    T evaluate(const Vec3T<T>& p) const {
        T min_dist = static_cast<T>(1e9);
        for(const auto& s : spheres) {
            if constexpr (std::is_same_v<T, Scalar>) {
                Point3 pt(p.x, p.y, p.z);
                min_dist = std::min(min_dist, s->eval(pt));
            } else if constexpr (std::is_same_v<T, DualScalar>) {
                Point3D pt(p.x, p.y, p.z);
                min_dist = min(min_dist, s->evalD(pt));
            } else if constexpr (std::is_same_v<T, Dual2Scalar>) {
                Point3D2 pt(p.x, p.y, p.z);
                min_dist = min(min_dist, s->evalD2(pt));
            }
        }
        return min_dist;
    }

    BoundingBox boundingBox() const override { return bounds; }
    size_t count() const { return spheres.size(); }
};

int main() {
    int gridSize = 3; 
    double spacing = 2.5;
    double radius = 1.0;
    
    std::cout << "Building Scene with " << (gridSize*gridSize*gridSize) << " spheres..." << std::endl;
    ComplexScene scene(gridSize, spacing, radius);
    BoundingBox bounds = scene.boundingBox();
    bounds.min = bounds.min - Vec3(1,1,1);
    bounds.max = bounds.max + Vec3(1,1,1);
    
    std::cout << "Bounds: " << bounds.min.x << "," << bounds.min.y << "," << bounds.min.z << " to " 
              << bounds.max.x << "," << bounds.max.y << "," << bounds.max.z << "\n\n";

    // Generate test points
    double res = 0.5;
    std::vector<Point3> points;
    for(double z = bounds.min.z; z <= bounds.max.z; z += res) {
        for(double y = bounds.min.y; y <= bounds.max.y; y += res) {
            for(double x = bounds.min.x; x <= bounds.max.x; x += res) {
                points.push_back(Point3(x, y, z));
            }
        }
    }
    
    std::cout << "Number of Evaluation Points: " << points.size() << "\n\n";
    
    // Disable Hessians initially
    Field::enableHessian = false;

    // 1. Eval Only
    auto start = std::chrono::high_resolution_clock::now();
    double sum = 0;
    for (const auto& p : points) {
        sum += scene.eval(p);
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> ms_eval = (end - start);
    std::cout << "1. Eval Only:          " << std::setw(8) << (ms_eval.count() * 1000.0) << " ms \t(Result: " << sum << ")\n";

    // 2. Eval + Gradient
    start = std::chrono::high_resolution_clock::now();
    double sum2 = 0;
    for (const auto& p : points) {
        sum2 += scene.gradient(p).x;
    }
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> ms_grad = (end - start);
    std::cout << "2. Eval + Gradient:    " << std::setw(8) << (ms_grad.count() * 1000.0) << " ms \t(Result: " << sum2 << ")\n";

    // 3. Eval + Hessian
    Field::enableHessian = true; // MUST enable for hessian
    start = std::chrono::high_resolution_clock::now();
    double sum3 = 0;
    for (const auto& p : points) {
        sum3 += scene.hessian(p)(0,0);
    }
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> ms_hess = (end - start);
    std::cout << "3. Eval + Hessian:     " << std::setw(8) << (ms_hess.count() * 1000.0) << " ms \t(Result: " << sum3 << ")\n";

    return 0;
}
