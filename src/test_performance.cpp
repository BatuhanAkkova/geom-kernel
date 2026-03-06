#include "Primitives.h"
#include "Field.h"
#include <iostream>
#include <vector>
#include <chrono>
#include <iomanip>
#include <memory>

using namespace Geom;

class BenchmarkScene : public SDFNode<BenchmarkScene> {
    std::vector<std::shared_ptr<Sphere>> spheres;
public:
    BenchmarkScene(int count) {
        for(int i=0; i<count; ++i) {
            spheres.push_back(std::make_shared<Sphere>(Point3(i*2.0, 0, 0), 1.0));
        }
    }

    template <typename T>
    T evaluate(const Vec3T<T>& p) const {
        using std::min;
        T d = static_cast<T>(1e9);
        for(const auto& s : spheres) {
            T sd = s->template evaluate<T>(p);
            if(sd < d) d = sd;
        }
        return d;
    }

    BoundingBox boundingBox() const override { return BoundingBox(); }
};

int main() {
    std::cout << "Running Performance Tests..." << std::endl;

    int sphereCount = 50;
    int evalCount = 1000; // Reduced for faster CI/tests
    BenchmarkScene scene(sphereCount);
    Point3 p(25, 0, 0);

    // 1. Eval
    auto start = std::chrono::high_resolution_clock::now();
    volatile Scalar sum = 0;
    for(int i=0; i<evalCount; ++i) {
        sum = sum + scene.eval(p); // explicit addition to avoid volatile warnings
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Eval x" << evalCount << ": " << std::chrono::duration<double, std::milli>(end-start).count() << " ms" << std::endl;

    // 2. Gradient
    start = std::chrono::high_resolution_clock::now();
    for(int i=0; i<evalCount; ++i) {
        Vec3 g = scene.gradient(p);
        sum = sum + g.x;
    }
    end = std::chrono::high_resolution_clock::now();
    std::cout << "Gradient x" << evalCount << ": " << std::chrono::duration<double, std::milli>(end-start).count() << " ms" << std::endl;

    // 3. Hessian
    Field::enableHessian = true;
    start = std::chrono::high_resolution_clock::now();
    for(int i=0; i<evalCount; ++i) {
        Mat3 h = scene.hessian(p);
        sum = sum + h(0,0);
    }
    end = std::chrono::high_resolution_clock::now();
    std::cout << "Hessian x" << evalCount << ": " << std::chrono::duration<double, std::milli>(end-start).count() << " ms" << std::endl;

    std::cout << "Performance tests completed. Final sum (discarded): " << sum << std::endl;
    return 0;
}
