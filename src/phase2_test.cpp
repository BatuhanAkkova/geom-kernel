#include "../include/SDF.h"
#include "../include/Primitives.h"
#include "../include/SmoothBooleans.h"
#include "../include/Modifiers.h"
#include <iostream>
#include <cassert>
#include <cmath>
#include <memory>

using namespace Geom;

void test_gradient() {
    std::cout << "Testing Gradient..." << std::endl;
    auto sphere = std::make_shared<Sphere>(Point3(0, 0, 0), 1.0);
    
    // Gradient on surface should be normal (normalized)
    // Point at (1, 0, 0) -> normal is (1, 0, 0)
    Point3 p(1.0, 0, 0);
    Vec3 g = sphere->gradient(p);
    
    // Finite difference might not be exact 1.0, but very close
    assert(std::abs(g.x - 1.0) < 0.001);
    assert(std::abs(g.y - 0.0) < 0.001);
    assert(std::abs(g.z - 0.0) < 0.001);

    // Inside checks
    Point3 p2(0.5, 0, 0);
    Vec3 g2 = sphere->gradient(p2); // Should point outwards (1,0,0)
    assert(std::abs(g2.x - 1.0) < 0.001);

    std::cout << "Gradient Test Passed" << std::endl;
}

void test_smooth_union() {
    std::cout << "Testing Smooth Union..." << std::endl;
    // Two spheres touching at origin
    auto s1 = std::make_shared<Sphere>(Point3(-0.9, 0, 0), 1.0);
    auto s2 = std::make_shared<Sphere>(Point3(0.9, 0, 0), 1.0);
    
    // Standard union eval at (0,0,0) should be > 0 (0.1 gap) if they were farther,
    // s1 dist at 0 is 0.9 - 1.0 = -0.1.
    // s2 dist at 0 is 0.9 - 1.0 = -0.1.
    // min(-0.1, -0.1) = -0.1.
    
    // Smooth union with k=0.5
    auto su = std::make_shared<SmoothUnion>(s1, s2, 0.5);
    Scalar d_smooth = su->eval(Point3(0, 0, 0));
    
    // smin(-0.1, -0.1) = -0.1 - h*h*k*0.25
    // h = max(k - 0, 0)/k = 1.
    // -0.1 - 1*1*0.5*0.25 = -0.1 - 0.125 = -0.225.
    // So it should be smaller (more negative)
    
    assert(d_smooth < -0.1);
    std::cout << "Smooth Union passed (Value: " << d_smooth << ")" << std::endl;
}

void test_offset() {
    std::cout << "Testing Offset..." << std::endl;
    auto box = std::make_shared<Box>(Point3(0,0,0), Vec3(1,1,1));
    // Box surface at x=1.
    // Offset by 0.1. Surface should be at x=1.1.
    auto expanded = std::make_shared<Offset>(box, 0.1);
    
    assert(std::abs(expanded->eval(Point3(1.1, 0, 0))) < 0.0001);
    // At 1.0, dist should be -0.1
    assert(std::abs(expanded->eval(Point3(1.0, 0, 0)) + 0.1) < 0.0001);

    std::cout << "Offset Test Passed" << std::endl;
}

void test_shell() {
    std::cout << "Testing Shell..." << std::endl;
    auto s = std::make_shared<Sphere>(Point3(0,0,0), 2.0);
    // Shell thickness 0.1.
    // Surface starts at dist 0 (radius 2).
    // Shell makes surface at |d| - 0.1 = 0 => |d| = 0.1.
    // So surfaces at r=1.9 and r=2.1.
    
    auto shell = std::make_shared<Shell>(s, 0.1);
    
    // At r=2.0 (original surface), SDF is 0.
    // Shell SDF = |0| - 0.1 = -0.1 (Inside wall center)
    assert(shell->eval(Point3(2,0,0)) < 0);
    assert(std::abs(shell->eval(Point3(2,0,0)) + 0.1) < 0.0001);
    
    // At r=2.1, original dist = 0.1. Shell dist = |0.1| - 0.1 = 0.
    assert(std::abs(shell->eval(Point3(2.1,0,0))) < 0.0001);
    
    // At r=1.9, original dist = -0.1. Shell dist = |-0.1| - 0.1 = 0.
    assert(std::abs(shell->eval(Point3(1.9,0,0))) < 0.0001);
    
    std::cout << "Shell Test Passed" << std::endl;
}

int main() {
    test_gradient();
    test_smooth_union();
    test_offset();
    test_shell();
    
    std::cout << "All Phase 2 tests passed!" << std::endl;
    return 0;
}
