#include <iostream>
#include "Geometry.h"
#include "SDF.h"
#include "BoundingBox.h"

// Mock SDF for testing
class Sphere : public Geom::SDF {
public:
    Geom::Point3 center;
    Geom::Scalar radius;

    Sphere(Geom::Point3 c, Geom::Scalar r) : center(c), radius(r) {}

    Geom::Scalar eval(const Geom::Point3& p) const override {
        return (p - center).length() - radius;
    }

    Geom::BoundingBox boundingBox() const override {
        Geom::Vec3 r(radius, radius, radius);
        return Geom::BoundingBox(center - r, center + r);
    }
};

int main() {
    Geom::Point3 p(10, 20, 30);
    std::cout << "Point: " << p.x << ", " << p.y << ", " << p.z << std::endl;

    Sphere s(Geom::Point3(0,0,0), 10.0);
    Geom::Scalar dist = s.eval(p);
    
    std::cout << "Distance to sphere: " << dist << std::endl;

    if (dist > 0) {
        std::cout << "Outside" << std::endl;
    } else {
        std::cout << "Inside" << std::endl;
    }

    return 0;
}
