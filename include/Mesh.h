#pragma once

#include "Geometry.h"
#include <vector>

namespace Geom {

    struct Triangle {
        Point3 v0, v1, v2;
        Vec3 normal;

        Triangle(Point3 a, Point3 b, Point3 c) : v0(a), v1(b), v2(c) {
            Vec3 edge1 = v1 - v0;
            Vec3 edge2 = v2 - v0;
            // Cross product
            normal = Vec3(
                edge1.y * edge2.z - edge1.z * edge2.y,
                edge1.z * edge2.x - edge1.x * edge2.z,
                edge1.x * edge2.y - edge1.y * edge2.x
            ).normalized();
        }
    };

    struct Mesh {
        std::vector<Triangle> triangles;

        void addTriangle(const Point3& a, const Point3& b, const Point3& c) {
            triangles.emplace_back(a, b, c);
        }

        void clear() {
            triangles.clear();
        }
    };

}
