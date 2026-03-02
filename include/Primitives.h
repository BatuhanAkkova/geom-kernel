#pragma once

#include "SDF.h"
#include <algorithm>

namespace Geom {

    class Sphere : public SDF {
    public:
        Point3 center;
        Scalar radius;

        Sphere(Point3 c, Scalar r) : center(c), radius(r) {}

        Scalar eval(const Point3& p) const override {
            return (p - center).length() - radius;
        }

        DualScalar evalD(const Point3D& p) const override {
            DualScalar dx = p.x - center.x;
            DualScalar dy = p.y - center.y;
            DualScalar dz = p.z - center.z;
            return sqrt(dx*dx + dy*dy + dz*dz) - radius;
        }

        BoundingBox boundingBox() const override {
            Vec3 r(radius, radius, radius);
            return BoundingBox(center - r, center + r);
        }
    };

    class Box : public SDF {
    public:
        Point3 center;
        Vec3 bounds; // Half-dimensions

        Box(Point3 c, Vec3 b) : center(c), bounds(b) {}

        Scalar eval(const Point3& p) const override {
            Vec3 q = Vec3(std::abs(p.x - center.x), std::abs(p.y - center.y), std::abs(p.z - center.z)) - bounds;
            
            // Length of max(q, 0.0) + min(max(q.x,max(q.y,q.z)),0.0)
            Vec3 q_positive = max(q, Vec3(0,0,0));
            Scalar outside_dist = q_positive.length();
            
            Scalar inside_dist = std::min(std::max(q.x, std::max(q.y, q.z)), 0.0);
            
            return outside_dist + inside_dist;
        }

        DualScalar evalD(const Point3D& p) const override {
            DualScalar qx = abs(p.x - center.x) - bounds.x;
            DualScalar qy = abs(p.y - center.y) - bounds.y;
            DualScalar qz = abs(p.z - center.z) - bounds.z;

            DualScalar u = max(qx, 0.0);
            DualScalar v = max(qy, 0.0);
            DualScalar w = max(qz, 0.0);
            DualScalar outside_dist = sqrt(u*u + v*v + w*w);

            DualScalar inside_dist = min(max(qx, max(qy, qz)), 0.0);
            return outside_dist + inside_dist;
        }

        BoundingBox boundingBox() const override {
            return BoundingBox(center - bounds, center + bounds);
        }
    };

    class Cylinder : public SDF {
    public:
        Point3 center; // Base center or mid-point? Let's assume infinite along Y for now, or simplify. 
                       // Standard SDF cylinder is usually infinite or capped. 
                       // Let's do a vertical cylinder (along Z) for MVP simplicity, centered at `center`.
        Scalar radius;
        Scalar height; // Total height

        Cylinder(Point3 c, Scalar r, Scalar h) : center(c), radius(r), height(h) {}

        Scalar eval(const Point3& p) const override {
            Point3 local = p - center;
            Scalar x = std::sqrt(local.x * local.x + local.y * local.y) - radius;
            Scalar y = std::abs(local.z) - height * 0.5;
            
            Scalar u = std::max(x, 0.0);
            Scalar v = std::max(y, 0.0);
            Scalar w = std::min(std::max(x, y), 0.0);
            
            return std::sqrt(u*u + v*v) + w;
        }

        DualScalar evalD(const Point3D& p) const override {
            DualScalar lx = p.x - center.x;
            DualScalar ly = p.y - center.y;
            DualScalar lz = p.z - center.z;

            DualScalar x = sqrt(lx*lx + ly*ly) - radius;
            DualScalar y = abs(lz) - height * 0.5;

            DualScalar u = max(x, 0.0);
            DualScalar v = max(y, 0.0);
            DualScalar w = min(max(x, y), 0.0);

            return sqrt(u*u + v*v) + w;
        }

        BoundingBox boundingBox() const override {
            Vec3 half_size(radius, radius, height * 0.5);
            return BoundingBox(center - half_size, center + half_size);
        }
    };
    
    class Plane : public SDF {
    public:
        Vec3 normal;
        Scalar d; // Distance from origin

        Plane(Vec3 n, Scalar dist) : normal(n.normalized()), d(dist) {}

        Scalar eval(const Point3& p) const override {
            return p.dot(normal) + d;
        }

        DualScalar evalD(const Point3D& p) const override {
            return p.x * normal.x + p.y * normal.y + p.z * normal.z + d;
        }

        BoundingBox boundingBox() const override {
            // Infinite bounding box
            return BoundingBox(); 
        }
    };

}
