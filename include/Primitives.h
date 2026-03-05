#pragma once

#include "SDF.h"
#include <algorithm>

namespace Geom {

    class Sphere : public SDFNode<Sphere> {
    public:
        Point3 center;
        Scalar radius;

        Sphere(Point3 c, Scalar r) : center(c), radius(r) {}

        template <typename T>
        T evaluate(const Vec3T<T>& p) const {
            Vec3T<T> tc(static_cast<T>(center.x), static_cast<T>(center.y), static_cast<T>(center.z));
            return (p - tc).length() - static_cast<T>(radius);
        }

        BoundingBox boundingBox() const override {
            Vec3 r_vec(radius, radius, radius);
            return BoundingBox(center - r_vec, center + r_vec);
        }

        size_t numParams() const override { return 4; } // cx, cy, cz, r
        
        Scalar getParam(size_t i) const override {
            if(i==0) return center.x;
            if(i==1) return center.y;
            if(i==2) return center.z;
            if(i==3) return radius;
            return 0.0;
        }

        void setParam(size_t i, Scalar val) override {
            if(i==0) center.x = val;
            else if(i==1) center.y = val;
            else if(i==2) center.z = val;
            else if(i==3) radius = val;
        }

        DualScalar evalParamD(const Point3& p, size_t paramIndex) const override {
            // center x,y,z (0, 1, 2)
            // radius (3)
            DualScalar cx(center.x, paramIndex == 0 ? 1.0 : 0.0);
            DualScalar cy(center.y, paramIndex == 1 ? 1.0 : 0.0);
            DualScalar cz(center.z, paramIndex == 2 ? 1.0 : 0.0);
            DualScalar r(radius, paramIndex == 3 ? 1.0 : 0.0);

            Vec3T<DualScalar> tc(cx, cy, cz);
            Vec3T<DualScalar> tp(DualScalar::constant(p.x), DualScalar::constant(p.y), DualScalar::constant(p.z));

            return (tp - tc).length() - r;
        }
    };

    class Box : public SDFNode<Box> {
    public:
        Point3 center;
        Vec3 bounds; // Half-dimensions

        Box(Point3 c, Vec3 b) : center(c), bounds(b) {}

        template <typename T>
        T evaluate(const Vec3T<T>& p) const {
            T tc_x = static_cast<T>(center.x);
            T tc_y = static_cast<T>(center.y);
            T tc_z = static_cast<T>(center.z);
            
            T tb_x = static_cast<T>(bounds.x);
            T tb_y = static_cast<T>(bounds.y);
            T tb_z = static_cast<T>(bounds.z);
            
            T qx, qy, qz;
            if constexpr (std::is_same_v<T, Scalar>) {
                qx = std::abs(p.x - tc_x) - tb_x;
                qy = std::abs(p.y - tc_y) - tb_y;
                qz = std::abs(p.z - tc_z) - tb_z;
            } else {
                qx = abs(p.x - tc_x) - tb_x;
                qy = abs(p.y - tc_y) - tb_y;
                qz = abs(p.z - tc_z) - tb_z;
            }

            T u, v, w, inside_dist, outside_dist;
            if constexpr (std::is_same_v<T, Scalar>) {
                u = std::max(qx, static_cast<T>(0.0));
                v = std::max(qy, static_cast<T>(0.0));
                w = std::max(qz, static_cast<T>(0.0));
                inside_dist = std::min(std::max(qx, std::max(qy, qz)), static_cast<T>(0.0));
            } else {
                u = max(qx, static_cast<T>(0.0));
                v = max(qy, static_cast<T>(0.0));
                w = max(qz, static_cast<T>(0.0));
                inside_dist = min(max(qx, max(qy, qz)), static_cast<T>(0.0));
            }
            
            if constexpr (std::is_same_v<T, Scalar>) {
                outside_dist = std::sqrt(u*u + v*v + w*w);
            } else {
                outside_dist = sqrt(u*u + v*v + w*w);
            }
            
            return outside_dist + inside_dist;
        }

        BoundingBox boundingBox() const override {
            return BoundingBox(center - bounds, center + bounds);
        }

        size_t numParams() const override { return 6; } // cx, cy, cz, bx, by, bz
        
        Scalar getParam(size_t i) const override {
            if(i==0) return center.x;
            if(i==1) return center.y;
            if(i==2) return center.z;
            if(i==3) return bounds.x;
            if(i==4) return bounds.y;
            if(i==5) return bounds.z;
            return 0.0;
        }

        void setParam(size_t i, Scalar val) override {
            if(i==0) center.x = val;
            else if(i==1) center.y = val;
            else if(i==2) center.z = val;
            else if(i==3) bounds.x = val;
            else if(i==4) bounds.y = val;
            else if(i==5) bounds.z = val;
        }
    };

    class Cylinder : public SDFNode<Cylinder> {
    public:
        Point3 center; 
        Scalar radius;
        Scalar height; 

        Cylinder(Point3 c, Scalar r, Scalar h) : center(c), radius(r), height(h) {}

        template <typename T>
        T evaluate(const Vec3T<T>& p) const {
            T tc_x = static_cast<T>(center.x);
            T tc_y = static_cast<T>(center.y);
            T tc_z = static_cast<T>(center.z);
            
            T lx = p.x - tc_x;
            T ly = p.y - tc_y;
            T lz = p.z - tc_z;

            T x, y;
            if constexpr (std::is_same_v<T, Scalar>) {
                x = std::sqrt(lx*lx + ly*ly) - static_cast<T>(radius);
                y = std::abs(lz) - static_cast<T>(height * 0.5);
            } else {
                x = sqrt(lx*lx + ly*ly) - static_cast<T>(radius);
                y = abs(lz) - static_cast<T>(height * 0.5);
            }

            T u, v, w;
            if constexpr (std::is_same_v<T, Scalar>) {
                u = std::max(x, static_cast<T>(0.0));
                v = std::max(y, static_cast<T>(0.0));
                w = std::min(std::max(x, y), static_cast<T>(0.0));
            } else {
                u = max(x, static_cast<T>(0.0));
                v = max(y, static_cast<T>(0.0));
                w = min(max(x, y), static_cast<T>(0.0));
            }

            if constexpr (std::is_same_v<T, Scalar>) {
                return std::sqrt(u*u + v*v) + w;
            } else {
                return sqrt(u*u + v*v) + w;
            }
        }

        BoundingBox boundingBox() const override {
            Vec3 half_size(radius, radius, height * 0.5);
            return BoundingBox(center - half_size, center + half_size);
        }

        size_t numParams() const override { return 5; } // cx, cy, cz, r, h
        
        Scalar getParam(size_t i) const override {
            if(i==0) return center.x;
            if(i==1) return center.y;
            if(i==2) return center.z;
            if(i==3) return radius;
            if(i==4) return height;
            return 0.0;
        }

        void setParam(size_t i, Scalar val) override {
            if(i==0) center.x = val;
            else if(i==1) center.y = val;
            else if(i==2) center.z = val;
            else if(i==3) radius = val;
            else if(i==4) height = val;
        }
    };
    
    class Torus : public SDFNode<Torus> {
    public:
        Point3 center;
        Scalar R; // Major radius
        Scalar r; // Minor radius

        Torus(Point3 c, Scalar majorR, Scalar minorR) : center(c), R(majorR), r(minorR) {}

        template <typename T>
        T evaluate(const Vec3T<T>& p) const {
            T tc_x = static_cast<T>(center.x);
            T tc_y = static_cast<T>(center.y);
            T tc_z = static_cast<T>(center.z);
            
            T lx = p.x - tc_x;
            T ly = p.y - tc_y;
            T lz = p.z - tc_z;

            if constexpr (std::is_same_v<T, Scalar>) {
                T xz_dist = std::sqrt(lx*lx + lz*lz) - static_cast<T>(R);
                return std::sqrt(xz_dist*xz_dist + ly*ly) - static_cast<T>(r);
            } else {
                T xz_dist = sqrt(lx*lx + lz*lz) - static_cast<T>(R);
                return sqrt(xz_dist*xz_dist + ly*ly) - static_cast<T>(r);
            }
        }

        BoundingBox boundingBox() const override {
            Vec3 range(R + r, r, R + r);
            return BoundingBox(center - range, center + range);
        }

        size_t numParams() const override { return 5; }
        
        Scalar getParam(size_t i) const override {
            if(i==0) return center.x;
            if(i==1) return center.y;
            if(i==2) return center.z;
            if(i==3) return R;
            if(i==4) return r;
            return 0.0;
        }

        void setParam(size_t i, Scalar val) override {
            if(i==0) center.x = val;
            else if(i==1) center.y = val;
            else if(i==2) center.z = val;
            else if(i==3) R = val;
            else if(i==4) r = val;
        }
    };

    class Plane : public SDFNode<Plane> {
    public:
        Vec3 normal;
        Scalar d; // Distance from origin

        Plane(Vec3 n, Scalar dist) : normal(n.normalized()), d(dist) {}

        template <typename T>
        T evaluate(const Vec3T<T>& p) const {
            return p.x * static_cast<T>(normal.x) + 
                   p.y * static_cast<T>(normal.y) + 
                   p.z * static_cast<T>(normal.z) + static_cast<T>(d);
        }

        BoundingBox boundingBox() const override {
            return BoundingBox(); // Infinite bounding box
        }

        size_t numParams() const override { return 4; } // nx, ny, nz, d
        
        Scalar getParam(size_t i) const override {
            if(i==0) return normal.x;
            if(i==1) return normal.y;
            if(i==2) return normal.z;
            if(i==3) return d;
            return 0.0;
        }

        void setParam(size_t i, Scalar val) override {
            if(i==0) normal.x = val; // Note: Setting parts of a normal might invalidate normalization...
            else if(i==1) normal.y = val;
            else if(i==2) normal.z = val;
            else if(i==3) d = val;
        }
    };
}
