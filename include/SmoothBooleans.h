#pragma once

#include "SDF.h"
#include <algorithm>
#include <cmath>

namespace Geom {

    inline Scalar smin(Scalar a, Scalar b, Scalar k) {
        Scalar h = std::max(k - std::abs(a - b), 0.0) / k;
        return std::min(a, b) - h * h * k * 0.25;
    }

    inline Scalar smax(Scalar a, Scalar b, Scalar k) {
        Scalar h = std::max(k - std::abs(a - b), 0.0) / k;
        return std::max(a, b) + h * h * k * 0.25;
    }

    template <typename T>
    inline Dual<T> smin(Dual<T> a, Dual<T> b, Scalar k) {
        Dual<T> h = max(Dual<T>(static_cast<T>(k)) - abs(a - b), static_cast<T>(0.0)) / static_cast<T>(k);
        return min(a, b) - h * h * static_cast<T>(k * 0.25);
    }

    template <typename T>
    inline Dual<T> smax(Dual<T> a, Dual<T> b, Scalar k) {
        Dual<T> h = max(Dual<T>(static_cast<T>(k)) - abs(a - b), static_cast<T>(0.0)) / static_cast<T>(k);
        return max(a, b) + h * h * static_cast<T>(k * 0.25);
    }

    class SmoothUnion : public SDF {
        SDFPtr a, b;
        Scalar k;
    public:
        SmoothUnion(SDFPtr a, SDFPtr b, Scalar k) : a(a), b(b), k(k) {}
        
        Scalar eval(const Point3& p) const override {
            return smin(a->eval(p), b->eval(p), k);
        }

        DualScalar evalD(const Point3D& p) const override {
            return smin(a->evalD(p), b->evalD(p), k);
        }

        Dual2Scalar evalD2(const Point3D2& p) const override {
            return smin(a->evalD2(p), b->evalD2(p), k);
        }
        
        BoundingBox boundingBox() const override {
            BoundingBox box = a->boundingBox();
            box.expand(b->boundingBox());
            box.min -= Vec3(k*0.25, k*0.25, k*0.25);
            box.max += Vec3(k*0.25, k*0.25, k*0.25);
            return box;
        }
    };

    class SmoothIntersection : public SDF {
        SDFPtr a, b;
        Scalar k;
    public:
        SmoothIntersection(SDFPtr a, SDFPtr b, Scalar k) : a(a), b(b), k(k) {}
        
        Scalar eval(const Point3& p) const override {
            return smax(a->eval(p), b->eval(p), k);
        }

        DualScalar evalD(const Point3D& p) const override {
            return smax(a->evalD(p), b->evalD(p), k);
        }

        Dual2Scalar evalD2(const Point3D2& p) const override {
            return smax(a->evalD2(p), b->evalD2(p), k);
        }
        
        BoundingBox boundingBox() const override {
            BoundingBox box_a = a->boundingBox();
            BoundingBox box_b = b->boundingBox();
            Point3 min_p = max(box_a.min, box_b.min);
            Point3 max_p = min(box_a.max, box_b.max);
            
            if (min_p.x > max_p.x || min_p.y > max_p.y || min_p.z > max_p.z) {
                 return BoundingBox(); // Empty or invalid
            }
            return BoundingBox(min_p, max_p);
        }
    };

    class SmoothDifference : public SDF {
        SDFPtr a, b; // a - b
        Scalar k;
    public:
        SmoothDifference(SDFPtr a, SDFPtr b, Scalar k) : a(a), b(b), k(k) {}
        
        Scalar eval(const Point3& p) const override {
            return smax(a->eval(p), -b->eval(p), k);
        }

        DualScalar evalD(const Point3D& p) const override {
            return smax(a->evalD(p), b->evalD(p) * -1.0, k);
        }

        Dual2Scalar evalD2(const Point3D2& p) const override {
            return smax(a->evalD2(p), b->evalD2(p) * -1.0, k);
        }
        
        BoundingBox boundingBox() const override {
            return a->boundingBox();
        }
    };
}
