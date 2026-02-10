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

    class SmoothUnion : public SDF {
        SDFPtr a, b;
        Scalar k;
    public:
        SmoothUnion(SDFPtr a, SDFPtr b, Scalar k) : a(a), b(b), k(k) {}
        
        Scalar eval(const Point3& p) const override {
            return smin(a->eval(p), b->eval(p), k);
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
        
        BoundingBox boundingBox() const override {
            return a->boundingBox();
        }
    };
}
