#pragma once

#include "SDF.h"
#include <algorithm>

namespace Geom {

    class Union : public SDF {
        SDFPtr a, b;
    public:
        Union(SDFPtr a, SDFPtr b) : a(a), b(b) {}
        
        Scalar eval(const Point3& p) const override {
            return std::min(a->eval(p), b->eval(p));
        }

        DualScalar evalD(const Point3D& p) const override {
            return min(a->evalD(p), b->evalD(p));
        }

        Dual2Scalar evalD2(const Point3D2& p) const override {
            return min(a->evalD2(p), b->evalD2(p));
        }
        
        BoundingBox boundingBox() const override {
            BoundingBox box = a->boundingBox();
            box.expand(b->boundingBox());
            return box;
        }
    };

    class Intersection : public SDF {
        SDFPtr a, b;
    public:
        Intersection(SDFPtr a, SDFPtr b) : a(a), b(b) {}
        
        Scalar eval(const Point3& p) const override {
            return std::max(a->eval(p), b->eval(p));
        }

        DualScalar evalD(const Point3D& p) const override {
            return max(a->evalD(p), b->evalD(p));
        }

        Dual2Scalar evalD2(const Point3D2& p) const override {
            return max(a->evalD2(p), b->evalD2(p));
        }
        
        BoundingBox boundingBox() const override {
            BoundingBox box_a = a->boundingBox();
            BoundingBox box_b = b->boundingBox();
            // Intersection box is the intersection of bounding boxes
            // Simplified: Min of maxes, Max of mins
            Point3 min_p = max(box_a.min, box_b.min);
            Point3 max_p = min(box_a.max, box_b.max);
            
            // Check if valid
            if (min_p.x > max_p.x || min_p.y > max_p.y || min_p.z > max_p.z) {
                 // Return empty box if no intersection roughly
                 // For now, return a naive intersection box
            }
            return BoundingBox(min_p, max_p);
        }
    };

    class Difference : public SDF {
        SDFPtr a, b; // a - b
    public:
        Difference(SDFPtr a, SDFPtr b) : a(a), b(b) {}
        
        Scalar eval(const Point3& p) const override {
            return std::max(a->eval(p), -b->eval(p));
        }

        DualScalar evalD(const Point3D& p) const override {
            return max(a->evalD(p), b->evalD(p) * -1.0);
        }

        Dual2Scalar evalD2(const Point3D2& p) const override {
            return max(a->evalD2(p), b->evalD2(p) * -1.0);
        }
        
        BoundingBox boundingBox() const override {
            // Subtracting doesn't expand the bounding box of A
            // It might shrink it, but calculating that is hard.
            // Safe upper bound is just A's box.
            return a->boundingBox();
        }
    };
}
