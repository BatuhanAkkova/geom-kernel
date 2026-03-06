#pragma once

#include "SDF.h"
#include <algorithm>
#include <cmath>

namespace Geom {

    // --- Smooth Minimum/Maximum Evaluators ---
    // --- Smooth Minimum/Maximum Overloads ---
    
    inline Scalar smin(Scalar a, Scalar b, Scalar k) {
        using std::max; using std::min; using std::abs;
        Scalar h = max(k - abs(a - b), 0.0) / k;
        return min(a, b) - h * h * k * 0.25;
    }

    inline Scalar smax(Scalar a, Scalar b, Scalar k) {
        using std::max; using std::min; using std::abs;
        Scalar h = max(k - abs(a - b), 0.0) / k;
        return max(a, b) + h * h * k * 0.25;
    }

    inline DualScalar smin(DualScalar a, DualScalar b, Scalar k) {
        DualScalar h = max(static_cast<DualScalar>(k) - abs(a - b), static_cast<DualScalar>(0.0)) / static_cast<DualScalar>(k);
        return min(a, b) - h * h * static_cast<DualScalar>(k * 0.25);
    }

    inline DualScalar smax(DualScalar a, DualScalar b, Scalar k) {
        DualScalar h = max(static_cast<DualScalar>(k) - abs(a - b), static_cast<DualScalar>(0.0)) / static_cast<DualScalar>(k);
        return max(a, b) + h * h * static_cast<DualScalar>(k * 0.25);
    }

    inline Dual2Scalar smin(Dual2Scalar a, Dual2Scalar b, Scalar k) {
        Dual2Scalar h = max(static_cast<Dual2Scalar>(k) - abs(a - b), static_cast<Dual2Scalar>(0.0)) / static_cast<Dual2Scalar>(k);
        return min(a, b) - h * h * static_cast<Dual2Scalar>(k * 0.25);
    }

    inline Dual2Scalar smax(Dual2Scalar a, Dual2Scalar b, Scalar k) {
        Dual2Scalar h = max(static_cast<Dual2Scalar>(k) - abs(a - b), static_cast<Dual2Scalar>(0.0)) / static_cast<Dual2Scalar>(k);
        return max(a, b) + h * h * static_cast<Dual2Scalar>(k * 0.25);
    }

    // --- Booleans ---

    class Union : public SDFNode<Union> {
        SDFPtr a, b;
    public:
        Union(SDFPtr a, SDFPtr b) : a(a), b(b) {}
        
        template <typename T>
        T evaluate(const Vec3T<T>& p) const {
            using std::max; using std::min; using std::abs; using std::sqrt; using std::pow; using std::sin; using std::cos;
             
            if constexpr (std::is_same_v<T, Scalar>) {
                Point3 pt(p.x, p.y, p.z);
                return std::min(a->eval(pt), b->eval(pt));
            } else if constexpr (std::is_same_v<T, DualScalar>) {
                Point3D pt(p.x, p.y, p.z);
                return min(a->evalD(pt), b->evalD(pt));
            } else if constexpr (std::is_same_v<T, Dual2Scalar>) {
                Point3D2 pt(p.x, p.y, p.z);
                return min(a->evalD2(pt), b->evalD2(pt));
            }
            return static_cast<T>(0);
        }
        
        BoundingBox boundingBox() const override {
            BoundingBox box = a->boundingBox();
            box.expand(b->boundingBox());
            return box;
        }
    };

    class Intersection : public SDFNode<Intersection> {
        SDFPtr a, b;
    public:
        Intersection(SDFPtr a, SDFPtr b) : a(a), b(b) {}
        
        template <typename T>
        T evaluate(const Vec3T<T>& p) const {
            using std::max; using std::min; using std::abs; using std::sqrt; using std::pow; using std::sin; using std::cos;
            if constexpr (std::is_same_v<T, Scalar>) {
                Point3 pt(p.x, p.y, p.z);
                return std::max(a->eval(pt), b->eval(pt));
            } else if constexpr (std::is_same_v<T, DualScalar>) {
                Point3D pt(p.x, p.y, p.z);
                return max(a->evalD(pt), b->evalD(pt));
            } else if constexpr (std::is_same_v<T, Dual2Scalar>) {
                Point3D2 pt(p.x, p.y, p.z);
                return max(a->evalD2(pt), b->evalD2(pt));
            }
            return static_cast<T>(0);
        }
        
        BoundingBox boundingBox() const override {
            BoundingBox box_a = a->boundingBox();
            BoundingBox box_b = b->boundingBox();
            Point3 min_p = max(box_a.min, box_b.min);
            Point3 max_p = min(box_a.max, box_b.max);
            
            if (min_p.x > max_p.x || min_p.y > max_p.y || min_p.z > max_p.z) {
                 // return naive intersection box
            }
            return BoundingBox(min_p, max_p);
        }
    };

    class Difference : public SDFNode<Difference> {
        SDFPtr a, b; // a - b
    public:
        Difference(SDFPtr a, SDFPtr b) : a(a), b(b) {}
        
        template <typename T>
        T evaluate(const Vec3T<T>& p) const {
            using std::max; using std::min; using std::abs; using std::sqrt; using std::pow; using std::sin; using std::cos;
            if constexpr (std::is_same_v<T, Scalar>) {
                Point3 pt(p.x, p.y, p.z);
                return std::max(a->eval(pt), -b->eval(pt));
            } else if constexpr (std::is_same_v<T, DualScalar>) {
                Point3D pt(p.x, p.y, p.z);
                return max(a->evalD(pt), b->evalD(pt) * -1.0);
            } else if constexpr (std::is_same_v<T, Dual2Scalar>) {
                Point3D2 pt(p.x, p.y, p.z);
                return max(a->evalD2(pt), b->evalD2(pt) * -1.0);
            }
            return static_cast<T>(0);
        }
        
        BoundingBox boundingBox() const override {
            return a->boundingBox();
        }
    };

    // --- Smooth Booleans ---

    class SmoothUnion : public SDFNode<SmoothUnion> {
        SDFPtr a, b;
        Scalar k;
    public:
        SmoothUnion(SDFPtr a, SDFPtr b, Scalar k) : a(a), b(b), k(k) {}
        
        template <typename T>
        T evaluate(const Vec3T<T>& p) const {
            using std::max; using std::min; using std::abs; using std::sqrt; using std::pow; using std::sin; using std::cos;
            if constexpr (std::is_same_v<T, Scalar>) {
                Point3 pt(p.x, p.y, p.z);
                return smin(a->eval(pt), b->eval(pt), k);
            } else if constexpr (std::is_same_v<T, DualScalar>) {
                Point3D pt(p.x, p.y, p.z);
                return smin(a->evalD(pt), b->evalD(pt), k);
            } else if constexpr (std::is_same_v<T, Dual2Scalar>) {
                Point3D2 pt(p.x, p.y, p.z);
                return smin(a->evalD2(pt), b->evalD2(pt), k);
            }
            return static_cast<T>(0);
        }
        
        BoundingBox boundingBox() const override {
            BoundingBox box = a->boundingBox();
            box.expand(b->boundingBox());
            box.min -= Vec3(k*0.25, k*0.25, k*0.25);
            box.max += Vec3(k*0.25, k*0.25, k*0.25);
            return box;
        }
    };

    class SmoothIntersection : public SDFNode<SmoothIntersection> {
        SDFPtr a, b;
        Scalar k;
    public:
        SmoothIntersection(SDFPtr a, SDFPtr b, Scalar k) : a(a), b(b), k(k) {}
        
        template <typename T>
        T evaluate(const Vec3T<T>& p) const {
            using std::max; using std::min; using std::abs; using std::sqrt; using std::pow; using std::sin; using std::cos;
            if constexpr (std::is_same_v<T, Scalar>) {
                Point3 pt(p.x, p.y, p.z);
                return smax(a->eval(pt), b->eval(pt), k);
            } else if constexpr (std::is_same_v<T, DualScalar>) {
                Point3D pt(p.x, p.y, p.z);
                return smax(a->evalD(pt), b->evalD(pt), k);
            } else if constexpr (std::is_same_v<T, Dual2Scalar>) {
                Point3D2 pt(p.x, p.y, p.z);
                return smax(a->evalD2(pt), b->evalD2(pt), k);
            }
            return static_cast<T>(0);
        }
        
        BoundingBox boundingBox() const override {
            BoundingBox box_a = a->boundingBox();
            BoundingBox box_b = b->boundingBox();
            Point3 min_p = max(box_a.min, box_b.min);
            Point3 max_p = min(box_a.max, box_b.max);
            
            if (min_p.x > max_p.x || min_p.y > max_p.y || min_p.z > max_p.z) {
                 return BoundingBox(); 
            }
            return BoundingBox(min_p, max_p);
        }
    };

    class SmoothDifference : public SDFNode<SmoothDifference> {
        SDFPtr a, b; 
        Scalar k;
    public:
        SmoothDifference(SDFPtr a, SDFPtr b, Scalar k) : a(a), b(b), k(k) {}
        
        template <typename T>
        T evaluate(const Vec3T<T>& p) const {
            using std::max; using std::min; using std::abs; using std::sqrt; using std::pow; using std::sin; using std::cos;
            if constexpr (std::is_same_v<T, Scalar>) {
                Point3 pt(p.x, p.y, p.z);
                return smax(a->eval(pt), -b->eval(pt), k);
            } else if constexpr (std::is_same_v<T, DualScalar>) {
                Point3D pt(p.x, p.y, p.z);
                return smax(a->evalD(pt), b->evalD(pt) * -1.0, k);
            } else if constexpr (std::is_same_v<T, Dual2Scalar>) {
                Point3D2 pt(p.x, p.y, p.z);
                return smax(a->evalD2(pt), b->evalD2(pt) * -1.0, k);
            }
            return static_cast<T>(0);
        }
        
        BoundingBox boundingBox() const override {
            return a->boundingBox();
        }
    };
}
