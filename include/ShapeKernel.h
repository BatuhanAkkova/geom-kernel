#pragma once

#include "SDF.h"
#include "Primitives.h"
#include "Booleans.h"
#include "Modifiers.h"
#include <vector>

namespace Geom {
namespace ShapeKernel {

    /**
     * @brief A Pipe is a line-segment based capsule.
     */
    class Pipe : public SDFNode<Pipe> {
        Point3 start, end;
        Scalar radius;
    public:
        Pipe(Point3 s, Point3 e, Scalar r) : start(s), end(e), radius(r) {}

        template <typename T>
        T evaluate(const Vec3T<T>& p) const {
            using std::max; using std::min; using std::abs; using std::sqrt; using std::pow; using std::sin; using std::cos;
            
            Vec3T<T> t_start(static_cast<T>(start.x), static_cast<T>(start.y), static_cast<T>(start.z));
            Vec3T<T> t_end(static_cast<T>(end.x), static_cast<T>(end.y), static_cast<T>(end.z));
            
            Vec3T<T> pa = p - t_start;
            Vec3T<T> ba = t_end - t_start;
            
            T h = pa.dot(ba) / ba.lengthSquared();
            h = max(static_cast<T>(0.0), min(static_cast<T>(1.0), h));
            
            return (pa - ba * h).length() - static_cast<T>(radius);
        }

        BoundingBox boundingBox() const override {
            Vec3 r(radius, radius, radius);
            BoundingBox box(min(start, end) - r, max(start, end) + r);
            return box;
        }
    };

    /**
     * @brief A SmoothJunction provides a filleted union between two SDFs.
     */
    class SmoothJunction : public SDFNode<SmoothJunction> {
        SDFPtr a, b;
        Scalar filletRadius;
    public:
        SmoothJunction(SDFPtr a, SDFPtr b, Scalar r) : a(a), b(b), filletRadius(r) {}

        template <typename T>
        T evaluate(const Vec3T<T>& p) const {
            using std::max; using std::min; using std::abs; using std::sqrt; using std::pow; using std::sin; using std::cos;
            
            if constexpr (std::is_same_v<T, Scalar>) return Geom::smin(a->eval(p), b->eval(p), filletRadius);
            else if constexpr (std::is_same_v<T, DualScalar>) return Geom::smin(a->evalD(p), b->evalD(p), filletRadius);
            else if constexpr (std::is_same_v<T, Dual2Scalar>) return Geom::smin(a->evalD2(p), b->evalD2(p), filletRadius);
            return static_cast<T>(0);
        }

        BoundingBox boundingBox() const override {
            BoundingBox box = a->boundingBox();
            box.expand(b->boundingBox());
            Vec3 r(filletRadius * 0.25, filletRadius * 0.25, filletRadius * 0.25);
            box.min -= r;
            box.max += r;
            return box;
        }
    };

    /**
     * @brief A StructuralRib adds a reinforcing blade to a base geometry.
     */
    class StructuralRib : public SDFNode<StructuralRib> {
        SDFPtr base;
        Plane midPlane;
        Scalar thickness;
        Scalar height;
    public:
        StructuralRib(SDFPtr base, Plane p, Scalar t, Scalar h) 
            : base(base), midPlane(p), thickness(t), height(h) {}

        template <typename T>
        T evaluate(const Vec3T<T>& p) const {
            using std::max; using std::min; using std::abs; using std::sqrt; using std::pow; using std::sin; using std::cos;
            
            T distToPlane = abs(midPlane.evaluate<T>(p));
            T ribDist = distToPlane - static_cast<T>(thickness * 0.5);
            
            T baseDist;
            if constexpr (std::is_same_v<T, Scalar>) baseDist = base->eval(p);
            else if constexpr (std::is_same_v<T, DualScalar>) baseDist = base->evalD(p);
            else if constexpr (std::is_same_v<T, Dual2Scalar>) baseDist = base->evalD2(p);
            
            T solidRib = max(ribDist, baseDist - static_cast<T>(height));
            return min(baseDist, solidRib);
        }

        BoundingBox boundingBox() const override {
            BoundingBox box = base->boundingBox();
            Vec3 h(height, height, height);
            box.min -= h;
            box.max += h;
            return box;
        }
    };

    /**
     * @brief A CoolingChannel is a subtracted internal path.
     */
    class CoolingChannel : public SDFNode<CoolingChannel> {
        SDFPtr body;
        SDFPtr path;
    public:
        CoolingChannel(SDFPtr body, SDFPtr path) : body(body), path(path) {}

        template <typename T>
        T evaluate(const Vec3T<T>& p) const {
            using std::max; using std::min; using std::abs; using std::sqrt; using std::pow; using std::sin; using std::cos;
            
            if constexpr (std::is_same_v<T, Scalar>) return std::max(body->eval(p), -path->eval(p));
            else if constexpr (std::is_same_v<T, DualScalar>) return max(body->evalD(p), path->evalD(p) * -1.0);
            else if constexpr (std::is_same_v<T, Dual2Scalar>) return max(body->evalD2(p), path->evalD2(p) * -1.0);
            return static_cast<T>(0);
        }

        BoundingBox boundingBox() const override {
            return body->boundingBox();
        }
    };

    /**
     * @brief EngineeringShell creates a hollow wall of specified thickness.
     */
    class EngineeringShell : public SDFNode<EngineeringShell> {
        SDFPtr body;
        Scalar thickness;
        bool internal;
    public:
        EngineeringShell(SDFPtr b, Scalar t, bool i = true) : body(b), thickness(t), internal(i) {}

        template <typename T>
        T evaluate(const Vec3T<T>& p) const {
            using std::max; using std::min; using std::abs; using std::sqrt; using std::pow; using std::sin; using std::cos;
            
            T d;
            if constexpr (std::is_same_v<T, Scalar>) d = body->eval(p);
            else if constexpr (std::is_same_v<T, DualScalar>) d = body->evalD(p);
            else if constexpr (std::is_same_v<T, Dual2Scalar>) d = body->evalD2(p);
            
            T t = static_cast<T>(thickness);
            
            if (internal) return max(d, (d + t) * -1.0);
            else return max(d * -1.0, d - t);
        }

        BoundingBox boundingBox() const override {
            BoundingBox box = body->boundingBox();
            if (!internal) {
                Vec3 t(thickness, thickness, thickness);
                box.min -= t;
                box.max += t;
            }
            return box;
        }
    };

} // namespace ShapeKernel
} // namespace Geom
