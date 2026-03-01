#pragma once

#include "SDF.h"
#include "Primitives.h"
#include "Booleans.h"
#include "SmoothBooleans.h"
#include "Modifiers.h"
#include <vector>

namespace Geom {
namespace ShapeKernel {

    /**
     * @brief A Pipe is a line-segment based capsule.
     * Useful for routing, framework, etc.
     */
    class Pipe : public SDF {
        Point3 start, end;
        Scalar radius;
    public:
        Pipe(Point3 s, Point3 e, Scalar r) : start(s), end(e), radius(r) {}

        Scalar eval(const Point3& p) const override {
            Vec3 pa = p - start;
            Vec3 ba = end - start;
            Scalar h = std::clamp(pa.dot(ba) / ba.lengthSquared(), 0.0, 1.0);
            return (pa - ba * h).length() - radius;
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
    class SmoothJunction : public SDF {
        SDFPtr a, b;
        Scalar filletRadius;
    public:
        SmoothJunction(SDFPtr a, SDFPtr b, Scalar r) : a(a), b(b), filletRadius(r) {}

        Scalar eval(const Point3& p) const override {
            return smin(a->eval(p), b->eval(p), filletRadius);
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
     * Defined by a mid-plane, thickness, and height limit.
     */
    class StructuralRib : public SDF {
        SDFPtr base;
        Plane midPlane;
        Scalar thickness;
        Scalar height;
    public:
        StructuralRib(SDFPtr base, Plane p, Scalar t, Scalar h) 
            : base(base), midPlane(p), thickness(t), height(h) {}

        Scalar eval(const Point3& p) const override {
            // Distance to plane (signed)
            Scalar distToPlane = std::abs(midPlane.eval(p));
            // Rib thickness constraint
            Scalar ribDist = distToPlane - thickness * 0.5;
            
            // Limit rib height relative to base surface? 
            // Simplified: Rib is an infinite slab intersected with an offset of the base
            // Or rib is only where base-SDF is within 'height' distance.
            Scalar baseDist = base->eval(p);
            
            // Rib exists where baseDist < height and |distToPlane| < thickness/2
            // We want it to be a smooth addition.
            Scalar heightConstraint = baseDist - height;
            
            // Intersection of thickness slab and height constraint
            Scalar solidRib = std::max(ribDist, heightConstraint);
            
            // Union with base
            return std::min(baseDist, solidRib);
        }

        BoundingBox boundingBox() const override {
            // Conservative: Box of base expanded by height
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
    class CoolingChannel : public SDF {
        SDFPtr body;
        SDFPtr path;
    public:
        CoolingChannel(SDFPtr body, SDFPtr path) : body(body), path(path) {}

        Scalar eval(const Point3& p) const override {
            // body - path
            return std::max(body->eval(p), -path->eval(p));
        }

        BoundingBox boundingBox() const override {
            return body->boundingBox();
        }
    };

    /**
     * @brief EngineeringShell creates a hollow wall of specified thickness.
     */
    class EngineeringShell : public SDF {
        SDFPtr body;
        Scalar thickness;
        bool internal;
    public:
        EngineeringShell(SDFPtr b, Scalar t, bool i = true) : body(b), thickness(t), internal(i) {}

        Scalar eval(const Point3& p) const override {
            Scalar d = body->eval(p);
            if (internal) {
                // Interior shell: surface is d=0 and d=-t
                return std::max(d, -(d + thickness));
            } else {
                // Exterior shell: surface is d=0 and d=t
                return std::max(-d, d - thickness);
            }
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
