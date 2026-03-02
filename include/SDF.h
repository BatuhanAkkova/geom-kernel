#pragma once

#include "Geometry.h"
#include "BoundingBox.h"
#include "Autodiff.h"
#include <memory>

namespace Geom {

    /**
     * @brief Abstract Base Class for all Signed Distance Functions.
     * 
     * An SDF represents a solid geometry where:
     * - eval(p) < 0 : Inside the object
     * - eval(p) > 0 : Outside the object
     * - eval(p) = 0 : On the surface
     */
    class SDF {
    public:
        virtual ~SDF() = default;

        /**
         * @brief Evaluate the signed distance at a given point in space.
         * @param p The query point.
         * @return Scalar Distance to the nearest surface. Negative if inside.
         */
        virtual Scalar eval(const Point3& p) const = 0;

        /**
         * @brief Templated evaluate for Differentiable points.
         * Default implementation uses central differences if not overridden, 
         * but we prefer to override for exact analytical derivatives.
         */
        virtual DualScalar evalD(const Point3D& p) const {
            // Default fallback: Finite differences or just evaluate val
            // Most primitives will override this for exact autodiff.
            return DualScalar(eval(Point3(p.x.val, p.y.val, p.z.val)));
        }

        /**
         * @brief Calculate the gradient (surface normal * length) of the SDF at p.
         * Default implementation uses central finite differences.
         * @param p The query point
         * @param h The step size for finite difference
         * @return Vec3 The gradient vector
         */
        virtual Vec3 gradient(const Point3& p, Scalar h = 0.0001) const {
            const Scalar dx = eval(p + Vec3(h, 0, 0)) - eval(p - Vec3(h, 0, 0));
            const Scalar dy = eval(p + Vec3(0, h, 0)) - eval(p - Vec3(0, h, 0));
            const Scalar dz = eval(p + Vec3(0, 0, h)) - eval(p - Vec3(0, 0, h));
            return Vec3(dx, dy, dz).normalized(); // Note: This is actually the normal if normalized
        }

        /**
         * @brief Get the axis-aligned bounding box of this SDF.
         * @return BoundingBox Guaranteed to contain the zero-isosurface.
         */
        virtual BoundingBox boundingBox() const = 0;
    };

    using SDFPtr = std::shared_ptr<SDF>;
}
