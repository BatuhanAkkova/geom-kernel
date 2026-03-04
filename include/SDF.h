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
            return DualScalar(eval(Point3(p.x.val, p.y.val, p.z.val)));
        }

        /**
         * @brief Templated evaluate for Second-Order Differentiable points.
         */
        virtual Dual2Scalar evalD2(const Point3D2& p) const {
            // Default fallback: evaluation at the base value
            return Dual2Scalar::constant(evalD(Point3D(p.x.val, p.y.val, p.z.val)));
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

        /**
         * @brief Evaluate the analytical gradient at point p.
         * Default implementation uses forward-mode AD (evalD).
         */
        virtual Vec3 computeAnalyticalGradient(const Point3& p) const {
            DualScalar dx = evalD(Point3D(DualScalar::variable(p.x), DualScalar::constant(p.y), DualScalar::constant(p.z)));
            DualScalar dy = evalD(Point3D(DualScalar::constant(p.x), DualScalar::variable(p.y), DualScalar::constant(p.z)));
            DualScalar dz = evalD(Point3D(DualScalar::constant(p.x), DualScalar::constant(p.y), DualScalar::variable(p.z)));
            return Vec3(dx.der, dy.der, dz.der);
        }

        /**
         * @brief Evaluate the analytical Hessian matrix at point p.
         * Uses nested forward-mode AD.
         */
        virtual Mat3 analyticalHessian(const Point3& p) const {
            Mat3 h;
            
            // We need to extract:
            // d2f/dx2, d2f/dy2, d2f/dz2
            // d2f/dxdy, d2f/dxdz, d2f/dydz
            
            auto seed = [&](int i, int j) -> Scalar {
                // To get d2f/di dj:
                // Seed outer-dual variable on axis i, inner-dual variable on axis j
                Dual2Scalar x(DualScalar::constant(p.x));
                Dual2Scalar y(DualScalar::constant(p.y));
                Dual2Scalar z(DualScalar::constant(p.z));
                
                // Outer derivative (der) corresponds to axis i
                // Inner derivative (val.der) corresponds to axis j
                // The mixed partial d2f/didj is the coefficient of epsilon_i * epsilon_j, 
                // which ends up in .der.der of Dual<Dual<Scalar>>
                
                if (i == 0) x.der = DualScalar::constant(1.0);
                else if (i == 1) y.der = DualScalar::constant(1.0);
                else if (i == 2) z.der = DualScalar::constant(1.0);
                
                if (j == 0) x.val.der = 1.0;
                else if (j == 1) y.val.der = 1.0;
                else if (j == 2) z.val.der = 1.0;
                
                // If i == j, it's a pure second derivative
                // If i != j, it's a mixed partial
                
                return evalD2(Point3D2(x, y, z)).der.der;
            };

            h(0, 0) = seed(0, 0);
            h(1, 1) = seed(1, 1);
            h(2, 2) = seed(2, 2);
            
            h(0, 1) = h(1, 0) = seed(0, 1);
            h(0, 2) = h(2, 0) = seed(0, 2);
            h(1, 2) = h(2, 1) = seed(1, 2);
            
            return h;
        }
    };

    using SDFPtr = std::shared_ptr<SDF>;
}
