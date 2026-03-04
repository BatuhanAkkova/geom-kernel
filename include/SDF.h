#pragma once

#include "Geometry.h"
#include "BoundingBox.h"
#include "Field.h"
#include <memory>

namespace Geom {

    /**
     * @brief Abstract Base Class for all Signed Distance Functions.
     * 
     * An SDF represents a solid geometry where:
     * - eval(p) < 0 : Inside the object
     * - eval(p) > 0 : Outside the object
     * - eval(p) = 0 : On the surface
     * 
     * It extends Field by providing BoundingBoxes and specific surface methods.
     */
    class SDF : public Field {
    public:
        virtual ~SDF() = default;

        /**
         * @brief Get the axis-aligned bounding box of this SDF.
         * @return BoundingBox Guaranteed to contain the zero-isosurface.
         */
        virtual BoundingBox boundingBox() const = 0;

        // The gradient and hessian default implementations use AD from Field evaluation.
        
        /**
         * @brief Calculate the gradient (surface normal * length) of the SDF at p.
         * Default implementation uses forward-mode AD (evalD).
         */
        virtual Vec3 gradient(const Point3& p) const {
            DualScalar dx = evalD(Point3D(DualScalar::variable(p.x), DualScalar::constant(p.y), DualScalar::constant(p.z)));
            DualScalar dy = evalD(Point3D(DualScalar::constant(p.x), DualScalar::variable(p.y), DualScalar::constant(p.z)));
            DualScalar dz = evalD(Point3D(DualScalar::constant(p.x), DualScalar::constant(p.y), DualScalar::variable(p.z)));
            return Vec3(dx.der, dy.der, dz.der);
        }

        /**
         * @brief Evaluate the analytical Hessian matrix at point p.
         * Uses nested second-order forward-mode AD (evalD2).
         */
        virtual Mat3 hessian(const Point3& p) const {
            if (!enableHessian) return Mat3();

            Mat3 h;
            auto seed = [&](int i, int j) -> Scalar {
                Dual2Scalar x(DualScalar::constant(p.x));
                Dual2Scalar y(DualScalar::constant(p.y));
                Dual2Scalar z(DualScalar::constant(p.z));
                
                if (i == 0) x.der = DualScalar::constant(1.0);
                else if (i == 1) y.der = DualScalar::constant(1.0);
                else if (i == 2) z.der = DualScalar::constant(1.0);
                
                if (j == 0) x.val.der = 1.0;
                else if (j == 1) y.val.der = 1.0;
                else if (j == 2) z.val.der = 1.0;
                
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

    /**
     * @brief Architecture equivalent of FieldNode for SDF nodes.
     * Inheriting from SDFNode automatically generates generic `eval`, `evalD`, `evalD2`.
     * Note: Some specific classes (like Booleans) might prefer extending FieldNode or SDFNode,
     * but usually geometric combinations should extent SDFNode to inherit BoundingBox handling.
     */
    template <typename Derived>
    class SDFNode : public SDF {
    public:
        Scalar eval(const Point3& p) const override {
            Vec3T<Scalar> pt(p.x, p.y, p.z);
            return static_cast<const Derived*>(this)->template evaluate<Scalar>(pt);
        }

        DualScalar evalD(const Point3D& p) const override {
            return static_cast<const Derived*>(this)->template evaluate<DualScalar>(p);
        }

        Dual2Scalar evalD2(const Point3D2& p) const override {
            if (!Field::enableHessian) {
                return Dual2Scalar::constant(evalD(Point3D(p.x.val.val, p.y.val.val, p.z.val.val)));
            }
            return static_cast<const Derived*>(this)->template evaluate<Dual2Scalar>(p);
        }
    };

    using SDFPtr = std::shared_ptr<SDF>;
}
