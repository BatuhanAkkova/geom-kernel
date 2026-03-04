#pragma once

#include "Geometry.h"
#include "Autodiff.h"
#include <memory>
#include <vector>
#include <cmath>
#include <algorithm>

namespace Geom {

    /**
     * @brief Abstract Base Class for all Scalar Fields.
     * 
     * A Field represents a scalar value that varies over space.
     * With refactoring, Fields also expose parameters native to their operation,
     * and automatic differentiation variables for gradients and Hessians.
     */
    class Field {
    public:
        virtual ~Field() = default;

        /**
         * @brief Evaluate the field value at a given point in space.
         * @param p The query point.
         * @return Scalar value of the field at p.
         */
        virtual Scalar eval(const Point3& p) const = 0;

        /**
         * @brief Evaluate with first-order dual numbers.
         */
        virtual DualScalar evalD(const Point3D& p) const = 0;

        /**
         * @brief Evaluate with second-order dual numbers.
         */
        virtual Dual2Scalar evalD2(const Point3D2& p) const = 0;

        // --- Parameter Exposure Options ---
        virtual size_t numParams() const { return 0; }
        virtual void setParam(size_t i, Scalar val) {}
        virtual Scalar getParam(size_t i) const { return 0.0; }
        
        /**
         * @brief Evaluates the derivative of the field at p with respect to a specific parameter.
         * @param p Evaluated position.
         * @param paramIndex Index of the parameter (defined per node).
         */
        virtual DualScalar evalParamD(const Point3& p, size_t paramIndex) const {
            return DualScalar::constant(eval(p));
        }

        // Global flag to disable Hessian computations across the entire tree if they are not needed.
        static bool enableHessian;
    };

    using FieldPtr = std::shared_ptr<Field>;

    // Inline definition of static flag for linking purposes.
    inline bool Field::enableHessian = true;

    /**
     * @brief CRTP Base Class for specific field nodes to auto-generate `eval`, `evalD`, and `evalD2`.
     * Derived classes only need to implement:
     * `template <typename T> T evaluate(const Vec3T<T>& p) const`
     */
    template <typename Derived>
    class FieldNode : public Field {
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
                // If Hessian is disabled, avoid deep template recursion and just return standard eval 
                // nested within dual constants.
                return Dual2Scalar::constant(evalD(Point3D(p.x.val.val, p.y.val.val, p.z.val.val)));
            }
            return static_cast<const Derived*>(this)->template evaluate<Dual2Scalar>(p);
        }
    };

    /**
     * @brief A field that returns a constant value everywhere.
     */
    class ConstantField : public FieldNode<ConstantField> {
        Scalar value;
    public:
        ConstantField(Scalar v) : value(v) {}
        
        template <typename T>
        T evaluate(const Vec3T<T>& p) const {
            return static_cast<T>(value);
        }

        size_t numParams() const override { return 1; }
        Scalar getParam(size_t i) const override { return value; }
        void setParam(size_t i, Scalar val) override { if(i==0) value = val; }
        
        DualScalar evalParamD(const Point3& p, size_t paramIndex) const override {
            if(paramIndex == 0) return DualScalar::variable(value);
            return DualScalar::constant(value);
        }
    };

    /**
     * @brief A linear ramp field along an axis.
     */
    class RampField : public FieldNode<RampField> {
        Point3 origin;
        Vec3 direction; 
        Scalar length;
        Scalar startVal;
        Scalar endVal;
        bool clamp_t;

    public:
        RampField(Point3 origin, Vec3 direction, Scalar length, Scalar startVal, Scalar endVal, bool clamp = true)
            : origin(origin), direction(direction.normalized()), length(length), startVal(startVal), endVal(endVal), clamp_t(clamp) {}

        template <typename T>
        T evaluate(const Vec3T<T>& p) const {
            Vec3T<T> t_origin(static_cast<T>(origin.x), static_cast<T>(origin.y), static_cast<T>(origin.z));
            Vec3T<T> p_rel = p - t_origin;
            
            Vec3T<T> t_dir(static_cast<T>(direction.x), static_cast<T>(direction.y), static_cast<T>(direction.z));
            
            T t = p_rel.dot(t_dir) / static_cast<T>(length);
            
            if (clamp_t) {
                if constexpr (std::is_same_v<T, Scalar>) {
                    t = std::max(static_cast<T>(0.0), std::min(static_cast<T>(1.0), t));
                } else {
                    t = max(static_cast<T>(0.0), min(static_cast<T>(1.0), t));
                }
            }
            return static_cast<T>(startVal) + (static_cast<T>(endVal) - static_cast<T>(startVal)) * t;
        }
    };

    /**
     * @brief A radial field that falls off from a center point.
     */
    class RadialField : public FieldNode<RadialField> {
        Point3 center;
        Scalar radius;
        Scalar startVal;
        Scalar endVal;
        bool clamp_t;

    public:
        RadialField(Point3 center, Scalar radius, Scalar startVal, Scalar endVal, bool clamp = true)
            : center(center), radius(radius), startVal(startVal), endVal(endVal), clamp_t(clamp) {}

        template <typename T>
        T evaluate(const Vec3T<T>& p) const {
            Vec3T<T> t_center(static_cast<T>(center.x), static_cast<T>(center.y), static_cast<T>(center.z));
            T dist;
            if constexpr (std::is_same_v<T, Scalar>) {
                dist = (p - t_center).length();
            } else {
                dist = (p - t_center).length();
            }
            T t = dist / static_cast<T>(radius);
            if (clamp_t) {
                if constexpr (std::is_same_v<T, Scalar>) {
                    t = std::max(static_cast<T>(0.0), std::min(static_cast<T>(1.0), t));
                } else {
                    t = max(static_cast<T>(0.0), min(static_cast<T>(1.0), t));
                }
            }
            return static_cast<T>(startVal) + (static_cast<T>(endVal) - static_cast<T>(startVal)) * t;
        }
    };

    class FieldAdd : public FieldNode<FieldAdd> {
        FieldPtr f1, f2;
    public:
        FieldAdd(FieldPtr f1, FieldPtr f2) : f1(f1), f2(f2) {}
        
        template <typename T>
        T evaluate(const Vec3T<T>& p) const {
            if constexpr (std::is_same_v<T, Scalar>) {
                Point3 pt(p.x, p.y, p.z);
                return f1->eval(pt) + f2->eval(pt);
            } else if constexpr (std::is_same_v<T, DualScalar>) {
                Point3D pt(p.x, p.y, p.z);
                return f1->evalD(pt) + f2->evalD(pt);
            } else if constexpr (std::is_same_v<T, Dual2Scalar>) {
                Point3D2 pt(p.x, p.y, p.z);
                return f1->evalD2(pt) + f2->evalD2(pt);
            }
            return static_cast<T>(0);
        }
    };

    class FieldMultiply : public FieldNode<FieldMultiply> {
        FieldPtr f1, f2;
    public:
        FieldMultiply(FieldPtr f1, FieldPtr f2) : f1(f1), f2(f2) {}
        
        template <typename T>
        T evaluate(const Vec3T<T>& p) const {
            if constexpr (std::is_same_v<T, Scalar>) {
                Point3 pt(p.x, p.y, p.z);
                return f1->eval(pt) * f2->eval(pt);
            } else if constexpr (std::is_same_v<T, DualScalar>) {
                Point3D pt(p.x, p.y, p.z);
                return f1->evalD(pt) * f2->evalD(pt);
            } else if constexpr (std::is_same_v<T, Dual2Scalar>) {
                Point3D2 pt(p.x, p.y, p.z);
                return f1->evalD2(pt) * f2->evalD2(pt);
            }
            return static_cast<T>(0);
        }
    };

    class FieldMix : public FieldNode<FieldMix> {
        FieldPtr f1, f2, weight;
    public:
        FieldMix(FieldPtr f1, FieldPtr f2, FieldPtr weight) : f1(f1), f2(f2), weight(weight) {}
        
        template <typename T>
        T evaluate(const Vec3T<T>& p) const {
            if constexpr (std::is_same_v<T, Scalar>) {
                Point3 pt(p.x, p.y, p.z);
                Scalar w = weight->eval(pt);
                return f1->eval(pt) * (1.0 - w) + f2->eval(pt) * w;
            } else if constexpr (std::is_same_v<T, DualScalar>) {
                Point3D pt(p.x, p.y, p.z);
                DualScalar w = weight->evalD(pt);
                return f1->evalD(pt) * (1.0 - w) + f2->evalD(pt) * w;
            } else if constexpr (std::is_same_v<T, Dual2Scalar>) {
                Point3D2 pt(p.x, p.y, p.z);
                Dual2Scalar w = weight->evalD2(pt);
                return f1->evalD2(pt) * (1.0 - w) + f2->evalD2(pt) * w;
            }
            return static_cast<T>(0);
        }
    };

}
