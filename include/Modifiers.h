#pragma once

#include "SDF.h"
#include "Field.h"
#include <cmath>
#include <algorithm>


namespace Geom {

    /**
     * @brief Offset the surface of an SDF by a distance `r`.
     */
    class Offset : public SDFNode<Offset> {
        SDFPtr sdf;
        Scalar r;
    public:
        Offset(SDFPtr sdf, Scalar r) : sdf(sdf), r(r) {}

        template <typename T>
        T evaluate(const Vec3T<T>& p) const {
            using std::max; using std::min; using std::abs; using std::sqrt; using std::pow; using std::sin; using std::cos;
            if constexpr (std::is_same_v<T, Scalar>) {
                Point3 pt(p.x, p.y, p.z);
                return sdf->eval(pt) - r;
            } else if constexpr (std::is_same_v<T, DualScalar>) {
                Point3D pt(p.x, p.y, p.z);
                return sdf->evalD(pt) - r;
            } else if constexpr (std::is_same_v<T, Dual2Scalar>) {
                Point3D2 pt(p.x, p.y, p.z);
                return sdf->evalD2(pt) - r;
            }
            return static_cast<T>(0);
        }

        BoundingBox boundingBox() const override {
            BoundingBox box = sdf->boundingBox();
            box.min -= Vec3(r, r, r);
            box.max += Vec3(r, r, r);
            return box;
        }
    };

    /**
     * @brief Create a shell (hollow object) from an SDF.
     * Effectively `abs(eval(p)) - thickness`.
     */
    class Shell : public SDFNode<Shell> {
        SDFPtr sdf;
        Scalar thickness;
    public:
        Shell(SDFPtr sdf, Scalar thickness) : sdf(sdf), thickness(thickness) {}

        template <typename T>
        T evaluate(const Vec3T<T>& p) const {
            using std::max; using std::min; using std::abs; using std::sqrt; using std::pow; using std::sin; using std::cos;
            if constexpr (std::is_same_v<T, Scalar>) {
                Point3 pt(p.x, p.y, p.z);
                return std::abs(sdf->eval(pt)) - thickness;
            } else if constexpr (std::is_same_v<T, DualScalar>) {
                Point3D pt(p.x, p.y, p.z);
                return abs(sdf->evalD(pt)) - thickness;
            } else if constexpr (std::is_same_v<T, Dual2Scalar>) {
                Point3D2 pt(p.x, p.y, p.z);
                return abs(sdf->evalD2(pt)) - thickness;
            }
            return static_cast<T>(0);
        }

        BoundingBox boundingBox() const override {
            BoundingBox box = sdf->boundingBox();
            box.min -= Vec3(thickness, thickness, thickness);
            box.max += Vec3(thickness, thickness, thickness);
            return box;
        }
    };

    /**
     * @brief Offset the surface of an SDF by a distance `r` defined by a Field.
     */
    class FieldOffset : public SDFNode<FieldOffset> {
        SDFPtr sdf;
        FieldPtr field;
    public:
        FieldOffset(SDFPtr sdf, FieldPtr field) : sdf(sdf), field(field) {}

        template <typename T>
        T evaluate(const Vec3T<T>& p) const {
            using std::max; using std::min; using std::abs; using std::sqrt; using std::pow; using std::sin; using std::cos;
            if constexpr (std::is_same_v<T, Scalar>) {
                Point3 pt(p.x, p.y, p.z);
                return sdf->eval(pt) - field->eval(pt);
            } else if constexpr (std::is_same_v<T, DualScalar>) {
                Point3D pt(p.x, p.y, p.z);
                return sdf->evalD(pt) - field->evalD(pt);
            } else if constexpr (std::is_same_v<T, Dual2Scalar>) {
                Point3D2 pt(p.x, p.y, p.z);
                return sdf->evalD2(pt) - field->evalD2(pt);
            }
            return static_cast<T>(0);
        }

        BoundingBox boundingBox() const override {
            BoundingBox box = sdf->boundingBox();
            Vec3 margin(5.0, 5.0, 5.0); 
            box.min -= margin;
            box.max += margin;
            return box;
        }
    };

    /**
     * @brief Twist the space around the Y axis.
     * angle = k * y
     */
    class Twist : public SDFNode<Twist> {
        SDFPtr sdf;
        Scalar k; // Twist amount
    public:
        Twist(SDFPtr sdf, Scalar k) : sdf(sdf), k(k) {}

        template <typename T>
        T evaluate(const Vec3T<T>& p) const {
            using std::max; using std::min; using std::abs; using std::sqrt; using std::pow; using std::sin; using std::cos;
            if constexpr (std::is_same_v<T, Scalar>) {
                Scalar c = std::cos(k * p.y);
                Scalar s = std::sin(k * p.y);
                Scalar qx = c * p.x - s * p.z;
                Scalar qz = s * p.x + c * p.z;
                Point3 q(qx, p.y, qz);
                return sdf->eval(q);
            } else if constexpr (std::is_same_v<T, DualScalar>) {
                using std::cos;
                using std::sin;
                DualScalar c = cos(p.y * k);
                DualScalar s = sin(p.y * k);
                DualScalar qx = c * p.x - s * p.z;
                DualScalar qz = s * p.x + c * p.z;
                Point3D q(qx, p.y, qz);
                return sdf->evalD(q);
            } else if constexpr (std::is_same_v<T, Dual2Scalar>) {
                using std::cos;
                using std::sin;
                Dual2Scalar c = cos(p.y * k);
                Dual2Scalar s = sin(p.y * k);
                Dual2Scalar qx = c * p.x - s * p.z;
                Dual2Scalar qz = s * p.x + c * p.z;
                Point3D2 q(qx, p.y, qz);
                return sdf->evalD2(q);
            }
            return static_cast<T>(0);
        }

        BoundingBox boundingBox() const override {
            BoundingBox box = sdf->boundingBox();
            return box; // TODO: Better bound
        }
    };

}
