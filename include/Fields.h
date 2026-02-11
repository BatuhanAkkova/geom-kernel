#pragma once

#include "Field.h"
#include <cmath>
#include <algorithm>

namespace Geom {

    /**
     * @brief A field that returns a constant value everywhere.
     */
    class ConstantField : public Field {
        Scalar value;
    public:
        ConstantField(Scalar v) : value(v) {}
        Scalar eval(const Point3& p) const override {
            return value;
        }
    };

    /**
     * @brief A linear ramp field along an axis.
     * Value = (p - origin) . direction / length
     * Clamped to [0, 1] if clamped is true, or just raw projection.
     * Let's make it more general: val = startVal + (endVal - startVal) * t
     * where t is the projection along the axis, normalized by length.
     */
    class RampField : public Field {
        Point3 origin;
        Vec3 direction; // Normalized direction
        Scalar length;
        Scalar startVal;
        Scalar endVal;
        bool clamp;

    public:
        RampField(Point3 origin, Vec3 direction, Scalar length, Scalar startVal, Scalar endVal, bool clamp = true)
            : origin(origin), direction(direction.normalized()), length(length), startVal(startVal), endVal(endVal), clamp(clamp) {}

        Scalar eval(const Point3& p) const override {
            Scalar t = (p - origin).dot(direction) / length;
            if (clamp) {
                t = std::max(0.0, std::min(1.0, t));
            }
            return startVal + (endVal - startVal) * t;
        }
    };

    /**
     * @brief A radial field that falls off from a center point.
     * Value = startVal at center, decreasing/increasing to endVal at radius.
     */
    class RadialField : public Field {
        Point3 center;
        Scalar radius;
        Scalar startVal;
        Scalar endVal;
        bool clamp;

    public:
        RadialField(Point3 center, Scalar radius, Scalar startVal, Scalar endVal, bool clamp = true)
            : center(center), radius(radius), startVal(startVal), endVal(endVal), clamp(clamp) {}

        Scalar eval(const Point3& p) const override {
            Scalar dist = (p - center).length();
            Scalar t = dist / radius;
            if (clamp) {
                t = std::max(0.0, std::min(1.0, t));
            }
            // Linear falloff for now. Could add easing.
            return startVal + (endVal - startVal) * t;
        }
    };

    class FieldAdd : public Field {
        FieldPtr f1, f2;
    public:
        FieldAdd(FieldPtr f1, FieldPtr f2) : f1(f1), f2(f2) {}
        Scalar eval(const Point3& p) const override {
            return f1->eval(p) + f2->eval(p);
        }
    };

    class FieldMultiply : public Field {
        FieldPtr f1, f2;
    public:
        FieldMultiply(FieldPtr f1, FieldPtr f2) : f1(f1), f2(f2) {}
        Scalar eval(const Point3& p) const override {
            return f1->eval(p) * f2->eval(p);
        }
    };

    class FieldMix : public Field {
        FieldPtr f1, f2, weight;
    public:
        FieldMix(FieldPtr f1, FieldPtr f2, FieldPtr weight) : f1(f1), f2(f2), weight(weight) {}
        Scalar eval(const Point3& p) const override {
            Scalar w = weight->eval(p);
            return f1->eval(p) * (1.0 - w) + f2->eval(p) * w;
        }
    };


}
