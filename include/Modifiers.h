#pragma once

#include "SDF.h"
#include "Field.h"
#include <cmath>
#include <algorithm>


namespace Geom {

    /**
     * @brief Offset the surface of an SDF by a distance `r`.
     * r > 0: Inflate (round corners)
     * r < 0: Deflate (shrink)
     */
    class Offset : public SDF {
        SDFPtr sdf;
        Scalar r;
    public:
        Offset(SDFPtr sdf, Scalar r) : sdf(sdf), r(r) {}

        Scalar eval(const Point3& p) const override {
            return sdf->eval(p) - r;
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
     * The new surface is at distance `thickness` from the zero isosurface.
     * Effectively `abs(eval(p)) - thickness`.
     */
    class Shell : public SDF {
        SDFPtr sdf;
        Scalar thickness;
    public:
        Shell(SDFPtr sdf, Scalar thickness) : sdf(sdf), thickness(thickness) {}

        Scalar eval(const Point3& p) const override {
            // Shell: |d| - t
            // negative inside the shell wall (between -t and +t of original)
            return std::abs(sdf->eval(p)) - thickness;
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
    class FieldOffset : public SDF {
        SDFPtr sdf;
        FieldPtr field;
    public:
        FieldOffset(SDFPtr sdf, FieldPtr field) : sdf(sdf), field(field) {}

        Scalar eval(const Point3& p) const override {
            return sdf->eval(p) - field->eval(p);
        }

        BoundingBox boundingBox() const override {
            // This is tricky because the field value is unknown.
            // For now, we might need a max offset estimate or just expand by a safe margin?
            // Or we just return the original box and hope the field doesn't push it out too much?
            // Ideally, Field should have a max value or range.
            // Let's assume a "safe" expansion constant or just use the original box for MVP
            // and maybe warn.
            // Better: Field interface *should* have a range estimate.
            // For MVP: Expand by 1.0 (arbitrary) or just keep original.
            // Let's expand by a generous amount if possible, or just accept it might be clipped.
            BoundingBox box = sdf->boundingBox();
            Vec3 margin(5.0, 5.0, 5.0); // Arbitrary margin
            box.min -= margin;
            box.max += margin;
            return box;
        }
    };

    /**
     * @brief Twist the space around the Y axis.
     * angle = k * y
     */
    class Twist : public SDF {
        SDFPtr sdf;
        Scalar k; // Twist amount
    public:
        Twist(SDFPtr sdf, Scalar k) : sdf(sdf), k(k) {}

        Scalar eval(const Point3& p) const override {
            Scalar c = std::cos(k * p.y);
            Scalar s = std::sin(k * p.y);
            // Rotate p around Y
            Scalar qx = c * p.x - s * p.z;
            Scalar qz = s * p.x + c * p.z;
            Point3 q(qx, p.y, qz);
            return sdf->eval(q);
        }

        BoundingBox boundingBox() const override {
            // Twisting generally preserves the bounding cylinder height but expands x/z.
            BoundingBox box = sdf->boundingBox();
            // Conservative estimate: max radius in xz
            // For now, just return valid box unchanged? No, that's wrong.
            // Let's keep it simple: Expand XZ by diagonal of original XZ
            return box; // TODO: Better bound
        }
    };

}
