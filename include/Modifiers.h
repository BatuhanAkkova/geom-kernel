#pragma once

#include "SDF.h"
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

}
