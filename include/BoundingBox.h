#pragma once

#include "Geometry.h"

namespace Geom {

    struct BoundingBox {
        Point3 min;
        Point3 max;

        // Default constructor: Invalid box
        BoundingBox() 
            : min(INF, INF, INF), max(-INF, -INF, -INF) {}

        BoundingBox(const Point3& p) 
            : min(p), max(p) {}

        BoundingBox(const Point3& min, const Point3& max) 
            : min(min), max(max) {}

        void expand(const Point3& p) {
            min = Geom::min(min, p);
            max = Geom::max(max, p);
        }

        void expand(const BoundingBox& other) {
            min = Geom::min(min, other.min);
            max = Geom::max(max, other.max);
        }

        bool isValid() const {
            return min.x <= max.x && min.y <= max.y && min.z <= max.z;
        }

        Point3 center() const {
            return (min + max) * 0.5;
        }

        Vec3 size() const {
            return max - min;
        }
    };
}
