#pragma once

#include "Geometry.h"
#include "BoundingBox.h"
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
         * @brief Get the axis-aligned bounding box of this SDF.
         * @return BoundingBox Guaranteed to contain the zero-isosurface.
         */
        virtual BoundingBox boundingBox() const = 0;
    };

    using SDFPtr = std::shared_ptr<SDF>;
}
