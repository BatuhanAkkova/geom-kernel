#pragma once

#include "Geometry.h"
#include <memory>

namespace Geom {

    /**
     * @brief Abstract Base Class for all Scalar Fields.
     * 
     * A Field represents a scalar value that varies over space.
     * Used for driving geometry parameters like thickness, offset, blend radius, etc.
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
    };

    using FieldPtr = std::shared_ptr<Field>;
}
