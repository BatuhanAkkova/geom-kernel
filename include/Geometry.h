#pragma once

#include <cmath>
#include <limits>
#include <concepts>

namespace Geom {

    // Precision choice: Double for engineering accuracy
    using Scalar = double;
    
    constexpr Scalar EPSILON = 1e-6;
    constexpr Scalar INF = std::numeric_limits<Scalar>::infinity();

    struct Vec3 {
        Scalar x, y, z;

        constexpr Vec3() : x(0), y(0), z(0) {}
        constexpr Vec3(Scalar x, Scalar y, Scalar z) : x(x), y(y), z(z) {}

        Vec3 operator+(const Vec3& other) const { return {x + other.x, y + other.y, z + other.z}; }
        Vec3 operator-(const Vec3& other) const { return {x - other.x, y - other.y, z - other.z}; }
        Vec3 operator*(Scalar s) const { return {x * s, y * s, z * s}; }
        Vec3 operator/(Scalar s) const { return {x / s, y / s, z / s}; }

        Scalar dot(const Vec3& other) const { return x * other.x + y * other.y + z * other.z; }
        
        Scalar lengthSquared() const { return x * x + y * y + z * z; }
        Scalar length() const { return std::sqrt(lengthSquared()); }

        Vec3 normalized() const {
            Scalar len = length();
            return (len > EPSILON) ? *this / len : Vec3(0, 0, 0);
        }
    };

    using Point3 = Vec3;

    inline Vec3 min(const Vec3& a, const Vec3& b) {
        return {std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z)};
    }

    inline Vec3 max(const Vec3& a, const Vec3& b) {
        return {std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z)};
    }
}
