#pragma once

#include "Geometry.h"
#include <cmath>

namespace Geom {

    /**
     * @brief Dual Number for Forward-Mode Automatic Differentiation.
     * Represents a value as (v + d*epsilon), where v is the value and d is the derivative.
     */
    template <typename T = Scalar>
    struct Dual {
        T val;
        T der;

        constexpr Dual(T v = 0) : val(v), der(0) {}
        constexpr Dual(T v, T d) : val(v), der(d) {}

        static Dual variable(T v) { return Dual(v, 1); }

        Dual operator+(const Dual& other) const { return {val + other.val, der + other.der}; }
        Dual operator-(const Dual& other) const { return {val - other.val, der - other.der}; }
        Dual operator*(const Dual& other) const { 
            return {val * other.val, val * other.der + der * other.val}; 
        }
        Dual operator/(const Dual& other) const {
            return {val / other.val, (der * other.val - val * other.der) / (other.val * other.val)};
        }

        Dual operator+(T s) const { return {val + s, der}; }
        Dual operator-(T s) const { return {val - s, der}; }
        Dual operator*(T s) const { return {val * s, der * s}; }
        Dual operator/(T s) const { return {val / s, der / s}; }

        friend Dual operator+(T s, const Dual& d) { return d + s; }
        friend Dual operator-(T s, const Dual& d) { return Dual(s - d.val, -d.der); }
        friend Dual operator*(T s, const Dual& d) { return d * s; }
        friend Dual operator/(T s, const Dual& d) { 
            return Dual(s / d.val, -s * d.der / (d.val * d.val)); 
        }

        bool operator<(const Dual& other) const { return val < other.val; }
        bool operator>(const Dual& other) const { return val > other.val; }
        bool operator<=(const Dual& other) const { return val <= other.val; }
        bool operator>=(const Dual& other) const { return val >= other.val; }
    };

    /**
     * @brief Differentiable 3D Vector.
     */
    template <typename T>
    struct Vec3T {
        T x, y, z;

        constexpr Vec3T() : x(0), y(0), z(0) {}
        constexpr Vec3T(T x, T y, T z) : x(x), y(y), z(z) {}

        Vec3T operator+(const Vec3T& other) const { return {x + other.x, y + other.y, z + other.z}; }
        Vec3T operator-(const Vec3T& other) const { return {x - other.x, y - other.y, z - other.z}; }
        
        template <typename S>
        Vec3T operator*(S s) const { return {x * s, y * s, z * s}; }
        
        template <typename S>
        Vec3T operator/(S s) const { return {x / s, y / s, z / s}; }

        T dot(const Vec3T& other) const { return x * other.x + y * other.y + z * other.z; }
        T lengthSquared() const { return x * x + y * y + z * z; }
        T length() const { using std::sqrt; return sqrt(lengthSquared()); }
    };

    using DualScalar = Dual<Scalar>;
    using Point3D = Vec3T<DualScalar>;
    using Vec3D = Vec3T<DualScalar>;

    // Math functions for Dual
    template <typename T>
    Dual<T> sqrt(const Dual<T>& d) {
        T s = std::sqrt(d.val);
        return {s, d.der / (static_cast<T>(2) * s)};
    }

    template <typename T>
    Dual<T> abs(const Dual<T>& d) {
        return {std::abs(d.val), d.val >= 0 ? d.der : -d.der};
    }

    template <typename T>
    Dual<T> max(const Dual<T>& a, const Dual<T>& b) {
        return a.val > b.val ? a : b;
    }

    template <typename T>
    Dual<T> min(const Dual<T>& a, const Dual<T>& b) {
        return a.val < b.val ? a : b;
    }

    template <typename T>
    Dual<T> max(const Dual<T>& a, T b) {
        return a.val > b ? a : Dual<T>(b, 0);
    }

    template <typename T>
    Dual<T> min(const Dual<T>& a, T b) {
        return a.val < b ? a : Dual<T>(b, 0);
    }
}
