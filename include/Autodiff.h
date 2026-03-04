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

        constexpr Dual(T v = 0) : val(v), der(static_cast<T>(0)) {}
        constexpr Dual(T v, T d) : val(v), der(d) {}

        static Dual variable(T v) { return Dual(v, static_cast<T>(1)); }
        static Dual constant(T v) { return Dual(v, static_cast<T>(0)); }

        // Binary operators with other Duals
        template <typename U>
        auto operator+(const Dual<U>& other) const { return Dual<decltype(val + other.val)>(val + other.val, der + other.der); }
        template <typename U>
        auto operator-(const Dual<U>& other) const { return Dual<decltype(val - other.val)>(val - other.val, der - other.der); }
        template <typename U>
        auto operator*(const Dual<U>& other) const { 
            return Dual<decltype(val * other.val)>(val * other.val, val * other.der + der * other.val); 
        }
        template <typename U>
        auto operator/(const Dual<U>& other) const {
            auto v = val / other.val;
            return Dual<decltype(v)>(v, (der * other.val - val * other.der) / (other.val * other.val));
        }

        // Binary operators with Scalars/Types
        template <typename U>
        auto operator+(U s) const { return Dual<T>(val + s, der); }
        template <typename U>
        auto operator-(U s) const { return Dual<T>(val - s, der); }
        template <typename U>
        auto operator*(U s) const { return Dual<T>(val * s, der * s); }
        template <typename U>
        auto operator/(U s) const { return Dual<T>(val / s, der / s); }

        template <typename U>
        friend auto operator+(U s, const Dual& d) { return d + s; }
        template <typename U>
        friend auto operator-(U s, const Dual& d) { return Dual<T>(s - d.val, -d.der); }
        template <typename U>
        friend auto operator*(U s, const Dual& d) { return d * s; }
        template <typename U>
        friend auto operator/(U s, const Dual& d) { 
            return Dual<T>(s / d.val, -static_cast<T>(s) * d.der / (d.val * d.val)); 
        }

        Dual operator-() const { return {-val, -der}; }

        // Comparison
        template <typename U>
        bool operator<(const Dual<U>& o) const { return val < o.val; }
        template <typename U>
        bool operator>(const Dual<U>& o) const { return val > o.val; }
        template <typename U>
        bool operator<=(const Dual<U>& o) const { return val <= o.val; }
        template <typename U>
        bool operator>=(const Dual<U>& o) const { return val >= o.val; }
        template <typename U>
        bool operator==(const Dual<U>& o) const { return val == o.val; }
        template <typename U>
        bool operator!=(const Dual<U>& o) const { return val != o.val; }

        template <typename U>
        bool operator<(U s) const { return val < s; }
        template <typename U>
        bool operator>(U s) const { return val > s; }
        template <typename U>
        bool operator<=(U s) const { return val <= s; }
        template <typename U>
        bool operator>=(U s) const { return val >= s; }
        template <typename U>
        bool operator==(U s) const { return val == s; }
        template <typename U>
        bool operator!=(U s) const { return val != s; }

        template <typename U>
        friend bool operator<(U s, const Dual& d) { return s < d.val; }
        template <typename U>
        friend bool operator>(U s, const Dual& d) { return s > d.val; }
        template <typename U>
        friend bool operator<=(U s, const Dual& d) { return s <= d.val; }
        template <typename U>
        friend bool operator>=(U s, const Dual& d) { return s >= d.val; }
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
        T length() const { 
            using std::sqrt; 
            T lsq = lengthSquared();
            // Protect against zero length causing NaN in gradients
            if constexpr (std::is_same_v<T, Scalar>) {
                return sqrt(lsq);
            } else {
                return lsq.val > 1e-12 ? sqrt(lsq) : T(0); 
            }
        }
    };

    using DualScalar = Dual<Scalar>;
    using Point3D = Vec3T<DualScalar>;
    using Vec3D = Vec3T<DualScalar>;

    // Second order AD types
    using Dual2Scalar = Dual<DualScalar>;
    using Point3D2 = Vec3T<Dual2Scalar>;

    // Math functions for Dual (Generic templates)
    template <typename T>
    auto sqrt(const Dual<T>& d) {
        using std::sqrt;
        using std::sqrt;
        auto s = sqrt(d.val);
        // Avoid division by zero when val is exactly 0
        if (s == 0) return Dual<decltype(s)>(0, 0); 
        return Dual<decltype(s)>(s, d.der / (static_cast<decltype(s)>(2) * s));
    }

    template <typename T>
    auto abs(const Dual<T>& d) {
        using std::abs;
        auto v = abs(d.val);
        return Dual<decltype(v)>(v, d.val >= 0 ? d.der : -d.der);
    }

    template <typename T, typename U>
    auto max(const Dual<T>& a, const Dual<U>& b) {
        return a.val > b.val ? a : b;
    }

    template <typename T, typename U>
    auto min(const Dual<T>& a, const Dual<U>& b) {
        return a.val < b.val ? a : b;
    }

    template <typename T, typename U>
    auto max(const Dual<T>& a, U b) {
        return a.val > b ? a : Dual<T>(static_cast<T>(b), static_cast<T>(0));
    }

    template <typename T, typename U>
    auto min(const Dual<T>& a, U b) {
        return a.val < b ? a : Dual<T>(static_cast<T>(b), static_cast<T>(0));
    }

    template <typename T>
    auto sin(const Dual<T>& d) {
        using std::sin; using std::cos;
        return Dual<decltype(sin(d.val))>(sin(d.val), d.der * cos(d.val));
    }

    template <typename T>
    auto cos(const Dual<T>& d) {
        using std::sin; using std::cos;
        return Dual<decltype(cos(d.val))>(cos(d.val), -d.der * sin(d.val));
    }

    template <typename T, typename S>
    auto pow(const Dual<T>& d, S p) {
        using std::pow;
        auto v = pow(d.val, p);
        return Dual<decltype(v)>(v, static_cast<decltype(v)>(p) * pow(d.val, p - 1) * d.der);
    }

    /**
     * @brief Enzyme AD Hooks Documentation
     * 
     * To use Enzyme (https://enzyme.mit.edu) for reverse-mode AD:
     * 1. Structure SDF::eval as a leaf function.
     * 2. Use __enzyme_autodiff<double>(SDF::eval, ...) to generate gradients.
     * 
     * Example:
     * extern double __enzyme_autodiff(void*, ...);
     * double grad_p[3];
     * __enzyme_autodiff((void*)my_sdf_eval, enzyme_dup, p, grad_p);
     */
    struct EnzymeHook {};
}
