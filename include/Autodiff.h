#pragma once

#include "Geometry.h"
#include <cmath>
#include <type_traits>

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

        // Basic arithmetic
        template <typename U>
        auto operator+(U s) const { return Dual<T>(val + s, der); }
        template <typename U>
        auto operator-(U s) const { return Dual<T>(val - s, der); }
        template <typename U>
        auto operator*(U s) const { return Dual<T>(val * s, der * s); }
        template <typename U>
        auto operator/(U s) const { return Dual<T>(val / s, der / s); }

        template <typename U>
        friend auto operator+(U s, const Dual& d) -> std::enable_if_t<std::is_arithmetic_v<U> || std::is_same_v<U, Scalar>, Dual<T>> { 
            return d + s; 
        }
        template <typename U>
        friend auto operator-(U s, const Dual& d) -> std::enable_if_t<std::is_arithmetic_v<U> || std::is_same_v<U, Scalar>, Dual<T>> { 
            return Dual<T>(s - d.val, -d.der); 
        }
        template <typename U>
        friend auto operator*(U s, const Dual& d) -> std::enable_if_t<std::is_arithmetic_v<U> || std::is_same_v<U, Scalar>, Dual<T>> { 
            return d * s; 
        }
        template <typename U>
        friend auto operator/(U s, const Dual& d) -> std::enable_if_t<std::is_arithmetic_v<U> || std::is_same_v<U, Scalar>, Dual<T>> { 
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


    using DualScalar = Dual<Scalar>;
    using Point3D = Vec3T<DualScalar>;
    using Vec3D = Vec3T<DualScalar>;

    // Second order AD types
    using Dual2Scalar = Dual<DualScalar>;
    using Point3D2 = Vec3T<Dual2Scalar>;

    // Math functions for DualScalar (First-order AD)
    inline DualScalar sqrt(const DualScalar& d) {
        using std::sqrt;
        auto s = sqrt(d.val);
        if (s == 0) return DualScalar(0, 0); 
        return DualScalar(s, d.der / (2.0 * s));
    }

    inline DualScalar abs(const DualScalar& d) {
        using std::abs;
        auto v = abs(d.val);
        return DualScalar(v, d.val >= 0 ? d.der : -d.der);
    }

    inline DualScalar sin(const DualScalar& d) {
        using std::sin; using std::cos;
        return DualScalar(sin(d.val), d.der * cos(d.val));
    }

    inline DualScalar cos(const DualScalar& d) {
        using std::sin; using std::cos;
        return DualScalar(cos(d.val), -d.der * sin(d.val));
    }

    inline DualScalar pow(const DualScalar& d, double p) {
        using std::pow;
        auto v = pow(d.val, p);
        return DualScalar(v, p * pow(d.val, p - 1.0) * d.der);
    }

    inline DualScalar max(const DualScalar& a, const DualScalar& b) {
        return a.val > b.val ? a : b;
    }

    inline DualScalar min(const DualScalar& a, const DualScalar& b) {
        return a.val < b.val ? a : b;
    }

    inline DualScalar max(const DualScalar& a, double b) {
        return a.val > b ? a : DualScalar(b, 0.0);
    }

    inline DualScalar max(double a, const DualScalar& b) {
        return a > b.val ? DualScalar(a, 0.0) : b;
    }

    inline DualScalar min(const DualScalar& a, double b) {
        return a.val < b ? a : DualScalar(b, 0.0);
    }

    inline DualScalar min(double a, const DualScalar& b) {
        return a < b.val ? DualScalar(a, 0.0) : b;
    }

    // Math functions for Dual2Scalar (Second-order AD)
    inline Dual2Scalar sqrt(const Dual2Scalar& d) {
        auto s = sqrt(d.val); // Uses DualScalar sqrt
        if (s.val == 0) return Dual2Scalar(DualScalar(0, 0), DualScalar(0, 0)); 
        return Dual2Scalar(s, d.der / (static_cast<DualScalar>(2.0) * s));
    }

    inline Dual2Scalar abs(const Dual2Scalar& d) {
        auto v = abs(d.val); // Uses DualScalar abs
        return Dual2Scalar(v, d.val.val >= 0 ? d.der : -d.der);
    }

    inline Dual2Scalar sin(const Dual2Scalar& d) {
        return Dual2Scalar(sin(d.val), d.der * cos(d.val));
    }

    inline Dual2Scalar cos(const Dual2Scalar& d) {
        return Dual2Scalar(cos(d.val), -d.der * sin(d.val));
    }

    inline Dual2Scalar pow(const Dual2Scalar& d, double p) {
        auto v = pow(d.val, p);
        return Dual2Scalar(v, static_cast<DualScalar>(p) * pow(d.val, p - 1.0) * d.der);
    }

    inline Dual2Scalar max(const Dual2Scalar& a, const Dual2Scalar& b) {
        return a.val.val > b.val.val ? a : b;
    }

    inline Dual2Scalar min(const Dual2Scalar& a, const Dual2Scalar& b) {
        return a.val.val < b.val.val ? a : b;
    }

    inline Dual2Scalar max(const Dual2Scalar& a, double b) {
        return a.val.val > b ? a : Dual2Scalar(DualScalar(b, 0.0), DualScalar(0.0, 0.0));
    }

    inline Dual2Scalar max(double a, const Dual2Scalar& b) {
        return a > b.val.val ? Dual2Scalar(DualScalar(a, 0.0), DualScalar(0.0, 0.0)) : b;
    }

    inline Dual2Scalar min(const Dual2Scalar& a, double b) {
        return a.val.val < b ? a : Dual2Scalar(DualScalar(b, 0.0), DualScalar(0.0, 0.0));
    }

    inline Dual2Scalar min(double a, const Dual2Scalar& b) {
        return a < b.val.val ? Dual2Scalar(DualScalar(a, 0.0), DualScalar(0.0, 0.0)) : b;
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
