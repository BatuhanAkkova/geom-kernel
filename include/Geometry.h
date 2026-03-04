#pragma once

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include <cmath>
#include <limits>
#include <concepts>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

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

        Vec3& operator+=(const Vec3& other) { x += other.x; y += other.y; z += other.z; return *this; }
        Vec3& operator-=(const Vec3& other) { x -= other.x; y -= other.y; z -= other.z; return *this; }
        Vec3& operator*=(Scalar s) { x *= s; y *= s; z *= s; return *this; }
        Vec3& operator/=(Scalar s) { x /= s; y /= s; z /= s; return *this; }

        Scalar dot(const Vec3& other) const { return x * other.x + y * other.y + z * other.z; }
        
        Scalar lengthSquared() const { return x * x + y * y + z * z; }
        Scalar length() const { return std::sqrt(lengthSquared()); }

        Vec3 normalized() const {
            Scalar len = length();
            return (len > EPSILON) ? *this / len : Vec3(0, 0, 0);
        }

        Scalar operator[](int i) const {
            if (i == 0) return x;
            if (i == 1) return y;
            return z;
        }

        Scalar& operator[](int i) {
            if (i == 0) return x;
            if (i == 1) return y;
            return z;
        }
    };

    /**
     * @brief 3x3 Matrix for Hessians and small transforms.
     */
    struct Mat3 {
        Scalar m[9];

        Mat3() {
            for (int i = 0; i < 9; ++i) m[i] = 0;
        }

        static Mat3 identity() {
            Mat3 mat;
            mat.m[0] = mat.m[4] = mat.m[8] = 1;
            return mat;
        }

        Scalar& operator()(int r, int c) { return m[c * 3 + r]; }
        Scalar operator()(int r, int c) const { return m[c * 3 + r]; }

        Mat3 operator+(const Mat3& other) const {
            Mat3 res;
            for (int i = 0; i < 9; ++i) res.m[i] = m[i] + other.m[i];
            return res;
        }

        Mat3 operator*(Scalar s) const {
            Mat3 res;
            for (int i = 0; i < 9; ++i) res.m[i] = m[i] * s;
            return res;
        }
    };

    using Point3 = Vec3;

    inline Vec3 min(const Vec3& a, const Vec3& b) {
        return {std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z)};
    }

    inline Vec3 max(const Vec3& a, const Vec3& b) {
        return {std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z)};
    }

    inline Vec3 operator*(Scalar s, const Vec3& v) {
        return v * s;
    }

    inline Vec3 cross(const Vec3& a, const Vec3& b) {
        return {
            a.y * b.z - a.z * b.y,
            a.z * b.x - a.x * b.z,
            a.x * b.y - a.y * b.x
        };
    }

    /**
     * @brief 4x4 Matrix for transformations.
     * Column-major order.
     */
    struct Mat4 {
        Scalar m[16];

        Mat4() {
            for (int i = 0; i < 16; ++i) m[i] = 0;
            m[0] = m[5] = m[10] = m[15] = 1; // Identity
        }

        static Mat4 translate(const Vec3& v) {
            Mat4 mat;
            mat.m[12] = v.x;
            mat.m[13] = v.y;
            mat.m[14] = v.z;
            return mat;
        }

        static Mat4 rotateX(Scalar angle) {
            Mat4 mat;
            Scalar c = std::cos(angle);
            Scalar s = std::sin(angle);
            mat.m[5] = c;  mat.m[6] = s;
            mat.m[9] = -s; mat.m[10] = c;
            return mat;
        }

        static Mat4 rotateY(Scalar angle) {
            Mat4 mat;
            Scalar c = std::cos(angle);
            Scalar s = std::sin(angle);
            mat.m[0] = c;  mat.m[2] = -s;
            mat.m[8] = s;  mat.m[10] = c;
            return mat;
        }

        static Mat4 rotateZ(Scalar angle) {
            Mat4 mat;
            Scalar c = std::cos(angle);
            Scalar s = std::sin(angle);
            mat.m[0] = c;  mat.m[1] = s;
            mat.m[4] = -s; mat.m[5] = c;
            return mat;
        }

        static Mat4 scale(const Vec3& s) {
            Mat4 mat;
            mat.m[0] = s.x;
            mat.m[5] = s.y;
            mat.m[10] = s.z;
            return mat;
        }

        Point3 transformPoint(const Point3& p) const {
            return {
                p.x * m[0] + p.y * m[4] + p.z * m[8]  + m[12],
                p.x * m[1] + p.y * m[5] + p.z * m[9]  + m[13],
                p.x * m[2] + p.y * m[6] + p.z * m[10] + m[14]
            };
        }

        Mat4 operator*(const Mat4& b) const {
            Mat4 out;
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    out.m[j * 4 + i] = 0;
                    for (int k = 0; k < 4; k++) {
                        out.m[j * 4 + i] += m[k * 4 + i] * b.m[j * 4 + k];
                    }
                }
            }
            return out;
        }

        // Simple inverse for orthogonal or simple TRS matrices
        // For full inverse, we use a standard cofactor-based approach
        Mat4 inverse() const {
            // Simplified version for the kernel (usually TRS matrices)
            // Let's implement a robust version
            Scalar inv[16];
            Scalar det;

            inv[0] = m[5]  * m[10] * m[15] - 
                     m[5]  * m[11] * m[14] - 
                     m[9]  * m[6]  * m[15] + 
                     m[9]  * m[7]  * m[14] +
                     m[13] * m[6]  * m[11] - 
                     m[13] * m[7]  * m[10];

            inv[4] = -m[4]  * m[10] * m[15] + 
                      m[4]  * m[11] * m[14] + 
                      m[8]  * m[6]  * m[15] - 
                      m[8]  * m[7]  * m[14] - 
                      m[12] * m[6]  * m[11] + 
                      m[12] * m[7]  * m[10];

            inv[8] = m[4]  * m[9] * m[15] - 
                     m[4]  * m[11] * m[13] - 
                     m[8]  * m[5] * m[15] + 
                     m[8]  * m[7] * m[13] + 
                     m[12] * m[5] * m[11] - 
                     m[12] * m[7] * m[9];

            inv[12] = -m[4]  * m[9] * m[14] + 
                       m[4]  * m[10] * m[13] +
                       m[8]  * m[5] * m[14] - 
                       m[8]  * m[6] * m[13] - 
                       m[12] * m[5] * m[10] + 
                       m[12] * m[6] * m[9];

            inv[1] = -m[1]  * m[10] * m[15] + 
                      m[1]  * m[11] * m[14] + 
                      m[9]  * m[2] * m[15] - 
                      m[9]  * m[3] * m[14] - 
                      m[13] * m[2] * m[11] + 
                      m[13] * m[3] * m[10];

            inv[5] = m[0]  * m[10] * m[15] - 
                     m[0]  * m[11] * m[14] - 
                     m[8]  * m[2] * m[15] + 
                     m[8]  * m[3] * m[14] + 
                     m[12] * m[2] * m[11] - 
                     m[12] * m[3] * m[10];

            inv[9] = -m[0]  * m[9] * m[15] + 
                      m[0]  * m[11] * m[13] + 
                      m[8]  * m[1] * m[15] - 
                      m[8]  * m[3] * m[13] - 
                      m[12] * m[1] * m[11] + 
                      m[12] * m[3] * m[9];

            inv[13] = m[0]  * m[9] * m[14] - 
                      m[0]  * m[10] * m[13] - 
                      m[8]  * m[1] * m[14] + 
                      m[8]  * m[2] * m[13] + 
                      m[12] * m[1] * m[10] - 
                      m[12] * m[2] * m[9];

            inv[2] = m[1]  * m[6] * m[15] - 
                     m[1]  * m[7] * m[14] - 
                     m[5]  * m[2] * m[15] + 
                     m[5]  * m[3] * m[14] + 
                     m[13] * m[2] * m[7] - 
                     m[13] * m[3] * m[6];

            inv[6] = -m[0]  * m[6] * m[15] + 
                      m[0]  * m[7] * m[14] + 
                      m[4]  * m[2] * m[15] - 
                      m[4]  * m[3] * m[14] - 
                      m[12] * m[2] * m[7] + 
                      m[12] * m[3] * m[6];

            inv[10] = m[0]  * m[5] * m[15] - 
                      m[0]  * m[7] * m[13] - 
                      m[4]  * m[1] * m[15] + 
                      m[4]  * m[3] * m[13] + 
                      m[12] * m[1] * m[7] - 
                      m[12] * m[3] * m[5];

            inv[14] = -m[0]  * m[5] * m[14] + 
                       m[0]  * m[6] * m[13] + 
                       m[4]  * m[1] * m[14] - 
                       m[4]  * m[2] * m[13] - 
                       m[12] * m[1] * m[6] + 
                       m[12] * m[2] * m[5];

            inv[3] = -m[1] * m[6] * m[11] + 
                      m[1] * m[7] * m[10] + 
                      m[5] * m[2] * m[11] - 
                      m[5] * m[3] * m[10] - 
                      m[9] * m[2] * m[7] + 
                      m[9] * m[3] * m[6];

            inv[7] = m[0] * m[6] * m[11] - 
                     m[0] * m[7] * m[10] - 
                     m[4] * m[2] * m[11] + 
                     m[4] * m[3] * m[10] + 
                     m[8] * m[2] * m[7] - 
                     m[8] * m[3] * m[6];

            inv[11] = -m[0] * m[5] * m[11] + 
                       m[0] * m[7] * m[9] + 
                       m[4] * m[1] * m[11] - 
                       m[4] * m[3] * m[9] - 
                       m[8] * m[1] * m[7] + 
                       m[8] * m[3] * m[5];

            inv[15] = m[0] * m[5] * m[10] - 
                      m[0] * m[6] * m[9] - 
                      m[4] * m[1] * m[10] + 
                      m[4] * m[2] * m[9] + 
                      m[8] * m[1] * m[6] - 
                      m[8] * m[2] * m[5];

            det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];

            Mat4 res;
            if (std::abs(det) < 1e-12) return res;

            Scalar invDet = 1.0 / det;
            for (int i = 0; i < 16; i++) res.m[i] = inv[i] * invDet;

            return res;
        }
    };
}
