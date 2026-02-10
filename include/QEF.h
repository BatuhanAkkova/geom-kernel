#pragma once

#include "Geometry.h"
#include <cmath>
#include <limits>

namespace Geom {

    class QEF {
    public:
        // A (3x3 symmetric) = sum(n * n^T)
        // b (3 vector) = sum((n . p) * n)
        // massPoint = sum(p) / count
        
        double ata[6]; // Symmetric Matrix: 00, 01, 02, 11, 12, 22
        Vec3 atb;
        Point3 massPointSum;
        int pointCount;

        QEF() : atb(0,0,0), massPointSum(0,0,0), pointCount(0) {
            for(int i=0; i<6; ++i) ata[i] = 0.0;
        }

        void add(const Point3& p, const Vec3& n) {
            // A += n * n^T
            ata[0] += n.x * n.x;
            ata[1] += n.x * n.y;
            ata[2] += n.x * n.z;
            ata[3] += n.y * n.y;
            ata[4] += n.y * n.z;
            ata[5] += n.z * n.z;

            // b += (n . p) * n
            double dot = n.dot(p);
            atb = atb + n * dot;
            
            massPointSum = massPointSum + p;
            pointCount++;
        }

        Point3 solve() const {
            if (pointCount == 0) return Point3(0,0,0);
            
            Point3 mp = massPointSum * (1.0 / pointCount);
            
            // Equation: A x = b
            // But usually we minimize E(x) = sum( (n.(x-p))^2 )
            // Gradient is 2 * sum( n * (n.(x-p)) ) = 0
            // sum( n * (n.x - n.p) ) = 0
            // sum( (n*n^T) x ) = sum( n * (n.p) )
            // A x = b
            
            // Try inverting A.
            // A is symmetric.
            // Matrix M:
            // 0 1 2
            // 1 3 4
            // 2 4 5
            
            double m[3][3];
            m[0][0] = ata[0]; m[0][1] = ata[1]; m[0][2] = ata[2];
            m[1][0] = ata[1]; m[1][1] = ata[3]; m[1][2] = ata[4];
            m[2][0] = ata[2]; m[2][1] = ata[4]; m[2][2] = ata[5];
            
            // Determinant
            double det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
                         m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
                         m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
                         
            if (std::abs(det) < 1e-6) {
                // Singular, fallback to mass point
                return mp;
            }
            
            // Inverse
            double invDet = 1.0 / det;
            double inv[3][3];
            
            inv[0][0] = (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * invDet;
            inv[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * invDet;
            inv[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * invDet;

            inv[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * invDet;
            inv[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * invDet;
            inv[1][2] = (m[1][0] * m[0][2] - m[0][0] * m[1][2]) * invDet;

            inv[2][0] = (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * invDet;
            inv[2][1] = (m[1][1] * m[0][2] - m[0][1] * m[2][1]) * invDet; // Error in minor sign?
             // Standard formula:
             // 11 = +(00*22 - 02*20) -> (m[0][0]m[2][2] - m[0][2]m[2][0]) Correct.
             // 12 = -(00*21 - 01*20) -> -(m[0][0]m[2][1] - m[0][1]m[2][0]) = (01*20 - 00*21).
             // My line above: (m10 m02 - m00 m12)... indices are mixed.
             
             // Let's rewrite carefully.
             // Adjugate matrix
             // C00 = +(m11m22 - m12m21)
             // C01 = -(m10m22 - m12m20) = m12m20 - m10m22
             // C02 = +(m10m21 - m11m20)
             
             // C10 = -(m01m22 - m02m21) = m02m21 - m01m22
             // C11 = +(m00m22 - m02m20)
             // C12 = -(m00m21 - m01m20) = m01m20 - m00m21
             
             // C20 = +(m01m12 - m02m11)
             // C21 = -(m00m12 - m02m10) = m02m10 - m00m12
             // C22 = +(m00m11 - m01m10)
             
             // Inv = Adj^T / Det. Since symmetric, Adj is symmetric.
             
             inv[0][0] = (m[1][1]*m[2][2] - m[1][2]*m[2][1]) * invDet;
             inv[0][1] = (m[0][2]*m[2][1] - m[0][1]*m[2][2]) * invDet;
             inv[0][2] = (m[0][1]*m[1][2] - m[0][2]*m[1][1]) * invDet; // symmetric with 20? 
             // C02 = m10m21 - m11m20.
             // m01m12 - m02m11. Same for symmetric.
             
             inv[1][1] = (m[0][0]*m[2][2] - m[0][2]*m[2][0]) * invDet;
             inv[1][2] = (m[0][2]*m[1][0] - m[0][0]*m[1][2]) * invDet;
             
             inv[2][2] = (m[0][0]*m[1][1] - m[0][1]*m[1][0]) * invDet;
             
             inv[1][0] = inv[0][1];
             inv[2][0] = inv[0][2];
             inv[2][1] = inv[1][2];

             // Result
             double x = inv[0][0]*atb.x + inv[0][1]*atb.y + inv[0][2]*atb.z;
             double y = inv[1][0]*atb.x + inv[1][1]*atb.y + inv[1][2]*atb.z;
             double z = inv[2][0]*atb.x + inv[2][1]*atb.y + inv[2][2]*atb.z;
             
             return Point3(x, y, z);
        }
    };
}
