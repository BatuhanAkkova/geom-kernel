#pragma once

#include "SDF.h"
#include "Geometry.h"
#include <algorithm>
#include <cmath>

namespace Geom {

    /**
     * @brief Transforms an SDF using a 4x4 matrix.
     * Stores the inverse transform to map world points to local space.
     * To maintain SDF bounds under non-uniform scale, we divide the local 
     * distance by the maximum singular value (or maximum scale factor) 
     * in the inverse transformation.
     */
    class Transform : public SDFNode<Transform> {
        SDFPtr sdf;
        Mat4 invMat;
        Scalar scaleFactor;

    public:
        Transform(SDFPtr s, const Mat4& m) : sdf(s) {
            invMat = m.inverse();
            
            // To guarantee the distance bound for an SDF after a transform (especially non-uniform scale),
            // the distance in world space is d_local / max_scale(invMat).
            // This is a lower bound, making it a valid sign-correct approximation.
            // max_scale(invMat) is the max singular value of the 3x3 rotational/scale part.
            // A simple conservative estimate is the maximum length of the mapped basis vectors.
            
            Vec3 vx(invMat.m[0], invMat.m[1], invMat.m[2]);
            Vec3 vy(invMat.m[4], invMat.m[5], invMat.m[6]);
            Vec3 vz(invMat.m[8], invMat.m[9], invMat.m[10]);
            
            scaleFactor = std::max({vx.length(), vy.length(), vz.length()});
            
            // Note: If scaleFactor is 0, the matrix is singular. We protect against div by zero.
            if(scaleFactor < EPSILON) scaleFactor = 1.0;
        }

        template <typename T>
        T evaluate(const Vec3T<T>& p) const {
            using std::max; using std::min; using std::abs; using std::sqrt; using std::pow; using std::sin; using std::cos;
            
            T lx = p.x * static_cast<T>(invMat.m[0]) + p.y * static_cast<T>(invMat.m[4]) + p.z * static_cast<T>(invMat.m[8])  + static_cast<T>(invMat.m[12]);
            T ly = p.x * static_cast<T>(invMat.m[1]) + p.y * static_cast<T>(invMat.m[5]) + p.z * static_cast<T>(invMat.m[9])  + static_cast<T>(invMat.m[13]);
            T lz = p.x * static_cast<T>(invMat.m[2]) + p.y * static_cast<T>(invMat.m[6]) + p.z * static_cast<T>(invMat.m[10]) + static_cast<T>(invMat.m[14]);

            if constexpr (std::is_same_v<T, Scalar>) {
                Point3 pt(lx, ly, lz);
                return sdf->eval(pt) / static_cast<T>(scaleFactor);
            } else if constexpr (std::is_same_v<T, DualScalar>) {
                Point3D pt(lx, ly, lz);
                return sdf->evalD(pt) / static_cast<T>(scaleFactor);
            } else if constexpr (std::is_same_v<T, Dual2Scalar>) {
                Point3D2 pt(lx, ly, lz);
                return sdf->evalD2(pt) / static_cast<T>(scaleFactor);
            }
            return static_cast<T>(0);
        }

        BoundingBox boundingBox() const override {
            BoundingBox local = sdf->boundingBox();
            Point3 corners[8];
            Point3 min = local.min;
            Point3 max = local.max;
            corners[0] = Point3(min.x, min.y, min.z);
            corners[1] = Point3(max.x, min.y, min.z);
            corners[2] = Point3(min.x, max.y, min.z);
            corners[3] = Point3(max.x, max.y, min.z);
            corners[4] = Point3(min.x, min.y, max.z);
            corners[5] = Point3(max.x, min.y, max.z);
            corners[6] = Point3(min.x, max.y, max.z);
            corners[7] = Point3(max.x, max.y, max.z);

            Mat4 forward = invMat.inverse(); 
            BoundingBox world;
            for (int i = 0; i < 8; ++i) {
                world.expand(forward.transformPoint(corners[i]));
            }
            return world;
        }
    };

    // Factory Helpers
    inline SDFPtr translate(SDFPtr sdf, const Vec3& v) {
        return std::make_shared<Transform>(sdf, Mat4::translate(v));
    }

    inline SDFPtr rotateX(SDFPtr sdf, Scalar angle) {
        return std::make_shared<Transform>(sdf, Mat4::rotateX(angle));
    }

    inline SDFPtr rotateY(SDFPtr sdf, Scalar angle) {
        return std::make_shared<Transform>(sdf, Mat4::rotateY(angle));
    }

    inline SDFPtr rotateZ(SDFPtr sdf, Scalar angle) {
        return std::make_shared<Transform>(sdf, Mat4::rotateZ(angle));
    }

    inline SDFPtr scale(SDFPtr sdf, const Vec3& s) {
        return std::make_shared<Transform>(sdf, Mat4::scale(s));
    }

    inline SDFPtr uniformScale(SDFPtr sdf, Scalar s) {
        return std::make_shared<Transform>(sdf, Mat4::scale(Vec3(s, s, s)));
    }

}
