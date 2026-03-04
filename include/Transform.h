#pragma once

#include "SDF.h"
#include "Geometry.h"

namespace Geom {

    /**
     * @brief Transforms an SDF using a 4x4 matrix.
     * Stores the inverse transform to map world points to local space.
     */
    class Transform : public SDF {
        SDFPtr sdf;
        Mat4 invMat;
        Scalar scaleFactor;

    public:
        Transform(SDFPtr s, const Mat4& m) : sdf(s) {
            invMat = m.inverse();
            // Determine average scale from the matrix (rough estimate for distance correction)
            // For uniform scale, it's exact.
            Vec3 v(1, 0, 0);
            Point3 p = m.transformPoint(Point3(0,0,0));
            Point3 p1 = m.transformPoint(Point3(1,0,0));
            scaleFactor = (p1 - p).length();
        }

        Scalar eval(const Point3& p) const override {
            Point3 local = invMat.transformPoint(p);
            return sdf->eval(local) * scaleFactor;
        }

        DualScalar evalD(const Point3D& p) const override {
            // Manual transformation of Dual coordinates
            DualScalar lx = p.x * invMat.m[0] + p.y * invMat.m[4] + p.z * invMat.m[8]  + invMat.m[12];
            DualScalar ly = p.x * invMat.m[1] + p.y * invMat.m[5] + p.z * invMat.m[9]  + invMat.m[13];
            DualScalar lz = p.x * invMat.m[2] + p.y * invMat.m[6] + p.z * invMat.m[10] + invMat.m[14];
            
            return sdf->evalD(Point3D(lx, ly, lz)) * scaleFactor;
        }

        Dual2Scalar evalD2(const Point3D2& p) const override {
            // Manual transformation of Dual2 coordinates
            Dual2Scalar lx = p.x * invMat.m[0] + p.y * invMat.m[4] + p.z * invMat.m[8]  + invMat.m[12];
            Dual2Scalar ly = p.x * invMat.m[1] + p.y * invMat.m[5] + p.z * invMat.m[9]  + invMat.m[13];
            Dual2Scalar lz = p.x * invMat.m[2] + p.y * invMat.m[6] + p.z * invMat.m[10] + invMat.m[14];

            return sdf->evalD2(Point3D2(lx, ly, lz)) * scaleFactor;
        }

        BoundingBox boundingBox() const override {
            BoundingBox local = sdf->boundingBox();
            // Transform 8 corners and find new min/max
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

            // Need to transform back to world (not using invMat)
            // Let's store the forward matrix or just invert invMat (which might be costly)
            // Better to pass both or re-invert.
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
