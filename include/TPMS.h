#pragma once

#include "SDF.h"
#include "Autodiff.h"
#include <cmath>

namespace Geom {

    /**
     * @brief Base class for Triply Periodic Minimal Surfaces (TPMS).
     * 
     * TPMS are used to create high-strength, lightweight lattice microstructures.
     */
    class TPMS : public SDF {
    public:
        Scalar period;
        Scalar isovalue;

        TPMS(Scalar p, Scalar iso = 0) : period(p), isovalue(iso) {}

        BoundingBox boundingBox() const override {
            // TPMS is theoretically infinite. return an infinite BB or a very large one.
            return BoundingBox(Point3(-INF, -INF, -INF), Point3(INF, INF, INF));
        }

    protected:
        Scalar k() const { return 2.0 * M_PI / period; }
    };

    /**
     * @brief Gyroid surface: sin(kx)cos(ky) + sin(ky)cos(kz) + sin(kz)cos(kx) = iso
     */
    class Gyroid : public TPMS {
    public:
        using TPMS::TPMS;

        Scalar eval(const Point3& p) const override {
            Scalar s = k();
            return (std::sin(s * p.x) * std::cos(s * p.y) + 
                    std::sin(s * p.y) * std::cos(s * p.z) + 
                    std::sin(s * p.z) * std::cos(s * p.x)) - isovalue;
        }

        DualScalar evalD(const Point3D& p) const override {
            Scalar s = k();
            return (sin(p.x * s) * cos(p.y * s) + 
                    sin(p.y * s) * cos(p.z * s) + 
                    sin(p.z * s) * cos(p.x * s)) - isovalue;
        }

        Dual2Scalar evalD2(const Point3D2& p) const override {
            Scalar s = k();
            return (sin(p.x * s) * cos(p.y * s) + 
                    sin(p.y * s) * cos(p.z * s) + 
                    sin(p.z * s) * cos(p.x * s)) - isovalue;
        }
    };

    /**
     * @brief Diamond TPMS: sin(kx)sin(ky)sin(kz) + sin(kx)cos(ky)cos(kz) + cos(kx)sin(ky)cos(kz) + cos(kx)cos(ky)sin(kz) = iso
     */
    class DiamondTPMS : public TPMS {
    public:
        using TPMS::TPMS;

        Scalar eval(const Point3& p) const override {
            Scalar s = k();
            Scalar sx = std::sin(s * p.x), sy = std::sin(s * p.y), sz = std::sin(s * p.z);
            Scalar cx = std::cos(s * p.x), cy = std::cos(s * p.y), cz = std::cos(s * p.z);
            return (sx*sy*sz + sx*cy*cz + cx*sy*cz + cx*cy*sz) - isovalue;
        }

        DualScalar evalD(const Point3D& p) const override {
            Scalar s = k();
            DualScalar sx = sin(p.x * s), sy = sin(p.y * s), sz = sin(p.z * s);
            DualScalar cx = cos(p.x * s), cy = cos(p.y * s), cz = cos(p.z * s);
            return (sx*sy*sz + sx*cy*cz + cx*sy*cz + cx*cy*sz) - isovalue;
        }

        Dual2Scalar evalD2(const Point3D2& p) const override {
            Scalar s = k();
            Dual2Scalar sx = sin(p.x * s), sy = sin(p.y * s), sz = sin(p.z * s);
            Dual2Scalar cx = cos(p.x * s), cy = cos(p.y * s), cz = cos(p.z * s);
            return (sx*sy*sz + sx*cy*cz + cx*sy*cz + cx*cy*sz) - isovalue;
        }
    };

    /**
     * @brief Schwarz P surface: cos(kx) + cos(ky) + cos(kz) = iso
     */
    class SchwarzP : public TPMS {
    public:
        using TPMS::TPMS;

        Scalar eval(const Point3& p) const override {
            Scalar s = k();
            return (std::cos(s * p.x) + std::cos(s * p.y) + std::cos(s * p.z)) - isovalue;
        }

        DualScalar evalD(const Point3D& p) const override {
            Scalar s = k();
            return (cos(p.x * s) + cos(p.y * s) + cos(p.z * s)) - isovalue;
        }

        Dual2Scalar evalD2(const Point3D2& p) const override {
            Scalar s = k();
            return (cos(p.x * s) + cos(p.y * s) + cos(p.z * s)) - isovalue;
        }
    };

}
