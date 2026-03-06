#pragma once

#include "SDF.h"
#include <vector>
#include <cmath>
#include <algorithm>

namespace Geom {

    /**
     * @brief Bridges a 3D density grid to the SDF interface.
     * Uses trilinear interpolation to sample densities at arbitrary points.
     */
    class DensityField : public SDFNode<DensityField> {
        const std::vector<double>& rho;
        int nx, ny, nz;
        BoundingBox domain;
        double threshold;

    public:
        DensityField(const std::vector<double>& densities, int _nx, int _ny, int _nz, 
                    const BoundingBox& bbox, double _threshold = 0.5)
            : rho(densities), nx(_nx), ny(_ny), nz(_nz), domain(bbox), threshold(_threshold) {}

        template <typename T>
        T evaluate(const Vec3T<T>& p) const {
            using std::max; using std::min; using std::abs; using std::sqrt; using std::pow; using std::sin; using std::cos;
            // Note: Since the grid is discrete, gradients/hessians of the grid itself 
            // usually requires finite differences or knowing the analytical derivative of trilinear interp.
            // For MVP, we will cast point to Scalar, perform trilinear eval, and return constant Dual.
            // A full implementation would propagate dual derivatives through the trilinear interpolation weights.
            
            // To be technically correct for AD:
            Vec3 size = domain.max - domain.min;
            double voxelSizeX = size.x / nx;
            double voxelSizeY = size.y / ny;
            double voxelSizeZ = size.z / nz;

            // Use .val chain to get raw doubles if T is Dual, but we can do it generically.
            // Actually, best to just write it with T to get full AD if we use standard operators!
            
            using std::max; using std::min; using std::abs;
            T gx = (p.x - static_cast<T>(domain.min.x)) / static_cast<T>(std::max(1e-9, voxelSizeX)) - static_cast<T>(0.5);
            T gy = (p.y - static_cast<T>(domain.min.y)) / static_cast<T>(std::max(1e-9, voxelSizeY)) - static_cast<T>(0.5);
            T gz = (p.z - static_cast<T>(domain.min.z)) / static_cast<T>(std::max(1e-9, voxelSizeZ)) - static_cast<T>(0.5);

            gx = max(static_cast<T>(0.0), min(static_cast<T>((double)nx - 1.00000001), gx));
            gy = max(static_cast<T>(0.0), min(static_cast<T>((double)ny - 1.00000001), gy));
            gz = max(static_cast<T>(0.0), min(static_cast<T>((double)nz - 1.00000001), gz));

            // Integer parts for grid indexing.
            // We need to extract the raw value to cast to int.
            auto getVal = [](const T& v) -> double {
                if constexpr (std::is_same_v<T, Scalar>) return v;
                else if constexpr (std::is_same_v<T, DualScalar>) return v.val;
                else if constexpr (std::is_same_v<T, Dual2Scalar>) return v.val.val;
                return 0;
            };

            int i0 = (int)getVal(gx), i1 = i0 + 1;
            int j0 = (int)getVal(gy), j1 = j0 + 1;
            int k0 = (int)getVal(gz), k1 = k0 + 1;

            T dx = gx - static_cast<T>(i0);
            T dy = gy - static_cast<T>(j0);
            T dz = gz - static_cast<T>(k0);

            auto sample = [&](int i, int j, int k) {
                return static_cast<T>(rho[(k * ny + j) * nx + i]);
            };

            T v000 = sample(i0, j0, k0);
            T v100 = sample(i1, j0, k0);
            T v010 = sample(i0, j1, k0);
            T v110 = sample(i1, j1, k0);
            T v001 = sample(i0, j0, k1);
            T v101 = sample(i1, j0, k1);
            T v011 = sample(i0, j1, k1);
            T v111 = sample(i1, j1, k1);

            T t1 = static_cast<T>(1.0);
            
            T v00 = v000 * (t1 - dx) + v100 * dx;
            T v10 = v010 * (t1 - dx) + v110 * dx;
            T v01 = v001 * (t1 - dx) + v101 * dx;
            T v11 = v011 * (t1 - dx) + v111 * dx;

            T v0 = v00 * (t1 - dy) + v10 * dy;
            T v1 = v01 * (t1 - dy) + v11 * dy;

            T density = v0 * (t1 - dz) + v1 * dz;

            return static_cast<T>(threshold) - density;
        }

        BoundingBox boundingBox() const override {
            return domain;
        }
    };

} // namespace Geom
