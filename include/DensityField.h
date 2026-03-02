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
    class DensityField : public SDF {
        const std::vector<double>& rho;
        int nx, ny, nz;
        BoundingBox domain;
        double threshold;

    public:
        DensityField(const std::vector<double>& densities, int _nx, int _ny, int _nz, 
                    const BoundingBox& bbox, double _threshold = 0.5)
            : rho(densities), nx(_nx), ny(_ny), nz(_nz), domain(bbox), threshold(_threshold) {}

        /**
         * @brief Evaluate the density field at point p.
         * Returns (threshold - rho) so that positive rho > threshold is INSIDE (negative SDF).
         */
        Scalar eval(const Point3& p) const override {
            // Map p from domain to [0, nx-1] x [0, ny-1] x [0, nz-1]
            Vec3 size = domain.max - domain.min;
            double gx = ((p.x - domain.min.x) / std::max(1e-9, size.x)) * (nx - 1);
            double gy = ((p.y - domain.min.y) / std::max(1e-9, size.y)) * (ny - 1);
            double gz = ((p.z - domain.min.z) / std::max(1e-9, size.z)) * (nz - 1);

            // Clamp to grid bounds with very small epsilon
            gx = std::max(0.0, std::min((double)nx - 1.00000001, gx));
            gy = std::max(0.0, std::min((double)ny - 1.00000001, gy));
            gz = std::max(0.0, std::min((double)nz - 1.00000001, gz));

            int i0 = (int)gx, i1 = i0 + 1;
            int j0 = (int)gy, j1 = j0 + 1;
            int k0 = (int)gz, k1 = k0 + 1;

            double dx = gx - i0;
            double dy = gy - j0;
            double dz = gz - k0;

            auto sample = [&](int i, int j, int k) {
                return rho[(k * ny + j) * nx + i];
            };

            // Trilinear interpolation
            double v000 = sample(i0, j0, k0);
            double v100 = sample(i1, j0, k0);
            double v010 = sample(i0, j1, k0);
            double v110 = sample(i1, j1, k0);
            double v001 = sample(i0, j0, k1);
            double v101 = sample(i1, j0, k1);
            double v011 = sample(i0, j1, k1);
            double v111 = sample(i1, j1, k1);

            double v00 = v000 * (1 - dx) + v100 * dx;
            double v10 = v010 * (1 - dx) + v110 * dx;
            double v01 = v001 * (1 - dx) + v101 * dx;
            double v11 = v011 * (1 - dx) + v111 * dx;

            double v0 = v00 * (1 - dy) + v10 * dy;
            double v1 = v01 * (1 - dy) + v11 * dy;

            double density = v0 * (1 - dz) + v1 * dz;

            return threshold - density;
        }

        BoundingBox boundingBox() const override {
            return domain;
        }
    };

} // namespace Geom
