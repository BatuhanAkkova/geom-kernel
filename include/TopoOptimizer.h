#pragma once

#include "Geometry.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>

namespace Geom {

    /**
     * @brief SIMP (Solid Isotropic Material with Penalization) Topology Optimizer.
     * Operates on a structured 3D grid of densities.
     */
    class TopoOptimizer {
    public:
        struct Config {
            int nx, ny, nz;
            double volfrac;
            double penal;
            double rmin; // filter radius in voxels
        };

    protected:
        Config cfg;
        std::vector<double> x;      // densities [0, 1]
        std::vector<double> dc;     // sensitivities (dC/dx)
        std::vector<double> xnew;   // buffer for OC update

    public:
        TopoOptimizer(Config c) : cfg(c) {
            int total = cfg.nx * cfg.ny * cfg.nz;
            x.assign(total, cfg.volfrac);
            dc.assign(total, 0.0);
            xnew.assign(total, 0.0);
        }

        // Accessors
        const std::vector<double>& densities() const { return x; }
        int nx() const { return cfg.nx; }
        int ny() const { return cfg.ny; }
        int nz() const { return cfg.nz; }

        /**
         * @brief Mock compliance sensitivity calculation.
         * In a real system, this would come from a FEA solver (U^T * K * U).
         * For the MVP, we use a load field: sensitivity is high near load points.
         */
        void computeSensitivities(const std::vector<double>& strainEnergyDensity) {
            int total = cfg.nx * cfg.ny * cfg.nz;
            for (int i = 0; i < total; ++i) {
                // SIMP sensitivity: dC/dx = -p * x^(p-1) * energy
                dc[i] = -cfg.penal * std::pow(x[i], cfg.penal - 1.0) * strainEnergyDensity[i];
            }
        }

        /**
         * @brief Mesh-independent filtering of sensitivities (Hat filter).
         */
        void filterSensitivities() {
            int total = cfg.nx * cfg.ny * cfg.nz;
            std::vector<double> dcf(total, 0.0);
            
            for (int k = 0; k < cfg.nz; ++k) {
                for (int j = 0; j < cfg.ny; ++j) {
                    for (int i = 0; i < cfg.nx; ++i) {
                        double sum_w = 0.0;
                        int idx = (k * cfg.ny + j) * cfg.nx + i;
                        
                        // Local neighborhood search
                        int r_int = (int)std::ceil(cfg.rmin);
                        for (int nk = std::max(0, k-r_int); nk <= std::min(cfg.nz-1, k+r_int); ++nk) {
                            for (int nj = std::max(0, j-r_int); nj <= std::min(cfg.ny-1, j+r_int); ++nj) {
                                for (int ni = std::max(0, i-r_int); ni <= std::min(cfg.nx-1, i+r_int); ++ni) {
                                    double d = std::sqrt(std::pow(i-ni, 2) + std::pow(j-nj, 2) + std::pow(k-nk, 2));
                                    double w = std::max(0.0, cfg.rmin - d);
                                    if (w > 0) {
                                        int nidx = (nk * cfg.ny + nj) * cfg.nx + ni;
                                        dcf[idx] += w * x[nidx] * dc[nidx];
                                        sum_w += w;
                                    }
                                }
                            }
                        }
                        dcf[idx] /= (std::max(1e-9, x[idx] * sum_w));
                    }
                }
            }
            dc = dcf;
        }

        /**
         * @brief Optimality Criteria (OC) update for densities.
         */
        double updateDensities() {
            double l1 = 0, l2 = 1e9, move = 0.2;
            int total = cfg.nx * cfg.ny * cfg.nz;
            double max_change = 0;

            while ((l2 - l1) / (l1 + l2) > 1e-4) {
                double lmid = 0.5 * (l2 + l1);
                for (int i = 0; i < total; ++i) {
                    // OC heuristic: x_new = x * sqrt(-dc / lambda)
                    // with move limits and 0..1 bounds
                    double target = x[i] * std::sqrt(-dc[i] / std::max(1e-12, lmid));
                    xnew[i] = std::max(0.0, std::max(x[i] - move, std::min(1.0, std::min(x[i] + move, target))));
                }
                
                double current_vol = std::accumulate(xnew.begin(), xnew.end(), 0.0) / total;
                if (current_vol > cfg.volfrac) l1 = lmid;
                else l2 = lmid;
            }

            for (int i = 0; i < total; ++i) {
                double change = std::abs(xnew[i] - x[i]);
                if (change > max_change) max_change = change;
                x[i] = xnew[i];
            }

            return max_change;
        }
    };

} // namespace Geom
