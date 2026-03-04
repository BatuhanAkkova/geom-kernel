#pragma once

#include "SDF.h"
#include <vector>
#include <unordered_map>
#include <future>
#include <cmath>

namespace Geom {

    /**
     * @brief Dense 3D Voxel Grid for SDF storage and sampling.
     * 
     * This class provides a compact 3D array for storing pre-sampled SDF values,
     * which can be used for fast lookups or transfer to GPU textures.
     */
    class VoxelGrid {
    public:
        std::vector<float> data;
        int nx, ny, nz;
        BoundingBox bounds;

        VoxelGrid(int x, int y, int z, const BoundingBox& b) 
            : nx(x), ny(y), nz(z), bounds(b) {
            data.resize(nx * ny * nz, 0.0f);
        }

        float& at(int i, int j, int k) {
            return data[(k * ny + j) * nx + i];
        }

        const float& at(int i, int j, int k) const {
            return data[(k * ny + j) * nx + i];
        }

        Point3 indexToWorld(int i, int j, int k) const {
            Vec3 size = bounds.size();
            return bounds.min + Vec3(
                (i / (Scalar)(nx - 1)) * size.x,
                (j / (Scalar)(ny - 1)) * size.y,
                (k / (Scalar)(nz - 1)) * size.z
            );
        }

        /**
         * @brief Parallel SDF sampling into the grid.
         * @param sdf The source SDF to sample.
         * @param b The bounding box in world space to cover.
         * @param nx, ny, nz Grid resolution.
         * @return A shared pointer to the populated VoxelGrid.
         */
        static std::shared_ptr<VoxelGrid> sampleSDF(const SDF& sdf, const BoundingBox& b, int nx, int ny, int nz) {
            auto grid = std::make_shared<VoxelGrid>(nx, ny, nz, b);
            
            int num_threads = 8;
            int chunk_size = nz / num_threads;
            std::vector<std::future<void>> futures;

            for (int t = 0; t < num_threads; ++t) {
                int k_start = t * chunk_size;
                int k_end = (t == num_threads - 1) ? nz : (t + 1) * chunk_size;

                futures.push_back(std::async(std::launch::async, [&sdf, grid, k_start, k_end, nx, ny]() {
                    for (int k = k_start; k < k_end; ++k) {
                        for (int j = 0; j < ny; ++j) {
                            for (int i = 0; i < nx; ++i) {
                                grid->at(i, j, k) = (float)sdf.eval(grid->indexToWorld(i, j, k));
                            }
                        }
                    }
                }));
            }

            for (auto& f : futures) f.wait();
            return grid;
        }
    };

    /**
     * @brief Sparse Voxel Grid using a hash map (VDB-lite).
     * Keys are bit-packed i, j, k.
     */
    class SparseVoxelGrid {
    public:
        std::unordered_map<uint64_t, float> voxels;
        int nx, ny, nz;
        BoundingBox bounds;

        /**
         * @brief Constructs a SparseVoxelGrid with specified dimensions and bounding box.
         * @param x Resolution along the X-axis.
         * @param y Resolution along the Y-axis.
         * @param z Resolution along the Z-axis.
         * @param b The world-space bounding box covered by the grid.
         */
        SparseVoxelGrid(int x, int y, int z, const BoundingBox& b) 
            : nx(x), ny(y), nz(z), bounds(b) {}

        uint64_t packKey(int i, int j, int k) const {
            return ((uint64_t)i & 0xFFFFF) | (((uint64_t)j & 0xFFFFF) << 20) | (((uint64_t)k & 0xFFFFF) << 40);
        }

        void set(int i, int j, int k, float val) {
            voxels[packKey(i, j, k)] = val;
        }

        float get(int i, int j, int k) const {
            auto it = voxels.find(packKey(i, j, k));
            return (it != voxels.end()) ? it->second : 1e10f; // Return far away if empty
        }

        /**
         * @brief Only samples points within a distance bandwidth of the surface.
         */
        void populateNarrowBand(const SDF& sdf, float bandwidthVoxels) {
            Vec3 size = bounds.size();
            Scalar rx = size.x / (nx - 1);
            Scalar ry = size.y / (ny - 1);
            Scalar rz = size.z / (nz - 1);
            Scalar voxelSize = std::max({rx, ry, rz});
            Scalar limit = bandwidthVoxels * voxelSize;

            for (int k = 0; k < nz; ++k) {
                for (int j = 0; j < ny; ++j) {
                    for (int i = 0; i < nx; ++i) {
                        Point3 p = bounds.min + Vec3(i*rx, j*ry, k*rz);
                        Scalar dist = sdf.eval(p);
                        if (std::abs(dist) <= limit) {
                            set(i, j, k, (float)dist);
                        }
                    }
                }
            }
        }
    };

    /**
     * @brief Partial Volume Calculation.
     * Determines how much of a voxel is "solid" (volume fraction).
     * Uses a smooth Heaviside function (Sigmoid).
     */
    inline Scalar partialVolume(Scalar sdf_val, Scalar voxel_size) {
        // Simple sigmoid approximation
        // At 0, returns 0.5. At -voxel_size, return ~1.0. At +voxel_size, return ~0.0.
        // k controls sharpness. 
        Scalar k = 4.0 / voxel_size; 
        return 1.0 / (1.0 + std::exp(sdf_val * k));
    }

}
