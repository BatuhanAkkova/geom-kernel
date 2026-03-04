#pragma once

#include "SDF.h"
#include "VoxelGrid.h"
#include <memory>

namespace Geom {

    /**
     * @brief Abstract Base Class representing a discretized domain.
     * 
     * Different solvers might require different domain representations
     * (e.g. dense grids, sparse trees, unstructured meshes). 
     * This base class provides a common interface that discretizers can return.
     */
    class DiscreteDomain {
    public:
        virtual ~DiscreteDomain() = default;
    };

    /**
     * @brief Discretized domain represented as a dense VoxelGrid.
     * 
     * Useful for solvers like Matrix-Free FEA or LBM that operate on fixed grids.
     */
    class StructuredDiscreteDomain : public DiscreteDomain {
    public:
        std::shared_ptr<VoxelGrid> grid;

        explicit StructuredDiscreteDomain(std::shared_ptr<VoxelGrid> g) : grid(std::move(g)) {}
    };

    /**
     * @brief Abstract interface for surface discretization methods.
     * 
     * Translates a continuous SDF into a specific `DiscreteDomain` representation.
     */
    class ISurfaceDiscretizer {
    public:
        virtual ~ISurfaceDiscretizer() = default;

        /**
         * @brief Discretize the given SDF.
         * @param sdf Continuous signed distance field to discretize.
         * @return A unique pointer to the generated DiscreteDomain.
         */
        virtual std::unique_ptr<DiscreteDomain> discretize(const SDF& sdf) const = 0;
    };

    /**
     * @brief A simple discretizer that creates a structured dense grid from an SDF.
     * 
     * @note Future updates should include more advanced discretizers such as:
     * - AdaptiveOctreeDiscretizer: Generate an octree for varying levels of detail, reducing memory footprint.
     * - CutFEM/XFEMDiscretizer: Create boundary-conforming background meshes for advanced finite element formulations.
     */
    class StructuredGridSampler : public ISurfaceDiscretizer {
    public:
        int nx, ny, nz;
        BoundingBox bounds;

        StructuredGridSampler(int resolution_x, int resolution_y, int resolution_z, const BoundingBox& domain_bounds)
            : nx(resolution_x), ny(resolution_y), nz(resolution_z), bounds(domain_bounds) {}

        std::unique_ptr<DiscreteDomain> discretize(const SDF& sdf) const override {
            auto grid = VoxelGrid::sampleSDF(sdf, bounds, nx, ny, nz);
            return std::make_unique<StructuredDiscreteDomain>(grid);
        }
    };

} // namespace Geom
