#pragma once

#include "SDF.h"
#include <vector>
#include <functional>
#include <thread>
#include <future>
#include <algorithm>
#include <execution>

namespace Geom {

    /**
     * @brief Result of a design sweep.
     */
    struct OptimizationResult {
        Scalar parameterValue;
        Scalar objectiveValue;
        Scalar sensitivity;
    };

    /**
     * @brief Utilities for batch geometry generation and optimization.
     */
    class Optimization {
    public:
        /**
         * @brief Perform a parallel parameter sweep over an SDF-generating function.
         * @param paramRange Vector of parameter values to test.
         * @param generator Function that takes a parameter and returns an SDF.
         * @param objective Function that evaluates an SDF and returns a scalar objective.
         * @return Vector of OptimizationResults.
         */
        static std::vector<OptimizationResult> sweep(
            const std::vector<Scalar>& paramRange,
            std::function<SDFPtr(Scalar)> generator,
            std::function<Scalar(SDFPtr)> objective) 
        {
            std::vector<OptimizationResult> results(paramRange.size());
            
            // Parallel execution using C++17 execution policies
            std::vector<size_t> indices(paramRange.size());
            for(size_t i = 0; i < indices.size(); ++i) indices[i] = i;

            std::for_each(std::execution::par, indices.begin(), indices.end(), [&](size_t i) {
                Scalar p = paramRange[i];
                SDFPtr sdf = generator(p);
                results[i] = {p, objective(sdf), 0.0};
            });

            return results;
        }

        /**
         * @brief Compute the sensitivity of an SDF objective at a specific point p with respect to parameter v.
         * Uses Automatic Differentiation.
         */
        static Scalar computeSensitivity(
            const Point3& queryPoint,
            std::function<DualScalar(Point3D, DualScalar)> differentiableSDFModel,
            Scalar parameterValue) 
        {
            // Seed the parameter with a derivative of 1
            DualScalar param = DualScalar::variable(parameterValue);
            Point3D pD(DualScalar(queryPoint.x), DualScalar(queryPoint.y), DualScalar(queryPoint.z));
            
            DualScalar result = differentiableSDFModel(pD, param);
            return result.der;
        }
    };
}
