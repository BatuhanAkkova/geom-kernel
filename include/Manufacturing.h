#pragma once

#include "SDF.h"
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace Geom {

    enum class IssueType {
        Overhang,
        HighCurvature,
        ThinFeature,
        DisconnectedIsland
    };

    struct ValidationIssue {
        Point3 position;
        IssueType type;
        Scalar value; // The measured value (e.g. angle, radius, thickness)
        std::string description;
    };

    struct ManufacturingConfig {
        Vec3 buildDirection = Vec3(0, 0, 1);
        Scalar maxOverhangAngle = 45.0; // Degrees from horizontal. 0 = flat ceiling, 90 = vertical.
        Scalar minCurvatureRadius = 2.0; // SI units (e.g. mm)
        Scalar minBeadWidth = 4.0;      // SI units (e.g. mm)
        Scalar layerHeight = 1.0;       // SI units (e.g. mm)
    };

    class ManufacturingValidator {
    public:
        ManufacturingValidator(const ManufacturingConfig& config) : m_config(config) {}

        /**
         * @brief Check if a point on the surface is an overhang.
         * @param sdf The geometry
         * @param p The point on the surface (sdf.eval(p) approx 0)
         * @return true if it's a critical overhang.
         */
        bool isOverhang(const SDFPtr& sdf, const Point3& p, Scalar& angleOut) const {
            Vec3 normal = sdf->gradient(p).normalized();
            // Angle between normal and build direction
            // If normal points exactly along -buildDirection, angle is 180 deg (flat ceiling)
            // If normal is perpendicular to buildDirection, angle is 90 deg (vertical wall)
            Scalar dot = normal.dot(m_config.buildDirection);
            
            // Convert dot to angle from horizontal
            // horizontal: dot = 0 -> angle = 90
            // ceiling: dot = -1 -> angle = 0
            // floor: dot = 1 -> angle = 180 (not an overhang)
            
            // angle_from_horizontal = 90 + asin(dot) * 180 / PI
            // But simpler: if dot < -sin(maxOverhangAngle)
            Scalar radLimit = m_config.maxOverhangAngle * M_PI / 180.0;
            Scalar sinLimit = -std::cos(radLimit); // e.g. 45 deg -> -0.707
            
            angleOut = std::acos(std::clamp(dot, -1.0, 1.0)) * 180.0 / M_PI;
            
            // We consider it an overhang if it faces downwards more than the limit.
            // A 45 deg overhang means the normal is 135 deg from build direction (0 is up).
            // So if angle > (180 - maxOverhangAngle) or dot < -cos(maxOverhangAngle)
            return dot < -std::cos(radLimit);
        }

        /**
         * @brief Check for thin features using ray marching into the solid.
         */
        bool isTooThin(const SDFPtr& sdf, const Point3& p, Scalar& thicknessOut) const {
            Vec3 normal = sdf->gradient(p).normalized();
            Vec3 inward = normal * -1.0;
            
            // Simple check: sample SDF along the inward normal at minBeadWidth
            // If sdf(p + inward * minBeadWidth) > 0, it means we exited the solid.
            // For better accuracy, we could search for the zero crossing.
            Scalar t = 0.1 * m_config.minBeadWidth;
            Scalar step = 0.1 * m_config.minBeadWidth;
            thicknessOut = m_config.minBeadWidth;
            
            for (int i = 0; i < 15; ++i) {
                Point3 sample = p + inward * t;
                Scalar dist = sdf->eval(sample);
                if (dist > 0) {
                    thicknessOut = t; // Found the other side
                    return true; 
                }
                // If dist is very negative, we are deep inside.
                // But SDF might not be perfectly normalized.
                t += step;
                if (t > m_config.minBeadWidth * 1.5) break;
            }
            
            return false;
        }

        /**
         * @brief Check curvature radius. 
         * Uses finite difference to estimate Hessian and then principal curvatures.
         */
        bool hasHighCurvature(const SDFPtr& sdf, const Point3& p, Scalar& radiusOut) const {
            // Simplified: compute gradient at offset points to see how much normal changes
            Scalar h = 0.01; // Step size
            Vec3 n = sdf->gradient(p).normalized();
            
            // Estimate divergence of normal (Mean Curvature H = div(n)/2)
            // div(n) = dnx/dx + dny/dy + dnz/dz
            Vec3 nx = sdf->gradient(p + Vec3(h, 0, 0)).normalized();
            Vec3 ny = sdf->gradient(p + Vec3(0, h, 0)).normalized();
            Vec3 nz = sdf->gradient(p + Vec3(0, 0, h)).normalized();
            
            Scalar div = (nx.x - n.x) / h + (ny.y - n.y) / h + (nz.z - n.z) / h;
            Scalar H = std::abs(div) / 2.0;
            
            if (H < 1e-6) {
                radiusOut = 1e10; // Flat
                return false;
            }
            
            radiusOut = 1.0 / H;
            return radiusOut < m_config.minCurvatureRadius;
        }

        std::vector<ValidationIssue> validate(const SDFPtr& sdf, const BoundingBox& bbox, Scalar resolution) const {
            std::vector<ValidationIssue> issues;
            
            // Sample the bounding box
            for (Scalar z = bbox.min.z; z <= bbox.max.z; z += resolution) {
                for (Scalar y = bbox.min.y; y <= bbox.max.y; y += resolution) {
                    for (Scalar x = bbox.min.x; x <= bbox.max.x; x += resolution) {
                        Point3 p(x, y, z);
                        Scalar d = sdf->eval(p);
                        
                        // We only care about points near the surface
                        if (std::abs(d) < resolution * 0.5) {
                            Scalar val;
                            if (isOverhang(sdf, p, val)) {
                                issues.push_back({p, IssueType::Overhang, val, "Overhang angle exceeds limit"});
                            }
                            if (hasHighCurvature(sdf, p, val)) {
                                issues.push_back({p, IssueType::HighCurvature, val, "Curvature radius too small"});
                            }
                            if (isTooThin(sdf, p, val)) {
                                issues.push_back({p, IssueType::ThinFeature, val, "Feature thickness below bead width"});
                            }
                        }
                    }
                }
            }
            return issues;
        }

    private:
        ManufacturingConfig m_config;
    };
}
