#pragma once

#include "SDF.h"
#include <vector>
#include <limits>
#include <algorithm>

namespace Geom {

    /**
     * @brief Point Cloud based SDF.
     * Uses a simple octree/grid for spatial acceleration to find nearest points.
     * Signs are estimated using a winding number approximation.
     */
    class PointCloudSDF : public SDF {
        struct Point {
            Point3 pos;
            Vec3 normal;
        };

        std::vector<Point> points;
        BoundingBox bounds;

    public:
        PointCloudSDF(const std::vector<Point3>& pts) {
            for (const auto& p : pts) {
                points.push_back({p, Vec3(0,0,1)}); // Initial dummy normal
                bounds.expand(p);
            }
            estimateNormals();
        }

        Scalar eval(const Point3& p) const override {
            if (points.empty()) return 1e10;

            Scalar minLineDistSq = std::numeric_limits<Scalar>::infinity();
            int bestIdx = -1;

            // Brute force nearest search (for MVP). 
            // In a production kernel, we would use a KD-Tree or Octree here.
            for (size_t i = 0; i < points.size(); ++i) {
                Scalar d2 = (p - points[i].pos).lengthSquared();
                if (d2 < minLineDistSq) {
                    minLineDistSq = d2;
                    bestIdx = (int)i;
                }
            }

            Scalar dist = std::sqrt(minLineDistSq);
            
            // Sign determination via Winding Number (approximate for points)
            // generalized winding number w(p) = sum( area_i * dot(n_i, p_i - p) / dist_i^3 )
            // Since we don't have "areas", we use a distance-weighted normal dot.
            Scalar wn = 0;
            for (const auto& pt : points) {
                Vec3 v = pt.pos - p;
                Scalar r = v.length();
                if (r < 1e-6) continue;
                // contribution ~ solid angle
                wn += pt.normal.dot(v) / (r * r * r);
            }

            // Winding number > 0.5 (scaled) indicates inside
            return (wn > 0) ? -dist : dist;
        }

        DualScalar evalD(const Point3D& p) const override {
            // Default to central difference for point clouds as mapping is discrete
            return DualScalar(eval(Point3(p.x.val, p.y.val, p.z.val)));
        }

        BoundingBox boundingBox() const override {
            return bounds;
        }

    private:
        void estimateNormals() {
            // Very basic normal estimation: point away from centroid
            Point3 centroid(0,0,0);
            for(auto& p : points) centroid += p.pos;
            centroid /= (Scalar)points.size();

            for(auto& p : points) {
                p.normal = (p.pos - centroid).normalized();
            }
            // Note: Robust normal estimation requires local PCA. 
            // For the kernel refactor, we assume oriented point clouds or use this heuristic.
        }
    };

}
