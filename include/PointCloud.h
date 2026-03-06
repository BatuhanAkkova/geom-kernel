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
    class PointCloudSDF : public SDFNode<PointCloudSDF> {
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

        template <typename T>
        T evaluate(const Vec3T<T>& p) const {
            using std::max; using std::min; using std::abs; using std::sqrt; using std::pow; using std::sin; using std::cos;
            if (points.empty()) return static_cast<T>(1e10);

            T minLineDistSq = static_cast<T>(1e18);
            int bestIdx = -1;

            // Brute force nearest search (for MVP). 
            for (size_t i = 0; i < points.size(); ++i) {
                Vec3T<T> pi(static_cast<T>(points[i].pos.x), static_cast<T>(points[i].pos.y), static_cast<T>(points[i].pos.z));
                T d2 = (p - pi).lengthSquared();
                if (d2 < minLineDistSq) {
                    minLineDistSq = d2;
                    bestIdx = (int)i;
                }
            }

            T dist = sqrt(minLineDistSq);
            
            T wn = static_cast<T>(0);
            for (const auto& pt : points) {
                Vec3T<T> p_pos(static_cast<T>(pt.pos.x), static_cast<T>(pt.pos.y), static_cast<T>(pt.pos.z));
                Vec3T<T> p_normal(static_cast<T>(pt.normal.x), static_cast<T>(pt.normal.y), static_cast<T>(pt.normal.z));
                Vec3T<T> v = p_pos - p;
                T r = v.length();
                if (r < static_cast<T>(1e-6)) continue;
                wn = wn + p_normal.dot(v) / (r * r * r);
            }

            return (wn > static_cast<T>(0)) ? -dist : dist;
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
