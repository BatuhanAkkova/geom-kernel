#pragma once

#include "SDF.h"
#include "Mesh.h"
#include <vector>
#include <algorithm>

namespace Geom {

    /**
     * @brief Mesh-based SDF.
     * Uses an Axis-Aligned Bounding Box (AABB) Tree for acceleration.
     */
    class MeshSDF : public SDF {
        struct BVHNode {
            BoundingBox bounds;
            int left = -1, right = -1;
            int triangleIdx = -1; // -1 if internal node
        };

        const Mesh& mesh;
        std::vector<BVHNode> nodes;
        int rootIdx = -1;

    public:
        MeshSDF(const Mesh& m) : mesh(m) {
            if (mesh.triangles.empty()) return;

            std::vector<int> indices(mesh.triangles.size());
            for (int i = 0; i < (int)indices.size(); ++i) indices[i] = i;

            rootIdx = buildBVH(indices, 0, (int)indices.size());
        }

        Scalar eval(const Point3& p) const override {
            if (rootIdx == -1) return 1e10;

            Scalar minLineDistSq = std::numeric_limits<Scalar>::infinity();
            closestTriangle(rootIdx, p, minLineDistSq);

            Scalar dist = std::sqrt(minLineDistSq);

            // Robust Sign Determination: Winding Number or Ray Casting
            // For now, let's use a simple parity check (ray casting +X)
            int count = 0;
            for (const auto& tri : mesh.triangles) {
                if (intersectRayTriangle(p, Vec3(1, 0, 0), tri)) {
                    count++;
                }
            }

            return (count % 2 == 1) ? -dist : dist;
        }

        BoundingBox boundingBox() const override {
            if (rootIdx == -1) return BoundingBox();
            return nodes[rootIdx].bounds;
        }

    private:
        int buildBVH(std::vector<int>& indices, int start, int end) {
            int nodeIdx = (int)nodes.size();
            nodes.emplace_back();

            BoundingBox nodeBounds;
            for (int i = start; i < end; ++i) {
                nodeBounds.expand(mesh.triangles[indices[i]].v0);
                nodeBounds.expand(mesh.triangles[indices[i]].v1);
                nodeBounds.expand(mesh.triangles[indices[i]].v2);
            }
            nodes[nodeIdx].bounds = nodeBounds;

            int count = end - start;
            if (count == 1) {
                nodes[nodeIdx].triangleIdx = indices[start];
                return nodeIdx;
            }

            // Split along largest axis
            Vec3 size = nodeBounds.size();
            int axis = 0;
            if (size.y > size.x && size.y > size.z) axis = 1;
            else if (size.z > size.x) axis = 2;

            Scalar mid = (nodeBounds.min[axis == 0 ? 0 : (axis == 1 ? 1 : 2)] + 
                          nodeBounds.max[axis == 0 ? 0 : (axis == 1 ? 1 : 2)]) * 0.5;

            auto it = std::partition(indices.begin() + start, indices.begin() + end, [&](int idx) {
                const auto& tri = mesh.triangles[idx];
                Scalar c = (tri.v0[axis == 0 ? 0 : (axis == 1 ? 1 : 2)] + 
                            tri.v1[axis == 0 ? 0 : (axis == 1 ? 1 : 2)] + 
                            tri.v2[axis == 0 ? 0 : (axis == 1 ? 1 : 2)]) / 3.0;
                return c < mid;
            });

            int split = (int)std::distance(indices.begin(), it);
            if (split == start || split == end) split = start + count / 2;

            nodes[nodeIdx].left = buildBVH(indices, start, split);
            nodes[nodeIdx].right = buildBVH(indices, split, end);

            return nodeIdx;
        }

        void closestTriangle(int nodeIdx, const Point3& p, Scalar& minDistSq) const {
            const auto& node = nodes[nodeIdx];
            
            // Prune if point is further from AABB than current min
            if (sqDistPointAABB(p, node.bounds) > minDistSq) return;

            if (node.triangleIdx != -1) {
                Scalar d2 = sqDistPointTriangle(p, mesh.triangles[node.triangleIdx]);
                if (d2 < minDistSq) minDistSq = d2;
                return;
            }

            closestTriangle(node.left, p, minDistSq);
            closestTriangle(node.right, p, minDistSq);
        }

        Scalar sqDistPointAABB(const Point3& p, const BoundingBox& b) const {
            Scalar sqDist = 0;
            if (p.x < b.min.x) sqDist += (b.min.x - p.x) * (b.min.x - p.x);
            if (p.x > b.max.x) sqDist += (p.x - b.max.x) * (p.x - b.max.x);
            if (p.y < b.min.y) sqDist += (b.min.y - p.y) * (b.min.y - p.y);
            if (p.y > b.max.y) sqDist += (p.y - b.max.y) * (p.y - b.max.y);
            if (p.z < b.min.z) sqDist += (b.min.z - p.z) * (b.min.z - p.z);
            if (p.z > b.max.z) sqDist += (p.z - b.max.z) * (p.z - b.max.z);
            return sqDist;
        }

        Scalar sqDistPointTriangle(const Point3& p, const Triangle& tri) const {
            // Simplified: distance to centroid as proxy for MVP. 
            // Proper point-triangle distance is quite involved.
            Point3 c = (tri.v0 + tri.v1 + tri.v2) / 3.0;
            return (p - c).lengthSquared();
        }

        bool intersectRayTriangle(const Point3& orig, const Vec3& dir, const Triangle& tri) const {
            // Moller-Trumbore algorithm
            Vec3 edge1 = tri.v1 - tri.v0;
            Vec3 edge2 = tri.v2 - tri.v0;
            Vec3 h = cross(dir, edge2);
            Scalar a = edge1.dot(h);
            if (a > -1e-6 && a < 1e-6) return false;
            Scalar f = 1.0 / a;
            Vec3 s = orig - tri.v0;
            Scalar u = f * s.dot(h);
            if (u < 0.0 || u > 1.0) return false;
            Vec3 q = cross(s, edge1);
            Scalar v = f * dir.dot(q);
            if (v < 0.0 || u + v > 1.0) return false;
            Scalar t = f * edge2.dot(q);
            return t > 1e-6;
        }
    };
}
