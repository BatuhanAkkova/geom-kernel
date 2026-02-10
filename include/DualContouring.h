#pragma once

#include "Octree.h"
#include "Mesh.h"
#include "QEF.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <future>
#include <mutex>

namespace Geom {

    class DualContouring {
    public:
        static void generateMesh(Octree& octree, const SDF& sdf, Mesh& mesh) {
            if (octree.nodes.empty()) return;
            
            // Parallel Vertex Generation
            // 8 Tasks for root children
            if (!octree.nodes[0].isLeaf()) {
                 std::vector<std::future<void>> futures;
                 Vec3 half = octree.rootBounds.size() * 0.5;
                 Point3 min = octree.rootBounds.min;
                 Point3 center = min + half;
                 
                 for(int i=0; i<8; ++i) {
                    double x = (i & 1) ? center.x : min.x;
                    double y = (i & 2) ? center.y : min.y;
                    double z = (i & 4) ? center.z : min.z;
                    BoundingBox childBounds(Point3(x,y,z), Point3(x + half.x, y + half.y, z + half.z));
                    int32_t childIdx = octree.nodes[0].firstChildIndex + i;
                    
                    futures.push_back(std::async(std::launch::async, [&octree, childIdx, childBounds, &sdf]() {
                        computeVertices(octree, childIdx, childBounds, sdf);
                    }));
                 }
                 for(auto& f : futures) f.get();
            } else {
                computeVertices(octree, 0, octree.rootBounds, sdf);
            }

            // Parallel Face/Cell Processing
            // Similar strategy: 8 sub-tasks return partial meshes
            if (!octree.nodes[0].isLeaf()) {
                 std::vector<std::future<Mesh>> futures;
                 Vec3 half = octree.rootBounds.size() * 0.5;
                 Point3 min = octree.rootBounds.min;
                 Point3 center = min + half;

                 for(int i=0; i<8; ++i) {
                    double x = (i & 1) ? center.x : min.x;
                    double y = (i & 2) ? center.y : min.y;
                    double z = (i & 4) ? center.z : min.z;
                    BoundingBox childBounds(Point3(x,y,z), Point3(x + half.x, y + half.y, z + half.z));
                    
                    int32_t childIdx = octree.nodes[0].firstChildIndex + i;

                    futures.push_back(std::async(std::launch::async, [&octree, childIdx, childBounds, &sdf]() {
                        Mesh subMesh;
                        processCell(octree, childIdx, childBounds, sdf, subMesh);
                        return subMesh;
                    }));
                 }
                 
                 // Collect results
                 std::vector<Mesh> subMeshes;
                 for(auto& f : futures) {
                     subMeshes.push_back(f.get());
                 }
                 
                 // Process Root Internal Faces (Synchronous or Task based)
                 // The root itself calls processCell BUT we only ran children recursively.
                 // We missed the Inter-Child faces of the root!
                 // The standard recursion is:
                 // processCell(node) -> calls Recurse(children) -> calls processFace(internal).
                 // We replaced Recurse(children) with tasks.
                 // We still need to run the rest of processCell(root).
                 
                 // Manually run root's face processing
                 Mesh rootInternalMesh;
                 processRootFaces(octree, 0, octree.rootBounds, sdf, rootInternalMesh);
                 
                 // Merge all
                 for(const auto& m : subMeshes) {
                     mesh.triangles.insert(mesh.triangles.end(), m.triangles.begin(), m.triangles.end());
                 }
                 mesh.triangles.insert(mesh.triangles.end(), rootInternalMesh.triangles.begin(), rootInternalMesh.triangles.end());
                 
            } else {
                 processCell(octree, 0, octree.rootBounds, sdf, mesh);
            }
        }

    private:
        // Expose internal face processing logic separately from full recursion
        static void processRootFaces(Octree& octree, int32_t nodeIdx, const BoundingBox& bounds, const SDF& sdf, Mesh& mesh) {
            // This is the second half of processCell
            OctreeNode& node = octree.nodes[nodeIdx];
            if (node.isLeaf()) return;
            
            int32_t c = node.firstChildIndex;
            Vec3 half = bounds.size() * 0.5;
            Point3 min = bounds.min;
            Point3 center = min + half;

            auto getChildBounds = [&](int i) {
                double x = (i & 1) ? center.x : min.x;
                double y = (i & 2) ? center.y : min.y;
                double z = (i & 4) ? center.z : min.z;
                return BoundingBox(Point3(x,y,z), Point3(x + half.x, y + half.y, z + half.z));
            };

            processFace(octree, c+0, c+1, 0, getChildBounds(0), getChildBounds(1), sdf, mesh);
            processFace(octree, c+2, c+3, 0, getChildBounds(2), getChildBounds(3), sdf, mesh);
            processFace(octree, c+4, c+5, 0, getChildBounds(4), getChildBounds(5), sdf, mesh);
            processFace(octree, c+6, c+7, 0, getChildBounds(6), getChildBounds(7), sdf, mesh);

            processFace(octree, c+0, c+2, 1, getChildBounds(0), getChildBounds(2), sdf, mesh);
            processFace(octree, c+1, c+3, 1, getChildBounds(1), getChildBounds(3), sdf, mesh);
            processFace(octree, c+4, c+6, 1, getChildBounds(4), getChildBounds(6), sdf, mesh);
            processFace(octree, c+5, c+7, 1, getChildBounds(5), getChildBounds(7), sdf, mesh);

            processFace(octree, c+0, c+4, 2, getChildBounds(0), getChildBounds(4), sdf, mesh);
            processFace(octree, c+1, c+5, 2, getChildBounds(1), getChildBounds(5), sdf, mesh);
            processFace(octree, c+2, c+6, 2, getChildBounds(2), getChildBounds(6), sdf, mesh);
            processFace(octree, c+3, c+7, 2, getChildBounds(3), getChildBounds(7), sdf, mesh);
            
            processEdge(octree, c+0, c+2, c+6, c+4, 0, getChildBounds(0), getChildBounds(2), getChildBounds(6), getChildBounds(4), sdf, mesh);
            processEdge(octree, c+0, c+1, c+5, c+4, 1, getChildBounds(0), getChildBounds(1), getChildBounds(5), getChildBounds(4), sdf, mesh);
            processEdge(octree, c+0, c+2, c+3, c+1, 2, getChildBounds(0), getChildBounds(2), getChildBounds(3), getChildBounds(1), sdf, mesh);
        }

        static void computeVertices(Octree& octree, int32_t nodeIdx, const BoundingBox& bounds, const SDF& sdf) {
            OctreeNode& node = octree.nodes[nodeIdx];
            
            if (node.isLeaf()) {
                calculateLeafVertex(node, bounds, sdf);
                return;
            }

            Vec3 half = bounds.size() * 0.5;
            Point3 min = bounds.min;
            Point3 center = min + half;

            for (int i = 0; i < 8; ++i) {
                 double x = (i & 1) ? center.x : min.x;
                 double y = (i & 2) ? center.y : min.y;
                 double z = (i & 4) ? center.z : min.z;
                 BoundingBox childBounds(Point3(x,y,z), Point3(x + half.x, y + half.y, z + half.z));
                 
                 computeVertices(octree, node.firstChildIndex + i, childBounds, sdf);
            }
        }

        static void calculateLeafVertex(OctreeNode& node, const BoundingBox& bounds, const SDF& sdf) {
            Vec3 s = bounds.size();
            Point3 min = bounds.min;
            static const int binToMC[8] = {0, 1, 3, 2, 4, 5, 7, 6};
            static const int edges[12][2] = {
                {0,1}, {2,3}, {4,5}, {6,7}, // X
                {0,2}, {1,3}, {4,6}, {5,7}, // Y
                {0,4}, {1,5}, {2,6}, {3,7}  // Z
            };
            
            Point3 corners[8];
            for(int i=0; i<8; ++i) {
                double x = (i & 1) ? min.x + s.x : min.x;
                double y = (i & 2) ? min.y + s.y : min.y;
                double z = (i & 4) ? min.z + s.z : min.z;
                corners[i] = Point3(x,y,z);
            }

            QEF qef;
            Point3 massPoint(0,0,0);
            int intersectionCount = 0;

            for (int e = 0; e < 12; ++e) {
                int c1 = edges[e][0];
                int c2 = edges[e][1];
                double val1 = node.cornerValues[binToMC[c1]];
                double val2 = node.cornerValues[binToMC[c2]];
                
                if ((val1 > 0 && val2 < 0) || (val1 < 0 && val2 > 0)) {
                    double mu = (0 - val1) / (val2 - val1);
                    Point3 inter = corners[c1] + (corners[c2] - corners[c1]) * mu;
                    
                    // Add to QEF
                    Vec3 normal = sdf.gradient(inter);
                    qef.add(inter, normal);
                    
                    massPoint = massPoint + inter;
                    intersectionCount++;
                }
            }

            if (intersectionCount > 0) {
                // Solve QEF
                Point3 qefPoint = qef.solve();
                massPoint = massPoint * (1.0 / intersectionCount);
                
                // Clamp QEF point to bounds?
                // QEF can shoot outside if matrix is ill-conditioned (handled by fallback) 
                // or if geometry implies vertex outside (should restrict).
                // Simple clamp to bounds:
                if (qefPoint.x < min.x || qefPoint.x > min.x + s.x ||
                    qefPoint.y < min.y || qefPoint.y > min.y + s.y ||
                    qefPoint.z < min.z || qefPoint.z > min.z + s.z) {
                     // Fallback to mass point if outside
                     node.vertex = massPoint;
                } else {
                     node.vertex = qefPoint;
                }
                
                node.hasVertex = true;
            } else {
                node.hasVertex = false;
            }
        }

        static void processCell(Octree& octree, int32_t nodeIdx, const BoundingBox& bounds, const SDF& sdf, Mesh& mesh) {
            OctreeNode& node = octree.nodes[nodeIdx];
            if (node.isLeaf()) return;

            Vec3 half = bounds.size() * 0.5;
            Point3 min = bounds.min;
            Point3 center = min + half;
            
            // Recurse Children
            for (int i = 0; i < 8; ++i) {
                 double x = (i & 1) ? center.x : min.x;
                 double y = (i & 2) ? center.y : min.y;
                 double z = (i & 4) ? center.z : min.z;
                 BoundingBox childBounds(Point3(x,y,z), Point3(x + half.x, y + half.y, z + half.z));
                 processCell(octree, node.firstChildIndex + i, childBounds, sdf, mesh);
            }

            // Face neighbors (same as processRootFaces but inline)
            processRootFaces(octree, nodeIdx, bounds, sdf, mesh);
        }

        static void processFace(Octree& octree, int32_t n1Idx, int32_t n2Idx, int dir, 
                              const BoundingBox& b1, const BoundingBox& b2, 
                              const SDF& sdf, Mesh& mesh) {
            OctreeNode& n1 = octree.nodes[n1Idx];
            OctreeNode& n2 = octree.nodes[n2Idx];
            bool leaf1 = n1.isLeaf();
            bool leaf2 = n2.isLeaf();
            
            if (leaf1 && leaf2) return;
            
            auto get = [&](int32_t parentIdx, bool isLeaf, int childOffset) {
                return isLeaf ? parentIdx : octree.nodes[parentIdx].firstChildIndex + childOffset;
            };
            auto getBounds = [&](const BoundingBox& parentBounds, int i) {
                Vec3 half = parentBounds.size() * 0.5;
                Point3 min = parentBounds.min;
                Point3 center = min + half;
                double x = (i & 1) ? center.x : min.x;
                double y = (i & 2) ? center.y : min.y;
                double z = (i & 4) ? center.z : min.z;
                return BoundingBox(Point3(x,y,z), Point3(x + half.x, y + half.y, z + half.z));
            };

            if (dir == 0) { // X
                processFace(octree, get(n1Idx, leaf1, 1), get(n2Idx, leaf2, 0), 0, getBounds(b1, 1), getBounds(b2, 0), sdf, mesh);
                processFace(octree, get(n1Idx, leaf1, 3), get(n2Idx, leaf2, 2), 0, getBounds(b1, 3), getBounds(b2, 2), sdf, mesh);
                processFace(octree, get(n1Idx, leaf1, 5), get(n2Idx, leaf2, 4), 0, getBounds(b1, 5), getBounds(b2, 4), sdf, mesh);
                processFace(octree, get(n1Idx, leaf1, 7), get(n2Idx, leaf2, 6), 0, getBounds(b1, 7), getBounds(b2, 6), sdf, mesh);
                
                processEdge(octree, get(n1Idx, leaf1, 1), get(n2Idx, leaf2, 0), get(n2Idx, leaf2, 4), get(n1Idx, leaf1, 5), 1, 
                            getBounds(b1, 1), getBounds(b2, 0), getBounds(b2, 4), getBounds(b1, 5), sdf, mesh);
                processEdge(octree, get(n1Idx, leaf1, 1), get(n2Idx, leaf2, 0), get(n2Idx, leaf2, 2), get(n1Idx, leaf1, 3), 2, 
                            getBounds(b1, 1), getBounds(b2, 0), getBounds(b2, 2), getBounds(b1, 3), sdf, mesh);
            } else if (dir == 1) { // Y
                processFace(octree, get(n1Idx, leaf1, 2), get(n2Idx, leaf2, 0), 1, getBounds(b1, 2), getBounds(b2, 0), sdf, mesh);
                processFace(octree, get(n1Idx, leaf1, 3), get(n2Idx, leaf2, 1), 1, getBounds(b1, 3), getBounds(b2, 1), sdf, mesh);
                processFace(octree, get(n1Idx, leaf1, 6), get(n2Idx, leaf2, 4), 1, getBounds(b1, 6), getBounds(b2, 4), sdf, mesh);
                processFace(octree, get(n1Idx, leaf1, 7), get(n2Idx, leaf2, 5), 1, getBounds(b1, 7), getBounds(b2, 5), sdf, mesh);
                
                processEdge(octree, get(n1Idx, leaf1, 2), get(n2Idx, leaf2, 0), get(n2Idx, leaf2, 4), get(n1Idx, leaf1, 6), 0,
                            getBounds(b1, 2), getBounds(b2, 0), getBounds(b2, 4), getBounds(b1, 6), sdf, mesh);
                processEdge(octree, get(n1Idx, leaf1, 2), get(n2Idx, leaf2, 0), get(n2Idx, leaf2, 1), get(n1Idx, leaf1, 3), 2,
                            getBounds(b1, 2), getBounds(b2, 0), getBounds(b2, 1), getBounds(b1, 3), sdf, mesh);
            } else { // Z
                processFace(octree, get(n1Idx, leaf1, 4), get(n2Idx, leaf2, 0), 2, getBounds(b1, 4), getBounds(b2, 0), sdf, mesh);
                processFace(octree, get(n1Idx, leaf1, 5), get(n2Idx, leaf2, 1), 2, getBounds(b1, 5), getBounds(b2, 1), sdf, mesh);
                processFace(octree, get(n1Idx, leaf1, 6), get(n2Idx, leaf2, 2), 2, getBounds(b1, 6), getBounds(b2, 2), sdf, mesh);
                processFace(octree, get(n1Idx, leaf1, 7), get(n2Idx, leaf2, 3), 2, getBounds(b1, 7), getBounds(b2, 3), sdf, mesh);
                
                processEdge(octree, get(n1Idx, leaf1, 4), get(n2Idx, leaf2, 0), get(n2Idx, leaf2, 1), get(n1Idx, leaf1, 5), 0,
                            getBounds(b1, 4), getBounds(b2, 0), getBounds(b2, 1), getBounds(b1, 5), sdf, mesh);
                processEdge(octree, get(n1Idx, leaf1, 4), get(n2Idx, leaf2, 0), get(n2Idx, leaf2, 2), get(n1Idx, leaf1, 6), 1,
                            getBounds(b1, 4), getBounds(b2, 0), getBounds(b2, 2), getBounds(b1, 6), sdf, mesh);
            }
        }

        static void processEdge(Octree& octree, int32_t n1Idx, int32_t n2Idx, int32_t n3Idx, int32_t n4Idx, int dir, 
                                const BoundingBox& b1, const BoundingBox& b2, const BoundingBox& b3, const BoundingBox& b4, 
                                const SDF& sdf, Mesh& mesh) {
            OctreeNode& n1 = octree.nodes[n1Idx];
            OctreeNode& n2 = octree.nodes[n2Idx];
            OctreeNode& n3 = octree.nodes[n3Idx];
            OctreeNode& n4 = octree.nodes[n4Idx];
            
            bool allLeaves = n1.isLeaf() && n2.isLeaf() && n3.isLeaf() && n4.isLeaf();
            
            if (allLeaves) {
                if (n1.hasVertex && n2.hasVertex && n3.hasVertex && n4.hasVertex) {
                    Point3 pStart, pEnd;
                    if (dir == 0) { // X
                         double y = b1.max.y;
                         double z = b1.max.z;
                         double x_min = std::max({b1.min.x, b2.min.x, b3.min.x, b4.min.x});
                         double x_max = std::min({b1.max.x, b2.max.x, b3.max.x, b4.max.x});
                         pStart = Point3(x_min, y, z);
                         pEnd = Point3(x_max, y, z);
                    } else if (dir == 1) { // Y
                         double x = b1.max.x;
                         double z = b1.max.z;
                         double y_min = std::max({b1.min.y, b2.min.y, b3.min.y, b4.min.y});
                         double y_max = std::min({b1.max.y, b2.max.y, b3.max.y, b4.max.y});
                         pStart = Point3(x, y_min, z);
                         pEnd = Point3(x, y_max, z);
                    } else { // Z
                         double x = b1.max.x;
                         double y = b1.max.y;
                         double z_min = std::max({b1.min.z, b2.min.z, b3.min.z, b4.min.z});
                         double z_max = std::min({b1.max.z, b2.max.z, b3.max.z, b4.max.z});
                         pStart = Point3(x, y, z_min);
                         pEnd = Point3(x, y, z_max);
                    }
                    
                    double val1 = sdf.eval(pStart);
                    double val2 = sdf.eval(pEnd);
                    
                    if ((val1 > 0 && val2 < 0) || (val1 < 0 && val2 > 0)) {
                         if (val1 < 0) {
                             mesh.triangles.emplace_back(n1.vertex, n2.vertex, n3.vertex);
                             mesh.triangles.emplace_back(n1.vertex, n3.vertex, n4.vertex);
                         } else {
                             mesh.triangles.emplace_back(n1.vertex, n3.vertex, n2.vertex);
                             mesh.triangles.emplace_back(n1.vertex, n4.vertex, n3.vertex);
                         }
                    }
                }
                return;
            }
            
            auto get = [&](int32_t parentIdx, bool isLeaf, int childOffset) {
                return isLeaf ? parentIdx : octree.nodes[parentIdx].firstChildIndex + childOffset;
            };
            auto getBounds = [&](const BoundingBox& parentBounds, int i) {
                Vec3 half = parentBounds.size() * 0.5;
                Point3 min = parentBounds.min;
                Point3 center = min + half;
                double x = (i & 1) ? center.x : min.x;
                double y = (i & 2) ? center.y : min.y;
                double z = (i & 4) ? center.z : min.z;
                return BoundingBox(Point3(x,y,z), Point3(x + half.x, y + half.y, z + half.z));
            };
            
            bool l1 = n1.isLeaf(), l2 = n2.isLeaf(), l3 = n3.isLeaf(), l4 = n4.isLeaf();

            if (dir == 0) { // X
                processEdge(octree, get(n1Idx, l1, 6), get(n2Idx, l2, 4), get(n3Idx, l3, 0), get(n4Idx, l4, 2), 0,
                            getBounds(b1, 6), getBounds(b2, 4), getBounds(b3, 0), getBounds(b4, 2), sdf, mesh);
                processEdge(octree, get(n1Idx, l1, 7), get(n2Idx, l2, 5), get(n3Idx, l3, 1), get(n4Idx, l4, 3), 0,
                            getBounds(b1, 7), getBounds(b2, 5), getBounds(b3, 1), getBounds(b4, 3), sdf, mesh);
            } else if (dir == 1) { // Y
                processEdge(octree, get(n1Idx, l1, 5), get(n2Idx, l2, 4), get(n3Idx, l3, 0), get(n4Idx, l4, 1), 1,
                            getBounds(b1, 5), getBounds(b2, 4), getBounds(b3, 0), getBounds(b4, 1), sdf, mesh);
                processEdge(octree, get(n1Idx, l1, 7), get(n2Idx, l2, 6), get(n3Idx, l3, 2), get(n4Idx, l4, 3), 1,
                            getBounds(b1, 7), getBounds(b2, 6), getBounds(b3, 2), getBounds(b4, 3), sdf, mesh);
            } else { // Z
                processEdge(octree, get(n1Idx, l1, 3), get(n2Idx, l2, 1), get(n3Idx, l3, 0), get(n4Idx, l4, 2), 2,
                            getBounds(b1, 3), getBounds(b2, 1), getBounds(b3, 0), getBounds(b4, 2), sdf, mesh);
                processEdge(octree, get(n1Idx, l1, 7), get(n2Idx, l2, 5), get(n3Idx, l3, 4), get(n4Idx, l4, 6), 2,
                            getBounds(b1, 7), getBounds(b2, 5), getBounds(b3, 4), getBounds(b4, 6), sdf, mesh);
            }
        }
    };
}
