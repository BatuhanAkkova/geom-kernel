#pragma once

#include "SDF.h"
#include "BoundingBox.h"
#include <vector>
#include <memory>
#include <array>
#include <cmath>
#include <future>
#include <algorithm>

namespace Geom {

    struct OctreeNode {
        int32_t firstChildIndex = -1; 
        Point3 vertex; 
        bool hasVertex = false;
        std::array<double, 8> cornerValues;

        bool isLeaf() const { return firstChildIndex == -1; }
    };

    class Octree {
    public:
        std::vector<OctreeNode> nodes;
        BoundingBox rootBounds;
        double minSize;

        Octree(const SDF& sdf, const BoundingBox& bounds, double minNodeSize) 
            : rootBounds(bounds), minSize(minNodeSize) {
            
            // Parallel Build Strategy:
            // 1. Create Root.
            // 2. Evaluate root corners.
            // 3. Subdivide root (create 8 children placeholders).
            // 4. Launch 8 threads to build subtrees.
            // 5. Merge subtrees.
            
            nodes.emplace_back(); // Root at 0
            
            // Corner eval for root
            evaluateCorners(0, rootBounds, sdf);
            // Check subdivision (root usually subdivides)
            if (shouldSubdivide(0, rootBounds)) {
                
                // Add 8 immediate children to main vector (slots 1..8)
                nodes.resize(1 + 8);
                nodes[0].firstChildIndex = 1;
                
                // Compute bounds for 8 children
                Vec3 half = rootBounds.size() * 0.5;
                Point3 min = rootBounds.min;
                Point3 center = min + half;
                
                auto getChildBounds = [&](int i) {
                    double x = (i & 1) ? center.x : min.x;
                    double y = (i & 2) ? center.y : min.y;
                    double z = (i & 4) ? center.z : min.z;
                    return BoundingBox(Point3(x,y,z), Point3(x + half.x, y + half.y, z + half.z));
                };

                // Launch tasks
                std::vector<std::future<std::vector<OctreeNode>>> futures;
                for(int i=0; i<8; ++i) {
                    // Copy bounds and sdf reference (thread safe SDF assumed)
                    BoundingBox b = getChildBounds(i);
                    futures.push_back(std::async(std::launch::async, [this, &sdf, b]() {
                        std::vector<OctreeNode> subNodes;
                        // Child root is index 0 in subNodes
                        subNodes.emplace_back();
                        buildRecursiveLocal(subNodes, 0, b, sdf);
                        return subNodes;
                    }));
                }
                
                // Wait and Merge
                // Child 0 corresponds to main nodes[1]
                // Child i corresponds to main nodes[1+i]
                
                // We need to fill nodes[1..8] with the roots of subtrees,
                // and append the rest.
                
                // Oops, `buildRecursiveLocal` builds a tree where 0 is the root.
                // We want to put that 0 into nodes[1+i].
                // And the rest of the subtree at end of `nodes`.
                // And adjust indices.
                
                for(int i=0; i<8; ++i) {
                    std::vector<OctreeNode> subTree = futures[i].get();
                    
                    if (subTree.empty()) continue; // Should have at least one node
                    
                    // The root of subTree (index 0) goes to nodes[1+i].
                    nodes[1+i] = subTree[0];
                    
                    // The rest (indices 1..end) are appended to `nodes`.
                    if (subTree.size() > 1) {
                        int32_t offset = static_cast<int32_t>(nodes.size()) - 1; // -1 because subTree indices start at 1
                        
                        // Fix indices in subTree (excluding root logic which we handled)
                        // Actually, simpler:
                        // subTree[0].firstChildIndex needs to be shifted by `offset`.
                        // subTree[k].firstChildIndex needs to be shifted by `offset`.
                        
                        // Wait, subTree[0] is copied to nodes[1+i].
                        // Its children start at subTree[1].
                        // In `nodes`, subTree[1] will be at `nodes.size()`.
                        // So shift = nodes.size() - 1.
                        
                        // We need to shift indices of ALL nodes in subTree before appending.
                        // (Except if childIndex is -1).
                        
                        int32_t shift = static_cast<int32_t>(nodes.size()) - 1;
                        
                        // Remap root's child index
                        if (nodes[1+i].firstChildIndex != -1) {
                            nodes[1+i].firstChildIndex += shift;
                        }
                        
                        // Append and remap others
                        // Optimization: reserve
                        nodes.reserve(nodes.size() + subTree.size() - 1);
                        
                        for (size_t k = 1; k < subTree.size(); ++k) {
                            OctreeNode& n = subTree[k];
                            if (n.firstChildIndex != -1) {
                                n.firstChildIndex += shift;
                            }
                            nodes.push_back(n);
                        }
                    }
                }
            }
        }

    private:
        // Local build builds into a local vector
        void buildRecursiveLocal(std::vector<OctreeNode>& localNodes, int32_t nodeIdx, const BoundingBox& bounds, const SDF& sdf) {
             // 1. Evaluate
             evaluateCornersLocal(localNodes[nodeIdx], bounds, sdf);
             
             // 2. Check
             if (!shouldSubdivideLocal(localNodes[nodeIdx], bounds)) return;
             
             // 3. Subdivide
             int32_t oldSize = static_cast<int32_t>(localNodes.size());
             localNodes.resize(oldSize + 8);
             localNodes[nodeIdx].firstChildIndex = oldSize;
             
             Vec3 half = bounds.size() * 0.5;
             Point3 min = bounds.min;
             Point3 center = min + half;
             
             for (int i = 0; i < 8; ++i) {
                double x = (i & 1) ? center.x : min.x;
                double y = (i & 2) ? center.y : min.y;
                double z = (i & 4) ? center.z : min.z;
                BoundingBox childBounds(Point3(x,y,z), Point3(x + half.x, y + half.y, z + half.z));
                buildRecursiveLocal(localNodes, oldSize + i, childBounds, sdf);
             }
        }

        bool shouldSubdivideLocal(const OctreeNode& node, const BoundingBox& bounds) {
            bool hasSignChange = false;
            double firstSign = std::copysign(1.0, node.cornerValues[0]);
            for (size_t i = 1; i < 8; ++i) {
                if (std::copysign(1.0, node.cornerValues[i]) != firstSign) {
                    hasSignChange = true;
                    break;
                }
            }
            
            double minAbsVal = std::abs(node.cornerValues[0]);
            for(double v : node.cornerValues) minAbsVal = std::min(minAbsVal, std::abs(v));
            
            Vec3 size = bounds.size();
            double maxSize = std::max({size.x, size.y, size.z});
            double radius = maxSize * 0.866;

            if (!hasSignChange && minAbsVal > radius) return false;
            if (maxSize <= minSize) return false;
            return true;
        }

        bool shouldSubdivide(int32_t nodeIdx, const BoundingBox& bounds) {
            return shouldSubdivideLocal(nodes[nodeIdx], bounds);
        }

        void evaluateCorners(int32_t nodeIdx, const BoundingBox& bounds, const SDF& sdf) {
             evaluateCornersLocal(nodes[nodeIdx], bounds, sdf);
        }

        void evaluateCornersLocal(OctreeNode& node, const BoundingBox& bounds, const SDF& sdf) {
            Vec3 s = bounds.size();
            Point3 p0 = bounds.min;
            Point3 p[8];
            p[0] = p0;
            p[1] = p0 + Vec3(s.x, 0, 0);
            p[2] = p0 + Vec3(s.x, s.y, 0);
            p[3] = p0 + Vec3(0, s.y, 0);
            p[4] = p0 + Vec3(0, 0, s.z);
            p[5] = p0 + Vec3(s.x, 0, s.z);
            p[6] = p0 + Vec3(s.x, s.y, s.z);
            p[7] = p0 + Vec3(0, s.y, s.z);

            for(int i=0; i<8; ++i) {
                node.cornerValues[i] = sdf.eval(p[i]);
            }
        }
        
        // Old recursive method for non-parallel if needed, but we used local instead.
        void buildRecursive(int32_t nodeIdx, const BoundingBox& bounds, const SDF& sdf) {
            // Unused in parallel version
        }
    };
}
