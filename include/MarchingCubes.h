#pragma once

#include "SDF.h"
#include "Mesh.h"
#include "MCTables.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <future> // Added for std::async
#include <thread> // Added for std::thread::hardware_concurrency

namespace Geom {

    class MarchingCubes {
    public:
        static void march(const SDF& sdf, const BoundingBox& bounds, double resolution, Mesh& mesh) {
            if (resolution <= 1e-6) return; // Use small epsilon

            // Calculate number of steps
            Vec3 size = bounds.size();
            int nx = static_cast<int>(std::ceil(size.x / resolution));
            int ny = static_cast<int>(std::ceil(size.y / resolution));
            int nz = static_cast<int>(std::ceil(size.z / resolution));

            // Multithreading setup
            int numThreads = std::thread::hardware_concurrency();
            if (numThreads == 0) numThreads = 1;
            
            // If easy task, don't spin up threads
            if (nz < numThreads * 2) numThreads = 1;

            std::vector<std::future<std::vector<Triangle>>> futures;
            int zChunkSize = nz / numThreads;

            for (int t = 0; t < numThreads; ++t) {
                int zStart = t * zChunkSize;
                int zEnd = (t == numThreads - 1) ? nz : (t + 1) * zChunkSize;

                futures.push_back(std::async(std::launch::async, 
                    [zStart, zEnd, nx, ny, &sdf, &bounds, resolution]() {
                        std::vector<Triangle> localTriangles;
                        // Precompute diagonal for space skipping
                        // Distance from center of a voxel to corner is sqrt(3) * (res/2)
                        // Safety factor 1.01
                         double voxelDiagHalf = std::sqrt(3.0) * resolution * 0.5 * 1.01;

                        for (int z = zStart; z < zEnd; ++z) {
                            for (int y = 0; y < ny; ++y) {
                                // Cache for X-axis reuse
                                // We store the 4 values of the "right" face of the previous cell
                                // which become the "left" face of the current cell.
                                // Right face indices: 1, 2, 5, 6
                                // Left face indices: 0, 3, 4, 7
                                double prevRight[4]; 
                                bool hasPrev = false;

                                for (int x = 0; x < nx; ++x) {
                                    Point3 p0 = bounds.min + Vec3(x * resolution, y * resolution, z * resolution);
                                    
                                    // Space Skipping (Early Out)
                                    Point3 center = p0 + Vec3(resolution * 0.5, resolution * 0.5, resolution * 0.5);
                                    double centerDist = sdf.eval(center);
                                    if (std::abs(centerDist) > voxelDiagHalf) {
                                        hasPrev = false; // Break the chain
                                        continue;
                                    }

                                    // 8 corners
                                    Point3 p[8];
                                    double val[8];
                                    
                                    // Coordinates
                                    // We need coord for all 8 points to be sure, or just the ones we evaluate.
                                    // To be safe and simple, let's compute all P positions (cheap additions)
                                    p[0] = p0;
                                    p[1] = p0 + Vec3(resolution, 0, 0);
                                    p[2] = p0 + Vec3(resolution, resolution, 0);
                                    p[3] = p0 + Vec3(0, resolution, 0);
                                    p[4] = p0 + Vec3(0, 0, resolution);
                                    p[5] = p0 + Vec3(resolution, 0, resolution);
                                    p[6] = p0 + Vec3(resolution, resolution, resolution);
                                    p[7] = p0 + Vec3(0, resolution, resolution);
                                    
                                    // Evaluate or Reuse
                                    if (hasPrev) {
                                        // Reuse left face from previous right face
                                        val[0] = prevRight[0]; // prev 1
                                        val[3] = prevRight[1]; // prev 2
                                        val[4] = prevRight[2]; // prev 5
                                        val[7] = prevRight[3]; // prev 6
                                    } else {
                                        // Compute left face
                                        val[0] = sdf.eval(p[0]);
                                        val[3] = sdf.eval(p[3]);
                                        val[4] = sdf.eval(p[4]);
                                        val[7] = sdf.eval(p[7]);
                                    }

                                    // Always compute right face
                                    val[1] = sdf.eval(p[1]);
                                    val[2] = sdf.eval(p[2]);
                                    val[5] = sdf.eval(p[5]);
                                    val[6] = sdf.eval(p[6]);
                                    
                                    // Update cache for next iteration
                                    prevRight[0] = val[1];
                                    prevRight[1] = val[2];
                                    prevRight[2] = val[5];
                                    prevRight[3] = val[6];
                                    hasPrev = true;
                                    
                                    int cubeIndex = 0;
                                    for(int i=0; i<8; ++i) {
                                        if (val[i] < 0) cubeIndex |= (1 << i);
                                    }
                                    
                                    if (edgeTable[cubeIndex] == 0) continue;
                                    
                                    // Interpolate vertices
                                    Point3 vertlist[12];
                                    
                                    if (edgeTable[cubeIndex] & 1)    vertlist[0] = vertexInterp(p[0], p[1], val[0], val[1]);
                                    if (edgeTable[cubeIndex] & 2)    vertlist[1] = vertexInterp(p[1], p[2], val[1], val[2]);
                                    if (edgeTable[cubeIndex] & 4)    vertlist[2] = vertexInterp(p[2], p[3], val[2], val[3]);
                                    if (edgeTable[cubeIndex] & 8)    vertlist[3] = vertexInterp(p[3], p[0], val[3], val[0]);
                                    if (edgeTable[cubeIndex] & 16)   vertlist[4] = vertexInterp(p[4], p[5], val[4], val[5]);
                                    if (edgeTable[cubeIndex] & 32)   vertlist[5] = vertexInterp(p[5], p[6], val[5], val[6]);
                                    if (edgeTable[cubeIndex] & 64)   vertlist[6] = vertexInterp(p[6], p[7], val[6], val[7]);
                                    if (edgeTable[cubeIndex] & 128)  vertlist[7] = vertexInterp(p[7], p[4], val[7], val[4]);
                                    if (edgeTable[cubeIndex] & 256)  vertlist[8] = vertexInterp(p[0], p[4], val[0], val[4]);
                                    if (edgeTable[cubeIndex] & 512)  vertlist[9] = vertexInterp(p[1], p[5], val[1], val[5]);
                                    if (edgeTable[cubeIndex] & 1024) vertlist[10] = vertexInterp(p[2], p[6], val[2], val[6]);
                                    if (edgeTable[cubeIndex] & 2048) vertlist[11] = vertexInterp(p[3], p[7], val[3], val[7]);
                                    
                                    // Create triangles
                                    for (int i = 0; triTable[cubeIndex][i] != -1; i += 3) {
                                        localTriangles.emplace_back(
                                            vertlist[triTable[cubeIndex][i]],
                                            vertlist[triTable[cubeIndex][i+1]],
                                            vertlist[triTable[cubeIndex][i+2]]
                                        );
                                    }
                                }
                            }
                        }
                        return localTriangles;
                    }
                ));
            }

            // Collect results
            for (auto& f : futures) {
                std::vector<Triangle> part = f.get();
                mesh.triangles.insert(mesh.triangles.end(), part.begin(), part.end());
            }
        }
        
    private:
        static Point3 vertexInterp(const Point3& p1, const Point3& p2, double val1, double val2) {
             if (std::abs(val1 - val2) < 1e-6) return p1;
            double mu = (0 - val1) / (val2 - val1);
            return p1 + (p2 - p1) * mu;
        }
    };
}
