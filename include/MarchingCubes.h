#pragma once

#include "SDF.h"
#include "Mesh.h"
#include "MCTables.h"
#include <vector>
#include <cmath>
#include <iostream>

namespace Geom {

    class MarchingCubes {
    public:
        static void march(const SDF& sdf, const BoundingBox& bounds, double resolution, Mesh& mesh) {
            if (resolution <= EPSILON) return;

            // Calculate number of steps
            Vec3 size = bounds.size();
            int nx = static_cast<int>(std::ceil(size.x / resolution));
            int ny = static_cast<int>(std::ceil(size.y / resolution));
            int nz = static_cast<int>(std::ceil(size.z / resolution));

            // Iterate over the grid
            // Note: We need samples at corners. So we need (nx+1) * (ny+1) * (nz+1) samples if we cache.
            // For MVP, simple re-evaluation is fine (slower but simpler).
            
            for (int z = 0; z < nz; ++z) {
                for (int y = 0; y < ny; ++y) {
                    for (int x = 0; x < nx; ++x) {
                        Point3 p0 = bounds.min + Vec3(x * resolution, y * resolution, z * resolution);
                        
                        // 8 corners
                        Point3 p[8];
                        double val[8];
                        
                        // 0: 0,0,0
                        // 1: 1,0,0
                        // 2: 1,0,1 (!) -> check standard indexing
                        // Standard MC indexing:
                        // 0: 0,0,0
                        // 1: 1,0,0
                        // 2: 1,1,0
                        // 3: 0,1,0
                        // 4: 0,0,1
                        // 5: 1,0,1
                        // 6: 1,1,1
                        // 7: 0,1,1
                        
                        // Paul Bourke indexing:
                        // 0: 0,0,0
                        // 1: 1,0,0
                        // 2: 1,0,1
                        // 3: 0,0,1
                        // ... wait, check edge table consistency.
                        
                        // Standard vtk/Paul Bourke:
                        // v0: 0,0,0
                        // v1: 1,0,0
                        // v2: 1,1,0
                        // v3: 0,1,0
                        // v4: 0,0,1
                        // v5: 1,0,1
                        // v6: 1,1,1
                        // v7: 0,1,1
                        
                        p[0] = p0;
                        p[1] = p0 + Vec3(resolution, 0, 0);
                        p[2] = p0 + Vec3(resolution, resolution, 0);
                        p[3] = p0 + Vec3(0, resolution, 0);
                        p[4] = p0 + Vec3(0, 0, resolution);
                        p[5] = p0 + Vec3(resolution, 0, resolution);
                        p[6] = p0 + Vec3(resolution, resolution, resolution);
                        p[7] = p0 + Vec3(0, resolution, resolution);
                        
                        for(int i=0; i<8; ++i) val[i] = sdf.eval(p[i]);
                        
                        // Determine cube index
                        int cubeIndex = 0;
                        if (val[0] < 0) cubeIndex |= 1;
                        if (val[1] < 0) cubeIndex |= 2;
                        if (val[2] < 0) cubeIndex |= 4;
                        if (val[3] < 0) cubeIndex |= 8;
                        if (val[4] < 0) cubeIndex |= 16;
                        if (val[5] < 0) cubeIndex |= 32;
                        if (val[6] < 0) cubeIndex |= 64;
                        if (val[7] < 0) cubeIndex |= 128; // Standard typically 128 for corner 7?
                        // Let's verify standard: 
                        // bit 0 -> v0, bit 1 -> v1, etc.
                        
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
                            mesh.addTriangle(
                                vertlist[triTable[cubeIndex][i]],
                                vertlist[triTable[cubeIndex][i+1]],
                                vertlist[triTable[cubeIndex][i+2]]
                            );
                        }
                    }
                }
            }
        }
        
    private:
        static Point3 vertexInterp(const Point3& p1, const Point3& p2, double val1, double val2) {
            double mu = (0 - val1) / (val2 - val1);
            return p1 + (p2 - p1) * mu;
        }
    };
}
