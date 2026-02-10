#pragma once

#include "Mesh.h"
#include <fstream>
#include <string>
#include <iostream>

namespace Geom {

    class STLExporter {
    public:
        static bool save(const Mesh& mesh, const std::string& filename) {
            std::ofstream file(filename, std::ios::binary);
            if (!file) return false;

            // 80 byte header
            char header[80] = "SDF Geometry Kernel Export";
            file.write(header, 80);

            // Number of triangles (uint32)
            uint32_t numTriangles = static_cast<uint32_t>(mesh.triangles.size());
            file.write(reinterpret_cast<const char*>(&numTriangles), sizeof(uint32_t));

            // Triangles
            for (const auto& tri : mesh.triangles) {
                // Normal (3 floats)
                float n[3] = { (float)tri.normal.x, (float)tri.normal.y, (float)tri.normal.z };
                file.write(reinterpret_cast<const char*>(n), 12);

                // Vertices (3 * 3 floats)
                float v1[3] = { (float)tri.v0.x, (float)tri.v0.y, (float)tri.v0.z };
                float v2[3] = { (float)tri.v1.x, (float)tri.v1.y, (float)tri.v1.z };
                float v3[3] = { (float)tri.v2.x, (float)tri.v2.y, (float)tri.v2.z };

                file.write(reinterpret_cast<const char*>(v1), 12);
                file.write(reinterpret_cast<const char*>(v2), 12);
                file.write(reinterpret_cast<const char*>(v3), 12);

                // Attribute byte count (uint16) - usually 0
                uint16_t attr = 0;
                file.write(reinterpret_cast<const char*>(&attr), 2);
            }

            return true;
        }
    };
}
