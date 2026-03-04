#pragma once

#include "Mesh.h"
#include <string>
#include <fstream>
#include <vector>

namespace Geom {

    class STLImporter {
    public:
        static Mesh load(const std::string& path) {
            Mesh mesh;
            std::ifstream file(path, std::ios::binary);
            if (!file) return mesh;

            char header[80];
            file.read(header, 80);

            uint32_t numTriangles;
            file.read(reinterpret_cast<char*>(&numTriangles), 4);

            mesh.triangles.reserve(numTriangles);

            for (uint32_t i = 0; i < numTriangles; ++i) {
                float n[3], v1[3], v2[3], v3[3];
                uint16_t attr;

                file.read(reinterpret_cast<char*>(n), 12);
                file.read(reinterpret_cast<char*>(v1), 12);
                file.read(reinterpret_cast<char*>(v2), 12);
                file.read(reinterpret_cast<char*>(v3), 12);
                file.read(reinterpret_cast<char*>(&attr), 2);

                mesh.triangles.push_back({
                    Point3(v1[0], v1[1], v1[2]),
                    Point3(v2[0], v2[1], v2[2]),
                    Point3(v3[0], v3[1], v3[2])
                });
            }

            return mesh;
        }
    };

}
