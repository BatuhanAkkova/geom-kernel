#include "../include/Geometry.h"
#include "../include/SDF.h"
#include "../include/Primitives.h"

#include <iostream>
#include <iomanip>

using namespace Geom;

void print_issues(const std::string& name, const std::vector<ValidationIssue>& issues) {
    std::cout << "--- Validation Results: " << name << " ---" << std::endl;
    if (issues.empty()) {
        std::cout << "No issues found." << std::endl;
        return;
    }
    
    int overhangs = 0, curvatures = 0, thins = 0;
    for (const auto& issue : issues) {
        if (issue.type == IssueType::Overhang) overhangs++;
        else if (issue.type == IssueType::HighCurvature) curvatures++;
        else if (issue.type == IssueType::ThinFeature) thins++;
    }

    std::cout << "Found " << issues.size() << " issues:" << std::endl;
    std::cout << " - Overhangs: " << overhangs << std::endl;
    std::cout << " - Curvature violations: " << curvatures << std::endl;
    std::cout << " - Thin features: " << thins << std::endl;

    // Show a few samples
    for (size_t i = 0; i < std::min(issues.size(), size_t(5)); ++i) {
        const auto& iss = issues[i];
        std::cout << "  At (" << iss.position.x << ", " << iss.position.y << ", " << iss.position.z << "): " 
                  << iss.description << " (Value: " << std::fixed << std::setprecision(2) << iss.value << ")" << std::endl;
    }
    if (issues.size() > 5) std::cout << "  ... and " << (issues.size() - 5) << " more." << std::endl;
    std::cout << std::endl;
}

int main() {
    ManufacturingConfig config;
    config.maxOverhangAngle = 45.0; 
    config.minCurvatureRadius = 5.0; 
    config.minBeadWidth = 3.0;
    
    ManufacturingValidator validator(config);

    // 1. A Sphere of Radius 10
    // Expected: Overhang at the bottom. Curvature R=10 (OK). No thin features (OK).
    {
        auto sphere = std::make_shared<Sphere>(Point3(0, 0, 10), 10.0);
        auto issues = validator.validate(sphere, sphere->boundingBox(), 1.0);
        print_issues("Sphere R=10", issues);
    }

    // 2. A Tiny Sphere of Radius 2
    // Expected: High Curvature (R=2 < 5). Overhang at bottom.
    {
        auto tinySphere = std::make_shared<Sphere>(Point3(20, 0, 10), 2.0);
        auto issues = validator.validate(tinySphere, tinySphere->boundingBox(), 0.5);
        print_issues("Tiny Sphere R=2", issues);
    }

    // 3. A Thin Box (Width 1.0)
    // Expected: Thin Feature (1.0 < 3.0). Overhang at top/bottom depends on orientation.
    // Box center (40,0,10), half-extents (0.5, 5, 5) -> thickness 1.0 along X.
    {
        auto thinBox = std::make_shared<Box>(Point3(40, 0, 10), Vec3(0.5, 5, 5));
        auto issues = validator.validate(thinBox, thinBox->boundingBox(), 0.5);
        print_issues("Thin Box (thick=1)", issues);
    }

    // 4. A Horizontal Plane/Plate
    // Expected: Massive overhang at the bottom (angle approx 0).
    {
        auto plate = std::make_shared<Box>(Point3(0, 30, 10), Vec3(10, 10, 1));
        auto issues = validator.validate(plate, plate->boundingBox(), 1.0);
        print_issues("Thin Horizontal Plate", issues);
    }

    return 0;
}
