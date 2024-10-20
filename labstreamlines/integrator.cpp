/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Wednesday, September 20, 2017 - 12:04:15
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <labstreamlines/integrator.h>

namespace inviwo {

// TODO: Implement a single integration step here

dvec2 Integrator::Euler(const VectorField2& vectorField, const dvec2& position, double stepSize)
{
    // Access the vector field with vectorField.interpolate(...)
    dvec2 currentPosition = position;
    auto velocity = vectorField.interpolate(position);

    currentPosition.x += stepSize * velocity.x;
    currentPosition.y += stepSize * velocity.y;
        
    return currentPosition;
}

// Runge-Kutta 4th Order (RK4) integration
dvec2 Integrator::RK4(const VectorField2& vectorField, const dvec2& position, const double stepSize, bool propDirection, bool propnormalizeVectorField)
{
    // Step 1: Compute k1 (slope at the beginning of the interval)
    auto k1 = vectorField.interpolate(position);
    if (glm::length(k1) == 0)
    {
        // If the vector field at this position is zero, return the current position (no movement)
        return position;
    }

    if (!propDirection)
    {
        k1 = -k1;
    }
    if (propnormalizeVectorField)
    {
        k1 = glm::normalize(k1);
    }

    // Step 2: Compute k2 (slope at the midpoint, using k1)
    dvec2 midPoint1 = position + 0.5 * stepSize * k1;
    auto k2 = vectorField.interpolate(midPoint1);
    if (glm::length(k2) == 0)
    {
        // If the vector field at this midpoint is zero, return the current position (no movement)
        return position;
    }

    if (!propDirection)
    {
        k2 = -k2;
    }
    if (propnormalizeVectorField)
    {
        k2 = glm::normalize(k2);
    }

    // Step 3: Compute k3 (another slope at the midpoint, using k2)
    dvec2 midPoint2 = position + 0.5 * stepSize * k2;
    auto k3 = vectorField.interpolate(midPoint2);
    if (glm::length(k3) == 0)
    {
        // If the vector field at this midpoint is zero, return the current position (no movement)
        return position;
    }

    if (!propDirection)
    {
        k3 = -k3;
    }
    if (propnormalizeVectorField)
    {
        k3 = glm::normalize(k3);
    }

    // Step 4: Compute k4 (slope at the end of the interval, using k3)
    dvec2 endPoint = position + stepSize * k3;
    auto k4 = vectorField.interpolate(endPoint);
    if (glm::length(k4) == 0)
    {
        // If the vector field at the endpoint is zero, return the current position (no movement)
        return position;
    }

    if (!propDirection)
    {
        k4 = -k4;
    }
    if (propnormalizeVectorField)
    {
        k4 = glm::normalize(k4);
    }

    // Final step: Combine all the slopes to estimate the next position
    dvec2 nextPosition = position + (stepSize / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);

    return nextPosition;
}

void Integrator::drawPoint(const dvec2& p, const vec4& color, IndexBufferRAM* indexBuffer,
                           std::vector<BasicMesh::Vertex>& vertices) {
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    vertices.push_back({vec3(p[0], p[1], 0), vec3(0, 0, 1), vec3(p[0], p[1], 0), color});
}

// Alias for draw point
void Integrator::drawNextPointInPolyline(const dvec2& p, const vec4& color,
                                         IndexBufferRAM* indexBuffer,
                                         std::vector<BasicMesh::Vertex>& vertices) {
    Integrator::drawPoint(p, color, indexBuffer, vertices);
}

void Integrator::drawLineSegment(const dvec2& v1, const dvec2& v2, const vec4& color,
                                 IndexBufferRAM* indexBuffer,
                                 std::vector<BasicMesh::Vertex>& vertices) {
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    vertices.push_back({vec3(v1[0], v1[1], 0), vec3(0, 0, 1), vec3(v1[0], v1[1], 0), color});
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    vertices.push_back({vec3(v2[0], v2[1], 0), vec3(0, 0, 1), vec3(v2[0], v2[1], 0), color});
}


std::vector<dvec2> Integrator::StreamLines(const VectorField2& vectorField, const dvec2& position,
            const double stepSize, bool backwardDirection,
            bool integrateDirectionField, double speedThreshold, int kernelSize) {
        
        std::vector<dvec2> streamLinePoints;
        double speed = DBL_MAX;
        dvec2 pos = position;
        dvec2 new_pos;
        int count = 0;
        
        while (vectorField.isInside(pos) && speed > speedThreshold && count <= kernelSize) {
            streamLinePoints.push_back(pos);
            
            new_pos = Integrator::RK4(vectorField, pos, stepSize, backwardDirection, integrateDirectionField);

            speed = sqrt(pow(new_pos.x - pos.x, 2) + pow(new_pos.y - pos.y, 2)) / double(stepSize);
            pos = new_pos;
            
            count++;
        }
        return streamLinePoints;
     }


}  // namespace inviwo

