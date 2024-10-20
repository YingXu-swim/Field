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
    bool divide_by_0 = false;

    // Step 1: Compute k1 (slope at the beginning of the interval)
    auto k1 = vectorField.interpolate(position);
    if (!propDirection)
    {
        k1 = -k1;
    }
    if (!propnormalizeVectorField)
    {
        if (pow(k1.x, 2) + pow(k1.y, 2) != 0)
        {
            k1 = k1 / sqrt(pow(k1.x, 2) + pow(k1.y, 2));
        }
        else
        {
            divide_by_0 = true;
        }
    }

    // Step 2: Compute k2 (slope at the midpoint, using k1)
    dvec2 midPoint1 = position + 0.5 * stepSize * k1;
    auto k2 = vectorField.interpolate(midPoint1);
    if (!propDirection)
    {
        k2 = -k2;
    }
    if (!propnormalizeVectorField)
    {
        if (pow(k2.x, 2) + pow(k2.y, 2) != 0)
        {
            k2 = k2 / sqrt(pow(k2.x, 2) + pow(k2.y, 2));
        }
        else
        {
            divide_by_0 = true;
        }
    }
    // Step 3: Compute k3 (another slope at the midpoint, using k2)
    dvec2 midPoint2 = position + 0.5 * stepSize * k2;
    auto k3 = vectorField.interpolate(midPoint2);
    if (!propDirection)
    {
        k3 = -k3;
    }
    if (!propnormalizeVectorField)
    {
        if (pow(k3.x, 2) + pow(k3.y, 2) != 0)
        {
            k3 = k3 / sqrt(pow(k3.x, 2) + pow(k3.y, 2));
        }
        else
        {
            divide_by_0 = true;
        }
    }
    // Step 4: Compute k4 (slope at the end of the interval, using k3)
    dvec2 endPoint = position + stepSize * k3;
    auto k4 = vectorField.interpolate(endPoint);
    if (!propDirection)
    {
        k4 = -k4;
    }
    if (!propnormalizeVectorField)
    {
        if (pow(k4.x, 2) + pow(k4.y, 2) != 0)
        {
            k4 = k4 / sqrt(pow(k4.x, 2) + pow(k4.y, 2));
        }
        else
        {
            divide_by_0 = true;
        }
    }
    // Final step: Combine all the slopes to estimate the next position
    dvec2 nextPosition = position + (stepSize / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
    if (divide_by_0)
    {
        return position;
    }
    else
    {
        return nextPosition;
    }
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

/*
dvec2 Integrator::drawStreamLine(const vec2& startPoint, const VectorField2& vectorField, const float stepSize, const int direction, const bool directionField, const int nSteps,
                                         const float maxArcLength, const float minSpeed, const bool displayPoints, std::shared_ptr<BasicMesh>& mesh, std::vector<BasicMesh::Vertex>& vertices)
{
    int stepsTaken = 0;// ����ִ�еĲ���
    vec2 currentPosition = startPoint;
    vec2 previousPosition = startPoint;

    vec2 BBoxMin = vectorField.getBBoxMin();
    vec2 BBoxMax = vectorField.getBBoxMax();

    auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);
    auto indexBufferLines = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);

    double arcLength = 0.0;

    if (displayPoints)
    {
        Integrator::drawPoint(startPoint, vec4(0, 0, 0, 1), indexBufferPoints.get(), vertices);
    }

    for (int i = 0; i < nSteps; ++i)//ִ�����߻���
    {
        // ѡ������л���
        if (direction == 0)
        {
            currentPosition = Integrator::RK4(vectorField, currentPosition, stepSize, true, directionField);
        }
        else
        {
            currentPosition = Integrator::RK4(vectorField, currentPosition, stepSize, false, directionField);
        }

        // �߽���
        if (currentPosition.x < BBoxMin.x || currentPosition.x > BBoxMax.x || currentPosition.y < BBoxMin.y || currentPosition.y > BBoxMax.y)
        {
            break;
        }

        // ֹͣ�����������ٶ�
        dvec2 velocity = vectorField.interpolate(currentPosition);
        if (glm::length(velocity) == 0.0 || glm::length(velocity) < minSpeed)
        {
            break;
        }

        // �������
        double segmentLength = glm::length(currentPosition - previousPosition);
        arcLength += segmentLength;
        if (arcLength > maxArcLength)
        {
            break;
        }

        // �����߶�
        Integrator::drawLineSegment(currentPosition, previousPosition, vec4(0, 0, 0, 1), indexBufferLines.get(), vertices);

        if (displayPoints)
        {
            Integrator::drawPoint(currentPosition, vec4(0, 0, 0, 1), indexBufferPoints.get(), vertices);
        }

        previousPosition = currentPosition;
        stepsTaken++;
    }
    return currentPosition; //����ִ�в���
}
*/
}  // namespace inviwo
