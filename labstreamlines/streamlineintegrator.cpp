/*********************************************************************
 *  Author  : Himangshu Saikia, Jiahui Liu
 *  Init    : Tuesday, September 19, 2017 - 15:08:33
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <inviwo/core/interaction/events/mouseevent.h>
#include <inviwo/core/util/utilities.h>
#include <labstreamlines/integrator.h>
#include <labstreamlines/streamlineintegrator.h>
#include <labutils/scalarvectorfield.h>

namespace inviwo
{

// The Class Identifier has to be globally unique. Use a reverse DNS naming
// scheme
const ProcessorInfo StreamlineIntegrator::processorInfo_{
    "org.inviwo.StreamlineIntegrator",// Class identifier
    "Streamline Integrator",            // Display name
    "KTH Lab",                          // Category
    CodeState::Experimental,            // Code state
    Tags::None,                         // Tags
};

const ProcessorInfo StreamlineIntegrator::getProcessorInfo() const
{
    return processorInfo_;
}

StreamlineIntegrator::StreamlineIntegrator()
    : Processor()
    , inData("volIn")
    , meshOut("meshOut")
    , meshBBoxOut("meshBBoxOut")
    , propDisplayPoints("displayPoints", "Display Points", true)
    , propStartPoint("startPoint", "Start Point", vec2(0.5, 0.5), vec2(-1), vec2(1), vec2(0.1))
    , propSeedMode("seedMode", "Seeds")
    , propNumStepsTaken("numstepstaken", "Number of actual steps", 0, 0, 1000)
    // Add the new properties according to task requirements
    , propNumSteps("numSteps", "Number of Steps", 1, 1, 1000)// task 1: Initialize additional properties
    , propStepSize("stepSize", "Step Size", 0.001, 0.001, 1)
    , propDirection("direction", "Direction")
    , propnormalizeVectorField("normalizeVectorField", "Normalize Vector Field", false) 
    , propStopWhenZero("stopWhenZero", "Stop at zeros of the vector field", false)
    , propmaxArcLength("maxArcLength", "Max Arc Length", 0.01, 0.001, 10)
    , propminVelocity("minVelocity", "Velocity Threshold", 0.1, 0.01, 1) 
    , propLineSeeding("lineSeedingMode", "Multipe line seeding mode")
    , propNumberOfStreamLines("numStreamLines", "Number of Stream lines", 1, 0, 1000)
    , propSeedLinesGridX("numStreamLinesX", "X Dim", 1, 0, 50)
    , propSeedLinesGridY("numStreamLinesY", "Y Dim", 1, 0, 50)
    , mouseMoveStart(
          "mouseMoveStart", "Move Start", [this](Event* e) { eventMoveStart(e); }, MouseButton::Left, MouseState::Press | MouseState::Move)

// TODO: Initialize additional properties
// propertyName("propertyIdentifier", "Display Name of the Propery",
// default value (optional), minimum value (optional), maximum value (optional),
// increment (optional)); propertyIdentifier cannot have spaces
{
    // Register Ports
    addPort(inData);
    addPort(meshOut);
    addPort(meshBBoxOut);

    // Register Properties
    propSeedMode.addOption("one", "Single Start Point", 0);
    propSeedMode.addOption("multiple", "Multiple Seeds", 1);
    addProperty(propSeedMode);
    propLineSeeding.addOption("random", "Random", 0);
    propLineSeeding.addOption("uniform", "Uniform", 1);
    propLineSeeding.addOption("magnitude", "Magnitude", 2);
    addProperty(propLineSeeding);
    addProperty(propStartPoint);
    addProperty(propDisplayPoints);
    addProperty(propNumStepsTaken);
    propNumStepsTaken.setReadOnly(true);
    propNumStepsTaken.setSemantics(PropertySemantics::Text);
    addProperty(mouseMoveStart);

    // TODO: Register additional properties
    // addProperty(propertyName);
    // Properties to meet the requirements of Task 4.2
    addProperty(propNumSteps);
    addProperty(propStepSize);
    addProperty(propDirection);
    addProperty(propnormalizeVectorField);
    addProperty(propmaxArcLength);
    addProperty(propminVelocity);
    addProperty(propNumberOfStreamLines);
    addProperty(propSeedLinesGridX);
    addProperty(propSeedLinesGridY);

    addProperty(propStopWhenZero);
    propDirection.addOption("forward", "Forward", 0);
    propDirection.addOption("backward", "Backward", 1);
    propDirection.addOption("both", "Both", 2);
    

    // Show properties for a single seed and hide properties for multiple seeds
    // (TODO)
    util::hide(propNumberOfStreamLines, propLineSeeding, propSeedLinesGridX, propSeedLinesGridY);
    util::hide(propStartPoint, mouseMoveStart, propNumStepsTaken);

    propSeedMode.onChange([this]() {
        if (propSeedMode.get() == 0)
        {
            util::show(propStartPoint, mouseMoveStart, propNumStepsTaken);
            util::hide(propNumberOfStreamLines, propLineSeeding, propSeedLinesGridX, propSeedLinesGridY);
        }
        else
        {
            util::hide(propStartPoint, mouseMoveStart, propNumStepsTaken);
            util::show(propLineSeeding);
            if (propLineSeeding.get() == 0)
            {
                util::hide(propSeedLinesGridX, propSeedLinesGridY);
                util::show(propNumberOfStreamLines);
            }
            else if (propLineSeeding.get() == 1)
            {
                util::show(propSeedLinesGridX, propSeedLinesGridY);
                util::hide(propNumberOfStreamLines);
            }
            else
            {
                util::show(propSeedLinesGridX, propSeedLinesGridY);
                util::hide(propNumberOfStreamLines);
            }
        }
    });
    propLineSeeding.onChange([this]() {
        if (propSeedMode.get() != 0)
        {
            if (propLineSeeding.get() == 0)
            {
                util::hide(propSeedLinesGridX, propSeedLinesGridY);
                util::show(propNumberOfStreamLines);
            }
            else if (propLineSeeding.get() == 1)
            {
                util::show(propSeedLinesGridX, propSeedLinesGridY);
                util::hide(propNumberOfStreamLines);
            }
            else
            {
                util::show(propSeedLinesGridX, propSeedLinesGridY);
                util::hide(propNumberOfStreamLines);
            }
        }
    });
}

void StreamlineIntegrator::eventMoveStart(Event* event)
{
    if (!inData.hasData()) return;
    auto mouseEvent = static_cast<MouseEvent*>(event);
    vec2 mousePos = mouseEvent->posNormalized();

    // Map to bounding box range
    mousePos[0] *= static_cast<float>(BBoxMax_[0] - BBoxMin_[0]);
    mousePos[1] *= static_cast<float>(BBoxMax_[1] - BBoxMin_[1]);
    mousePos += static_cast<vec2>(BBoxMin_);

    // Update starting point
    propStartPoint.set(mousePos);
    event->markAsUsed();
}


int StreamlineIntegrator::drawStreamLine(const vec2& startPoint, const VectorField2& vectorField, const float stepSize, const int direction, const bool directionField, const int nSteps,
                                          const float maxArcLength, const float minSpeed, const bool displayPoints, std::shared_ptr<BasicMesh>& mesh, std::vector<BasicMesh::Vertex>& vertices)
{
    int stepsTaken = 0;// 计数执行的步数
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

    for (int i = 0; i < nSteps; ++i)  //执行流线积分
    {
        // 选择方向进行积分
        if (direction == 0)
        {
            currentPosition = Integrator::RK4(vectorField, currentPosition, stepSize, true, directionField);
        }
        else
        {
            currentPosition = Integrator::RK4(vectorField, currentPosition, stepSize, false, directionField);
        }

        // 边界检查
        if (currentPosition.x < BBoxMin.x || currentPosition.x > BBoxMax.x || currentPosition.y < BBoxMin.y || currentPosition.y > BBoxMax.y)
        {
            break;
        }

        // 停止条件：根据速度
        dvec2 velocity = vectorField.interpolate(currentPosition);
        if (glm::length(velocity) == 0.0 || glm::length(velocity) < minSpeed)
        {
            break;
        }

        // 弧长检查
        double segmentLength = glm::length(currentPosition - previousPosition);
        arcLength += segmentLength;
        if (arcLength > maxArcLength)
        {
            break;
        }

        // 绘制线段
        Integrator::drawLineSegment(currentPosition, previousPosition, vec4(0, 0, 0, 1), indexBufferLines.get(), vertices);

        if (displayPoints)
        {
            Integrator::drawPoint(currentPosition, vec4(0, 0, 0, 1), indexBufferPoints.get(), vertices);
        }

        previousPosition = currentPosition;
        stepsTaken++;
    }
    return stepsTaken;  //返回执行步数
}

    
void StreamlineIntegrator::process()
{
    // Get input
    if (!inData.hasData())
    {
        return;
    }
    auto vol = inData.getData();
    if (!vol)
    {
        std::cerr << "Error: inData has no data" << std::endl;
        return;
    }

    // Retreive data in a form that we can access it
    auto vectorField = VectorField2::createFieldFromVolume(vol);
    BBoxMin_ = vectorField.getBBoxMin();
    BBoxMax_ = vectorField.getBBoxMax();

    // 检查BBox是否有效
    if (BBoxMin_.x >= BBoxMax_.x || BBoxMin_.y >= BBoxMax_.y)
    {
        std::cerr << "Error: Invalid bounding box dimensions" << std::endl;
        return;
    }

    // The start point should be inside the volume (set maximum to the upper right corner)
    propStartPoint.setMinValue(BBoxMin_ - dvec2(1, 1));
    propStartPoint.setMaxValue(BBoxMax_ + dvec2(1, 1));
    
    // Initialize mesh, vertices and index buffers for the two streamlines and the points
    auto bboxMesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> bboxVertices;


    // Make bounding box without vertex duplication, instead of line segments which duplicate
    // vertices, create line segments between each added points with connectivity type of the index
    // buffer
    auto indexBufferBBox = bboxMesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
    // Bounding Box vertex 0
    vec4 black = vec4(0, 0, 0, 1);
    vec4 blue = vec4(0, 0, 1, 1);
    Integrator::drawNextPointInPolyline(BBoxMin_, black, indexBufferBBox.get(), bboxVertices);
    Integrator::drawNextPointInPolyline(vec2(BBoxMin_[0], BBoxMax_[1]), black,
                                        indexBufferBBox.get(), bboxVertices);
    Integrator::drawNextPointInPolyline(BBoxMax_, black, indexBufferBBox.get(), bboxVertices);
    Integrator::drawNextPointInPolyline(vec2(BBoxMax_[0], BBoxMin_[1]), black,
                                        indexBufferBBox.get(), bboxVertices);
    // Connect back to the first point, to make a full rectangle
    indexBufferBBox->add(static_cast<std::uint32_t>(0));
    bboxMesh->addVertices(bboxVertices);
    meshBBoxOut.setData(bboxMesh);

    auto mesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> vertices;

    auto indexBufferRK = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
    auto indexBufferRK2 = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
    auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);
    auto indexBufferPoints2 = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);

    int totalStepsTaken = 0;// 初始化实际执行的步数

    // get input of propNumSteps and propStepSize
    auto numSteps = propNumSteps.get();
    auto stepSize = propStepSize.get();
    auto direction = propDirection.get();
    auto if_norm = propnormalizeVectorField.get();


    // If seed mode is 0, we have a single stream line
    if (propSeedMode.get() == 0)
    {
        vec2 startPoint = propStartPoint.get();
        // Draw start point
        if (propDisplayPoints.get() != 0)
            Integrator::drawPoint(startPoint, vec4(0, 0, 0, 1), indexBufferPoints.get(), vertices);

        // 前向和后向积分
        if (direction == 2)
        {
            totalStepsTaken += drawStreamLine(startPoint, vectorField, stepSize, 0, if_norm, numSteps, propmaxArcLength.get(), propminVelocity.get(), propDisplayPoints.get(), mesh, vertices);
            totalStepsTaken += drawStreamLine(startPoint, vectorField, stepSize, 1, if_norm, numSteps, propmaxArcLength.get(), propminVelocity.get(), propDisplayPoints.get(), mesh, vertices);
        }
        else
        {
            totalStepsTaken += drawStreamLine(startPoint, vectorField, stepSize, direction, if_norm, numSteps, propmaxArcLength.get(), propminVelocity.get(), propDisplayPoints.get(), mesh, vertices);
        }
    } 
  
    else if (propSeedMode.get() == 1)// Multiple seed mode
        {
        if (propLineSeeding.get() == 0)// Random seeding
        {
            for (int i = 0; i < propNumberOfStreamLines.get(); ++i)
            {
                // 得到一个在0到1之间的浮点数
                double randomX = BBoxMin_.x + static_cast<double>(rand()) / RAND_MAX * (BBoxMax_.x - BBoxMin_.x);
                double randomY = BBoxMin_.y + static_cast<double>(rand()) / RAND_MAX * (BBoxMax_.y - BBoxMin_.y);
                vec2 seedPoint = vec2(randomX, randomY);

                 // 前向和后向积分
                if (propDirection.get() == 2)
                {
                    totalStepsTaken += drawStreamLine(seedPoint, vectorField, propStepSize.get(), 0, propnormalizeVectorField.get(), propNumSteps.get(), propmaxArcLength.get(), propminVelocity.get(),
                                                        propDisplayPoints.get(), mesh, vertices);
                    totalStepsTaken += drawStreamLine(seedPoint, vectorField, propStepSize.get(), 1, propnormalizeVectorField.get(), propNumSteps.get(), propmaxArcLength.get(), propminVelocity.get(),
                                                        propDisplayPoints.get(), mesh, vertices);
                }
                else
                {
                    totalStepsTaken += drawStreamLine(seedPoint, vectorField, propStepSize.get(), propDirection.get(), propnormalizeVectorField.get(), propNumSteps.get(), propmaxArcLength.get(),
                                                        propminVelocity.get(), propDisplayPoints.get(), mesh, vertices);
                }
            }
        }
        else if (propLineSeeding.get() == 1)// Uniform grid seeding
        {
            int gridX = propSeedLinesGridX.get();
            int gridY = propSeedLinesGridY.get();

            if (gridX <= 0 || gridY <= 0)
            {
                std::cerr << "Error: Invalid grid dimensions" << std::endl;
                return;
            }

            //计算步长，计算X和Y方向上的步长。步长是边界框宽度和高度除以网格的种子点数量，用于确定种子点在网格上的位置。
            double stepX = (BBoxMax_.x - BBoxMin_.x) / static_cast<double>(gridX);
            double stepY = (BBoxMax_.y - BBoxMin_.y) / static_cast<double>(gridY);

            // 遍历网格的每个点，计算种子点的坐标
            for (int i = 0; i <= gridX; ++i)
            {
                for (int j = 0; j <= gridY; ++j)
                {
                    vec2 seedPoint = vec2(BBoxMin_.x + i * stepX, BBoxMin_.y + j * stepY);

                    // 前向和后向积分
                    if (direction == 2)
                    {
                        totalStepsTaken += drawStreamLine(seedPoint, vectorField, stepSize, 0, if_norm, numSteps, propmaxArcLength.get(), propminVelocity.get(), propDisplayPoints.get(), mesh, vertices);
                        totalStepsTaken += drawStreamLine(seedPoint, vectorField, stepSize, 1, if_norm, numSteps, propmaxArcLength.get(), propminVelocity.get(), propDisplayPoints.get(), mesh, vertices);
                    }
                    else
                    {
                        totalStepsTaken += drawStreamLine(seedPoint, vectorField, stepSize, direction, if_norm, numSteps, propmaxArcLength.get(), propminVelocity.get(), propDisplayPoints.get(), mesh, vertices);
                    }
                }
            }
        }
        else if (propLineSeeding.get() == 2)
        {// Magnitude-based seeding
            int numSeeds = propNumberOfStreamLines.get();
            std::vector<std::pair<vec2, double>> seeds;

            int gridX = propSeedLinesGridX.get();
            int gridY = propSeedLinesGridY.get();

            if (gridX <= 0 || gridY <= 0)
            {
                std::cerr << "Error: Invalid grid dimensions" << std::endl;
                return;
            }

            double stepX = (BBoxMax_.x - BBoxMin_.x) / static_cast<double>(gridX);
            double stepY = (BBoxMax_.y - BBoxMin_.y) / static_cast<double>(gridY);

            double maxMagnitude = 0.0;

            // 遍历整个向量场，计算幅值，并找到最大幅值
            for (double x = BBoxMin_.x; x <= BBoxMax_.x; x += stepX)
            {
                for (double y = BBoxMin_.y; y <= BBoxMax_.y; y += stepY)
                {
                    vec2 point(x, y);
                    double magnitude = glm::length(vectorField.interpolate(point));
                    seeds.push_back({point, magnitude});
                    maxMagnitude = std::max(maxMagnitude, magnitude);// 记录最大幅值
                }
            }

            // 根据幅值对种子点进行概率性采样
            for (int i = 0; i < seeds.size(); ++i)
            {
                vec2 seedPoint = seeds[i].first;
                double magnitude = seeds[i].second;

                // 根据幅值的归一化值来进行概率采样
                double probability = magnitude / maxMagnitude;              // 幅值归一化为 [0, 1]
                double randomValue = static_cast<double>(rand()) / RAND_MAX;// 生成 [0, 1] 之间的随机数

                // 仅当随机值小于幅值归一化值时才生成种子点
                if (randomValue <= probability)
                {
                    // 使用 drawStreamLine 进行流线积分和绘制
                    if (propDirection.get() == 2)
                    {// 双向积分
                        totalStepsTaken += drawStreamLine(seedPoint, vectorField, propStepSize.get(), 0, propnormalizeVectorField.get(), propNumSteps.get(), propmaxArcLength.get(),
                                                          propminVelocity.get(), propDisplayPoints.get(), mesh, vertices);
                        totalStepsTaken += drawStreamLine(seedPoint, vectorField, propStepSize.get(), 1, propnormalizeVectorField.get(), propNumSteps.get(), propmaxArcLength.get(),
                                                          propminVelocity.get(), propDisplayPoints.get(), mesh, vertices);
                    }
                    else
                    {// 单向积分
                        totalStepsTaken += drawStreamLine(seedPoint, vectorField, propStepSize.get(), propDirection.get(), propnormalizeVectorField.get(), propNumSteps.get(), propmaxArcLength.get(),
                                                          propminVelocity.get(), propDisplayPoints.get(), mesh, vertices);
                    }
                }
            }
        }
    }

    //设置实际执行步数
    propNumStepsTaken.set(totalStepsTaken);

    //添加顶点并输出网络
    mesh->addVertices(vertices);
    meshOut.setData(mesh);
}

}// namespace inviwo
