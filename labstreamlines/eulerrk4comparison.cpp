/*********************************************************************
 *  Author  : Himangshu Saikia, Wiebke Koepp
 *  Init    : Tuesday, September 19, 2017 - 15:08:24
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <inviwo/core/datastructures/geometry/basicmesh.h>
#include <inviwo/core/algorithm/boundingbox.h>
#include <inviwo/core/interaction/events/mouseevent.h>
#include <labstreamlines/eulerrk4comparison.h>
#include <labstreamlines/integrator.h>

namespace inviwo {

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo EulerRK4Comparison::processorInfo_{
    "org.inviwo.EulerRK4Comparison",  // Class identifier
    "Euler RK4 Comparison",           // Display name
    "KTH Lab",                        // Category
    CodeState::Experimental,          // Code state
    Tags::None,                       // Tags
};

const ProcessorInfo EulerRK4Comparison::getProcessorInfo() const { return processorInfo_; }

EulerRK4Comparison::EulerRK4Comparison() // 构造函数
    : Processor()
    , inData("inData")
    , meshOut("meshOut")
    , meshBBoxOut("meshBBoxOut")
    , propStartPoint("startPoint", "Start Point", vec2(0.5, 0.5), vec2(0), vec2(1024), vec2(0.5))
    , mouseMoveStart(
          "mouseMoveStart", "Move Start", [this](Event* e) { eventMoveStart(e); },
          MouseButton::Left, MouseState::Press | MouseState::Move)
    ,propNumSteps("numSteps", "Number of Steps", 1000, 1, 100000) // task 1: Initialize additional properties
    , propStepSize("stepSize", "Step Size", 0.01f, 0.001f, 100.0f)// task 1: Initialize additional properties
    , propMehtod("Euler", "Runge-Kutta (4th order)")
// propertyName("propertyIdentifier", "Display Name of the Property",
// default value (optional), minimum value (optional), maximum value (optional), increment
// (optional)); propertyIdentifier cannot have spaces

{
    // Register Ports
    addPort(meshOut);
    addPort(meshBBoxOut);
    addPort(inData);

    addProperties(propMehtod);
    propMehtod.addOption("euler", "Euler", 0);
    propMehtod.addOption("runge_Kutta_4order", "Runge_Kutta_4order", 1);
    
    // Register Properties
    addProperty(propStartPoint);
    addProperty(mouseMoveStart);

    // TODO: Register additional 
    // addProperty(propertyName);
    // task 1
    addProperty(propNumSteps);
    addProperty(propStepSize);


}
// eventMoveStart函数处理鼠标移动事件.它将鼠标的位置映射到数据的边界框范围内，并更新起始点的位置。
void EulerRK4Comparison::eventMoveStart(Event* event) {
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

void EulerRK4Comparison::process() {
    // Get input
    if (!inData.hasData()) {
        return;
    }
    auto vol = inData.getData();

    // Retrieve data in a form that we can access it 
    const VectorField2 vectorField = VectorField2::createFieldFromVolume(vol);
    BBoxMin_ = vectorField.getBBoxMin();
    BBoxMax_ = vectorField.getBBoxMax();

    // The start point should be inside the volume (set maximum to the upper right corner)
    propStartPoint.setMinValue(BBoxMin_ - dvec2(1, 1));
    propStartPoint.setMaxValue(BBoxMax_ + dvec2(1, 1));

    // Initialize mesh, vertices and index buffers for the two streamlines and the points
    auto mesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> vertices;

    auto indexBufferEuler = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
    auto indexBufferRK = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
    auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);

    auto bboxMesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> bboxVertices;

    // Make bounding box without vertex duplication, instead of line segments which duplicate
    // vertices, create line segments between each added points with connectivity type of the index
    // buffer
    auto indexBufferBBox = bboxMesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
    // Bounding Box vertex 0
    // task 1: define color
    // note: 
    // R,G,B,A: A is transparent
    vec4 black = vec4(0, 0, 0, 1);
    vec4 red = vec4(1, 0, 0, 1);
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

    // Draw start point
    dvec2 startPoint = propStartPoint.get();
    Integrator::drawPoint(startPoint, black, indexBufferPoints.get(), vertices);

    // TODO: Implement the Euler and Runge-Kutta of 4th order integration schemes
    // and then integrate forward for a specified number of integration steps and a given stepsize
    // (these should be additional properties of the processor)

    // get input of propNumSteps and propStepSize
    auto numSteps = propNumSteps.get();
    auto stepSize = propStepSize.get();

    // task 1
    dvec2 lastPosition = startPoint;
    if (propMehtod.get() == 0)// Euler
    {
        for (size_t i = 0; i < numSteps; ++i)
        {
            dvec2 currentPosition = Integrator::Euler(vectorField, lastPosition, stepSize);
            // 检查新位置是否在向量场的边界内
            if (currentPosition.x < BBoxMin_.x || currentPosition.x > BBoxMax_.x || currentPosition.y < BBoxMin_.y || currentPosition.y > BBoxMax_.y)
            {
                break;
            }
            // 绘制
             Integrator::drawPoint(lastPosition, red, indexBufferPoints.get(), vertices);
             Integrator::drawPoint(currentPosition, red, indexBufferPoints.get(), vertices);
             Integrator::drawPoint(lastPosition, red, indexBufferEuler.get(), vertices);
             Integrator::drawPoint(currentPosition, red, indexBufferEuler.get(), vertices);

            lastPosition = currentPosition;
        }    
    }
    else
    {
        for (size_t i = 0; i < numSteps; ++i)
        {
            dvec2 currentPosition = Integrator::RK4(vectorField, lastPosition, stepSize, false, false);
            // 检查新位置是否在向量场的边界内
            if (currentPosition.x < BBoxMin_.x || currentPosition.x > BBoxMax_.x || currentPosition.y < BBoxMin_.y || currentPosition.y > BBoxMax_.y)
            {
                break;
            }
            // 绘制
            Integrator::drawPoint(lastPosition, blue, indexBufferPoints.get(), vertices);
            Integrator::drawPoint(currentPosition, blue, indexBufferPoints.get(), vertices);
            Integrator::drawPoint(lastPosition, blue, indexBufferRK.get(), vertices);
            Integrator::drawPoint(currentPosition, blue, indexBufferRK.get(), vertices);

            lastPosition = currentPosition;
        }
    }


    // Integrator::Rk4(vectorField, dims, startPoint, ...);

    mesh->addVertices(vertices);
    meshOut.setData(mesh);
}

}  // namespace inviwo
