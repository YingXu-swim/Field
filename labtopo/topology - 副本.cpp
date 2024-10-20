/*********************************************************************
 *  Author  : Anke Friederici
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 **********************************************************************/

#include <inviwo/core/datastructures/geometry/basicmesh.h>
#include <inviwo/core/datastructures/volume/volumeram.h>
#include <labstreamlines/integrator.h>
#include <labutils/scalarvectorfield.h>
#include <labtopo/topology.h>
#include <labtopo/utils/gradients.h>
#include <queue>

namespace inviwo
{

const vec4 Topology::ColorsCP[6] = {
    vec4(1, 0.95, 0.5, 1), // Saddle - Yellow
    vec4(0.5, 0.5, 0.9, 1),// AttractingNode - Sink - Blue
    vec4(0.9, 0.5, 0.5, 1),// RepellingNode - Source - Red
    vec4(0.5, 0, 0.9, 1),  // AttractingFocus - Purple
    vec4(0.9, 0.5, 0.0, 1),// RepellingFocus - Orange
    vec4(0.3, 0.6, 0.3, 1) // Center - Green
};

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo Topology::processorInfo_{
    "org.inviwo.Topology",  // Class identifier
    "Vector Field Topology",// Display name
    "KTH Lab",              // Category
    CodeState::Experimental,// Code state
    Tags::None,             // Tags
};

const ProcessorInfo Topology::getProcessorInfo() const
{
    return processorInfo_;
}

Topology::Topology()
    : Processor()
    , inData("inData")
    , outMesh("meshOut")
    , meshBBoxOut("meshBBoxOut")
    // TODO: Initialize additional properties
    // propertyName("propertyIdentifier", "Display Name of the Propery",
    // default value (optional), minimum value (optional), maximum value (optional), increment
    // (optional)); propertyIdentifier cannot have spaces
    , propBoundarySwitchPoint("BoundarySwitchPoint", "Boundary Switch Point")

{
    // Register Ports
    addPort(outMesh);
    addPort(inData);
    addPort(meshBBoxOut);

    // TODO: Register additional properties
    // addProperty(propertyName);
    addProperty(propBoundarySwitchPoint);
}

void Topology::process()
{
    // Get input
    if (!inData.hasData())
    {
        return;
    }
    auto vol = inData.getData();

    // Retrieve data in a form that we can access it
    const VectorField2 vectorField = VectorField2::createFieldFromVolume(vol);

    // Add a bounding box to the mesh
    const dvec2& BBoxMin = vectorField.getBBoxMin();
    const dvec2& BBoxMax = vectorField.getBBoxMax();
    auto bboxMesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> bboxVertices;
    auto indexBufferBBox = bboxMesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
    // Bounding Box vertex 0
    vec4 black = vec4(0, 0, 0, 1);
    Integrator::drawNextPointInPolyline(BBoxMin, black, indexBufferBBox.get(), bboxVertices);
    Integrator::drawNextPointInPolyline(vec2(BBoxMin[0], BBoxMax[1]), black, indexBufferBBox.get(), bboxVertices);
    Integrator::drawNextPointInPolyline(BBoxMax, black, indexBufferBBox.get(), bboxVertices);
    Integrator::drawNextPointInPolyline(vec2(BBoxMax[0], BBoxMin[1]), black, indexBufferBBox.get(), bboxVertices);
    // Connect back to the first point, to make a full rectangle
    indexBufferBBox->add(static_cast<std::uint32_t>(0));
    bboxMesh->addVertices(bboxVertices);
    meshBBoxOut.setData(bboxMesh);

    // Initialize mesh, vertices and index buffers for separatrices
    auto mesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> vertices;
    // Either add all line segments to this index buffer (one large buffer, two consecutive points
    // make up one line), or use several index buffers with connectivity type strip.
    auto indexBufferSeparatrices = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);
    // auto indexBufferSeparatrices = mesh->addIndexBuffer(DrawType::Lines,
    // ConnectivityType::Strip);

    auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);

    // TODO: Compute the topological skeleton of the input vector field.
    // Find the critical points and color them according to their type.
    // Integrate all separatrices.

    double err = 1e-5;
    std::vector<dvec2> criticalPoints = extractCriticalPoints(vectorField, err);
    size2_t dims = vectorField.getNumVerticesPerDim();
    std::vector<std::vector<dvec2>> separatrices;

    //展示critical points，且不同的points进行六分类
    for (int i = 0; i < criticalPoints.size(); i++)
    {
        std::cout << "critical point " << i << std::endl;
        //dvec4 color = vec4(1, 0, 0, 1);
        dvec4 color = classifyCriticalPoints(criticalPoints[i], vectorField, separatrices);
        Integrator::drawPoint(criticalPoints[i], color, indexBufferPoints.get(), vertices);
    }
    //画separatrices
    vec4 white = vec4(1, 1, 1, 1);
    for (int i = 0; i < separatrices.size(); i++)
    {
        //draw each separatrice
        auto indexBufferSeparatrices = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
        for (int i_sep = 0; i_sep < separatrices[i].size(); i_sep++)
        {
            Integrator::drawNextPointInPolyline(separatrices[i][i_sep], white, indexBufferSeparatrices.get(), vertices);
        }
    }

    // 展示包含边缘的所有点（灰色显示）
    if (propBoundarySwitchPoint)
    {
        vec4 grey = vec4(0.5, 0.5, 0.5, 1);
        std::vector<dvec2> switchPointsLeft = extractBoundarySwitchPoints(vectorField, BBoxMin, vec2(BBoxMin[0], BBoxMax[1])); //vertical left
        std::vector<dvec2> switchPointsRight = extractBoundarySwitchPoints(vectorField, vec2(BBoxMax[0], BBoxMin[1]), BBoxMax);//vertical right
        std::vector<dvec2> switchPointsUp = extractBoundarySwitchPoints(vectorField, vec2(BBoxMin[0], BBoxMax[1]), BBoxMax);   //horizontal top
        std::vector<dvec2> switchPointsDown = extractBoundarySwitchPoints(vectorField, BBoxMin, vec2(BBoxMax[0], BBoxMin[1])); //horizontal top

        std::vector<std::vector<dvec2>> switchPointsVectors = {switchPointsLeft, switchPointsRight, switchPointsUp, switchPointsDown};

        for (std::vector<dvec2> switchPointsVector : switchPointsVectors)
        {
            //draw switch point
            for (dvec2 switchPoint : switchPointsVector)
            {
                Integrator::drawPoint(switchPoint, grey, indexBufferPoints.get(), vertices);
                //integrate its seperatrice and draw them
                mat2 eigenvectorsJacobian = util::eigenAnalysis(vectorField.derive(switchPoint)).eigenvectors;
                std::vector<std::vector<dvec2>> separatrices = computeSeparatrices(switchPoint, vectorField, eigenvectorsJacobian);
                for (int i = 0; i < separatrices.size(); i++)
                {
                    //draw each separatrice
                    auto indexBufferSeparatrices = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
                    for (int i_sep = 0; i_sep < separatrices[i].size(); i_sep++)
                    {
                        Integrator::drawNextPointInPolyline(separatrices[i][i_sep], grey, indexBufferSeparatrices.get(), vertices);
                    }
                }
            }
        }
    }

    // Other helpful functions
    // dvec2 pos = vectorField.getPositionAtVertex(size2_t(i, j));
    // Computing the jacobian at a position
    // dmat2 jacobian = vectorField.derive(pos);
    // Doing the eigen analysis
    // auto eigenResult = util::eigenAnalysis(jacobian);
    // The result of the eigen analysis has attributed eigenvaluesRe eigenvaluesIm and
    // eigenvectors


    // Accessing the colors
    vec4 colorCenter = ColorsCP[static_cast<int>(TypeCP::Center)];

    mesh->addVertices(vertices);
    outMesh.setData(mesh);
}


void Topology::drawLineSegment(const dvec2& v1, const dvec2& v2, const vec4& color, IndexBufferRAM* indexBuffer, std::vector<BasicMesh::Vertex>& vertices)
{
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    vertices.push_back({vec3(v1[0], v1[1], 0), vec3(0, 0, 1), vec3(v1[0], v1[1], 0), color});
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    vertices.push_back({vec3(v2[0], v2[1], 0), vec3(0, 0, 1), vec3(v2[0], v2[1], 0), color});
}


dvec2 Topology::cellPosition(double i, double j, const VectorField2 vectorField)
{
    const dvec2 bboxMin = vectorField.getBBoxMin();
    const dvec2 bboxMax = vectorField.getBBoxMax();
    //将边界框的x或y方向长度按照dims[0]-1等分，每一小段代表一个单元格的x或y方向跨度。最终返回二维空间中某个单元格的实际坐标。
    const size2_t dims = vectorField.getNumVerticesPerDim();
    const dvec2 bboxSize = bboxMax - bboxMin;
    // 计算并返回单元格的实际坐标
    return dvec2(bboxMin[0] + bboxSize[0] * i / (dims[0] - 1), bboxMin[1] + bboxSize[1] * j / (dims[1] - 1));
}


std::vector<dvec2> Topology::extractCriticalPoints(const VectorField2 vectorField, const double eps)
{
    std::vector<dvec2> criticalPoints;
    size2_t dims = vectorField.getNumVerticesPerDim();

    for (size_t j = 0; j < dims[1]; ++j)
    {
        for (size_t i = 0; i < dims[0]; ++i)
        {
            //initial_corners是左下角和右上角的点
            std::vector<dvec2> initial_corners = {cellPosition(i, j, vectorField), cellPosition(i + 1, j + 1, vectorField)};
            std::queue<std::vector<dvec2>> cornersQueue;
            cornersQueue.push(initial_corners);
            while (!cornersQueue.empty())
            {

                std::vector<dvec2> bottomLeftAndTopRight = cornersQueue.front();
                dvec2 bottomLeft = bottomLeftAndTopRight[0];
                dvec2 topRight = bottomLeftAndTopRight[1];
                cornersQueue.pop();

                std::vector<dvec2> corners;
                //把四个点都push到corners里
                corners.push_back(bottomLeft);                      //bottom left
                corners.push_back(vec2(topRight[0], bottomLeft[1]));//bottom right
                corners.push_back(vec2(bottomLeft[0], topRight[1]));//top left
                corners.push_back(topRight);                        //top right

                bool zeroFound = false;
                ivec2 signChanges = vec2(0, 0);
                findZerosAndSignChanges(criticalPoints, zeroFound, signChanges, corners, vectorField, eps);

                //找到0的话就break出循环
                if (zeroFound)
                {
                    break;// If zero is found, exit the loop
                }
                //domain decomposition
                //必须x和y都有正负时候，才把四个子区域的坐标也压入队列cornersQueue
                if (signChanges[0] == 1 && signChanges[1] == 1)
                {
                    dvec2 center = 0.5 * (bottomLeft + topRight);
                    cornersQueue.push({bottomLeft, center});                                          // bottom left cell
                    cornersQueue.push({vec2(center[0], bottomLeft[1]), vec2(topRight[0], center[1])});// bottom right cell
                    cornersQueue.push({vec2(bottomLeft[0], center[1]), vec2(center[0], topRight[1])});// top left cell
                    cornersQueue.push({center, topRight});                                            // top right cell
                }
            }
        }
    }
    return criticalPoints;
}

//change-of-sign test
void Topology::findZerosAndSignChanges(std::vector<dvec2>& criticalPoints, bool& zeroFound, ivec2& signChanges, const std::vector<dvec2> corners, const VectorField2 vectorField, const double eps)
{
    std::array<int, 2> lastSign = {0, 0};// 初始化为0

    for (int i_corner = 0; i_corner < 4; i_corner++)
    {
        // 线性插值
        vec2 valueAtCorner = vectorField.interpolate(corners[i_corner]);

        // eps是用于判断停止的误差，来决定现在选的值是否足够靠近0
        if (glm::length(valueAtCorner) <= eps)
        {// 找到零点
            criticalPoints.push_back(corners[i_corner]);
            zeroFound = true;
            break;
        }
        else
        {
            // 计算当前符号
            std::array<int, 2> currentSign;
            for (int i_component = 0; i_component < 2; ++i_component)
            {
                if (valueAtCorner[i_component] > 0)
                {
                    currentSign[i_component] = 1;
                }
                else if (valueAtCorner[i_component] < 0)
                {
                    currentSign[i_component] = -1;
                }
                else
                {
                    currentSign[i_component] = 0;
                }
            }

            // 检查符号变化
            for (int i_component = 0; i_component < 2; ++i_component)
            {
                if (currentSign[i_component] * lastSign[i_component] < 0)
                {
                    signChanges[i_component] = true;// 记录符号变化
                }
                // 仅在当前符号非零时更新lastSign
                if (currentSign[i_component] != 0)
                {
                    lastSign[i_component] = currentSign[i_component];
                }
            }
        }
    }
}

vec4 Topology::classifyCriticalPoints(const vec2 pos, const VectorField2 vectorField, std::vector<std::vector<dvec2>>& separatrices)
{
    dmat2 jacobian = vectorField.derive(pos);
    auto eigenResult = util::eigenAnalysis(jacobian);

    double R1 = eigenResult.eigenvaluesRe[0];
    double R2 = eigenResult.eigenvaluesRe[1];
    double I1 = eigenResult.eigenvaluesIm[0];
    double I2 = eigenResult.eigenvaluesIm[1];

    double eps = 1e-7;

    // 判断特征值的实部和虚部来确定关键点的类型
    //vec4(1, 0.95, 0.5, 1),  // Saddle - Yellow
    //vec4(0.5, 0.5, 0.9, 1),  // AttractingNode - Sink - Blue
    //vec4(0.9, 0.5, 0.5, 1),  // RepellingNode - Source - Red
    //vec4(0.5, 0, 0.9, 1),// AttractingFocus - Purple
    //vec4(0.9, 0.5, 0.0, 1),// RepellingFocus - Orange
    //vec4(0.3, 0.6, 0.3, 1)   // Center - Green
    if ((abs(I2) <= eps) && (abs(I1) <= eps) && (R1 * R2 < 0))
    {
        // Saddle Point
        std::vector<std::vector<dvec2>> separatricesSaddle = computeSeparatrices(pos, vectorField, eigenResult.eigenvectors);
        separatrices.insert(separatrices.end(), separatricesSaddle.begin(), separatricesSaddle.end());
        return ColorsCP[0];// 'Saddlepoint'
    }
    if ((abs(I2) <= eps) && (abs(I1) <= eps) && (R1 < 0) && (R2 < 0))
    {
        return ColorsCP[1];// 'Attracting Node'
    }
    if ((abs(I2) <= eps) && (abs(I1) <= eps) && (R1 > 0) && (R2 > 0))
    {
        return ColorsCP[2];// 'Repelling Node'
    }
    if ((abs(R1 - R2) <= eps) && (R1 < 0) && (R2 < 0) && (abs(I1 + I2) <= eps) && (I1 != 0) && (I2 != 0))
    {
        return ColorsCP[3];// 'Attracting Focus'
    }
    if ((abs(R1 - R2) <= eps) && (R1 > 0) && (R2 > 0) && (abs(I1 + I2) <= eps) && (I1 != 0) && (I2 != 0))
    {
        return ColorsCP[4];// 'Repelling Focus'
    }
    if ((abs(R2) <= eps) && (abs(R1) <= eps) && (abs(I1 + I2) <= eps) && (I2 != 0))
    {
        return ColorsCP[5];// 'Center'
    }

    return vec4(0, 0, 0, 0);// 返回透明颜色，表示异常情况
}


std::vector<std::vector<dvec2>> Topology::computeSeparatrices(dvec2 saddle, const VectorField2& vectorField, mat2 eigenvectors)
{
    std::vector<std::vector<dvec2>> separatrices;

    double stepSizeInitial = 0.01;
    double stepSize = 0.01;
    double minVelocity = 0.01;

    // Define directions based on eigenvectors
    std::vector<dvec2> directions = {eigenvectors[0], eigenvectors[1], -eigenvectors[0], -eigenvectors[1]};

    for (const auto& direction : directions)
    {
        // Take a step in the direction
        dvec2 position = saddle + stepSizeInitial * glm::normalize(direction);// Euler step

        // Integrate the streamline from the new position
        std::vector<dvec2> streamline = Integrator::StreamLines(vectorField, position, stepSize,
                                                                false,      // backwardDirection, set to true or false based on your needs
                                                                true,       // integrateDirectionField, set to true or false based on your needs
                                                                minVelocity,// speedThreshold
                                                                300         // kernelSize
        );
        std::vector<dvec2> streamline1 = Integrator::StreamLines(vectorField, position, stepSize,
                                                                 true,       // backwardDirection, set to true or false based on your needs
                                                                 true,       // integrateDirectionField, set to true or false based on your needs
                                                                 minVelocity,// speedThreshold
                                                                 300         // kernelSize
        );

        // Store the computed streamline in separatrices
        separatrices.push_back(streamline);
        separatrices.push_back(streamline1);
    }


    return separatrices;
}

std::vector<dvec2> Topology::extractBoundarySwitchPoints(const VectorField2 vectorField, const vec2 startPoint, const vec2 endPoint)
{
    std::vector<dvec2> boundarySwitchPoints;
    size2_t dims = vectorField.getNumVerticesPerDim();
    //we will look for a zero in the dimension (0 or 1) of the vector field that is perpendicular to the edge
    int dim = (int)(startPoint[0] != endPoint[0]);//dim = 1 for a horizontal edge and dim = 0 for a vertical egde

    // Looping through all values on the boundaries in the vector field
    for (size_t i = 0; i < dims[dim] - 1; ++i)
    {
        //for each line segment i, we do domain decomposition to check for a boundary switch point
        //we stop when the vector field component is close enough to zero
        vec2 left = startPoint + i / (double)(dims[dim] - 1) * (endPoint - startPoint);
        vec2 right = startPoint + (i + 1) / (double)(dims[dim] - 1) * (endPoint - startPoint);
        double valueLeft = vectorField.interpolate(left)[dim];
        double valueRight = vectorField.interpolate(right)[dim];

        if ((valueLeft - valueRight != 0) && (valueLeft * valueRight <= 0))
        {//we know this is linear interpolation
            vec2 root = left - valueLeft / (valueRight - valueLeft) * (right - left);
            boundarySwitchPoints.push_back(root);
        }
    }
    return boundarySwitchPoints;
}

}// namespace inviwo
