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
    // 1. 输入数据的处理
        // 首先检查输入数据 inData 是否有效，如果没有数据就退出函数。
    if (!inData.hasData())
    {
        return;
    }
    auto vol = inData.getData();
        // Retrieve data in a form that we can access it
        // 然后将体数据 vol 转换为一个二维矢量场 vectorField，这个矢量场定义了每个点的向量信息，供后续计算使用。
    const VectorField2 vectorField = VectorField2::createFieldFromVolume(vol);
    
    // 2. 绘制矢量场的包围盒 // Add a bounding box to the mesh    
        // BBoxMin 和 BBoxMax 分别是矢量场的包围盒的最小和最大点坐标，用来定义场的边界。        
    const dvec2& BBoxMin = vectorField.getBBoxMin();
    const dvec2& BBoxMax = vectorField.getBBoxMax();
        // 创建了一个 BasicMesh 对象 bboxMesh，用于存储包围盒的网格信息。接下来将包围盒的四个顶点绘制成一个矩形，并添加到 indexBufferBBox 中。
    auto bboxMesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> bboxVertices;
    auto indexBufferBBox = bboxMesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
        // 使用 Integrator::drawNextPointInPolyline 函数将每个顶点添加到包围盒的线段中，最后将第一个顶点与最后一个顶点连接，构成完整的矩形。
    vec4 black = vec4(0, 0, 0, 1);// Bounding Box vertex 0
    Integrator::drawNextPointInPolyline(BBoxMin, black, indexBufferBBox.get(), bboxVertices);
    Integrator::drawNextPointInPolyline(vec2(BBoxMin[0], BBoxMax[1]), black, indexBufferBBox.get(), bboxVertices);
    Integrator::drawNextPointInPolyline(BBoxMax, black, indexBufferBBox.get(), bboxVertices);
    Integrator::drawNextPointInPolyline(vec2(BBoxMax[0], BBoxMin[1]), black, indexBufferBBox.get(), bboxVertices);
        // 回连到第一个点，构成完整的矩形
        // Connect back to the first point, to make a full rectangle
    indexBufferBBox->add(static_cast<std::uint32_t>(0));
    bboxMesh->addVertices(bboxVertices);
    meshBBoxOut.setData(bboxMesh);
    
    
    // 3. 初始化分离线和关键点的网格  Initialize mesh, vertices and index buffers for separatrices
        // 为分离线和关键点的绘制初始化网格 mesh，
        // 并为分离线和点分别创建了索引缓冲区 indexBufferSeparatrices 和 indexBufferPoints
    auto mesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> vertices;
    auto indexBufferSeparatrices = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);
    auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);

    // TODO: Compute the topological skeleton of the input vector field.
    // Find the critical points and color them according to their type.
    // Integrate all separatrices.

    // 4. 提取关键点和分离线
        // 使用 extractCriticalPoints 函数从矢量场中提取关键点，误差容忍度设置为 1e-5。关键点是矢量场中的重要位置（例如鞍点、中心点等）。
        // 使用二维矢量场的尺寸信息 dims，初始化 separatrices 来存储分离线的信息。
    double err = 1e-5;
    std::vector<dvec2> criticalPoints = extractCriticalPoints(vectorField, err);
    size2_t dims = vectorField.getNumVerticesPerDim();
    std::vector<std::vector<dvec2>> separatrices;

    // 5. 展示critical points，且不同的points进行六分类
        // 对每个关键点的分离线进行绘制，使用白色 vec4(1, 1, 1, 1) 作为分离线的颜色。分离线将矢量场分割为不同的区域。
        // 使用 Integrator::drawNextPointInPolyline 逐个绘制分离线的顶点，并将其连接。
    for (int i = 0; i < criticalPoints.size(); i++)
    {
        // std::cout << "critical point " << i << std::endl;
        //dvec4 color = vec4(1, 0, 0, 1);
        dvec4 color = classifyCriticalPoints(criticalPoints[i], vectorField, separatrices);
        Integrator::drawPoint(criticalPoints[i], color, indexBufferPoints.get(), vertices);
    }
    // 6. 绘制分隔线
    vec4 white = vec4(1, 1, 1, 1);
    for (int i = 0; i < separatrices.size(); i++)
    {
        //draw each separatrice
        // 绘制每条分离线
        auto indexBufferSeparatrices = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
        for (int i_sep = 0; i_sep < separatrices[i].size(); i_sep++)
        {
            Integrator::drawNextPointInPolyline(separatrices[i][i_sep], white, indexBufferSeparatrices.get(), vertices);
        }
    }

    // 7. 处理边界开关点  展示包含边缘的所有点（灰色显示）
        // extractBoundarySwitchPoints 函数用于提取包围盒边界上的开关点
        // 接下来这些点也会像关键点一样绘制，并使用雅可比矩阵的特征向量计算分离线。
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
            // 绘制开关点
            for (dvec2 switchPoint : switchPointsVector)
            {
                Integrator::drawPoint(switchPoint, grey, indexBufferPoints.get(), vertices);
                // 计算并绘制开关点的分离线
                mat2 eigenvectorsJacobian = util::eigenAnalysis(vectorField.derive(switchPoint)).eigenvectors;
                std::vector<std::vector<dvec2>> separatrices = computeSeparatrices(switchPoint, vectorField, eigenvectorsJacobian);
                for (int i = 0; i < separatrices.size(); i++)
                {
                    // 添加顶点并输出网格数据
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

// 画出一个线段
void Topology::drawLineSegment(const dvec2& v1, const dvec2& v2, const vec4& color, IndexBufferRAM* indexBuffer, std::vector<BasicMesh::Vertex>& vertices)
{
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    vertices.push_back({vec3(v1[0], v1[1], 0), vec3(0, 0, 1), vec3(v1[0], v1[1], 0), color});
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    vertices.push_back({vec3(v2[0], v2[1], 0), vec3(0, 0, 1), vec3(v2[0], v2[1], 0), color});
}

// 计算单元格的实际坐标
// 该函数计算矢量场中网格单元的实际坐标，基于给定的 i 和 j 值计算出单元格在二维空间中的位置。矢量场的范围通过 bboxMin 和 bboxMax 确定，函数会按网格的维度 dims 将范围划分为若干等分，从而得出指定点的真实坐标。
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


// 提取关键点的函数，根据误差阈值 `eps` 找到矢量场中的关键点
/*
该函数遍历整个矢量场并寻找关键点，关键点是矢量场中速度为零的位置（即矢量为零的点）。
步骤：

遍历矢量场中的每一个单元格。
使用 cellPosition 获取单元格的边界角点坐标。
使用四分法细分单元格，如果在单元格的角点之间发现矢量符号发生变化，则说明可能存在零点，进一步细化查找。
如果找到速度为零的点（即零点），则记录该点为关键点。
*/
std::vector<dvec2> Topology::extractCriticalPoints(const VectorField2 vectorField, const double eps)
{
    std::vector<dvec2> criticalPoints;
    size2_t dims = vectorField.getNumVerticesPerDim();
    // 遍历矢量场的每个单元格
    for (size_t j = 0; j < dims[1]; ++j)
    {
        for (size_t i = 0; i < dims[0]; ++i)
        {
            // initial_corners 是单元格的左下角和右上角的坐标
            std::vector<dvec2> initial_corners = {cellPosition(i, j, vectorField), cellPosition(i + 1, j + 1, vectorField)};
            std::queue<std::vector<dvec2>> cornersQueue;
            cornersQueue.push(initial_corners);
            while (!cornersQueue.empty())
            {
                // 将四个顶点依次加入 corners
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
// 符号变化和零点检测的函数 
/*
这个函数用来检测给定单元格的四个角点是否存在矢量符号变化，从而判断该单元格是否包含关键点。

通过线性插值获取角点处的矢量值。
通过矢量的模与误差阈值 eps 比较，判断该点是否接近零点。
如果当前角点的矢量值在 x 或 y 方向上出现符号变化，则记录下符号变化，并可能进一步细分区域。
*/
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
// 分类关键点的函数 classifyCriticalPoints
/*
该函数通过分析关键点处的雅可比矩阵的特征值和特征向量，来判断关键点的类型（如鞍点、吸引子、排斥子等）。具体根据特征值的实部和虚部来判断关键点属于哪种类型，并返回相应的颜色用于后续的可视化。

Saddle 鞍点：特征值的实部异号，且虚部接近零。
Attracting Node / Sink	 吸引子节点：特征值的实部均为负数，虚部接近零。
Repelling Node / Source	排斥子节点：特征值的实部均为正数，虚部接近零。
Attracting Focus	吸引焦点：特征值的实部相等且为负数，虚部不为零。
Repelling Focus	排斥焦点：特征值的实部相等且为正数，虚部不为零。
Center 中心点：特征值的实部为零，虚部不为零。
每种关键点类型通过不同的颜色进行标识（例如鞍点为黄色，吸引节点为蓝色，排斥节点为红色等），最后返回分类结果供可视化使用。
*/
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


// 分割线是从关键点发出的曲线或线段，它们将矢量场划分为不同的区域。在这些区域中，流场的方向和行为模式会有所不同。
// 分割线通常从鞍点（saddle points）开始延伸，沿着流场的方向一直延伸到场的边界或者另一个关键点。
std::vector<std::vector<dvec2>> Topology::computeSeparatrices(dvec2 saddle, const VectorField2& vectorField, mat2 eigenvectors)
{
    /*
    初始化：

    创建一个 separatrices 用于存储最终计算的分割线。
    设置初始步长 stepSizeInitial 和 stepSize，用于欧拉积分法（Euler Step）时的步长大小。
    设置最小速度阈值 minVelocity，当速度低于这个值时，积分将停止。
    */
    std::vector<std::vector<dvec2>> separatrices;

    double stepSizeInitial = 0.01;
    double stepSize = 0.01;
    double minVelocity = 0.01;

    // Define directions based on eigenvectors
    // 使用特征向量矩阵 eigenvectors 定义分割线的方向。特征向量决定了从鞍点出发的初始方向。将特征向量以及它们的负向量作为可能的分割线方向。
    std::vector<dvec2> directions = {eigenvectors[0], eigenvectors[1], -eigenvectors[0], -eigenvectors[1]};

    /*
    计算分割线：

    对于每个方向 direction，首先使用欧拉法从鞍点出发沿着方向向前移动一个小步长 stepSizeInitial。
    使用 Integrator::StreamLines 方法进行流线积分（Streamline Integration），沿着流场的方向生成一条流线：
    false 表示不向后积分，即顺着流场方向前进。
    true 表示使用方向场进行积分。
    minVelocity 用于控制速度阈值，当速度小于这个值时，积分停止。
    300 为内核大小，用于控制流线的积分长度。
    */
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
        // 将每次计算的流线 streamline 和反向流线 streamline1 存储到 separatrices 中。
        separatrices.push_back(streamline);
        separatrices.push_back(streamline1);
    }


    return separatrices;
}

// 该函数用于提取边界上的交换点（Boundary Switch Points）
// 即流场中某个方向的矢量分量从正数变为负数（或反之）的点。
// 通过在流场的边界上寻找零交叉点来确定交换点，这些点在分析流场时具有重要意义。
std::vector<dvec2> Topology::extractBoundarySwitchPoints(const VectorField2 vectorField, const vec2 startPoint, const vec2 endPoint)
{
    // 1. 获取维度：获取向量场的维度信息，表示在每个维度上有多少个顶点。
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