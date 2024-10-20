/*********************************************************************
 *  Author  : Himangshu Saikia, Wiebke Koepp, Anke Friederici
 *  Init    : Monday, September 11, 2017 - 12:58:42
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <labmarchingsquares/marchingsquares.h>
#include <inviwo/core/util/utilities.h>

namespace inviwo
{

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo MarchingSquares::processorInfo_{
    "org.inviwo.MarchingSquares",// Class identifier
    "Marching Squares",          // Display name
    "KTH Lab",                   // Category
    CodeState::Experimental,     // Code state
    Tags::None,                  // Tags
};

const ProcessorInfo MarchingSquares::getProcessorInfo() const
{
    return processorInfo_;
}

MarchingSquares::MarchingSquares()
    : Processor()
    , inData("volumeIn")
    , meshIsoOut("meshIsoOut")
    , meshGridOut("meshGridOut")
    , propShowGrid("showGrid", "Show Grid")
    , propGridColor("gridColor", "Grid Lines Color", vec4(0.0f, 0.0f, 0.0f, 1.0f), vec4(0.0f), vec4(1.0f), vec4(0.1f), InvalidationLevel::InvalidOutput, PropertySemantics::Color)
    , propDeciderType("deciderType", "Decider Type")
    , propRandomSeed("seed", "Random Seed", 0, 0, std::mt19937::max())
    , propMultiple("multiple", "Iso Levels")
    , propIsoValue("isovalue", "Iso Value")
    , propIsoColor("isoColor", "Color", vec4(0.0f, 0.0f, 1.0f, 1.0f), vec4(0.0f), vec4(1.0f), vec4(0.1f), InvalidationLevel::InvalidOutput, PropertySemantics::Color)
    , propNumContours("numContours", "Number of Contours", 1, 1, 50, 1)
    , propIsoTransferFunc("isoTransferFunc", "Colors", &inData)


{
    // Register ports
    addPort(inData);
    addPort(meshIsoOut);
    addPort(meshGridOut);

    // Register properties
    addProperty(propShowGrid);
    addProperty(propGridColor);

    addProperty(propApplyGaussianFilter);

    addProperty(propDeciderType);
    propDeciderType.addOption("asymptotic", "Asymptotic", 0);
    propDeciderType.addOption("random", "Random", 1);

    addProperty(propRandomSeed);
    propRandomSeed.setSemantics(PropertySemantics::Text);

    addProperty(propMultiple);

    propMultiple.addOption("single", "Single", 0);
    addProperty(propIsoValue);
    addProperty(propIsoColor);

    propMultiple.addOption("multiple", "Multiple", 1);
    addProperty(propNumContours);
    addProperty(propIsoTransferFunc);

    //**Task 3.2(b)
    // The default transfer function has just two blue points
    propIsoTransferFunc.get().clear();
    propIsoTransferFunc.get().add(0.0f, vec4(0.0f, 0.0f, 1.0f, 1.0f)); //blue
    //propIsoTransferFunc.get().add(0.3f, vec4(0.0f, 1.0f, 0.0f, 1.0f));// green
    //propIsoTransferFunc.get().add(0.6f, vec4(1.0f, 1.0f, 0.0f, 1.0f));// yellow
    propIsoTransferFunc.get().add(1.0f, vec4(1.0f, 0.0f, 0.0f, 1.0f)); //red 
    propIsoTransferFunc.setCurrentStateAsDefault();


    util::hide(propGridColor, propRandomSeed, propNumContours, propIsoTransferFunc);

    propDeciderType.onChange([this]() {
        if (propDeciderType.get() == 1)
        {
            util::show(propRandomSeed);
        }
        else
        {
            util::hide(propRandomSeed);
        }
    });

    // Show the grid color property only if grid is actually displayed
    propShowGrid.onChange([this]() {
        if (propShowGrid.get())
        {
            util::show(propGridColor);
        }
        else
        {
            util::hide(propGridColor);
        }
    });

    // Show options based on display of one or multiple iso contours
    propMultiple.onChange([this]() {
        if (propMultiple.get() == 0)
        {
            util::show(propIsoValue, propIsoColor);
            util::hide(propNumContours, propIsoTransferFunc);
        }
        else
        {
            //util::hide(propIsoValue);
            //util::show(propIsoColor, propNumContours);

            // TODO (Bonus): Comment out above if you are using the transfer function
            // and comment in below instead
             util::hide(propIsoValue, propIsoColor);
             util::show(propNumContours, propIsoTransferFunc);
        }
    });
}

void applyGaussianFilter(const ScalarField2& inputField, ScalarField2& smoothedField, int nVertPerDim, double sigma)
{
    // Define Gaussian kernel size and standard deviation:
    const int kernelRadius = 3;// radius=3, kernal size = 7x7 
    const double pi = 3.14159265358979323846;

    // Create and calculate the Gaussian kernel:
    double kernel[2 * kernelRadius + 1][2 * kernelRadius + 1];
    double sum = 0.0;

    for (int i = -kernelRadius; i <= kernelRadius; i++)
    {
        for (int j = -kernelRadius; j <= kernelRadius; j++)
        {
            kernel[i + kernelRadius][j + kernelRadius] = exp(-(i * i + j * j) / (2 * sigma * sigma)) / (2 * pi * sigma * sigma);
            sum += kernel[i + kernelRadius][j + kernelRadius];
        }
    }

    // Normalize the Gaussian kernel:
    for (int i = 0; i <= 2 * kernelRadius; i++)
    {
        for (int j = 0; j <= 2 * kernelRadius; j++)
        {
            kernel[i][j] /= sum;
        }
    }

    // Apply the Gaussian filter:
    for (int i = 0; i < nVertPerDim; i++)
    {
        for (int j = 0; j < nVertPerDim; j++)
        {
            double smoothedValue = 0.0;
            // Apply Gaussian Kernal
            for (int ki = -kernelRadius; ki <= kernelRadius; ki++)
            {
                for (int kj = -kernelRadius; kj <= kernelRadius; kj++)
                {
                    // ensure that the index stays within the grid bounds
                    int ni = std::min(std::max(i + ki, 0), nVertPerDim - 1);
                    int nj = std::min(std::max(j + kj, 0), nVertPerDim - 1);
                    smoothedValue += inputField.getValueAtVertex({ni, nj}) * kernel[ki + kernelRadius][kj + kernelRadius];
                }
            }
            smoothedField.setValueAtVertex({i, j}, smoothedValue);
        }
    }
}

void MarchingSquares::process()
{
    if (!inData.hasData())
    {
        return;
    }

    // Create a structured grid from the input volume
    auto vol = inData.getData();
    auto grid = ScalarField2::createFieldFromVolume(vol);
    const ivec2 nVertPerDim = grid.getNumVerticesPerDim();

    // Check if Gaussian filtering is enabled
    if (propApplyGaussianFilter.get()==1)
    {
        // Create a new smoothedField to store the smoothed data
        ScalarField2 smoothedField = ScalarField2(nVertPerDim, grid.getBBoxMin(), grid.getBBoxMax() - grid.getBBoxMin());

        // Apply Gaussian filter
        applyGaussianFilter(grid, smoothedField, nVertPerDim.x, 1.0);// Assum sigma = 1.0

         // Continue processing with the smoothedField
        grid = smoothedField;// Use the smoothed data for further contour drawing
    }


    // Extract the minimum and maximum value from the input data
    const double minValue = grid.getMinValue();
    const double maxValue = grid.getMaxValue();

    // Set the range for the isovalue to that minimum and maximum
    propIsoValue.setMinValue(minValue);
    propIsoValue.setMaxValue(maxValue);

    // You can print to the Inviwo console with Log-commands:
    LogProcessorInfo("This scalar field contains values between " << minValue << " and " << maxValue << ".");
    // You can also inform about errors and warnings:
    // LogProcessorWarn("I am warning about something"); // Will print warning message in yellow
    // LogProcessorError("I am letting you know about an error"); // Will print error message in red
    // (There is also LogNetwork...() and just Log...(), these display a different source,
    // LogProcessor...() for example displays the name of the processor in the workspace while
    // Log...() displays the identifier of the processor (thus with multiple processors of the
    // same kind you would not know which one the information is coming from

    // Get the definition of our structured grid with
    // - number of vertices in each dimension {nx, ny}
    //const ivec2 nVertPerDim = grid.getNumVerticesPerDim();
    // - bounding box {xmin, ymin} - {xmax, ymax}
    const dvec2 bBoxMin = grid.getBBoxMin();
    const dvec2 bBoxMax = grid.getBBoxMax();
    float boundingBoxWidth = bBoxMax.x - bBoxMin.x;
    float boundingBoxHeight = bBoxMax.y - bBoxMin.y;
    // - cell size {dx, dy}
    const dvec2 cellSize = grid.getCellSize();


    // Values at the vertex positions can be accessed by the indices of the vertex
    // with index i ranging between [0, nx-1] and j in [0, ny-1]
    ivec2 ij = {0, 0};
    double valueAt00 = grid.getValueAtVertex(ij);
    LogProcessorInfo("The value at (0,0) is: " << valueAt00 << ".");

    // Initialize the output: mesh and vertices for the grid and bounding box
    auto gridmesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> gridvertices;

    // Set the random seed to the one selected in the interface
    randGenerator.seed(static_cast<std::mt19937::result_type>(propRandomSeed.get()));
    // You can create a random sample between min and max with
    float minRand = 0.0;
    float maxRand = 1.0;
    float rand = randomValue(minRand, maxRand);
    LogProcessorInfo("The first random sample for seed " << propRandomSeed.get() << " between " << minRand << " and " << maxRand << " is " << rand << ".");


        // Properties are accessed with propertyName.get()
    if (propShowGrid.get())
    {
        // TODO: Add grid lines of the given color

        // The function drawLineSegments creates two vertices at the specified positions,
        // that are placed into the Vertex vector defining our mesh.
        // An index buffer specifies which of those vertices should be grouped into to make up lines/trianges/quads.
        // Here two vertices make up a line segment.
        auto indexBufferGrid = gridmesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);

        // Draw a line segment from v1 to v2 with a color, the coordinates in the final
        // image range from 0 to 1 for both x and y
        //vec2 v1 = vec2(0.5, 0.5);
        //vec2 v2 = vec2(0.7, 0.7);
        //drawLineSegment(v1, v2, propGridColor.get(), indexBufferGrid, vertices);
        
        // task 3.1 
        // vertical line  
        for (size_t i = 0; i < nVertPerDim.x; i++)
        {
            float ix = 1.0 / (nVertPerDim.x - 1) * i;
            drawLineSegment(vec2(ix, 0), vec2(ix, 1), propGridColor.get(), indexBufferGrid.get(), gridvertices);
        }

        // Horizontal line 
        for (size_t j = 0; j < nVertPerDim.y; j++)
        {
            float iy = 1.0 / (nVertPerDim.y - 1) * j;
            drawLineSegment(vec2(0, iy), vec2(1, iy), propGridColor.get(), indexBufferGrid.get(), gridvertices);
        }
    }

    // Set the created grid mesh as output
    gridmesh->addVertices(gridvertices);
    meshGridOut.setData(gridmesh);

    // TODO (Bonus) Gaussian filter
    // Our input is const (i.e. cannot be altered), but you need to compute smoothed data and write
    // it somewhere
    // Create an editable structured grid with ScalarField2 smoothedField =
    // ScalarField2(nVertPerDim, bBoxMin, bBoxMax - bBoxMin); Values can be set with
    // smoothedField.setValueAtVertex({0, 0}, 4.2);
    // and read again in the same way as before
    // smoothedField.getValueAtVertex(ij);
    // Initialize the output: mesh and vertices


    auto mesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> vertices;

    auto myVolume = inData.getData();                           // Get volume data from the input port
    auto dims = myVolume->getDimensions();                      // Get the dimensions of the volume data
    auto myVolumeRAM = myVolume->getRepresentation<VolumeRAM>();// Get VolumeRAM representation


    if (propMultiple.get() == 0)
    {
        // TODO: Draw a single isoline at the specified isovalue (propIsoValue)
        // and color it with the specified color (propIsoColor)
        auto isoBufferGrid = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);
        drawIsolineSingleValue(propIsoValue, propIsoColor.get(), myVolumeRAM, dims, isoBufferGrid.get(), vertices);

    }

    else
    {
        // TODO: Draw the given number (propNumContours) of isolines between
        // the minimum and maximum value
        
        // Calculate the contour line interval
        // The interval between contour lines is calculated by subtracting the minimum value from the maximum value,
        // then dividing by the number of contours plus 1.
        float w = (propIsoValue.getMaxValue() - propIsoValue.getMinValue()) / (propNumContours.get() + 1);
        
        // Loop through and draw each contour line
        for (size_t iso = 0; iso < propNumContours.get(); iso++)
        {
            // Calculate the value of the current contour line and its normalized value
            const double isoVal = propIsoValue.getMinValue() + (iso + 1) * w;
            // Normalize the current contour value to the range [0, 1] for color mapping
            double isoValNormalized = (isoVal - propIsoValue.getMinValue()) / (propIsoValue.getMaxValue() - propIsoValue.getMinValue());
            // Get the corresponding color for the contour line
            const vec4& color = propIsoTransferFunc.get().sample(isoValNormalized);
            // Draw a single contour line
            drawIsolineSingleValue(isoVal, color, myVolumeRAM, dims, mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None).get(), vertices);
        }

        // TODO (Bonus): Use the transfer function property to assign a color
        // The transfer function normalizes the input data and sampling colors
        // from the transfer function assumes normalized input, that means
        // vec4 color = propIsoTransferFunc.get().sample(0.0f);
        // is the color for the minimum value in the data
        // vec4 color = propIsoTransferFunc.get().sample(1.0f);
        // is the color for the maximum value in the data
    }

    // Note: It is possible to add multiple index buffers to the same mesh,
    // thus you could for example add one for the grid lines and one for
    // each isoline
    // Also, consider to write helper functions to avoid code duplication
    // e.g. for the computation of a single iso contour

    mesh->addVertices(vertices);
    meshIsoOut.setData(mesh);
}

float MarchingSquares::randomValue(const float min, const float max) const
{
    return min + uniformReal(randGenerator) * (max - min);
}

void MarchingSquares::drawIsolineSingleValue(const double c, const vec4& color, const VolumeRAM* vr, const size3_t dims, IndexBufferRAM* isoBufferGrid, std::vector<BasicMesh::Vertex>& vertices)
{
    // Loop over the x dimension of the grid, skipping the last column to avoid out-of-bounds access
    for (float ix = 0; ix < (dims.x - 1); ix++)
    {
        // Loop over the y dimension of the grid, skipping the last row to avoid out-of-bounds access
        for (float iy = 0; iy < (dims.y - 1); iy++)
        {
            // Get the scalar values at the four corners of the current cell
            float f00 = getInputValue(vr, dims, ix, iy);         // bottom-left
            float f01 = getInputValue(vr, dims, ix, iy + 1);     // top-left
            float f11 = getInputValue(vr, dims, ix + 1, iy + 1); // top-right
            float f10 = getInputValue(vr, dims, ix + 1, iy);     // bottom-right

            // Find the minimum and maximum scalar values among the four corners of the cell
            float fmin = std::min({f00, f01, f10, f11});
            float fmax = std::max({f00, f01, f10, f11});

            // Check if the contour value 'c' lies between the min and max values in this cell
            // This means there is an isoline passing through the cell
            if (fmin < c && fmax > c)    // If the contour lies within this cell
            {
                // If the condition is true, draw the isoline for this cell
                drawSingleIsoline(ix, iy, c, color, vr, dims, isoBufferGrid, vertices);
            }
        }
    }
}

void MarchingSquares::drawSingleIsoline(float ix, float iy, const double c, const vec4& color, const VolumeRAM* vr, const size3_t dims, IndexBufferRAM* isoBufferGrid,
                                        std::vector<BasicMesh::Vertex>& vertices)
{
    // Get the scalar values at the four corners of the current cell
    float f00 = getInputValue(vr, dims, ix, iy);
    float f01 = getInputValue(vr, dims, ix, iy + 1);
    float f11 = getInputValue(vr, dims, ix + 1, iy + 1);
    float f10 = getInputValue(vr, dims, ix + 1, iy);

    // Define the coordinates of the four corners
    vec2 coord[] = {vec2(ix, iy), vec2(ix, iy + 1), vec2(ix + 1, iy + 1), vec2(ix + 1, iy)};

    // Define directions for each edge of the cell
    // These indicate the direction in which the contour intersects
    vec2 dir[] = {vec2(0, 1), vec2(1, 0), vec2(0, -1), vec2(-1, 0)};

    float f[] = {f00, f01, f11, f10};

    std::vector<vec2> p; // A vector to store the points where the contour intersects the cell edges
    float f1, f2;

    // Loop over each edge of the cell
    for (size_t it = 0; it < 4; it++)
    {
        // Get scalar values for the two vertices of the current edge
        if (it < 3)// For edges 0, 1, and 2, use adjacent corners
        {
            f1 = f[it];
            f2 = f[it + 1];
        }
        else  // For edge 3, use corners 3 and 0
        {
            f1 = f[it];
            f2 = f[0];
        }
        
        // Check if the contour line crosses the current edge
        if ((f1 < c && f2 >= c) || (f2 < c && f1 >= c))
        {
            // Compute the intersection point of the contour line along the edge
            float x = (c - f1) / (f2 - f1);// Linear interpolation factor
            vec2 pIso = vec2((coord[it].x + x * dir[it].x) / (dims.x - 1), (coord[it].y + x * dir[it].y) / (dims.y - 1));
            p.push_back(pIso);// Store the intersection point
        }
    }

    // Draw isoline segment(s)
    int count = p.size(); //Number of intersection points (should be 2 or 4)
    if (count == 2) // If there are two intersection points, simply draw a line between them
    {
        drawLineSegment(p[0], p[1], color, isoBufferGrid, vertices);
    }
    else if (count == 4) // If there are four intersection points, there is ambiguity
    {
        // Handle ambiguity: decide how to connect the points using a "random" or "asymptotic" decider
        if (propDeciderType.get() == 1)    // Random decider: randomly connect the points
        {
            float randomVal = this->randomValue(0.0f, 1.0f);// generate a random value

            if (this->randomValue(0.0f, 1.0f) < 0.5f) 
            {
                drawLineSegment(p[0], p[1], color, isoBufferGrid, vertices);
                drawLineSegment(p[2], p[3], color, isoBufferGrid, vertices);
            }
            else
            {
                drawLineSegment(p[0], p[3], color, isoBufferGrid, vertices);
                drawLineSegment(p[1], p[2], color, isoBufferGrid, vertices);
            }
        }
        else  // Asymptotic decider: resolve the ambiguity using a calculated value
        {
            // Calculate a value 'fab' that helps decide how to connect the points
            float fab = (f00 * f11 - f10 * f01) / (f11 + f00 - f01 - f10);

            // Based on 'fab', decide how to connect the points
            if (((c >= fab) && (f00 >= c)) || ((c < fab) && (f00 < c)))
            {
                drawLineSegment(p[0], p[3], color, isoBufferGrid, vertices);
                drawLineSegment(p[1], p[2], color, isoBufferGrid, vertices);
            }
            else
            {
                drawLineSegment(p[0], p[1], color, isoBufferGrid, vertices);
                drawLineSegment(p[2], p[3], color, isoBufferGrid, vertices);
            }
        }
    }
}

double MarchingSquares::getInputValue(const VolumeRAM* data, const size3_t dims, const size_t i, const size_t j)
{
    // Check if the indices are withing the dimensions of the volume
    if (i < dims.x && j < dims.y)
    {
        return data->getAsDouble(size3_t(i, j, 0));
    }
    else
    {
        LogProcessorError("Attempting to access data outside the boundaries of the volume, value is set to 0");
        return 0;
    }
}


void MarchingSquares::drawLineSegment(const vec2& v1, const vec2& v2, const vec4& color, IndexBufferRAM* indexBuffer, std::vector<BasicMesh::Vertex>& vertices)
{
    // Add first vertex
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    // A vertex has a position, a normal, a texture coordinate and a color
    // we do not use normal or texture coordinate, but still have to specify them
    vertices.push_back({vec3(v1[0], v1[1], 0), vec3(0, 0, 1), vec3(v1[0], v1[1], 0), color});
    // Add second vertex
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    vertices.push_back({vec3(v2[0], v2[1], 0), vec3(0, 0, 1), vec3(v2[0], v2[1], 0), color});
}


}// namespace inviwo
