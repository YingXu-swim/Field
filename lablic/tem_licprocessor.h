/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Monday, October 02, 2017 - 13:31:17
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#pragma once

#include <inviwo/core/common/inviwo.h>
#include <inviwo/core/datastructures/image/imageram.h>
#include <inviwo/core/ports/imageport.h>
#include <inviwo/core/ports/volumeport.h>
#include <inviwo/core/processors/processor.h>
#include <inviwo/core/properties/boolproperty.h>
#include <inviwo/core/properties/ordinalproperty.h>
#include <inviwo/core/properties/transferfunctionproperty.h>
#include <lablic/lablicmoduledefine.h>
#include <labutils/scalarvectorfield.h>
#include <labutils/rgbaimage.h>

namespace inviwo {

/** \docpage{org.inviwo.LICProcessor, LICProcessor}
    ![](org.inviwo.LICProcessor.png?classIdentifier=org.inviwo.LICProcessor)

    Line Integral Convolution with a box kernel.

    ### Inports
      * __vectorField__ 2-dimensional vector field (with vectors of
      two components thus two values within each voxel)
      This processor deals with 2-dimensional data only, therefore it is assumed
      the z-dimension will have size 1 otherwise the 0th slice of the volume
      will be processed.
      * __texture__ Texture to be convolved along the streamlines.

    ### Outports
      * __image__ The image resulting from smearing the given texture
      the streamlines of the given vector field.
*/
class IVW_MODULE_LABLIC_API LICProcessor : public Processor {
    // Friends
    // Types
public:
    // Construction / Deconstruction
public:
    LICProcessor();
    virtual ~LICProcessor() = default;

    // Methods
public:
    virtual const ProcessorInfo getProcessorInfo() const override;
    static const ProcessorInfo processorInfo_;

protected:
    /// Our main computation function
    virtual void process() override;

    // (TODO: Helper functions can be defined here and then implemented in the .cpp)
    // e.g. something like a function for standardLIC, fastLIC, autoContrast, ...
    std::vector<dvec2>& streamLineLIC(const VectorField2& vectorField, int kernelSize, const dvec2& position, double step, int skipNPoints = 0);
    void LIC(RGBAImage& licImage, const RGBAImage& texture, const VectorField2& vectorField, int kernelSize, double step, int skipNPoints);

    void textureColoring(RGBAImage& outImage, std::vector<std::vector<double>>& texture, std::vector<std::vector<double>>& magnitudeAtPoints, const double maxMagnitude);


    double LICProcessor::maxMagnitude(const VectorField2& vectorField, dvec2& scaleSlow);

    // Ports
public:
    // Input vector field
    VolumeInport volumeIn_;

    // Input texture
    ImageInport noiseTexIn_;

    // Output image
    ImageOutport licOut_;

    // Properties
public:
    // TODO: Declare properties
    // IntProperty prop1;
    // BoolProperty prop2;
    IntProperty propKernelSize;
    IntProperty propSkipNPoints;
    FloatProperty propStepSize;
    BoolProperty proptextureColor;
    TransferFunctionProperty propTransferFunc;
    /*
    IntProperty propNumSteps;
    FloatProperty propStepSize;// For controlling step size
    BoolProperty propnormalizeVectorField;// Normalize the vector field for integration
    FloatProperty propmaxArcLength;       // Max arc length for integration
    FloatProperty propminVelocity;        // Threshold to stop if velocity becomes too slow
    BoolProperty propStopWhenZero;
    */
    // Attributes
private:
    size3_t vectorFieldDims_;
    size2_t texDims_;
};

}  // namespace inviwo
