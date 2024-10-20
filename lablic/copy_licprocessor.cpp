/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Monday, October 02, 2017 - 13:31:17
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <inviwo/core/datastructures/volume/volumeram.h>
#include <lablic/licprocessor.h>
#include <labstreamlines/integrator.h>

namespace inviwo {

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo LICProcessor::processorInfo_{
    "org.inviwo.LICProcessor",  // Class identifier
    "LICProcessor",             // Display name
    "KTH Labs",                 // Category
    CodeState::Experimental,    // Code state
    Tags::None,                 // Tags
};

const ProcessorInfo LICProcessor::getProcessorInfo() const { return processorInfo_; }

LICProcessor::LICProcessor()
    : Processor()
    , volumeIn_("volIn")
    , noiseTexIn_("noiseTexIn")
    , licOut_("licOut")
// TODO: Register additional properties
    , kernelSize("kernelSize", "Kernel Size", 20, 1, 200, 1)
    , fastLIC("fastLIC", "Use FastLIC")
    , textureColor("textureColor", "Use Magnitude of VF for Texture Color")
    , contrastEnhancement("contrastEnhancement", "Contrast Enhancement")
    , my("my", "Mean", 0.0, 0.0, 1.0)
    , sigma("sigma", "Standard Deviation", 0.0, 0.0, 1.0)
    , propTransferFunc("isoTransferFunc", "Colors", &volumeIn_)
{
    // Register ports
    addPort(volumeIn_);
    addPort(noiseTexIn_);
    addPort(licOut_);

    // Register properties
    // TODO: Register additional properties
    addProperty(kernelSize);
    addProperty(fastLIC);
    addProperty(textureColor);
    fastLIC.set(true);
    addProperty(contrastEnhancement);
    addProperty(my);
    addProperty(sigma);
    
    addProperty(propTransferFunc);
    propTransferFunc.get().clear();
    propTransferFunc.get().add(0.0f, vec4(0.0f, 0.0f, 1.0f, 1.0f));
    propTransferFunc.get().add(0.5f, vec4(0.0f, 1.0f, 0.0f, 1.0f));
    propTransferFunc.get().add(1.0f, vec4(1.0f, 0.0f, 0.0f, 1.0f));
    propTransferFunc.setCurrentStateAsDefault();
    
}

void LICProcessor::process() {
    // Get input
    if (!volumeIn_.hasData()) {
        return;
    }

    if (!noiseTexIn_.hasData()) {
        return;
    }

    // 从端口获取体积数据 vol 和噪声纹理 tex
    auto vol = volumeIn_.getData();
    const VectorField2 vectorField = VectorField2::createFieldFromVolume(vol);
    vectorFieldDims_ = vol->getDimensions();

    auto tex = noiseTexIn_.getData();
    const RGBAImage texture = RGBAImage::createFromImage(tex);
    texDims_ = tex->getDimensions();


    // Prepare the output, it has the same dimensions as the texture and rgba values in [0,255]
    auto outImage = std::make_shared<Image>(texDims_, DataVec4UInt8::get());
    RGBAImage licImage(outImage);

        // 创建二维向量 licTexture，用于存储每个像素的 LIC 值。
    std::vector<std::vector<double>> licTexture(texDims_.x, std::vector<double>(texDims_.y));
        // 创建访问状态矩阵 visited，用于记录在算法中哪些像素被访问过 [used in fastlic]
    std::vector<std::vector<bool>> visited(texDims_.x, std::vector<bool>(texDims_.y, false));

    // TODO: Implement LIC and FastLIC
    // This code instead just creates a black image

    std::vector<dvec2> listForward;
    std::vector<dvec2> listBackward;

    // 计算缩放因子
        // 计算每个像素对应的向量场宽度 pxWidth 和高度 pxHeight
    double pxWidth = (double)(vectorField.getBBoxMax().x - vectorField.getBBoxMin().x) / texDims_.x,
           pxHeight =
               (double)(vectorField.getBBoxMax().y - vectorField.getBBoxMin().y) / texDims_.y;
        // 计算步长 step，用于确定流线的步长
    double step = pxWidth > pxHeight ? pxWidth : pxHeight;
        // 获取向量场的最小坐标 vf_bboxmin
    dvec2 vf_bboxmin = vectorField.getBBoxMin();
        // 计算缩放因子 scaleFactors 和 scaleSlow，用于在向量场和纹理坐标之间转换
    dvec2 scaleFactors = {
        (texDims_.x - 1) / (vectorField.getBBoxMax().x - vectorField.getBBoxMin().x),
        (texDims_.y - 1) / (vectorField.getBBoxMax().y - vectorField.getBBoxMin().y)};
    dvec2 scaleSlow = {
        (vectorField.getBBoxMax().x - vectorField.getBBoxMin().x) / (texDims_.x - 1),
        (vectorField.getBBoxMax().y - vectorField.getBBoxMin().y) / (texDims_.y - 1)};

    // 映射纹理到向量场
        // 创建映射矩阵 vectorFieldMapX 和 vectorFieldMapY，用于存储纹理坐标在向量场中的对应位置。
    std::vector<std::vector<int>> interVarForward, interVarBackward;
        // 映射纹理到向量场

    std::vector<std::vector<double>> vectorFieldMapX(texDims_.x, std::vector<double>(texDims_.y)),
        vectorFieldMapY(texDims_.x, std::vector<double>(texDims_.y));

    int pixelX, pixelY;
    double kernelValue = 1.0 / (2 * kernelSize + 1);
    double conv = 0.0;
        // magnitudeAtPoints 矩阵，用于存储每个纹理坐标点的向量场幅度
    std::vector<std::vector<double>> magnitudeAtPoints(texDims_.x, std::vector<double>(texDims_.y));
        //  mapTextureToVectorField 函数，将纹理坐标映射到向量场，并计算最大幅度 maxMagnitude。
    double maxMagnitude =
        mapTextureToVectorField(vectorField, scaleSlow, vf_bboxmin, texDims_, vectorFieldMapX,
                                vectorFieldMapY, magnitudeAtPoints);

    //Values for contrast enhancement
    double mean = 0.0, P = 0.0;
    double non_black_cnt = 0.0;  // keep count of all non-black pixels

    for (auto i = 0; i < (int)texDims_.x; i++) {
        for (auto j = 0; j < (int)texDims_.y; j++) { 

            if (fastLIC == false) { 
                // 1调用 Integrator::StreamLines 函数计算前向和后向流线
                // 2 调用 convolution 函数计算当前像素点的 LIC 值。 
                // 将计算结果设置到 licImage 中
                listForward = Integrator::StreamLines(
                    vectorField, {vectorFieldMapX[i][j], vectorFieldMapY[i][j]}, step, 0, true, -1,
                    kernelSize);
                listBackward = Integrator::StreamLines(
                    vectorField, {vectorFieldMapX[i][j], vectorFieldMapY[i][j]}, step, 1, true, -1,
                    kernelSize);
                

                licTexture[i][j] = convolution(vectorField, listForward, listBackward, texture, 1);

                int val = int(licTexture[i][j]);
                licImage.setPixel(size2_t(i, j), dvec4(val, val, val, 255));

                // If enabled, start computing mean and P in the for
                // loop, to not do it outside
                if (contrastEnhancement &&
                    licTexture[i][j] > 1.0) {
                    mean += licTexture[i][j];
                    P += licTexture[i][j] * licTexture[i][j];
                    non_black_cnt++;
                }

            } else { // fast LIC
                if (visited[i][j] == false) {
                    // 前向和后向流线
                    listForward = Integrator::StreamLines(
                        vectorField, {vectorFieldMapX[i][j], vectorFieldMapY[i][j]}, step, 0, true,
                        1e-6, 2000);
                    listBackward = Integrator::StreamLines(
                        vectorField, {vectorFieldMapX[i][j], vectorFieldMapY[i][j]}, step, 1, true,
                        1e-6, 2000);
                    // 将流线映射回纹理坐标。
                    interVarForward =
                        mapVectorToTextureField(listForward, vf_bboxmin, scaleFactors);
                    interVarBackward =
                        mapVectorToTextureField(listBackward, vf_bboxmin, scaleFactors);
                    // 使用卷积核值 kernelValue 和纹理图像 texture 计算卷积和
                    conv = 0.0;

                    int k = -(int)listBackward.size() + 1 - kernelSize;
                    while (k < (int)listForward.size()) {

                        if (k + kernelSize < 0 && -k - kernelSize < (int)listBackward.size())
                        {
                            conv += kernelValue * texture.readPixelGrayScale({interVarBackward[-k - kernelSize][0], interVarBackward[-k - kernelSize][1]});
                        }
                        else if (k + kernelSize >= 0 && k + kernelSize < (int)listForward.size())
                        {
                            conv += kernelValue * texture.readPixelGrayScale({interVarForward[k + kernelSize][0], interVarForward[k + kernelSize][1]});
                        }

                        if (k - kernelSize < 0 && -k + kernelSize < (int)listBackward.size())
                        {
                            conv -= kernelValue * texture.readPixelGrayScale({interVarBackward[-k + kernelSize][0], interVarBackward[-k + kernelSize][1]});
                        }
                        else if (k - kernelSize >= 0 && k - kernelSize < (int)listForward.size())
                        {
                            conv -= kernelValue * texture.readPixelGrayScale({interVarForward[k - kernelSize][0], interVarForward[k - kernelSize][1]});
                        }

                        if (k < 0 && -k < (int)interVarBackward.size() /*&& -k + kernelSize < (int)interVarBackward.size()
                                && kernelSize + k < (int)interVarForward.size()*/) {
                            pixelX = interVarBackward[-k][0];
                            pixelY = interVarBackward[-k][1];
                        } else if (k >= 0 && k < (int)interVarForward.size() /*&& k + kernelSize < (int)interVarForward.size()
                                       && kernelSize - k < (int)interVarBackward.size()*/) {
                            pixelX = interVarForward[k][0];
                            pixelY = interVarForward[k][1];
                        } else {
                            k++;
                            continue;
                        }

                        if (!visited[pixelX][pixelY]) {
                            licTexture[pixelX][pixelY] = conv;
                            int val = int(licTexture[pixelX][pixelY]);
                            licImage.setPixel(size2_t(pixelX, pixelY), dvec4(val, val, val, 255));
                            visited[pixelX][pixelY] = true;

                            // If enabled, start computing mean and P
                            // in the for loop, to not do it outside
                            if (contrastEnhancement &&
                                licTexture[i][j] > 1.0) {
                                mean += licTexture[i][j];
                                P += licTexture[i][j] * licTexture[i][j];
                                non_black_cnt++;
                            }
                        }

                        k++;
                    }
                }
            }
        }
    }
    //Call to contrast enhancement
    if (contrastEnhancement) {  // If enabled, compute the contrast thingy
        mean /= non_black_cnt;
        double sigma = sqrt((P - non_black_cnt * mean * mean) / (non_black_cnt - 1));
        contrastEnhance(licImage, licTexture, mean, sigma, texDims_);
    }

    //Call to map gray-scale to RGB or something, after contrast has been enhanced.
    if (textureColor) {
        textureColoring(licImage, licTexture, magnitudeAtPoints, maxMagnitude);

    }

    licOut_.setData(outImage);
}
// 将纹理坐标映射到向量场，并计算最大幅度。
double LICProcessor::mapTextureToVectorField(const VectorField2& vectorField, const dvec2 scaleSlow,
                                             const dvec2 vf_bboxmin, dvec2 texDims_,
                                             std::vector<std::vector<double>>& vectorListX,
                                             std::vector<std::vector<double>>& vectorListY,
                                             std::vector<std::vector<double>>& magnitudeVector) {

    double mapX, mapY;
    double maxMagnitude = 0.0;
    double magnitude;
    for (auto i = 0; i < (int)texDims_.x; i++) {
        for (auto j = 0; j < (int)texDims_.y; j++) {
            mapX = ((double)i) * scaleSlow.x + vf_bboxmin.x;
            mapY = ((double)j) * scaleSlow.y + vf_bboxmin.y;
            vectorListX[i][j] = mapX;
            vectorListY[i][j] = mapY;
            dvec2 pos = vectorField.interpolate({mapX, mapY});
            magnitude = sqrt(pow(pos.x, 2) + pow(pos.y, 2));
            magnitudeVector[i][j] = magnitude;
            if (magnitude > maxMagnitude) {
                maxMagnitude = magnitude;
            }
        }
    }
    return maxMagnitude;
}

std::vector<std::vector<int>> LICProcessor::mapVectorToTextureField(
    const std::vector<dvec2>& vectorList, const dvec2& vf_minbbox, const dvec2& scaleFactors) {
    std::vector<std::vector<int>> pixelList((int)vectorList.size(), std::vector<int>(2));
    double mapX, mapY;
    for (int i = 0; i < (int)vectorList.size(); i++) {
        mapX = (vectorList[i].x - vf_minbbox.x) * scaleFactors.x;
        mapY = (vectorList[i].y - vf_minbbox.y) * scaleFactors.y;

        pixelList[i][0] = (int)round(mapX);
        pixelList[i][1] = (int)round(mapY);
    }
    return pixelList;
}

double LICProcessor::LinearIntegralConvolution(const RGBAImage& inputImage,
                                               const std::vector<std::vector<int>>& forwardList,
                                               const std::vector<std::vector<int>>& backwardList,
                                               const int& filterStartIndex) {

    double conv = 0.0;

    for (int i = filterStartIndex - kernelSize; i < filterStartIndex + kernelSize; i++) {
        if (i < 0 && -i < (int)backwardList.size()) {
            conv += boxKernel(i - filterStartIndex) *
                    inputImage.readPixelGrayScale({backwardList[-i][0], backwardList[-i][1]});
        } else if (i >= 0 && i < (int)forwardList.size()) {
            conv += boxKernel(i - filterStartIndex) *
                    inputImage.readPixelGrayScale({forwardList[i][0], forwardList[i][1]});
        }
    }
    return conv;
}

double LICProcessor::convolution(const VectorField2& vectorField,
                                 const std::vector<dvec2>& forwardList,
                                 const std::vector<dvec2>& backwardList,
                                 const RGBAImage& inputImage, const double magnitudeRatio) {
    double conv = 0.0;
    double mapX;
    double mapY;

    for (int i = 0; i < (int)forwardList.size(); i++) {
        double kernelValue = boxKernel(i);
        mapX = (forwardList[i].x - vectorField.getBBoxMin().x) * (double)(texDims_.x - 1) /
               (vectorField.getBBoxMax().x - vectorField.getBBoxMin().x);
        mapY = (forwardList[i].y - vectorField.getBBoxMin().y) * (double)(texDims_.y - 1) /
               (vectorField.getBBoxMax().y - vectorField.getBBoxMin().y);
        if (textureColor) {
            conv += kernelValue * magnitudeRatio *
                    inputImage.readPixelGrayScale({(int)mapX, (int)mapY});
        } else {
            conv += kernelValue * inputImage.readPixelGrayScale({(int)mapX, (int)mapY});
        }
    }

    for (int i = 1; i < (int)backwardList.size(); i++) {
        double kernelValue = boxKernel(i);
        mapX = (backwardList[i].x - vectorField.getBBoxMin().x) * (double)(texDims_.x - 1) /
               (vectorField.getBBoxMax().x - vectorField.getBBoxMin().x);
        mapY = (backwardList[i].y - vectorField.getBBoxMin().y) * (double)(texDims_.y - 1) /
               (vectorField.getBBoxMax().y - vectorField.getBBoxMin().y);
        conv +=
            kernelValue * magnitudeRatio * inputImage.readPixelGrayScale({(int)mapX, (int)mapY});
    }

    return conv;
}

double LICProcessor::boxKernel(const int index) {
    if (index > kernelSize)
        return 0.0;
    else
        return 1.0 / (2 * kernelSize + 1);  // Maybe 2*kernelSize + 1.
}

void LICProcessor::contrastEnhance(RGBAImage& outImage, std::vector<std::vector<double>>& texture,
                                   const double mean, const double sig, const dvec2& texDims_) {
    /*
    Enhances the contrast of a texture by computing mean and standard deviation
        input:
                -> texture
                -> mean
                -> sig, standard deviation
        output:
                -> None, change texture.
    */
    double tmp_sigma = sigma * 255,
           tmp_my = my * 255;  // transform STD and mean to pixel values from normalized
    double s_f = tmp_sigma / sig;
    int val;
    for (int x = 0; x < texDims_.x; x++) {
        for (int y = 0; y < texDims_.y; y++) {
            val = (tmp_my + s_f * (texture[x][y] - mean));
            texture[x][y] = val;
            outImage.setPixel(size2_t(x, y), dvec4(val, val, val, 255));
        }
    }
    // return
}

void LICProcessor::textureColoring(RGBAImage& outImage,
                                   std::vector<std::vector<double>>& texture,
                                   std::vector<std::vector<double>> &magnitudeAtPoints,
    const double maxMagnitude) {

    //Go through all pixels (hit me; it's slow, I know)
    int val;
    for (int x = 0; x < texDims_.x; x++) {
        for (int y = 0; y < texDims_.y; y++) {
            val = texture[x][y] * magnitudeAtPoints[x][y] / maxMagnitude; //value is this even needed anymore...?
            texture[x][y] = val;
            double tmp = magnitudeAtPoints[x][y] / maxMagnitude;
            dvec4 color = 255 * propTransferFunc.get().sample(tmp);
            //outImage.setPixel(size2_t(x, y), dvec4(val, val, val, 255));
            outImage.setPixel(size2_t(x, y), (color + dvec4(val, val, val, 255)) * 0.5); //Take half of the lines and half of the colors to get best of both worlds
        }
    }

}

}  // namespace inviwo

