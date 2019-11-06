#include "BinaryCurveCollectionIO.h"
#include <Ziran/CS/Util/Logging.h>
#include <Ziran/CS/Util/Debug.h>

namespace ZIRAN {

template <class T, int dim>
void readBinaryCurveCollectionFile(const std::string& filename, StdVector<Vector<T, dim>>& X, StdVector<Vector<int, 2>>& segments)
{
    ZIRAN_ASSERT(dim == 3, "BCC file only supports 3d");
    struct BCCHeader {
        char sign[3];
        unsigned char byteCount;
        char curveType[2];
        char dimensions;
        char upDimension;
        uint64_t curveCount;
        uint64_t totalControlPointCount;
        char fileInfo[40];
    };

    struct p3f {
        float data[3];
    };

    BCCHeader header;
    FILE* pFile = fopen(filename.c_str(), "rb");
    ZIRAN_ASSERT(pFile);
    fread(&header, sizeof(header), 1, pFile);
    ZIRAN_ASSERT(header.sign[0] == 'B' && header.sign[1] == 'C' && header.sign[2] == 'C' && header.byteCount == 0x44);

    std::vector<p3f> controlPoints(header.totalControlPointCount);
    std::vector<int> firstControlPoint(header.curveCount + 1);
    std::vector<char> isCurveLoop(header.curveCount);

    std::vector<float> cc;
    cc.resize(controlPoints.size() * 3);
    float* cp = &cc[0];

    int prevCP = 0;
    ZIRAN_INFO("control points size:", controlPoints.size());
    ZIRAN_INFO("curveCount: ", header.curveCount);

    for (uint64_t i = 0; i < header.curveCount; i++) {
        int curveControlPointCount;
        fread(&curveControlPointCount, sizeof(int), 1, pFile);
        ZIRAN_INFO(i, " curveControlPointCount: ", curveControlPointCount);

        isCurveLoop[i] = curveControlPointCount < 0;
        if (curveControlPointCount < 0)
            curveControlPointCount = -curveControlPointCount;
        fread(cp, sizeof(float), curveControlPointCount * 3, pFile);
        cp += curveControlPointCount * 3;
        firstControlPoint[i] = prevCP;
        prevCP += curveControlPointCount;
    }
    firstControlPoint[header.curveCount] = prevCP;

    for (size_t k = 0; k < controlPoints.size(); k++) {
        controlPoints[k].data[0] = cc[3 * k];
        controlPoints[k].data[1] = cc[3 * k + 1];
        controlPoints[k].data[2] = cc[3 * k + 2];
    }

    segments.clear();
    X.clear();

    for (uint64_t i = 0; i < header.curveCount; i++) {

        bool loop = isCurveLoop[i];

        int curve_point_count = firstControlPoint[i + 1] - firstControlPoint[i];

        for (int k = 0; k < curve_point_count; ++k) {
            Vector<float, dim> this_point;
            this_point << controlPoints[firstControlPoint[i] + k].data[0], controlPoints[firstControlPoint[i] + k].data[1], controlPoints[firstControlPoint[i] + k].data[2];
            Vector<T, dim> P;
            P << this_point[0], this_point[1], this_point[2];
            X.emplace_back(P);
            if (k > 0) {
                Vector<int, 2> seg;
                seg << firstControlPoint[i] + k - 1, firstControlPoint[i] + k;
                segments.emplace_back(seg);
            }
        }
        if (loop) {
            Vector<int, 2> seg;
            seg << (int)(firstControlPoint[i] + curve_point_count - 1), firstControlPoint[i];
            segments.emplace_back(seg);
        }
    }

    fclose(pFile);
}
} // namespace ZIRAN
