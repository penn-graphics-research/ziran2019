#include "VoronoiNoise.h"
#include <Ziran/Math/MathTools.h>

namespace ZIRAN {

// Generate voronoi noise
float voronoiDistance(const Vector<float, 3>& x, const Vector<int, 3>& seed)
{
    using std::min;
    using namespace ZIRAN::MATH_TOOLS;
    Vector<int, 3> p = int_floor(x);
    Vector<float, 3> f = x - p.cast<float>();
    p += seed;

    Vector<int, 3> mb = Vector<int, 3>::Zero();
    Vector<float, 3> mr = Vector<float, 3>::Zero();

    float res = 8.0;
    for (int k = -1; k <= 1; k++)
        for (int j = -1; j <= 1; j++)
            for (int i = -1; i <= 1; i++) {
                Vector<int, 3> b(i, j, k);
                Vector<float, 3> r = b.template cast<float>() + pseudorand_3f(p + b) - f;
                float d = r.squaredNorm();

                if (d < res) {
                    res = d;
                    mr = r;
                    mb = b;
                }
            }

    res = 8.0;
    for (int k = -2; k <= 2; k++)
        for (int j = -2; j <= 2; j++)
            for (int i = -2; i <= 2; i++) {
                Vector<int, 3> b = mb + Vector<int, 3>(i, j, k);
                Vector<float, 3> r = b.template cast<float>() + pseudorand_3f(p + b) - f;
                float d = 0.5 * (mr + r).dot((r - mr).normalized());

                res = min(res, d);
            }

    return res;
}

// Generate voronoi noise
float voronoiDistance(const Vector<float, 2>& x, const Vector<int, 2>& seed)
{
    using std::min;
    using namespace ZIRAN::MATH_TOOLS;
    Vector<int, 2> p = int_floor(x);
    Vector<float, 2> f = x - p.cast<float>();
    p += seed;

    Vector<int, 2> mb = Vector<int, 2>::Zero();
    Vector<float, 2> mr = Vector<float, 2>::Zero();

    float res = 8.0;
    for (int j = -1; j <= 1; j++)
        for (int i = -1; i <= 1; i++) {
            Vector<int, 2> b(i, j);
            Vector<float, 2> r = b.template cast<float>() + pseudorand_2f(p + b) - f;
            float d = r.squaredNorm();

            if (d < res) {
                res = d;
                mr = r;
                mb = b;
            }
        }

    res = 8.0;
    for (int j = -2; j <= 2; j++)
        for (int i = -2; i <= 2; i++) {
            Vector<int, 2> b = mb + Vector<int, 2>(i, j);
            Vector<float, 2> r = b.template cast<float>() + pseudorand_2f(p + b) - f;
            float d = 0.5 * (mr + r).dot((r - mr).normalized());
            res = min(res, d);
        }

    return res;
}
} // namespace ZIRAN
