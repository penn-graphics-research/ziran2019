#include <Ziran/CS/Util/RandomNumber.h>
#include <Eigen/Geometry>
#include <chrono>
#include <math.h>

namespace ZIRAN {

template <class T>
RandomNumber<T>::RandomNumber(unsigned s)
    : generator(s)
{
}
template <class T>
RandomNumber<T>::~RandomNumber()
{
}

template <class T>
void RandomNumber<T>::resetSeed(T s)
{
    generator.seed(s);
}

template <class T>
void RandomNumber<T>::resetSeedUsingTime()
{
    auto s = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    generator.seed(s);
}

template <class T>
T RandomNumber<T>::randReal()
{
    std::uniform_real_distribution<T> distribution((T)0, (T)1);
    return distribution(generator);
}

template <class T>
T RandomNumber<T>::randReal(T a, T b)
{
    std::uniform_real_distribution<T> distribution(a, b);
    return distribution(generator);
}

template <class T>
int RandomNumber<T>::randInt(int a, int b)
{
    std::uniform_int_distribution<> distribution(a, b);
    return distribution(generator);
}

template <class T>
template <int d, int flags>
Vector<T, d> RandomNumber<T>::randInBox(const Eigen::Matrix<T, d, 1, flags, d, 1>& min_corner, const Eigen::Matrix<T, d, 1, flags, d, 1>& max_corner)
{
    Vector<T, d> r;
    for (int i = 0; i < d; i++)
        r(i) = randReal(min_corner(i), max_corner(i));
    return r;
}

template <class T>
template <int d>
Vector<T, d> RandomNumber<T>::randomBarycentricWeights()
{
    Vector<T, d> r;
    T sum;
    do {
        sum = 0;
        for (int i = 0; i < d - 1; i++) {
            r(i) = randReal();
            sum += r(i);
        }
    } while (sum > 1);
    r(d - 1) = 1 - sum;
    return r;
}

template <class T>
template <int d, int flags>
Vector<T, d> RandomNumber<T>::randInBall(const Eigen::Matrix<T, d, 1, flags, d, 1>& center, T radius)
{
    Vector<T, d> min_corner = center.array() - radius;
    Vector<T, d> max_corner = center.array() + radius;

    Vector<T, d> r;
    do {
        r = randInBox(min_corner, max_corner);
    } while ((r - center).squaredNorm() > radius * radius);
    return r;
}

template <class T>
void RandomNumber<T>::randRotation(Matrix<T, 2, 2>& R)
{
    T theta = randReal(0, 2 * M_PI);
    T c = std::cos(theta);
    T s = std::sin(theta);
    R << c, -s, s, c;
}

template <class T>
void RandomNumber<T>::randRotation(Matrix<T, 3, 3>& R)
{
    std::normal_distribution<T> n;
    Eigen::Quaternion<T> q(n(generator), n(generator), n(generator), n(generator));
    q.normalize();
    R = q.toRotationMatrix();
}

template <class T, class Derived>
static void fillHelper(RandomNumber<T>& rand, Eigen::DenseBase<Derived>& x, T a, T b)
{
    for (typename Derived::Index i = 0; i < x.size(); i++)
        x(i) = rand.randReal(a, b);
}

template <class T, class T2>
static std::enable_if_t<std::is_arithmetic<T2>::value> fillHelper(RandomNumber<T>& rand, T2& x, T a, T b)
{
    x = rand.randReal(a, b);
}

template <class T>
template <class T2>
void RandomNumber<T>::fill(T2& x, T a, T b)
{
    fillHelper(*this, x, a, b);
}

template <class T, class Derived>
static void fillIntHelper(RandomNumber<T>& rand, Eigen::DenseBase<Derived>& x, int a, int b)
{
    for (typename Derived::Index i = 0; i < x.size(); i++)
        x(i) = rand.randInt(a, b);
}

template <class T, class T2>
static std::enable_if_t<std::is_arithmetic<T2>::value> fillIntHelper(RandomNumber<T>& rand, T2& x, int a, int b)
{
    x = rand.randInt(a, b);
}

template <class T>
template <class T2>
void RandomNumber<T>::fillInt(T2& x, int a, int b)
{
    fillIntHelper(*this, x, a, b);
}

template Eigen::Matrix<double, 2, 1, 0, 2, 1> RandomNumber<double>::randInBall<2, 0>(Eigen::Matrix<double, 2, 1, 0, 2, 1> const&, double);
template Eigen::Matrix<double, 2, 1, 0, 2, 1> RandomNumber<double>::randInBall<2, 2>(Eigen::Matrix<double, 2, 1, 2, 2, 1> const&, double);
template Eigen::Matrix<double, 3, 1, 0, 3, 1> RandomNumber<double>::randInBox<3, 0>(Eigen::Matrix<double, 3, 1, 0, 3, 1> const&, Eigen::Matrix<double, 3, 1, 0, 3, 1> const&);
template Eigen::Matrix<double, 3, 1, 0, 3, 1> RandomNumber<double>::randInBall<3, 0>(Eigen::Matrix<double, 3, 1, 0, 3, 1> const&, double);
template Eigen::Matrix<double, 3, 1, 0, 3, 1> RandomNumber<double>::randInBall<3, 2>(Eigen::Matrix<double, 3, 1, 2, 3, 1> const&, double);
template Eigen::Matrix<double, 3, 1, 0, 3, 1> RandomNumber<double>::randomBarycentricWeights<3>();
template Eigen::Matrix<float, 2, 1, 0, 2, 1> RandomNumber<float>::randInBall<2, 0>(Eigen::Matrix<float, 2, 1, 0, 2, 1> const&, float);
template Eigen::Matrix<float, 2, 1, 0, 2, 1> RandomNumber<float>::randInBall<2, 2>(Eigen::Matrix<float, 2, 1, 2, 2, 1> const&, float);
template Eigen::Matrix<float, 2, 1, 0, 2, 1> RandomNumber<float>::randInBox<2, 0>(Eigen::Matrix<float, 2, 1, 0, 2, 1> const&, Eigen::Matrix<float, 2, 1, 0, 2, 1> const&);
template Eigen::Matrix<float, 3, 1, 0, 3, 1> RandomNumber<float>::randInBall<3, 0>(Eigen::Matrix<float, 3, 1, 0, 3, 1> const&, float);
template Eigen::Matrix<float, 3, 1, 0, 3, 1> RandomNumber<float>::randInBall<3, 2>(Eigen::Matrix<float, 3, 1, 2, 3, 1> const&, float);
template Eigen::Matrix<float, 3, 1, 0, 3, 1> RandomNumber<float>::randInBox<3, 0>(Eigen::Matrix<float, 3, 1, 0, 3, 1> const&, Eigen::Matrix<float, 3, 1, 0, 3, 1> const&);
template Eigen::Matrix<float, 3, 1, 0, 3, 1> RandomNumber<float>::randomBarycentricWeights<3>();
template Vector<double, 2> RandomNumber<double>::randInBox<2, 0>(Eigen::Matrix<double, 2, 1, 0, 2, 1> const&, Eigen::Matrix<double, 2, 1, 0, 2, 1> const&);
template class RandomNumber<double>;
template class RandomNumber<float>;
template void RandomNumber<double>::fill<double>(double&, double, double);
template void RandomNumber<float>::fill<float>(float&, float, float);
template void RandomNumber<double>::fill<Eigen::Matrix<double, 2, 1, 0, 2, 1>>(Eigen::Matrix<double, 2, 1, 0, 2, 1>&, double, double);
template void RandomNumber<double>::fill<Eigen::Matrix<double, 3, 1, 0, 3, 1>>(Eigen::Matrix<double, 3, 1, 0, 3, 1>&, double, double);
template void RandomNumber<double>::fill<Eigen::Matrix<double, -1, -1, 0, -1, -1>>(Eigen::Matrix<double, -1, -1, 0, -1, -1>&, double, double);
template void RandomNumber<double>::fill<Eigen::Matrix<double, 2, -1, 0, 2, -1>>(Eigen::Matrix<double, 2, -1, 0, 2, -1>&, double, double);
template void RandomNumber<double>::fill<Eigen::Matrix<double, 2, 2, 0, 2, 2>>(Eigen::Matrix<double, 2, 2, 0, 2, 2>&, double, double);
template void RandomNumber<double>::fill<Eigen::Matrix<double, 3, -1, 0, 3, -1>>(Eigen::Matrix<double, 3, -1, 0, 3, -1>&, double, double);
template void RandomNumber<double>::fill<Eigen::Matrix<double, 3, 3, 0, 3, 3>>(Eigen::Matrix<double, 3, 3, 0, 3, 3>&, double, double);
template void RandomNumber<double>::fill<Eigen::Matrix<double, 3, 2, 0, 3, 2>>(Eigen::Matrix<double, 3, 2, 0, 3, 2>&, double, double);
template void RandomNumber<double>::fill<Eigen::Matrix<double, 4, 2, 0, 4, 2>>(Eigen::Matrix<double, 4, 2, 0, 4, 2>&, double, double);
template void RandomNumber<double>::fill<Eigen::Matrix<double, 4, 3, 0, 4, 3>>(Eigen::Matrix<double, 4, 3, 0, 4, 3>&, double, double);
template void RandomNumber<double>::fill<Eigen::Matrix<double, 8, 2, 0, 8, 2>>(Eigen::Matrix<double, 8, 2, 0, 8, 2>&, double, double);
template void RandomNumber<double>::fill<Eigen::Matrix<double, 7, 3, 0, 7, 3>>(Eigen::Matrix<double, 7, 3, 0, 7, 3>&, double, double);
template void RandomNumber<double>::fill<Eigen::Matrix<double, 8, 3, 0, 8, 3>>(Eigen::Matrix<double, 8, 3, 0, 8, 3>&, double, double);
template void RandomNumber<double>::fill<Eigen::Matrix<double, 26, 3, 0, 26, 3>>(Eigen::Matrix<double, 26, 3, 0, 26, 3>&, double, double);
template void RandomNumber<double>::fillInt<Eigen::Matrix<double, 2, 2, 0, 2, 2>>(Eigen::Matrix<double, 2, 2, 0, 2, 2>&, int, int);
template void RandomNumber<double>::fillInt<Eigen::Matrix<double, 3, 3, 0, 3, 3>>(Eigen::Matrix<double, 3, 3, 0, 3, 3>&, int, int);

template void RandomNumber<float>::fill<Eigen::Array<float, 10000, 1, 0, 10000, 1>>(Eigen::Array<float, 10000, 1, 0, 10000, 1>&, float, float);
template void RandomNumber<float>::fill<Eigen::Matrix<float, 3, 1, 0, 3, 1>>(Eigen::Matrix<float, 3, 1, 0, 3, 1>&, float, float);
template void RandomNumber<float>::fill<Eigen::Matrix<float, 3, 2, 0, 3, 2>>(Eigen::Matrix<float, 3, 2, 0, 3, 2>&, float, float);
template void RandomNumber<float>::fill<Eigen::Matrix<float, -1, 1, 0, -1, 1>>(Eigen::Matrix<float, -1, 1, 0, -1, 1>&, float, float);
template void RandomNumber<float>::fill<Eigen::Matrix<float, 1, 1, 0, 1, 1>>(Eigen::Matrix<float, 1, 1, 0, 1, 1>&, float, float);
template void RandomNumber<float>::fill<Eigen::Matrix<float, 2, 1, 0, 2, 1>>(Eigen::Matrix<float, 2, 1, 0, 2, 1>&, float, float);
template void RandomNumber<float>::fill<Eigen::Matrix<float, 2, 2, 0, 2, 2>>(Eigen::Matrix<float, 2, 2, 0, 2, 2>&, float, float);
template void RandomNumber<float>::fill<Eigen::Matrix<float, 3, 3, 0, 3, 3>>(Eigen::Matrix<float, 3, 3, 0, 3, 3>&, float, float);
template void RandomNumber<float>::fill<Eigen::Matrix<float, 4, 2, 0, 4, 2>>(Eigen::Matrix<float, 4, 2, 0, 4, 2>&, float, float);
template void RandomNumber<float>::fill<Eigen::Matrix<float, 4, 3, 0, 4, 3>>(Eigen::Matrix<float, 4, 3, 0, 4, 3>&, float, float);
template void RandomNumber<float>::fill<Eigen::Matrix<float, 8, 2, 0, 8, 2>>(Eigen::Matrix<float, 8, 2, 0, 8, 2>&, float, float);
template void RandomNumber<float>::fill<Eigen::Matrix<float, 8, 3, 0, 8, 3>>(Eigen::Matrix<float, 8, 3, 0, 8, 3>&, float, float);
template void RandomNumber<float>::fill<Eigen::Matrix<float, 7, 3, 0, 7, 3>>(Eigen::Matrix<float, 7, 3, 0, 7, 3>&, float, float);
template void RandomNumber<float>::fill<Eigen::Matrix<float, 26, 3, 0, 26, 3>>(Eigen::Matrix<float, 26, 3, 0, 26, 3>&, float, float);
template void RandomNumber<float>::fillInt<Eigen::Matrix<float, 2, 2, 0, 2, 2>>(Eigen::Matrix<float, 2, 2, 0, 2, 2>&, int, int);
template void RandomNumber<float>::fillInt<Eigen::Matrix<float, 3, 3, 0, 3, 3>>(Eigen::Matrix<float, 3, 3, 0, 3, 3>&, int, int);
} // namespace ZIRAN
