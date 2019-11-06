#ifndef RANDOM_NUMBER_H
#define RANDOM_NUMBER_H

#include <random>
#include <Ziran/CS/Util/Forward.h>

namespace ZIRAN {

/**
    Random number generator.
*/
template <class T>
class RandomNumber {
public:
    std::mt19937 generator;

    RandomNumber(unsigned s = 123);

    ~RandomNumber();

    /**
       Reset the random seed.
    */
    void resetSeed(T s);

    /**
       Reset seed using time
    */
    void resetSeedUsingTime();

    /**
       Random real number from 0 to 1
    */
    T randReal();

    /**
       Random real number from an interval
    */
    T randReal(T a, T b);

    /**
       Random integer number from a to b inclusive, i.e. in the both closed interval [a,b]
    */
    int randInt(int a, int b);

    /**
      Random vector in box
      */
    template <int d, int flags>
    Vector<T, d> randInBox(const Eigen::Matrix<T, d, 1, flags, d, 1>& min_corner, const Eigen::Matrix<T, d, 1, flags, d, 1>& max_corner);

    /**
       Random barycentric weights
    */
    template <int d>
    Vector<T, d> randomBarycentricWeights();

    /**
      Random vector in ball
      */
    template <int d, int flags>
    Vector<T, d> randInBall(const Eigen::Matrix<T, d, 1, flags, d, 1>& center, T radius);

    /**
       Random rotation matrix
    */
    void randRotation(Matrix<T, 2, 2>& R);

    /**
      Random rotation matrix
      */
    void randRotation(Matrix<T, 3, 3>& R);

    /**
      Fill with uniform random numbers
    */
    template <class T2>
    void fill(T2& x, T a = 0, T b = 1);

    /**
      Fill with random integers
    */
    template <class T2>
    void fillInt(T2& x, int a = 0, int b = 1);
};
} // namespace ZIRAN
#endif
