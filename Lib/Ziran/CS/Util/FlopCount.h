#ifndef FLOP_COUNT_H
#define FLOP_COUNT_H
#include <Ziran/Math/Linear/DenseExt.h>

namespace ZIRAN {
class Double {
public:
    double val;
    Double()
        : val((double)(0))
    {
    }
    Double(const Double& x)
        : val(x.val)
    {
    }
    Double(double v)
        : val(v)
    {
    }
    static unsigned int mul_count;
    static unsigned int div_count;
    static unsigned int add_count;
    static unsigned int sub_count;

    //reset
    void reset()
    {
        mul_count = 0;
        div_count = 0;
        add_count = 0;
        sub_count = 0;
    }

    //print report
    static void printReport()
    {
        ZIRAN_INFO("Number of multiplications=", mul_count);
        ZIRAN_INFO("Number of divisions=", div_count);
        ZIRAN_INFO("Number of additions=", add_count);
        ZIRAN_INFO("Number of subtractions=", sub_count);
    }

    //assignment operator
    //allows operator chaining
    Double& operator=(double x)
    {
        val = x;
        return *this;
    }

    Double& operator=(const Double& x)
    {
        if (!(this == &x))
            val = x.val;
        return *this;
    }

    Double& operator*=(const Double& x)
    {
        mul_count++;
        val *= x.val;
        return *this;
    }

    Double& operator/=(const Double& x)
    {
        div_count++;
        val /= x.val;
        return *this;
    }

    Double& operator+=(const Double& x)
    {
        add_count++;
        val += x.val;
        return *this;
    }

    Double& operator-=(const Double& x)
    {
        sub_count++;
        val -= x.val;
        return *this;
    }

    const Double operator*(const Double& x) const
    {
        Double result = *this;
        result *= x;
        return result;
    }

    const Double operator/(const Double& x) const
    {
        Double result = *this;
        result /= x;
        return result;
    }

    const Double operator+(const Double& x) const
    {
        Double result = *this;
        result += x;
        return result;
    }

    const Double operator-(const Double& x) const
    {
        Double result = *this;
        result -= x;
        return result;
    }

    const Double operator-() const
    {
        Double result(-val);
        return result;
    }

    //comparison operators
    bool operator==(const Double& x) const
    {
        return val == x.val;
    }

    bool operator!=(const Double& x) const
    {
        return !(*this == x);
    }

    bool operator<(const Double& x) const
    {
        return val < x.val;
    }

    bool operator<=(const Double& x) const
    {
        return val <= x.val;
    }

    bool operator>(const Double& x) const
    {
        return val > x.val;
    }

    bool operator>=(const Double& x) const
    {
        return val >= x.val;
    }

    friend std::ostream& operator<<(std::ostream& os, const Double& x);
};

unsigned int Double::mul_count = 0;
unsigned int Double::div_count = 0;
unsigned int Double::add_count = 0;
unsigned int Double::sub_count = 0;

std::ostream& operator<<(std::ostream& os, const Double& x)
{
    os << x.val;
    return os;
}

inline Double operator+(const double& x, const Double& y)
{
    return y + x;
}

inline Double operator-(const double& x, const Double& y)
{
    return -(y - x);
}

inline Double operator*(const double& x, const Double& y)
{
    return y * x;
}

inline Double operator/(const double& x, const Double& y)
{
    y.div_count++;
    return Double(x / y.val);
}

inline Double sqrt(const Double& x)
{
    return Double(std::sqrt(x.val));
}

inline Double pow(const Double& x, const double& y)
{
    return Double(pow(x.val, y));
}

inline Double sin(const Double& x)
{
    return Double(sin(x.val));
}

inline Double cos(const Double& x)
{
    return Double(cos(x.val));
}

inline Double tan(const Double& x)
{
    return Double(tan(x.val));
}

inline const Double& conj(const Double& x)
{
    return x;
}

inline const Double& real(const Double& x)
{
    return x;
}

inline Double imag(const Double& x)
{
    return 0;
}

inline Double abs(const Double& x)
{
    return Double(std::abs(x.val));
}

inline Double abs2(const Double& x)
{
    return x * x;
}

inline Double fabs(const Double& x)
{
    return Double(std::fabs(x.val));
}

inline Double copysign(const Double& x, const Double& y)
{
    return Double(std::copysign(x.val, y.val));
}
} // namespace ZIRAN

#endif
