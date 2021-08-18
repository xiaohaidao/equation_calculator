
#include "EquationCalculator.h"

#include <algorithm>
#include <cmath>
#include <string>

#include "Math.h"

class Complex {
public:
    Complex() : r_(0), i_(0) {}
    Complex(double real, double imag) : r_(real), i_(imag) {}
    Complex(double a) : r_(a), i_(0) {}

    double Abs2() const {
        return r_ * r_ + i_ * i_;
    }

    double Abs() const {
        return sqrt(Abs2());
    }

    double Arg() const {
        return atan2(i_, r_);
    }

    Complex Conj() const {
        return Complex(r_, -i_);
    }

    Complex Sqrt() const {
        return Complex(cos(this->Arg() / 2), sin(this->Arg() / 2)) *
            sqrt(this->Abs());
    }

    Complex Cbrt() const {
        double t = (IsEqual(this->i_, 0) && IsLesser(this->r_, 0)) ?
            M_PI : this->Arg() / 3;

        return Complex(cos(t), sin(t)) * pow(this->Abs(), 1.0 / 3);
    }

    Complex operator+(const Complex &other) const {
        return Complex(this->r_ + other.r_, this->i_ + other.i_);
    }

    Complex operator-(const Complex &other) const {
        return Complex(this->r_ - other.r_, this->i_ - other.i_);
    }

    Complex operator-() const {
        return Complex(-this->r_, -this->i_);
    }

    Complex operator*(const Complex &other) const {
        return Complex(this->r_ * other.r_ - this->i_ * other.i_,
                this->r_ * other.i_ + this->i_ * other.r_);
    }

    Complex operator/(const Complex &other) const {
        return *this * other.Conj() * (1 / other.Abs2());
    }

    std::string ToString() const {
        double a = r_;
        double b = i_;
        if (IsEqual(b, 0)) {
            return std::to_string(a);
        }
        std::string re;
        if (!IsEqual(a, 0)) {
            re = std::to_string(a);
        }
        if (IsEqual(b, 1)) {
            re += "+i";
        } else if (IsEqual(b, -1)) {
            re += "-i";
        } else {
            if (IsGreater(b, 0)) {
                re += "+";
            }
            re += (std::to_string(b) + "i");
        }
        return re;
    }

    double Real() const {
        return r_;
    }

    double Imag() const {
        return i_;
    }

private:
    double r_;
    double i_;

}; // class Complex

bool IsEqual(const Complex &left, const Complex &right) {
    Complex t = left - right;
    return IsEqual(t.Real(), 0) && IsEqual(t.Imag(), 0);
}

// ax+b=0
std::vector<Complex> Fun(const Complex &a, const Complex &b) {
    std::vector<Complex> re;
    re.push_back( -b / a);
    return re;
}

// ax^2+bx+c=0
std::vector<Complex> Fun2(const Complex &a, const Complex &b, const Complex &c) {
    Complex t = (b * b - a * c * 4).Sqrt();
    std::vector<Complex> re;
    re.push_back((-t -b) / (a * 2));
    re.push_back((t -b) / (a * 2));
    return re;
}

double Fun3Root(double v, double a, double b, double c) {
    a /= v;
    b /= v;
    c /= v;
    double p = b - a * a / 3;
    double q = (2 * a * a - 9 * b) * a / 27 + c;
    Complex t = Complex(q * q / 4 + p * p * p / 27).Sqrt();
    return ((t - q / 2).Cbrt() + (-t - q / 2).Cbrt()).Real() - a / 3;
}

// vx^3+ax^2+bx+c=0
std::vector<Complex> Fun3(double v, double a, double b, double c) {
    a /= v;
    b /= v;
    c /= v;

    double x = Fun3Root(1, a, b, c);
    std::vector<Complex> re;
    re.push_back(x);
    std::vector<Complex> re2 = Fun2(1, x + a, x * (x + a) +b);

    re.insert(re.end(), re2.begin(), re2.end());
    return re;
}

// vx^4+ax^3+bx^2+cx+d=0
std::vector<Complex> Fun4Roots(double v, double a, double b, double c, double d) {
    a /= v;
    b /= v;
    c /= v;
    d /= v;
    double p = b - 3 * a * a / 8;
    double q = a * (a * a - 4 * b) / 8 + c;
    double r = d - a * (a * (3 * a * a - 16 * b) + 64 * c) / 256;
    Complex u, v1, v2, t;
    if (IsEqual(q, 0)) {
        u = 0;
        v1 = (Complex(p * p - 4 * r).Sqrt() + p) / 2;
        v2 = Complex(p) - v1;
    } else {
        u = Complex(Fun3Root(1, 2 * p, p * p - 4 * r, -q * q)).Sqrt();
        v1 = (Complex(p) + u * u - Complex(q) / u) / 2;
        v2 = v1 + Complex(q) / u;
    }

    std::vector<Complex> list;
    t = (u * u - Complex(4) * v1).Sqrt();
    list.push_back((-u + t) / 2 - a / 4);
    list.push_back((-u - t) / 2 - a / 4);
    t = (u * u - Complex(4) * v2).Sqrt();
    list.push_back((u + t) / 2 - a / 4);
    list.push_back((u - t) / 2 - a / 4);
    return list;
}

// x^5+ax^4+bx^3+cx^2+dx+e
double Func5Calc(double a, double b, double c, double d, double e, double x) {
    return x * (x * (x * (x * (x + a) + b) + c) + d) + e;
}

// 5x^4+4ax^3+3bx^2+2cx+d
double DFun5(double a, double b, double c, double d, double x) {
    return x * (x * (x * (5 * x + 4 * a) + 3 * b) + 2 * c) + d;
}

double Start(double a, double b, double c, double d, double e,
        const std::vector<double> &roots) {

    std::vector<double> list = roots;
    std::sort(list.begin(), list.end());
    return Func5Calc(a, b, c, d, e, list[0]) > 0 ? list[0] - 0.5 :
        Func5Calc(a, b, c, d, e, list[list.size() - 1]) < 0 ?
        list[list.size() - 1] + 0.5 : (list[1] + list[2]) / 2;
}

// x^5+ax^4+bx^3+cx^2+dx+e=0
std::vector<Complex> Func5(double a, double b, double c, double d, double e) {
    std::vector<Complex> list = Fun4Roots(5, 4 * a, 3 * b, 2 * c, d);
    std::vector<double> roots;
    for (auto const &z: list) {
        if (IsZero(z.Imag())) {
            roots.push_back(z.Real());
        }
    }
    double t = 0;
    double x = 0;
    int i = 0;

    for (auto const &root: roots) {
        if (IsZero(Func5Calc(a, b, c, d, e, root))) {
            x = root;
            i = -1;
            break;
        }
    }

    if (i == 0) {
        if (roots.size() > 0) {
            x - Start(a, b, c, d, e, roots);
        }
        if (!IsZero(Func5Calc(a, b, c, d, e, x))) {
            do {
                t = x;
                x -= Func5Calc(a, b, c, d, e, x) / DFun5(x, b, c, d, x);
                ++i;
            } while (i < 6000 && !IsEqual(x, t));
        }
    }
    if (IsZero(Func5Calc(a, b, c, d, e, t))) {
        x = t;
    }
    std::vector<Complex> re;
    re.push_back(x);
    std::vector<Complex> fun4 = Fun4Roots(1, x + a, x * (x + a) + b,
            x * (x * (x + a) + b) + c,
            x * (x * (x * (x + a) + b) + c) + d);

    re.insert(re.end(), fun4.begin(), fun4.end());
    return re;
}

std::vector<Complex> Fun5(double v, double a, double b, double c, double d,
        double e) {

    return Func5(a / v, b / v, c / v, d / v, e / v);
}

std::pair<double, double> ComplexToDouble(const Complex &number) {
    return std::make_pair(number.Real(), number.Imag());
}

std::vector<std::pair<double, double>> ComplexToDouble(
        const std::vector<Complex> &number) {

    std::vector<std::pair<double, double>> re;
    for (auto iter = number.begin(); iter != number.end(); ++iter) {
        re.push_back(ComplexToDouble(*iter));
    }
    return re;
}

std::vector<std::pair<double, double> > Func5(double a, double b, double c,
        double d, double e, double f) {

    return ComplexToDouble(Fun5(a, b, c, d, e, f));
}

std::vector<std::pair<double, double> > Func4(double a, double b, double c,
        double d, double e) {

    return ComplexToDouble(Fun4Roots(a, b, c, d, e));
}

std::vector<std::pair<double, double> > Func3(double a, double b, double c,
        double d) {

    return ComplexToDouble(Fun3(a, b, c, d));
}

std::vector<std::pair<double, double> > Func2(double a, double b, double c) {
    return ComplexToDouble(Fun2(a, b, c));
}

std::vector<std::pair<double, double> > Func(double a, double b) {
    return ComplexToDouble(Fun(a, b));
}

double Func5Calc(double x, double a, double b, double c, double d, double e,
        double f) {

    return (x * (x * (x * (x * (x * a + b)) + c) + d) + e) + f;
}

double Func4Calc(double x, double a, double b, double c, double d, double e) {
    return (x * (x * (x * (x * a + b)) + c) + d) + e;
}

double Func3Calc(double x, double a, double b, double c, double d) {
    return (x * (x * (x * a + b)) + c) + d;
}

double Func2Calc(double x, double a, double b, double c) {
    return (x * (x * a + b)) + c;
}

double FuncCalc(double x, double a, double b) {
    return x * a + b;
}

