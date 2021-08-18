
#ifndef EQUATIONCALCULATOR_H
#define EQUATIONCALCULATOR_H

#include <vector>

/// reference: https://github.com/CappuccinoZ/Equation-solver
/// @return return the complex number, the value type is `<Real, Imaginary>[]`
/// ax^5 + bx^4 + cx^3 + dx^2 + ex + f = 0
std::vector<std::pair<double, double> > Func5(double a, double b, double c, double d, double e, double f);
/// ax^4 + bx^3 + cx^2 + dx + e = 0
std::vector<std::pair<double, double> > Func4(double a, double b, double c, double d, double e);
/// ax^3 + bx^2 + cx + d = 0
std::vector<std::pair<double, double> > Func3(double a, double b, double c, double d);
/// ax^2 + bx + c = 0
std::vector<std::pair<double, double> > Func2(double a, double b, double c);
/// ax + b = 0
std::vector<std::pair<double, double> > Func(double a, double b);


/// @return get the calculated value
/// ax^5 + bx^4 + cx^3 + dx^2 + ex + f
double Func5Calc(double x, double a, double b, double c, double d, double e, double f);
/// ax^4 + bx^3 + cx^2 + dx + e
double Func4Calc(double x, double a, double b, double c, double d, double e);
/// ax^3 + bx^2 + cx + d
double Func3Calc(double x, double a, double b, double c, double d);
/// ax^2 + bx + c
double Func2Calc(double x, double a, double b, double c);
/// ax + b
double FuncCalc(double x, double a, double b);

#endif // EQUATIONCALCULATOR_H
