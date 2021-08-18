
#include <iostream>

#include "EquationCalculator.h"
#include "Math.h"

void PrintComplex(const std::vector<std::pair<double, double>> &coms) {
    auto complex_to_str = [](const std::pair<double, double> &comlex) {
        double a = comlex.first;
        double b = comlex.second;
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
    };

    int index = 1;
    for (auto iter = coms.begin(); iter != coms.end(); ++iter) {
        std::cout << "x" << index << ": " << complex_to_str(*iter) << "\n";
        ++index;
    }
}

int main(int args, char **argv) {
    if (args == 3) {
        auto re = Func(atof(argv[1]), atof(argv[2]));
        PrintComplex(re);
        std::cout << "value(x=1):" << FuncCalc(1, atof(argv[1]), atof(argv[2])
                ) << "\n";
    } else if (args == 4) {
        auto re = Func2(atof(argv[1]), atof(argv[2]), atof(argv[3]));
        PrintComplex(re);
        std::cout << "value(x=1):" << Func2Calc(1, atof(argv[1]), atof(argv[2]),
                atof(argv[3])) << "\n";
    } else if (args == 5) {
        auto re = Func3(atof(argv[1]), atof(argv[2]), atof(argv[3]),
                atof(argv[4]));

        PrintComplex(re);
        std::cout << "value(x=1):" << Func3Calc(1, atof(argv[1]), atof(argv[2]),
                atof(argv[3]), atof(argv[4])) << "\n";
    } else if (args == 6) {
        auto re = Func4(atof(argv[1]), atof(argv[2]), atof(argv[3]),
                atof(argv[4]), atof(argv[5]));

        PrintComplex(re);
        std::cout << "value(x=1):" << Func4Calc(1, atof(argv[1]), atof(argv[2]),
                atof(argv[3]), atof(argv[4]), atof(argv[5])) << "\n";
    } else if (args == 7) {
        auto re = Func5(atof(argv[1]), atof(argv[2]), atof(argv[3]),
                atof(argv[4]), atof(argv[5]), atof(argv[6]));

        PrintComplex(re);
        std::cout << "value(x=1):" << Func5Calc(1, atof(argv[1]), atof(argv[2]),
                atof(argv[3]), atof(argv[4]), atof(argv[5]), atof(argv[6])
                ) << "\n";
    }
}

