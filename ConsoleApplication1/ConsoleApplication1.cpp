#include <clocale>
#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>

// произв.
double function(double x)
{
    return std::sin(x) * x * x;
}

// используя метод Гаусса аппроксимируем производную
double derivGauss(double (*func)(double), double x0, double step)
{
    return (func(x0 + step) - func(x0 - step)) / (2 * step);
}

// вычисление полинома
double polLeg(int n, double x)
{
    if (n == 0)
    {
        return 1.0;
    }
    else if (n == 1)
    {
        return x;
    }
    else
    {
        return ((2.0 * n - 1.0) * x * polLeg(n - 1, x) - (n - 1) * polLeg(n - 2, x)) / n;
    }
}

// вычисление производной полинома
double derivPol_leg(int n, double x)
{
    if (n == 0)
    {
        return 0.0;
    }
    else if (n == 1)
    {
        return 1.0;
    }
    else
    {
        return (n / (x * x - 1.0)) * (x * polLeg(n, x) - polLeg(n - 1, x));
    }
}

// Вычислени интеграла при помощи метода Гаусса-Лежандра
double gaussLeg(double (*f)(double), double a, double b, int n)
{
    double* x = new double[n];
    double* w = new double[n];

    int    m = (n + 1) / 2;
    double xm = 0.5 * (b + a);
    double xl = 0.5 * (b - a);

    for (int i = 1; i <= m; i++)
    {
        double z = cos(M_PI * (i - 0.25) / (n + 0.5));
        double p0 = 1.0;
        double p1 = z;

        for (int j = 2; j <= n; j++)
        {
            double pj = ((2.0 * j - 1.0) * z * p1 - (j - 1.0) * p0) / j;
            p0 = p1;
            p1 = pj;
        }

        x[i - 1] = xm - xl * z;
        x[n - i] = xm + xl * z;
        w[i - 1] = 2.0 * xl / ((1.0 - z * z) * pow(derivPol_leg(n, z), 2));
        w[n - i] = w[i - 1];
    }

    double integral = 0.0;
    for (int i = 0; i < n; i++)
    {
        integral += w[i] * (*f)(x[i]);
    }

    delete[] x;
    delete[] w;

    return integral;
}

int main(void)
{
    std::setlocale(LC_ALL, "");

    std::string function_str = "x^2 * sin(x)";

    double step, x;

    std::cout << "Введите Число для нахождения производной от x: ";
    std::cin >> x;

    std::cout << "Введите размер шага (чем меньше, тем точнее): ";
    std::cin >> step;

    std::cout << std::endl
        << "Функция - " << function_str << "; Производная x: " << derivGauss(function, x, step) << std::endl
        << std::endl;

    double a, b;
    int    n;

    std::cout << "Введите нижнее значение производной (a): ";
    std::cin >> a;

    std::cout << "Введите верхнее значение производной (b): ";
    std::cin >> b;

    std::cout << "Введите желаемое кол-во итераций вычисления интеграла (чем больше, тем точнее)): ";
    std::cin >> n;

    std::cout << std::endl
        << "Функция - " << function_str << "; Интеграл от a до b: " << gaussLeg(function, a, b, n)
        << std::endl;

    return EXIT_SUCCESS;
}