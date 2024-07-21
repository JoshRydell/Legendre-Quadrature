#include <stdio.h>
#include <stdint.h>

const double PI = 3.141592654;

double power(double x, int index);

double exp(double x);

double sin(double x);

int factorial(int n);

double absolute(double x);

double integrate(double (*func)(double),double a, double b);

double log(double x);

double func(double t);

struct Matrix;