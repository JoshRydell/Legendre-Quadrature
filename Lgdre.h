#include <stdlib.h>
double nLegendre(double x, int N);

double nDerivLegendre(double x, int N);

double* findRoots(int N);

void printRootsAndWeights(int N);

double integrateGL(double (*func)(double),double a, double b, int N);