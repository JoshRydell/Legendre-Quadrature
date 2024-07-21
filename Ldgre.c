#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
double nLegendreRecurse(double x, int N, int n, double Pnm1, double Pnm2)
{
    double Pn = 1.0/n*((2*n-1)*x*Pnm1-(n-1)*Pnm2);
    if(n == N) return Pn;
    return nLegendreRecurse(x, N, n+1, Pn, Pnm1);
    
}

double nLegendre(double x, int N)
{
    if(N < 0) return -1;
    if(N == 0) return 1;
    if(N == 1) return x;

    return nLegendreRecurse(x, N, 2, x, 1);
}

double nDerivLegendre(double x, int N)
{
    if(N < 0) return -1;
    if(N == 0) return 0;
    if(N == 1) return 1;
    if(x == 1) return (double)(N*(N+1))/2.0;
    if(x == -1) return (N % 2 ==0 ) ? -(double)(N*(N+1))/2.0 : (double)(N*(N+1))/2.0;

    return ((double)N)/(x*x-1)*(x*nLegendre(x,N)-nLegendre(x,N-1));
}

double* findRoots(int N)
{
    double *roots = calloc(N,sizeof(double));
    double *adjustment = calloc(N,sizeof(double));
    for(int i = 0; i < N; i++) *(roots+i) = (i % 2 == 0) ? -1.0 + i*2.0/N - (double)i/(N*N) : -1.0 + i*2.0/N + (double)i/(N*N);
    for(int i =0; i < 30; i++)
    {
        for(int j = 0; j < N; j++)
        {
            
            double ratio = nLegendre(*(roots+j), N)/nDerivLegendre(*(roots+j),N); 
            double recipSum = 0.0;
            for(int k = 0; k < N; k++)
            {
                if(k != j)
                {
                    recipSum += 1.0/(*(roots+j)-*(roots+k));
                }
            }
            *(adjustment+j) = ratio/(1-ratio*recipSum);
        }
        
        bool check = true;
        double acrcy = 0.0001;
        for(int j = 0; j < N; j++)
        {
            *(roots+j) -= *(adjustment +j);
            if(check &&  (*(adjustment+j) > acrcy || *(adjustment+j) < -acrcy)) check = false;
        }
        
        if(check)
        {
            free(adjustment);
            return roots;
        }
    }

    free(adjustment);
    return roots;
}

double* getWeights(double *roots,int N)
{
    double * weights = calloc(N, sizeof(double));
    for(int i =0; i < N; i++)
    {
        *(weights+i) = 2.0/((1-(*(roots+i))*(*(roots+i)))*nDerivLegendre(*(roots+i),N)*nDerivLegendre(*(roots+i),N));
    }

    return weights;
}

void printRootsAndWeights(int N)
{
    double *roots = findRoots(N);
    double *weights = getWeights(roots,N);
    printf("ROOTS:\n");
    for(int i = 0; i < N; i++)
    {
        printf("%lf,", *(roots+i));
    }
    printf("\nWEIGHTS:\n");
    for(int i = 0; i < N; i++)
    {
        printf("%lf,", *(weights+i));
    }
    
    free(roots); 
    free(weights);

}

double integrateGL(double (*func)(double),double a, double b, int N)
{
    double * roots = findRoots(N);
    double * weights = getWeights(roots, N);

    double dx = ((double)(b-a)/2);
    double val = 0.0;
    double trans = +(double)(a+b)/2;
    for(int i = 0; i < N; i++)
    {
        double root = *(roots+i);
        double scaledValue = dx*root + trans;
        double w = *(weights+i);
        double funcVal = (*func)(scaledValue);
        val += funcVal*w;
    }

    free(roots);
    free(weights);
    return dx*val;
}
