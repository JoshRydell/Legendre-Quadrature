#include <stdio.h>
#include "mathLib.h"
#include "Lgdre.h"
int main()
{
    printf("\n");
    int N = 50;
    
    //printRootsAndWeights(N);
    double t = integrateGL(&log,1,10,N);
    printf("\n%lf", t);
    return 0;
}


  