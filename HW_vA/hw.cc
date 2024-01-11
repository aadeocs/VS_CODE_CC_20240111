#include "stdio.h"
#include <atomic>
#include <cmath>

//std::atomic<int> abc(0);

int main(int argc, char **argv)
{
    printf(" +-----------------------+\n");
    printf(" | Hello world, Maggie!! |\n");

    int a(3);
    int b(4);
    int c = sqrt(a*a+b*b); 
    printf(" | sqrt( %d^2 + %d^2 ) = %d |\n",a,b,c);
    printf(" +-----------------------+\n");
    printf("\n");

}

