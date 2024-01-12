#include "stdio.h"
#include "math.h"

const int ROWS=3;
const int COLS=3;
float MatA[ROWS][COLS];
float ColB[ROWS];
float ColX[ROWS];

void FillABX()
{
    int val_a=2;
    int val_b=1;
    for(int r=0; r<ROWS; r++)
    {
        for(int c=0; c<COLS; c++)
        {
            MatA[r][c]=val_a++;
        }
        ColB[r]=val_b++;
        ColX[r]=0;
    }
}

void PrintAB()
{
    for(int r=0; r<ROWS; r++)
    {
        printf("  |");
        for(int c=0; c<COLS; c++)
        {
            printf(" %4.1f ",MatA[r][c]);
        }
        printf(" %4.1f",ColB[r]);
        printf(" |\n");
    }
    printf("\n");
}

bool GaussJordan()
{
    for(int c = 0; c<COLS; c++)                 // For each col
    {
        for(int r=0;r<ROWS;r++)                 // For each row in that col
        {
            if(r==c)
            {
                int dv = MatA[r][r];
                for(int cp=0; cp<COLS; cp++)        // For each col and pivot
                    MatA[r][cp] /= dv;
                ColB[r]/=dv;
            }
            else
            {
                int dv = MatA[r][r];
                for(int cp=0; cp<COLS; cp++)        // For each col and pivot
                    MatA[r][cp] /= dv;
                ColB[r]/=dv;
            }
        }
        printf(" Col=%d:\n",c);
        PrintAB();
        return true;
    }
}

int main(int argc, char *argv[])
{
    FillABX();
    printf("MatA:\n");
    PrintAB();
    GaussJordan();
}