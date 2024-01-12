#include "stdio.h"
#include "math.h"

const int ROWS=3;               // matrix rows
const int COLS=3;               // matrix cols

float MatA[ROWS][COLS];         // Matrix A, original
float Mata[ROWS][COLS];         // Matrix a, working version of A
float ColB[ROWS];               // Column vector B, original
float Colb[ROWS];               // Column vector b, working version of B
float ColX[ROWS];               // Column vector X

/*--------------------------------------------------------------
| Prepare matrices and column vectors
--------------------------------------------------------------*/
void FillAXB()
{
    Mata[0][0]=2; Mata[0][1]=1; Mata[0][2]=1;
    Mata[1][0]=1; Mata[1][1]=2; Mata[1][2]=1;
    Mata[2][0]=1; Mata[2][1]=1; Mata[2][2]=2;
    int val_b=1;
    for(int r=0; r<ROWS; r++)
    {
        MatA[r][0]=Mata[r][0]; MatA[r][1]=Mata[r][1]; MatA[r][2]=Mata[r][2];
        ColB[r]=Colb[r]=val_b++;
        ColX[r]=0;
    }
}

/*--------------------------------------------------------------
| Print Ax=b in matrix with column vectors format
| - print_x, false print symbol, else print X values
--------------------------------------------------------------*/
void PrintAXB(bool print_x=false)
{
    for(int r=0; r<ROWS; r++)
    {
        printf("  |");
        for(int c=0; c<COLS; c++)
        {
            printf(" %5.2f ",MatA[r][c]);
        }
        if(print_x==false)
        {
            if(r==1)
                printf(" |   x%d  | =",r);
            else
                printf(" |   x%d  |  ",r); 
        }
        else
                {
            if(r==1)
                printf(" | %5.2f | =",ColX[r]);
            else
                printf(" | %5.2f |  ",ColX[r]); 
        }
        printf(" |%5.2f | ",ColB[r]);
        printf("\n");
    }
    printf("\n");
}

/*--------------------------------------------------------------
| Gaussian Elimination solver for Ax=b
| - Works with Mata, ColX, and Colb
--------------------------------------------------------------*/
bool GaussianElimination()
{
    for(int c = 0; c<COLS; c++)                 // For each col
    {
        float dv = Mata[c][c];
        for(int cp=0; cp<COLS; cp++)        // For each col and pivot
            Mata[c][cp] /= dv;
        Colb[c]/=dv;
        //PrintAB();
        for(int r=c;r<ROWS;r++)                 // For each row in that col
        {
            if(r!=c)
            {
                float dv = Mata[r][c];
                for(int cp=c; cp<COLS; cp++)        // For each col and pivot
                    Mata[r][cp] /= dv;
                Colb[r]/=dv;
                //PrintAB();
                for(int cp=c; cp<COLS; cp++)        // For each col and pivot
                    Mata[r][cp] -= Mata[c][cp];
                Colb[r]-=Colb[c];
                //PrintAB();
            }
        }
    }
    for(int r=ROWS-1;r>=0;r--)          // r: 2, 1, 0
    {
        ColX[r]=Colb[r];
        for(int c=r+1; c<COLS; c++)      // c: 2, 2 1, 2 1 0
        {
            ColX[r]-=Mata[r][c]*ColX[c];
        }
    }
     return true;
}

/*--------------------------------------------------------------
| Solve Ax=b system with Gaussian Elimination
--------------------------------------------------------------*/
int main(int argc, char *argv[])
{
    FillAXB();                          // fill Ax=b
    printf("Original system:\n");
    PrintAXB(false);                    // print original
    GaussianElimination();              // solve
    printf("Solved system:\n");
    PrintAXB(true);                     // print solution
}
