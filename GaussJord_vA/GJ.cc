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
| - Should detect singular, dv is very small
| - Add partial pivoting
--------------------------------------------------------------*/
bool GaussianElimination()
{
    // elimination
    for(int c = 0; c<COLS; c++)                     // For all cols
    {
        float dv = Mata[c][c];
        for(int cp=0; cp<COLS; cp++)                // Normalize the cth row
            Mata[c][cp] /= dv;
        Colb[c]/=dv;

        for(int r=c;r<ROWS;r++)                     // For all rows in that col
        {
            if(r!=c)                                // only process non-cth row
            {
                float dv = Mata[r][c];
                for(int cp=0; cp<COLS; cp++)        // for all cols (could start at c)
                    Mata[r][cp] /= dv;              // normalize row
                Colb[r]/=dv;                        // normalize b
                for(int cp=0; cp<COLS; cp++)        // for all cols (could start at c)
                    Mata[r][cp] -= Mata[c][cp];     // subtract ref row from current row
                Colb[r]-=Colb[c];                   // subtract ref b from current b
            }
        }
    }
    // back substitute
    for(int r=ROWS-1;r>=0;r--)          // r: 2, 1, 0, then c: 2, 2 1, 2 1 0
    {
        ColX[r]=Colb[r];                    // start with b
        for(int c=r+1; c<COLS; c++)         // c: 2, 2 1, 2 1 0
        {
            ColX[r]-=Mata[r][c]*ColX[c];    // subtract precompute terms
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
