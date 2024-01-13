/*###############################################################
|
+----------------------------------------------------------------
| Compile [Ctrl]+[Shift]+[B]
| Run:    [F5]
###############################################################*/

#include "stdio.h"
#include <stdarg.h>
#include "math.h"

const int ROWS=3;               // matrix rows
const int COLS=3;               // matrix cols

float MatA[ROWS][COLS];         // Matrix A, original
float Mata[ROWS][COLS];         // Matrix a, working version of A
float ColB[ROWS];               // Column vector B, original
float Colb[ROWS];               // Column vector b, working version of B
float ColX[ROWS];               // Column vector X
float Colx[ROWS];               // Column vector X, working version of X
int   Pvt[ROWS];

const int MAX_VARG_LEN = 1024;
char VargBuffer[MAX_VARG_LEN];

/*--------------------------------------------------------------
| Prepare matrices and column vectors
--------------------------------------------------------------*/
void FillAXB()
{
    float dv = 0.1;
    bool near_sing = false;

    Mata[0][0]=1; Mata[0][1]=2; Mata[0][2]=3;
    Mata[1][0]=2; Mata[1][1]=3; Mata[1][2]=1;
    Mata[2][0]=3; Mata[2][1]=1; Mata[2][2]=2;

    float x0=1;
    float x1=2;
    float x2=4;
    Colb[0]=ColB[0]=x0*Mata[0][0]+x1*Mata[0][1]+x2*Mata[0][2];
    Colb[1]=ColB[1]=x0*Mata[1][0]+x1*Mata[1][1]+x2*Mata[1][2];
    Colb[2]=ColB[2]=x0*Mata[2][0]+x1*Mata[2][1]+x2*Mata[2][2];
    for(int r=0; r<ROWS; r++)
    {
        MatA[r][0]=Mata[r][0]; MatA[r][1]=Mata[r][1]; MatA[r][2]=Mata[r][2];
        ColX[r]=0;
    }

    Pvt[0]=0;
    Pvt[1]=1;
    Pvt[2]=2;

    Pvt[0]=2;
    Pvt[1]=1;
    Pvt[2]=0;
}

/*--------------------------------------------------------------
| Print Ax=b in matrix with column vectors format
| - print_x, false print symbol, else print X values
--------------------------------------------------------------*/
void PrintAXB(bool print_x, const char *format, ...)
{
    va_list arg;
    int32_t done;
    va_start (arg, format);
        done = (int32_t)vsprintf(VargBuffer,format, arg);
    va_end (arg);

    printf("%s\n",VargBuffer);
    float b[COLS];
    if(print_x)
    {
        for(int r=0; r<ROWS; r++)
        {
            b[r]=0;
            for(int c=0; c<COLS; c++)
            {
                b[r]+=MatA[r][c]*ColX[c];
            }
        }
    }
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
                printf(" | %5.3f | =",ColX[r]);
            else
                printf(" | %5.3f |  ",ColX[r]); 
        }
        printf(" |%6.2f | ",ColB[r]);
        if(print_x)
            printf(" < %6.3f > ",b[r]);

        printf(" ( %d ) ",Pvt[r]);
        printf("\n");
    }
    printf("\n"); 
}

/*--------------------------------------------------------------
| //Print Ax=b in matrix with column vectors format
| //- print_x, false print symbol, else print X values
--------------------------------------------------------------*/
void Printab(int row, const char *format, ...)
{
    va_list arg;
    int32_t done;
    va_start (arg, format);
        done = (int32_t)vsprintf(VargBuffer,format, arg);
    va_end (arg);

    printf("%s\n",VargBuffer);
    for(int r=0; r<ROWS; r++)
    {
        printf("  |");
        for(int c=0; c<COLS; c++)
        {
            printf(" %5.2f ",Mata[r][c]);
        }
        printf(" :%6.2f |",Colb[r]);
        if(r==Pvt[row])
            printf("*");
        printf("\n");
    }
    printf("\n"); 
}

/*--------------------------------------------------------------
| Gauss-Jordan Elimination solver for Ax=b
| - Works with Mata, ColX, and Colb
| - Should detect singular, dv is very small
| - Add partial pivoting
--------------------------------------------------------------*/
bool GaussJordan(float mat_a[ROWS][COLS], float col_x[COLS],float col_b[COLS])
{

    // elimination
    for(int c = 0; c<COLS; c++)                     // For all cols
    {
        float piv = mat_a[ Pvt[c]][c];
        for(int cp=0; cp<COLS; cp++)                // Normalize the cth row
            mat_a[ Pvt[c]][cp] /= piv;
        col_b[ Pvt[c]]/=piv;
        
        Printab(c,"Apply Pivot, Pvt[%d]=Row[%d]",c,Pvt[c]);
        for(int r=0;r<ROWS;r++)                     // For all rows in that col
        {
            if(r!=c)                                // only process non-cth row
            {
                int zero_c=c;
                float mpy = mat_a[ Pvt[r]][c];
                for(int cp=0; cp<COLS; cp++)        // for all cols (could start at c)
                    mat_a[ Pvt[r]][cp] -= mat_a[ Pvt[c]][cp] *mpy;     // subtract ref row from current row
                col_b[ Pvt[r]]-=col_b[ Pvt[c]]*mpy;                   // subtract ref b from current b
                Printab(r,"Zero Entry, Mpv[%d][%d]=Mat[%d][%d]",r,zero_c,Pvt[r],zero_c);
            }
        }
    }
    // back substitute
    for(int r=ROWS-1;r>=0;r--)          // r: 2, 1, 0, then c: 2, 2 1, 2 1 0
    {
        col_x[ Pvt[r]]=col_b[ Pvt[r]];                    // start with b
        //for(int c=r+1; c<COLS; c++)         // c: 2, 2 1, 2 1 0
        //{
        //    col_x[ Pvt[r]]-=mat_a[ Pvt[r]][c]*col_x[ Pvt[c]];    // subtract precompute terms
        //}
    }
    return true;
}

/*--------------------------------------------------------------
| Solve Ax=b system with Gaussian Elimination
--------------------------------------------------------------*/
int main(int argc, char *argv[])
{
    FillAXB();                          // fill Ax=b
    //printf("Original system:\n");
    PrintAXB(false,"Orignal System:");                    // print original
    GaussJordan(Mata,Colx,Colb);              // solve
    for(int r=0; r<ROWS; r++)
        ColX[r]=Colx[Pvt[r]];
    //printf("Solved system:\n");
    PrintAXB(true,"Solved System:");                     // print solution
}
