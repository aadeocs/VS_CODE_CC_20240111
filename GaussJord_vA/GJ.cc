#include "stdio.h"
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
    
    Colb[0]=ColB[0]=6; 
    Colb[1]=ColB[1]=6; 
    Colb[2]=ColB[2]=6;

    for(int r=0; r<ROWS; r++)
    {
        MatA[r][0]=Mata[r][0]; MatA[r][1]=Mata[r][1]; MatA[r][2]=Mata[r][2];
        ColX[r]=0;
    }
}

/*--------------------------------------------------------------
| Print Ax=b in matrix with column vectors format
| - print_x, false print symbol, else print X values
--------------------------------------------------------------*/
void PrintAXB(bool print_x=false)
{
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
                printf(" | %6.3f | =",ColX[r]);
            else
                printf(" | %6.3f |  ",ColX[r]); 
        }
        printf(" |%5.2f | ",ColB[r]);
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
void Printab(char *msg, int row, bool print_x=false)
{
    printf("%s [%d] = %d\n",msg,row,Pvt[row]);
    for(int r=0; r<ROWS; r++)
    {
        printf("  |");
        for(int c=0; c<COLS; c++)
        {
            printf(" %5.2f ",Mata[r][c]);
        }
        printf(" :%5.2f | \n",Colb[r]);
    }
    printf("\n"); 
}

/*--------------------------------------------------------------
| Gaussian Elimination solver for Ax=b
| - Works with Mata, ColX, and Colb
| - Should detect singular, dv is very small
| - Add partial pivoting
--------------------------------------------------------------*/
bool GaussJordan(float mat_a[ROWS][COLS], float col_x[COLS],float col_b[COLS])
{
    Pvt[0]=0;
    Pvt[1]=1;
    Pvt[2]=2;

    //Pvt[0]=2;
    //Pvt[1]=1;
    //Pvt[2]=0;

    // elimination
    for(int c = 0; c<COLS; c++)                     // For all cols
    {
        //for(int r=0; r<ROWS; r++)
        //{
        //    
        //}
        float piv = mat_a[ Pvt[c]][c];
        for(int cp=0; cp<COLS; cp++)                // Normalize the cth row
            mat_a[ Pvt[c]][cp] /= piv;
        col_b[ Pvt[c]]/=piv;
        
        Printab("Apply Pivot",c);
        for(int r=0;r<ROWS;r++)                     // For all rows in that col
        {
            if(r!=c)                                // only process non-cth row
            {
                float mpy = mat_a[ Pvt[r]][c];
                for(int cp=0; cp<COLS; cp++)        // for all cols (could start at c)
                    mat_a[ Pvt[r]][cp] -= mat_a[ Pvt[c]][cp] *mpy;     // subtract ref row from current row
                col_b[ Pvt[r]]-=col_b[ Pvt[c]]*mpy;                   // subtract ref b from current b
                Printab("Eliminate Entry",r);
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
    printf("Original system:\n");
    PrintAXB(false);                    // print original
    GaussJordan(Mata,Colx,Colb);              // solve
    for(int r=0; r<ROWS; r++)
        ColX[r]=Colx[Pvt[r]];
    printf("Solved system:\n");
    PrintAXB(true);                     // print solution
}
