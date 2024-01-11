#include "stdio.h"
#include<iostream>

//using namespace std;

void LUdecomposition(float A[10][10], float L[10][10], float U[10][10], int n) 
{
   for (int i = 0; i < n; i++)
   {
      for (int j = 0; j < n; j++)
      {
         if (j < i)
            L[j][i] = 0;
         else
         {
            L[j][i] = A[j][i];
            for (int k = 0; k < i; k++)
            {
               L[j][i] = L[j][i] - L[j][k] * U[k][i];
            }
         }
      }
      for (int j = 0; j < n; j++)
      {
         if (j < i)
            U[i][j] = 0;
         else if (j == i)
            U[i][j] = 1;
         else
         {
            U[i][j] = A[i][j] / L[i][i];
            for (int k = 0; k < i; k++) 
            {
               U[i][j] = U[i][j] - ((L[i][k] * U[k][j]) / L[i][i]);
            }
         }
      }
   }
}

int main()
{
   float matA[10][10], matL[10][10], matU[10][10];
   int n = 0;

   printf("Enter size of square matrix: ");
   std::cin >> n;

   printf("Enter matrix values: ");
   for (int i = 0; i < n; i++)
   {
      for (int j = 0; j < n; j++)
      {
         std::cin >> matA[i][j];
      }
   }

   printf("A Matrix:\n");
   for (int i = 0; i < n; i++)
   {
      for (int j = 0; j < n; j++)
      {
         printf(" %3.0f ",matA[i][j]);
         //cout<<a[i][j]<<" ";
      }
      printf("\n");
      //cout << endl;
   }

   LUdecomposition(matA, matL, matU, n);

   printf("L Decomposition:\n");
   for (int i = 0; i < n; i++)
   {
      for (int j = 0; j < n; j++)
      {
         printf(" %3.0f ",matL[i][j]);
      }
      printf("\n");
   }

   printf("U Decomposition:\n");
   for (int i = 0; i < n; i++)
   {
      for (int j = 0; j < n; j++)
      {
         printf(" %3.0f ",matU[i][j]);
      }
      printf("\n");
   }
   return 0;
}
