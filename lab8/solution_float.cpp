#include <iostream>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <stdlib.h>
#include <thread>

float calculate_max_error(float *x, float *result, int n) {
    float max_error=0;
    for(int i=0;i<n;i++)
    {
       
        if(max_error<std::fabs(x[i]-result[i]))
        {
            max_error=std::fabs(x[i]-result[i]);
            
        }
    }
    return max_error;
}

float calculate_mean_error(float *x, float *result, int n) {
    float mean_error=0;
    for(int i=0;i<n;i++)
    {
        mean_error += std::fabs(x[i] - result[i]) / n;
        
    }
    return mean_error;
}

void print_matrix(float **A, int n) {
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            std::cout<<A[i][j]<<" ";
        }
        std::cout<<std::endl;
    }
}

void print_vector(float *A, int n) {
    for(int i=0;i<n;i++)
    {
        std::cout<<A[i]<<" ";
    }
    std::cout<<std::endl;
}

void create_matrix1(float **A, int n) {
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            if(i==0)
            {
                 A[i][j]=1;
            }
            else
            {
                A[i][j]=(float)1/((float)i+1+(float)j);
            }
        }
    }
}


void create_matrix2(float **A, int n) {
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            if(j>=i)
            {
                 A[i][j]=2*((float)i+1)/((float)j+1);
            }
           
        }
    }

    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            if(j<i)
            {
                 A[i][j]=A[j][i];
            }
           
        }
    }
}

void create_matrix3(float **A, int n) {
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            if(i==j)
            {
                A[i][j]=-2*(i+1)-4;
            }
            else if(j==i+1)
            {
                A[i][j]=i+1;
            }
            else if(j==i-1)
            {
                A[i][j]=2/(i+1);
            }
            else
            {
                A[i][j]=0;
            }
        }
    }
}

void create_three_diagonal_matrix(float **A, int n, float *a, float *b, float *c) {
    for(int i=0;i<n;i++)
    {
        if(i-1>=0)
        {
            a[i]=A[i][i-1];
        }
        b[i]=A[i][i];   
        if(i+1<n)
        {
            c[i]=A[i][i+1];
        }
        
    }
    
}

void float_thomas(float *a, float *b, float *c, float *d, float *result, int n) {
    for(int i=1;i<n;i++)
    {
       float m=a[i]/b[i-1];
       b[i]-=m*c[i-1];
       d[i]-=m*d[i-1];
    }
    
    result[n-1]=d[n-1]/b[n-1];
    for(int i=n-2;i>=0;i--)
    {

         result[i] = (d[i] - c[i] * result[i + 1]) / b[i];
    }
}




void float_gauss_elimination(float **A, float *B, float *result, int n) {
    for (int i = 0; i < n - 1; i++) {
        // Partial pivoting
        int maxRow = i;
        for (int k = i + 1; k < n; k++) {
            if (fabs(A[k][i]) > fabs(A[maxRow][i])) {
                maxRow = k;
            }
        }

        // Swap rows in A
        if (maxRow != i) {
            float *tempRow = A[i];
            A[i] = A[maxRow];
            A[maxRow] = tempRow;

            // Swap corresponding B values
            float tempB = B[i];
            B[i] = B[maxRow];
            B[maxRow] = tempB;
        }

        // Elimination
        for (int j = i + 1; j < n; j++) {
            if (fabs(A[i][i]) < 1e-12) {
                printf("Error: Division by zero or nearly zero pivot.\n");
                return;
            }

            float multiplier = A[j][i] / A[i][i];
            for (int k = i; k < n; k++) {
                A[j][k] -= multiplier * A[i][k];
            }
            B[j] -= multiplier * B[i];
        }
    }

    // Back substitution
    result[n - 1] = B[n - 1] / A[n - 1][n - 1];
    for (int i = n - 2; i >= 0; i--) {
        float sum = 0;
        for (int j = i + 1; j < n; j++) {
            sum += A[i][j] * result[j];
        }
        result[i] = (B[i] - sum) / A[i][i];
    }
}




int main()
{
    const int n=15000;
    float **A=new float*[n];
    for(int i=0;i<n;i++)
    {
        A[i]=new float[n];
    }
    float Bt[n];
    float Bg[n];
    float result[n];
    float x[n];
    float a[n];
    float b[n];
    float c[n];


    create_matrix3(A,n);



    for(int i=0;i<n;i++)
    {
        if(i%2==0)
        {
            x[i]=1;
        }
        else
        {
            x[i]=-1;
        }
    }

    for(int i=0;i<n;i++)
    {
        Bt[i]=0;
        Bg[i]=0;
        for(int j=0;j<n;j++)
        {
            Bt[i]+=A[i][j]*x[j];
            Bg[i]+=A[i][j]*x[j];
        }
    }

    

    create_three_diagonal_matrix(A,n,a,b,c);
    

    auto start = std::chrono::high_resolution_clock::now();


  
    float_thomas(a,b,c,Bt,result,n);
    
    

    
     auto end = std::chrono::high_resolution_clock::now();

    
    std::chrono::duration<long double> elapsed = end - start;

    
    std::cout << "Czas wykonania: " << elapsed.count() << " sekund" << std::endl;
    

    std::cout << "Max error: " << calculate_max_error(x, result, n) << std::endl;
    
    
   
    auto start2 = std::chrono::high_resolution_clock::now();

    float_gauss_elimination(A,Bg,result,n);

    auto end2 = std::chrono::high_resolution_clock::now();

    std::chrono::duration<long double> elapsed2 = end2 - start2;

    std::cout << "Czas wykonania: " << elapsed2.count() << " sekund" << std::endl;


    
    std::cout << "Max error: " << calculate_max_error(x, result, n) << std::endl;
   

    for(int i=0;i<n;i++)
    {
        delete[] A[i];
    }
    delete[] A;

}

