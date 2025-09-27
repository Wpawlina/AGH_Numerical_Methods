#include <iostream>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <stdlib.h>
#include <thread>

double calculate_max_error(double *x, double *result, int n) {
    double max_error=0;
    for(int i=0;i<n;i++)
    {
       
        if(max_error<std::fabs(x[i]-result[i]))
        {
            max_error=std::fabs(x[i]-result[i]);
            
        }
    }
    return max_error;
}

double calculate_mean_error(double *x, double *result, int n) {
    double mean_error=0;
    for(int i=0;i<n;i++)
    {
        mean_error += std::fabs(x[i] - result[i]) / n;
        
    }
    return mean_error;
}

void print_matrix(double **A, int n) {
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            std::cout<<A[i][j]<<" ";
        }
        std::cout<<std::endl;
    }
}

void print_vector(double *A, int n) {
    for(int i=0;i<n;i++)
    {
        std::cout<<A[i]<<" ";
    }
    std::cout<<std::endl;
}

void create_matrix1(double **A, int n) {
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
                A[i][j]=(double)1/((double)i+1+(double)j);
            }
        }
    }
}




void create_matrix2(double **A, int n) {
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            if(j>=i)
            {
                 A[i][j]=2*((double)i+1)/((double)j+1);
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

void create_matrix3(double **A, int n) {
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

void create_three_diagonal_matrix(double **A, int n, double *a, double *b, double *c) {
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

void double_thomas(double *a, double *b, double *c, double *d, double *result, int n) {
    for(int i=1;i<n;i++)
    {
       double m=a[i]/b[i-1];
       b[i]-=m*c[i-1];
       d[i]-=m*d[i-1];
    }
    
    result[n-1]=d[n-1]/b[n-1];
    for(int i=n-2;i>=0;i--)
    {

         result[i] = (d[i] - c[i] * result[i + 1]) / b[i];
    }
}



void double_gauss_elimination(double **A, double *B, double *result, int n) {
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
            double *tempRow = A[i];
            A[i] = A[maxRow];
            A[maxRow] = tempRow;

            // Swap corresponding B values
            double tempB = B[i];
            B[i] = B[maxRow];
            B[maxRow] = tempB;
        }

        // Elimination
        for (int j = i + 1; j < n; j++) {
            

            double multiplier = A[j][i] / A[i][i];
            for (int k = i; k < n; k++) {
                A[j][k] -= multiplier * A[i][k];
            }
            B[j] -= multiplier * B[i];
        }
    }

    // Back substitution
    result[n - 1] = B[n - 1] / A[n - 1][n - 1];
    for (int i = n - 2; i >= 0; i--) {
        double sum = 0;
        for (int j = i + 1; j < n; j++) {
            sum += A[i][j] * result[j];
        }
        result[i] = (B[i] - sum) / A[i][i];
    }
}





int main()
{
    const int n=10000;
    double **A=new double*[n];
    for(int i=0;i<n;i++)
    {
        A[i]=new double[n];
    }
    double Bt[n];
    double Bg[n];
    double result[n];
    double x[n];
    double a[n];
    double b[n];
    double c[n];


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


  
    double_thomas(a,b,c,Bt,result,n);
    
    

    
     auto end = std::chrono::high_resolution_clock::now();

    
    std::chrono::duration<long double> elapsed = end - start;

    
    std::cout << "Czas wykonania: " << elapsed.count() << " sekund" << std::endl;
    

    std::cout << "Max error: " << calculate_max_error(x, result, n) << std::endl;
    
    
   
    auto start2 = std::chrono::high_resolution_clock::now();

    double_gauss_elimination(A,Bg,result,n);

    auto end2 = std::chrono::high_resolution_clock::now();

    std::chrono::duration< long double> elapsed2 = end2 - start2;

    std::cout << "Czas wykonania: " << elapsed2.count() << " sekund" << std::endl;


    
    std::cout << "Max error: " << calculate_max_error(x, result, n) << std::endl;
   

    for(int i=0;i<n;i++)
    {
        delete[] A[i];
    }
    delete[] A;

}


