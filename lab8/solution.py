import numpy as np
from matplotlib import pyplot as plt
import math
import random
import time


def spectral_radius(matrix):
    n= len(matrix)
    B=[[0 if i==j else -1/A[i][i]*A[i][j] for j in range(n) ] for i in range(n)]
    return max(abs(np.linalg.eigvals(B)))

def gaussian_elimination(A,B):
    n=len(A)
    for i in range(n-1):
        for j in range(i+1,n):
            mulitplier=A[j][i]/A[i][i]
            for k in range(i,n):
                A[j][k]-=mulitplier*A[i][k]
            B[j]-=mulitplier*B[i]
    solution=[0 for _ in range(n)]
    solution[-1]=B[-1]/A[-1][-1]
    for i in range(n-2,-1,-1):
        sum=0
        for j in range(i+1,n):
            sum+=A[i][j]*solution[j]
        solution[i]=(B[i]-sum)/A[i][i]
    return solution    

def create_three_diagonal_matrix(A):
    n = len(A)
    a = [0] * n  
    b = [0] * n  
    c = [0] * n 
    for i in range(n):
        if(i-1>=0):
            a[i]=A[i][i-1]
        
        b[i]=A[i][i];   
        if(i+1<n):
            c[i]=A[i][i+1]
    return a, b, c        
    
        
        
   
    




def thomas(a, b,c, d,) :
    n = len(a)
    for i in range(1,n):
       m=a[i]/b[i-1]
       b[i]-=m*c[i-1]
       d[i]-=m*d[i-1]
    
    result=[0 for _ in range(n)]
    result[n-1]=d[n-1]/b[n-1]
    for i in range(n-2,-1,-1) :
         result[i] = (d[i] - c[i] * result[i + 1]) / b[i]
    return result
    


    
    
def create_matrix1(n):
    A=[[0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i==0:
                A[i][j]=1
            else:
                A[i][j]=1/(i+1+j+1-1)
    return A

def create_matrix2(n):
    A=[[0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if j>=i:
                A[i][j]=2*(i+1)/(j+1)
    
    for i in range(n):
        for j in range(n):
            if j<i:
                A[i][j]=A[j][i]
            
    return A   

def create_matrix3(n):
    A=[[0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j:
                A[i][j] = -2 * (i + 1) - 4
            elif j == i + 1:
                A[i][j] = i + 1
            elif j == i - 1:
                A[i][j] = 2 / (i + 1)
            else:
                A[i][j] = 0
    return A 
   
def calculate_max_error(solution, results):
    return max(abs(solution[i] - results[i]) for i in range(len(solution)))

    


def create_results(n):
    results=[0 for _ in range(n)]
    for i in range(n):
        if random.random() < 0.5:
            results[i] = 1
        else:
            results[i] = -1
    return results



def create_B_matrix(n,A,results):
    B=[0 for _ in range(n)]
    for i in range(n):
        B[i]=sum([A[i][j]*results[j] for j in range(n)])
    return B

def jacoby_solution(A,B,x0,eps=1e-10, max_iter=1000,residual=False):
    n= len(A)
    for iter in range(max_iter):
        x0_new = [0 for _ in range(n)]
        for i in range(n):
            sum1=0
            for j in range(n):
                if i != j:
                    sum1 += A[i][j] * x0[j]
            x0_new[i] = (B[i] - sum1) / A[i][i]
                
        if residual:
            if math.sqrt(sum(sum(A[i][j]*x0_new[j] for j in range(n))- B[i] for i in range(n))**2) < eps:
                
                return (x0_new,iter+1)
        else :
            if math.sqrt(sum([(x0_new[i] - x0[i]) ** 2 for i in range(n)])) < eps:
                return (x0_new,iter+1)
        x0 = x0_new
    return (x0, max_iter)
    








ns = [3,5,10,15,20,30,40,50,75,100,150,200,250,300,350,400,450,500,600,700,800,900,1000]
eps=[1e-1, 1e-3, 1e-5, 1e-10, 1e-15]

for n in ns:
    A= create_matrix3(n)
    result = create_results(n)
    B = create_B_matrix(n, A, result)
    a, b, c = create_three_diagonal_matrix(A)
    x0 = [0 for _ in range(n)]
    print(str(n)+";"+str(np.linalg.cond(A,2)))
   
            
            
            
            
        
       
           
        
        
          
        



    
    
    
    
  


    











           
            
            
            
            
        
    