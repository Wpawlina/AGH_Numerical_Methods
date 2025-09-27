import random
import numpy as np
import math
import time

def spectral_radius(matrix):
    n= len(matrix)
    B=[[0 if i==j else -1/A[i][i]*A[i][j] for j in range(n) ] for i in range(n)]
    return max(abs(np.linalg.eigvals(B)))

def calculate_max_error(solution, results):
    return max(abs(solution[i] - results[i]) for i in range(len(solution)))

def calculate_mean_error(solution, results):
    return sum(abs(solution[i] - results[i]) for i in range(len(solution))) / len(solution)
    


def create_results(n):
    results=[0 for _ in range(n)]
    for i in range(n):
        if random.random() < 0.5:
            results[i] = 1
        else:
            results[i] = -1
    return results

def create_matrix(n):
    matrix=[[0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i==j:
                matrix[i][j] = 8
            elif j>i:
                matrix[i][j] = (-1)**(j+1)*4.5/(j+1)
            elif j==i-1:
                matrix[i][j] = 4.5/(i+1)
    return matrix

def create_B_matrix(n,A,results):
    B=[0 for _ in range(n)]
    for i in range(n):
        B[i]=sum([A[i][j]*results[j] for j in range(n)])
    return B

def jacoby_solution(A,B,x0,eps=1e-10, max_iter=1000,residual=True):
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

ns=[25,50,75,100,125,150,175,200,250,300,350,400,450,500,600,700,800,900,1000,1200,1400,1600,1800,2000]
eps=[1e-1,1e-2,1e-3,1e-5,1e-10,1e-15]

A=create_matrix(2000)
results=create_results(2000)
B=create_B_matrix(2000,A,results)
x0 = [0 for _ in range(2000)]
start_time= time.time()
solution, iterations = jacoby_solution(A, B, x0, eps=1e-15, residual=True)
end_time= time.time()
time_seconds = end_time - start_time
print("Time taken for 2000x2000 matrix:", time_seconds, "seconds",iterations,calculate_max_error(results,solution))











            
            
        
        