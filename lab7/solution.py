import math

def calculate_zero_newton(f, df,start, residual=False, eps=1e-6, max_iter=1000):
    xn = start
    for i in range(max_iter):
        dfxn = df(xn)
        if abs(dfxn) < 1e-20:
            raise ZeroDivisionError("Pochodna bliska zeru")
        xn1 = xn - f(xn) / dfxn
        if residual:
            if abs(f(xn1)) < eps:
                print(f"[Newton] Iteracja {i+1}: x = {xn1}")
                return (xn1,i+1)
        else:
            if abs(xn1 - xn) < eps:
                print(f"[Newton] Iteracja {i+1}: x = {xn1}")
                return (xn1,i+1)
        xn = xn1
    raise RuntimeError("Nie osiągnięto zbieżności (Newton)")

def calculate_zero_secant(f, start1,start2,residual=False, eps=1e-6, max_iter=1000):
    x1, x2 = start1, start2
    for i in range(max_iter):
        denom = f(x2) - f(x1)
        if abs(denom) < 1e-20:
            raise ZeroDivisionError("Dzielnik bliski zeru")
        x3 = x2 - f(x2) * (x2 - x1) / denom
        if residual:
            if abs(f(x3)) < eps:
                print(f"[Sieczne] Iteracja {i+1}: x = {x3}")
                return (x3,i+1)
        else:
            if abs(x3 - x2) < eps:
                print(f"[Sieczne] Iteracja {i+1}: x = {x3}")
                return (x3,i+1)
        x1, x2 = x2, x3
    raise RuntimeError("Nie osiągnięto zbieżności (Sieczne)")

f=lambda x: (x-1)*math.exp(-13*x)+x**15
df=lambda x: (14-13*x)*math.exp(-13*x)+15*x**14

a=-0.1


epss=[1e-2, 1e-3, 1e-4, 1e-5, 1e-7,1e-10,1e-15]


    



    
        
    
    