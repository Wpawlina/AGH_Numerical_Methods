from matplotlib import pyplot as plt
import math
import numpy as np

def mean_error(sulution,f):
    x_calculated = [x for x, _ in sulution]
    y_calculated = [y for _, y in sulution]
    y_orginal = [f(x) for x in x_calculated]
    return sum(abs(y_calculated[i] - y_orginal[i]) for i in range(len(y_calculated))) / len(y_calculated)

def Euler(x0,y0,x1,h,df):
    values=[(x0,y0)]
    n=int((x1-x0)/h)
    for _ in range(n):
        x_next=values[-1][0] + h
        y_next=values[-1][1] + h * df(values[-1][0], values[-1][1])
        values.append((x_next, y_next))
    return values

def RungeKutta(x0,y0,x1,h,df):
    n=int((x1-x0)/h)
    values=[(x0,y0)]
    for _ in range(n):
        k1= df(values[-1][0], values[-1][1])
        k2= df(values[-1][0] + h, values[-1][1] + h * k1)
        k3= df(values[-1][0] + h/2, values[-1][1] + h * k1/2)
        k4= df(values[-1][0] + h/2, values[-1][1] + h * k2/2)
        x_next=values[-1][0] + h
        y_next= values[-1][1] + (h/6) * (k1 + 2*k2 + 2*k3 + k4)
        values.append((x_next, y_next))
    return values


def compare_to_orginal(values,f):
    x_calculated = [x for x, _ in values]
    y_calculated = [y for _, y in values]
    x_orginal = [x for x in np.arange(x_calculated[0], x_calculated[-1], 10e-5)]
    y_orginal = [f(x) for x in x_orginal]
    plt.plot(x_calculated, y_calculated, label='Calculated', color='blue')
    plt.plot(x_orginal, y_orginal, label='Original', color='red')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Comparison of Calculated and Original Functions')
    plt.legend()
    plt.grid()
    plt.show()

df= lambda x,y: 3*math.sin(3*x) * math.cos(3*x) + 3*y*math.sin(3*x) 
    
f=lambda x: math.e**(-math.cos(3*x))-math.cos(3*x) + 1

x0 = math.pi / 6
x1 = 2 * math.pi
ns=[0.5,0.4,0.2,0.1,10e-3,10e-4,10e-5]

for n in ns:
    values = Euler(x0, f(x0), x1, n, df)
    mean_errorE= mean_error(values, f)
    values2 = RungeKutta(x0, f(x0), x1, n, df)
    mean_errorR= mean_error(values2, f)
    with open('lab10.txt', 'a') as file:
        file.write(f"{n}; {mean_errorE}; {mean_errorR}\n")
    


   
