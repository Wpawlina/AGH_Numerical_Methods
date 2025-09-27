import matplotlib.pyplot as plt
import math

class CubicSpline:
    def __init__(self, func, start, end, num_points, boundary_condition=False, clamped_values=(0, 0)):
        self.x = [start + i * (end - start) / (num_points - 1) for i in range(num_points)]
        self.y = [func(xi) for xi in self.x]
        self.n = len(self.x) - 1
        self.boundary_condition = boundary_condition
        self.clamped_values = clamped_values
        self.h = [self.x[i+1] - self.x[i] for i in range(self.n)]
        self.M = self._compute_M()
        self.f=func
    
    def _solve_tridiagonal(self, a, b, c, d):
        
        n = len(d)
        c_ = [0] * (n-1)
        d_ = [0] * n
        x = [0] * n
        
        
        c_[0] = c[0] / b[0]
        d_[0] = d[0] / b[0]
        
        for i in range(1, n-1):
            denom = b[i] - a[i] * c_[i-1]
            if denom == 0:
                raise ValueError("Dzielenie przez zero w eliminacji")
            c_[i] = c[i] / denom
            d_[i] = (d[i] - a[i] * d_[i-1]) / denom
        
        if n > 1: 
            denom = b[n-1] - a[n-2] * c_[n-2] if n > 2 else b[n-1]
            d_[n-1] = (d[n-1] - a[n-2] * d_[n-2]) / denom
        
       
        x[n-1] = d_[n-1]
        for i in range(n-2, -1, -1):
            x[i] = d_[i] - c_[i] * x[i+1]
        
        return x
    
    def _compute_M(self):
       
        a = self.h[:-1]
        b = [2 * (self.h[i-1] + self.h[i]) for i in range(1, self.n)]
        c = self.h[1:]
        d = [6 * ((self.y[i+1] - self.y[i]) / self.h[i] - (self.y[i] - self.y[i-1]) / self.h[i-1]) for i in range(1, self.n)]
        
        if self.boundary_condition:
            b.insert(0, 2 * self.h[0])
            d.insert(0, 6 * ((self.y[1] - self.y[0]) / self.h[0] - self.clamped_values[0]))
            b.append(2 * self.h[-1])
            d.append(6 * (self.clamped_values[1] - (self.y[-1] - self.y[-2]) / self.h[-1]))
        else:
            b.insert(0, 1)
            d.insert(0, 0)
            b.append(1)
            d.append(0)
        
        a.insert(0, 0)
        c.append(0)
        
        return self._solve_tridiagonal(a, b, c, d) 
    
    def calculate(self, x_val):
        i = 0
        while i < self.n and x_val > self.x[i+1]:
            i += 1
        
        h_i = self.x[i+1] - self.x[i]
        A_i = (self.x[i+1] - x_val) / h_i
        B_i = (x_val - self.x[i]) / h_i
        
        return (A_i * self.y[i] + B_i * self.y[i+1] +
                ((A_i**3 - A_i) * self.M[i] + (B_i**3 - B_i) * self.M[i+1]) * (h_i**2) / 6)
    
    def draw(self):
        
        start, end = self.x[0], self.x[-1]
        sx = [start + i * (end - start) / 100000 for i in range(100001)]
        y_interp = [self.calculate(i) for i in sx]
        y_actual = [ self.f(i) for i in sx] 
        
        plt.plot(sx, y_interp, '-', label='Interpolacja sklejana')
        plt.plot(sx, y_actual, '-', label='Oryginalna funkcja')
        plt.scatter(self.x, self.y, color='red', label='Węzły')
        plt.legend()
        plt.show()
    
    def calculate_error(self):
        sx = [self.x[0] + i * (self.x[-1] - self.x[0]) / 1000 for i in range(1000)]
        y_interp = [self.calculate(i) for i in sx]
        y_actual = [self.f(i) for i in sx] 
        
        errors = [abs(y_actual[i] - y_interp[i]) for i in range(1000)]
        max_error = max(errors)
        mean_error = sum(errors) / len(errors)
        
        return max_error, mean_error
    
    
class QuadraticSpline:
    def __init__(self, func, start, end, num_points, boundary_condition=False, clamped_values=(0, 0)):
       
        self.x = [start + i*(end-start)/(num_points-1) for i in range(num_points)]
        self.y = [func(xi) for xi in self.x]
        self.n = len(self.x) - 1 
        self.h = [self.x[i+1] - self.x[i] for i in range(self.n)]
        self.boundary_condition = boundary_condition
        self.clamped_values = clamped_values  
        self.f=func
      
        self.delta = [(self.y[i+1]-self.y[i]) / self.h[i] for i in range(self.n)]
 
        self.m = self._compute_node_derivatives()
     
        self.c = [ (self.y[i+1] - self.y[i] - self.m[i]*self.h[i]) / (self.h[i]**2)
                   for i in range(self.n) ]
        
    def _compute_node_derivatives(self):
        m = [0]*(self.n+1)
        if self.boundary_condition:
            
            m[0] = self.clamped_values[0]
            m[self.n] = self.clamped_values[1]
            
            if self.n < 2:
              
                m[1] = 2*self.delta[0] - m[0]
                return m
            
            N = self.n - 1  
    
            r = [0] * N
            r[0] = 2*(self.delta[0] + self.delta[1]) - m[0]
            for i in range(1, N-1):
                r[i] = 2*(self.delta[i] + self.delta[i+1])
            r[N-1] = 2*(self.delta[self.n-2] + self.delta[self.n-1]) - m[self.n]
            
           
            diag = [2] * N
            lower = [1] * (N - 1)
            upper = [1] * (N - 1)
            
            M_internal = self._solve_tridiagonal(lower, diag, upper, r)
         
            for i in range(1, self.n):
                m[i] = M_internal[i-1]
        else:
          
            m[0] = self.delta[0]
            m[self.n] = self.delta[-1]
        
            for i in range(0, self.n):
                m[i+1] = 2*self.delta[i] - m[i]
        return m
    
    def _solve_tridiagonal(self, lower, diag, upper, r):
        
        n = len(diag)
        c_prime = [0]*n
        d_prime = [0]*n
        c_prime[0] = upper[0] / diag[0]
        d_prime[0] = r[0] / diag[0]
        for i in range(1, n):
            denom = diag[i] - lower[i-1] * c_prime[i-1]
            if i < n:
                if i < n - 1:
                    c_prime[i] = upper[i] / denom
                d_prime[i] = (r[i] - lower[i-1] * d_prime[i-1]) / denom
        x = [0]*n
        x[-1] = d_prime[-1]
        for i in range(n-2, -1, -1):
            x[i] = d_prime[i] - c_prime[i] * x[i+1]
        return x
    
    def calculate(self, x_val):
    
        i = 0
        while i < self.n and x_val > self.x[i+1]:
            i += 1
        dx = x_val - self.x[i]
        return self.y[i] + self.m[i]*dx + self.c[i]*(dx**2)

    def draw(self):
      
        start, end = self.x[0], self.x[-1]
        sx = [start + i * (end - start) / 1000 for i in range(1001)]
        y_interp = [self.calculate(i) for i in sx]
        y_actual = [self.f(i) for i in sx]

        plt.plot(sx, y_interp, '-', label='Interpolacja kwadratowa')
        plt.plot(sx, y_actual, '-', label='Oryginalna funkcja')
        plt.scatter(self.x, self.y, color='red', label='Węzły')
        plt.legend()
        plt.show()

    def calculate_error(self):
        sx = [self.x[0] + i * (self.x[-1] - self.x[0]) / 1000 for i in range(1000)]
        y_interp = [self.calculate(i) for i in sx]
        y_actual = [self.f(i) for i in sx]

        errors = [abs(y_actual[i] - y_interp[i]) for i in range(1000)]
        max_error = max(errors)
        mean_error = sum(errors) / len(errors)

        return max_error, mean_error

f=lambda x: math.e**(-math.sin(2*x))+math.cos(2*x)
spline =CubicSpline(f, -2*math.pi, 3*math.pi,45, boundary_condition=False,clamped_values=(-2, -2))
# spline.draw()
print(spline.calculate_error())
