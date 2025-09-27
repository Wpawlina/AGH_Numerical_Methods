import math
import matplotlib.pyplot as plt

class LeastSquaresPolynomial:
    def __init__(self, f, start, end):
        self.f = f
        self.start = start
        self.end = end
        self.nodes = []         
        self.coefficients = [] 
        self.degree = None    
           
    def set_line_space(self, n, degree):
        self.degree = degree
        h = (self.end - self.start) / (n - 1)
        self.nodes = []
        for i in range(n):
            x_val = self.start + i * h
            self.nodes.append( (x_val, self.f(x_val)) )
        self.compute_coefficients()
        
    def set_czebyszew_space(self, n, degree):
       
        self.degree = degree
        self.nodes = []
        for i in range(n):
           
            x_val = (self.start + self.end)/2 + (self.end - self.start)/2 * math.cos((2*i+1)/(2*n)*math.pi)
            self.nodes.append((x_val, self.f(x_val)))
       
        self.nodes.sort(key=lambda pair: pair[0])
        self.compute_coefficients()



    def compute_coefficients(self):
     
        m = len(self.nodes)
        n = self.degree
       
        xs = [pt[0] for pt in self.nodes]
        ys = [pt[1] for pt in self.nodes]
        
        A = []
        b = []
        for k in range(n+1):
            row = []
            for j in range(n+1):
                s = 0.0
                for i in range(m):
                    s += xs[i]**(k+j)
                row.append(s)
            A.append(row)
            s = 0.0
            for i in range(m):
                s += ys[i] * (xs[i]**k)
            b.append(s)
       
        self.coefficients = self.__gaussian_elimination(A, b)

    def __gaussian_elimination(self, A, b):
    
        n = len(b)
       
        for i in range(n):
            A[i].append(b[i])
       
        for i in range(n):
           
            max_row = i
            for j in range(i+1, n):
                if abs(A[j][i]) > abs(A[max_row][i]):
                    max_row = j
            if max_row != i:
                A[i], A[max_row] = A[max_row], A[i]
            pivot = A[i][i]
            if abs(pivot) < 1e-12:
                raise Exception("Macierz osobliwa lub bliska osobliwości!")
        
            for j in range(i, n+1):
                A[i][j] /= pivot
          
            for k in range(i+1, n):
                factor = A[k][i]
                for j in range(i, n+1):
                    A[k][j] -= factor * A[i][j]
       
        x = [0.0] * n
        for i in range(n-1, -1, -1):
            x[i] = A[i][n]
            for j in range(i+1, n):
                x[i] -= A[i][j] * x[j]
        return x

    def evaluate(self, x_val):
       
        result = 0.0
        power = 1.0
        for a in self.coefficients:
            result += a * power
            power *= x_val
        return result

    def draw(self):
        
        num_points = 1000
        xs_plot = []
        ys_approx = []
        ys_orig = []
        for i in range(num_points+1):
            t = self.start + i * (self.end - self.start) / num_points
            xs_plot.append(t)
            ys_approx.append(self.evaluate(t))
            ys_orig.append(self.f(t))
       
        nodes_x = [pt[0] for pt in self.nodes]
        nodes_y = [pt[1] for pt in self.nodes]

        plt.plot(xs_plot, ys_approx, label="Aproksymacja")
        plt.plot(xs_plot, ys_orig, label="Funkcja")
        plt.scatter(nodes_x, nodes_y, color="red", label="Węzły")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.title("Aproksymacja średniokwadratowa")
        plt.legend()
        plt.show()

    def calculate_error(self):
        
        num_points = 1000
        max_err = 0.0
        sum_err = 0.0
        for i in range(num_points+1):
            t = self.start + i * (self.end - self.start) / num_points
            error = abs(self.f(t) - self.evaluate(t))
            sum_err += error
            if error > max_err:
                max_err = error
        mean_err = sum_err / (num_points+1)
        return max_err, mean_err


def f(x):
    return math.e**(-math.sin(2*x))+math.cos(2*x)


approx = LeastSquaresPolynomial(f, -3*math.pi,2*math.pi)


approx.set_line_space(50,17)


approx.draw()

print(approx.calculate_error())