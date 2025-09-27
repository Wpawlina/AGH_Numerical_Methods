import math
import matplotlib.pyplot as plt

class Lagrange():
    def __init__(self,start,end,function ):
        self.start = start
        self.end = end
        self.f=function
        self.denominator=[]
        
    def __create_denominator(self):
        for i in range(self.numpoints):
            self.denominator.append(1)
            for j in range(self.numpoints):
                if i!=j:
                    self.denominator[i]*=(self.x[i]-self.x[j])
                    
         
    def set_line_space(self,numpoints):
        self.numpoints = numpoints
        dist=(self.end-self.start)/(self.numpoints-1)
        self.x = [self.start+i*dist for i in range(self.numpoints)]
        self.y = [self.f(i) for i in self.x]
        self.denominator=[]
        self.__create_denominator()
        
        
    def set_czybyszew_space(self,numpoints):
        self.numpoints = numpoints
        self.x = [(self.start+self.end)/2+(self.end-self.start)/2*math.cos((2*i+1)/(2*self.numpoints)*math.pi) for i in range(self.numpoints)]
        self.y = [self.f(i) for i in self.x]
        self.denominator=[]
        self.__create_denominator()
        
        
    def calculate(self,x):
        if len(self.denominator)==0:
            raise ValueError("You must set the space first")
        result=0
        numerator=[]
        for i in range(self.numpoints):
            numerator.append(self.y[i])
            for j in range(self.numpoints):
                if i!=j:
                    numerator[i]*=(x-self.x[j])
            result+=numerator[i]/self.denominator[i]
        return result
            
    def draw(self):
        if len(self.denominator)==0:
            raise ValueError("You must set the space first")

        sx=[self.start+i*(self.end-self.start)/1000 for i in range(1001)]
        y1=[self.f(i) for i in sx]
        y2=[self.calculate(i) for i in sx]
        plt.plot(sx,y1,label="Funkcja")
        plt.plot(sx,y2,label="Interpolacja")
        plt.scatter(self.x,self.y,label="Orginalne punkty",color="red")
        plt.legend()
        plt.title("Porównanie")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.show()
    
    
    def calculate_error(self):
        if len(self.denominator)==0:
            raise ValueError("You must set the space first")
        sx=[self.start+i*(self.end-self.start)/1000 for i in range(1001)]
        y1=[self.f(i) for i in sx]
        y2=[self.calculate(i) for i in sx]
        max_error=0
        mean_error=0
        for i in range(1001):
            max_error=max(max_error,abs(y1[i]-y2[i]))
            mean_error+=abs(y1[i]-y2[i])
        mean_error/=1001
        return max_error,mean_error
            
   
        
        

class Newton():
    def __init__(self,start,end,function ):
        self.start = start
        self.end = end
        self.f=function
        self.a=[]
        
    def __create_a(self): 
        self.a=[i for i in self.y]
        for i in range(1, self.numpoints):
            for j in range(self.numpoints - 1, i - 1, -1):
                self.a[j] = (self.a[j] - self.a[j - 1]) / (self.x[j] - self.x[j - i])
            
        
                    
         
    def set_line_space(self,numpoints):
        self.numpoints = numpoints
        if numpoints<2:
            raise ValueError("Number of points must be greater than 1")
        
        dist=(self.end-self.start)/(self.numpoints-1)
        self.x = [self.start+i*dist for i in range(self.numpoints)]
        self.y = [self.f(i) for i in self.x]
        self.a=[]
        self.__create_a()
       
    
        
        
    def set_czybyszew_space(self,numpoints):
        self.numpoints = numpoints
        self.x = [(self.start+self.end)/2+(self.end-self.start)/2*math.cos((2*i+1)/(2*self.numpoints)*math.pi) for i in range(self.numpoints)]
        self.y = [self.f(i) for i in self.x]
        self.a=[]
        self.__create_a()
        
        
        
    def calculate(self,x):
        if len(self.a)==0:
            raise ValueError("You must set the space first")
        result=0
        for i in range(self.numpoints):
            next=self.a[i]
            for j in range(i):
                next*=(x-self.x[j])
            result+=next
        return result
        
            
    def draw(self):
        if len(self.a)==0:
            raise ValueError("You must set the space first")

        sx=[self.start+i*(self.end-self.start)/1000 for i in range(1001)]
        y1=[self.f(i) for i in sx]
        y2=[self.calculate(i) for i in sx]
        plt.plot(sx,y1,label="Funkcja")
        plt.plot(sx,y2,label="Interpolacja")
        plt.scatter(self.x,self.y,label="Orginalne punkty",color="red")
        plt.legend()
        plt.title("Porównanie")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.show()
    
    def calculate_error(self):
        if len(self.a)==0:
            raise ValueError("You must set the space first")
        sx=[self.start+i*(self.end-self.start)/1000 for i in range(1001)]
        y1=[self.f(i) for i in sx]
        y2=[self.calculate(i) for i in sx]
        max_error=0
        mean_error=0
        for i in range(1001):
            max_error=max(max_error,abs(y1[i]-y2[i]))
            mean_error+=abs(y1[i]-y2[i])
        mean_error/=1001
        return max_error,mean_error
    
   
   
        
        
        
        
        
def draw_function(start,end,f):
    x=[start+i*(end-start)/1000 for i in range(1001)]
    y=[f(i) for i in x]
    plt.plot(x,y)
    plt.title("Funkcja")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.show()
    
def draw_points(points):
    x=[i[0] for i in points]
    y=[i[1] for i in points]
    plt.plot(x,y)
    plt.title("Punkty")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.show()
        
def draw_compare(start,end,f1,f2):
    x=[start+i*(end-start)/1000 for i in range(1001)]
    y1=[f1(i) for i in x]
    y2=[f2(i) for i in x]
    plt.plot(x,y1,label="Metoda Newton'a")
    plt.plot(x,y2,label="Metoda Lagrange'a")
    plt.legend()
    plt.title("Porównanie")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.show()
        
        

def f(x):
    return math.e**(-math.sin(2*x))+math.cos(2*x)



# lagrange =Lagrange(-2*math.pi,3*math.pi,f)

# lagrange.set_czybyszew_space(420)
# print(lagrange.calculate_error())
n=22
newton=Newton(-2*math.pi,3*math.pi,f)
lagrange=Lagrange(-2*math.pi,3*math.pi,f)

lagrange.set_czybyszew_space(n)
lagrange.draw()
print(lagrange.calculate_error())

# lagrange.set_czybyszew_space(n)
# lagrange.draw()


    
    
# newton.set_czybyszew_space(9)
# lagrange.set_czybyszew_space(9)
# draw_compare(-2*math.pi,3*math.pi,newton.calculate,lagrange.calculate)





# newton.draw()








# newton = Newton(-2*math.pi,3*math.pi,f)
# newton.set_line_space(n)

# draw_function(-2*math.pi,3*math.pi,lagrange.calculate)

# draw_function(-2*math.pi,3*math.pi,newton.calculate)

# print(newton.calculate_error())
# newton.draw()



#linspace Lagrange  najlepsze przyblizenie 8 efekt Rungego przy 9 i nie ustepuje,Newton najlepsze dla 8 efekt Rungego przy 9 do 66 wtedy blad numeryczny   , Lagrange i Newton rozne od 66

#czybyszew  Lagrange  najlepsze przyblizenie okolo 400 potem blad numeryczny efekt Rungego nie wystepuje ,Newton najlepsze dla 34 efekt Rungego nie wystepuje   , Lagrange i Newton rozne od 35



        
        
        
    
        
       
        
