import math
from matplotlib import pyplot as plt

class FourierApprox:
    def __init__(self,n,m,f,start,end):
        self.f=f
        self.n=n
        self.m=m
        self.start=start
        self.end=end
        self.__line_space()
        self.__a0()
        self.__an()
        self.__bn()
    def __transform(self,x):
        return ((x-self.start)/(self.end-self.start))*2*math.pi-math.pi
    
    def __line_space(self):
        dist=(self.end-self.start)/(self.n-1)
        self.x = [self.start+i*dist for i in range(self.n)]
        self.xp=[ self.__transform(i) for i in self.x]
    
    def __a0(self):
        self.a0= (2/self.n)*sum([self.f(self.x[i]) for i in range(self.n)])
    
    def __an(self):
        self.an=[0]*(self.m-1)
        for i in range(1,self.m):
            self.an[i-1]= (2/self.n)*sum([self.f(self.x[j])*math.cos(i*self.xp[j]) for j in range(self.n)])
    def __bn(self):
        self.bn=[0]*(self.m-1)
        for i in range(1,self.m):
            self.bn[i-1]= (2/self.n)*sum([self.f(self.x[j])*math.sin(i*self.xp[j]) for j in range(self.n)])
    
    def calculate(self,x):
        xp=self.__transform(x)
        result=self.a0/2
        for i in range(1,self.m):
            result+=self.an[i-1]*math.cos(i*xp)+self.bn[i-1]*math.sin(i*xp)
        return result
    
    def draw(self):

        
        start=self.start
        end=self.end
        sx=[start+i*(end-start)/1000 for i in range(1001)]
        y2=[self.calculate(i) for i in sx]
        y3=[self.f(i) for i in sx]

        plt.plot(sx,y3,label="Funkcja")
        plt.plot(sx,y2,label="Aproksymacja")
        y=[self.f(i) for i in self.x]
        
        plt.scatter(self.x,y,label="Orginalne punkty",color="red")
        plt.legend()
        plt.title("Porównanie")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.show()
        
    def calculate_error(self):
    
        sx=[self.start+i*(self.end-self.start)/1000 for i in range(1000)]
        y1=[self.f(i) for i in sx]
        y2=[self.calculate(i) for i in sx]
        max_error=0
        mean_error=0
        for i in range(1000):
            max_error=max(max_error,abs(y1[i]-y2[i]))
            mean_error+=abs(y1[i]-y2[i])
        mean_error/=1001
        return max_error,mean_error
                  
      
      

def f(x):
    return math.e**(-math.sin(2*x))+math.cos(2*x)
  

fa=FourierApprox(7,3,f,-3*math.pi,3*math.pi)
fa.draw()
print(fa.calculate_error())

        