from functools import reduce
from matplotlib import pyplot as plt
import math

class Hermit():
    def __init__(self,f,df,start,end):
        self.start=start
        self.end=end
        self.df=df   
        self.fun=f
        
    def set_line_space(self,n):
        self.n=n
        h=(self.end-self.start)/(n-1)
        self.nodes=[[self.start+i*h,self.fun(self.start+i*h),self.df(self.start+i*h)] for i in range(n)]
        self.e=reduce(lambda x,acc: x+acc  ,map(lambda x:len(x)-1,self.nodes))
        self.b=[0 for _ in range(self.e)]
        self.F={}
        self.__calcualte_b()
    
    
    def set_czybyszew_space(self,n):
        self.n=n
        h=(self.end-self.start)/(n-1)
        x=[(self.start+self.end)/2+(self.end-self.start)/2*math.cos((2*i+1)/(2*n)*math.pi) for i in range(n)]
        self.nodes=[[x[i],self.fun(x[i]),self.df(x[i])] for i in range(n)]
        self.e=reduce(lambda x,acc: x+acc  ,map(lambda x:len(x)-1,self.nodes))
        self.F={}
        self.b=[0 for _ in range(self.e)]
        self.__calcualte_b()
           
    

    def __calculate_F(self,index):
        if self.F.get(index)!=None: 
            return self.F[index]
        else:
            if len(index)==1:
                self.F[index]=self.nodes[index[0]][1]
            elif len(index)==2 and  index[0]==index[1]:
                self.F[index]=self.nodes[index[0]][2]
            elif len(index)==3 and index[0]==index[1] and index[1]==index[2]:
                self.F[index]=self.nodes[index[0]][3]/2
            else:
                self.F[index]=(self.__calculate_F(index[1:])-self.__calculate_F(index[:-1]))/(self.nodes[index[-1]][0]-self.nodes[index[0]][0])
            return self.F[index]
                
            
    
    def __calcualte_b(self):
        j=0
        k=0
        index=[]
        for node in self.nodes:
            for i in range(1,len(node)):
                index.append(j)                
                self.b[k]=self.__calculate_F(tuple(index))
                k+=1
            j+=1
            
    def calculate(self,x):
        if len(self.F)==0:
            raise ValueError("You must set the space first")
        result=0
        j=0
        multiplier=1
        for node in self.nodes:
            for i in range(1,len(node)):
                  result+=self.b[j]*multiplier
                  multiplier*=(x-node[0])
                  j+=1
        return result
                  
                  
    def draw(self):
        if len(self.F)==0:
            raise ValueError("You must set the space first")
        
        start=self.nodes[0][0]
        end=self.nodes[-1][0]
        sx=[start+i*(end-start)/1000 for i in range(1001)]
        y2=[self.calculate(i) for i in sx]
        y3=[self.fun(i) for i in sx]
        x=list(map(lambda x:x[0],self.nodes))
        y=list(map(lambda x:x[1],self.nodes))
        plt.plot(sx,y3,label="Funkcja")
        plt.plot(sx,y2,label="Interpolacja")
        
        plt.scatter(x,y,label="Orginalne punkty",color="red")
        plt.legend()
        plt.title("Porównanie")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.show()
        
    def calculate_error(self):
        if len(self.F)==0:
            raise ValueError("You must set the space first")
        sx=[self.start+i*(self.end-self.start)/1000 for i in range(1000)]
        y1=[self.fun(i) for i in sx]
        y2=[self.calculate(i) for i in sx]
        max_error=0
        mean_error=0
        for i in range(1000):
            max_error=max(max_error,abs(y1[i]-y2[i]))
            mean_error+=abs(y1[i]-y2[i])
        mean_error/=1001
        return max_error,mean_error
                  
              

f=lambda x: math.e**(-math.sin(2*x))+math.cos(2*x)
df=lambda x: math.e**(-math.sin(2*x))*(-2*math.cos(2*x))-2*math.sin(2*x)
hermit=Hermit(f,df,-2*math.pi,3*math.pi)     
hermit.set_czybyszew_space(8)
hermit.draw()
print(hermit.calculate_error())



  
                
               
            
        
    
            
            
        
        
        
            
        