"""
Assignment-7 

Using RK4 method to solve system of differential equations

"""

import numpy as np
import matplotlib.pyplot as plt
import math

 #Initial Conditions
 
x0 = 0
y0 = 6786 #enter the value of H here
vx0 = 7.66
vy0 = 0

x = np.array([x0]) #to store x coordinates
y = np.array([y0]) #To store y coordinates
vx = np.array([vx0]) #To store velocity along the x direction
vy = np.array([vy0]) #To store velocity along the y direction

def f(x,y):
    global k,M
    
    #we can switch the order of x and y to give f1 and f2 as given in the question.
    
    return (-398588)*(x/math.pow(x**2 + y**2, 1.5))    


def rk4(h):
    
    global x,y,vx,vy
    
    #i = 0
    
    for i in range(0,100):
        
        """
        kxi = is the ith constant for x
        kyi = is the ith constant for y
        kvxi = is the ith constant for vx
        kvyi = is the ith constant for vy
        
        where ki have the usual definition as the variables used in the increment function
        
        the double differential equation 
        
        (d^2/dt^2)x = f(x,y)
        
        can be broken into a system of first order linear differential equations
        
        (d/dt)v = f(x,y)
        (d/dt)x = v
        
        same is repeated for the y variable
        
        Orbit is found by finding all the x and y variables
        
        """
        
        kvx1 = f(x[i], y[i])
        kvy1 = f(y[i], x[i])
        kx1 = vx[i]
        ky1 = vy[i]
        
        kvx2 = f(x[i] + kx1*h/2, y[i] + ky1*h/2)
        kvy2 = f(y[i] + ky1*h/2, x[i] + kx1*h/2)
        kx2 = vx[i] + kvx1*h/2
        ky2 = vy[i] + kvy1*h/2
        
        kvx3 = f(x[i] + kx2*h/2, y[i] + ky2*h/2)
        kvy3 = f(y[i] + ky2*h/2, x[i] + kx2*h/2)
        kx3 = vx[i] + kvx2*h/2
        ky3 = vy[i] + kvy2*h/2
        
        kvx4 = f(x[i] + kx3*h, y[i] + ky3*h)
        kvy4 = f(y[i] + ky3*h, x[i] + kx3*h)
        kx4 = vx[i] + kvx3*h
        ky4 = vy[i] + kvy3*h
        
        #Using all the k's to calculate the increment function and thereby find (i+1)th variable
        
        vx = np.append(vx, [vx[i] + h*(kvx1 + 2*kvx2 + 2*kvx3 + kvx4)/6])
        vy = np.append(vy, [vy[i] + h*(kvy1 + 2*kvy2 + 2*kvy3 + kvy4)/6])
        x = np.append(x, [x[i] + h*(kx1 + 2*kx2 + 2*kx3 + kx4)/6])
        y = np.append(y, [y[i] + h*(ky1 + 2*ky2 + 2*ky3 + ky4)/6])
                
        #i+=1
        
    print(x)
    print(y)
       
#function to plot the found orbit
       
def plot():
    
    plt.plot(x, y, label = 'Step size: 60\n#of steps:100')
    
    plt.xlabel('X axis')
    plt.ylabel('Y axis')
    plt.legend(loc = 1)
    plt.title('Orbit of planets')
    
    plt.show()
    
rk4(60)

plot()




    
    
        
        
        
        
        
        
        
        
        
        
        
       
            
            
        
        
            
            
            
    
            
            
            
            
        
        
    
    
    
    
    
    
    
    
    
    
