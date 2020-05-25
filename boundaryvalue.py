# -*- coding: utf-8 -*-
"""
Created on Sun Apr 12 21:35:18 2020

Assignment 8:
    Finite difference and shooting method to solve the boundary value problem
    of heat flow

@author: Siddharth
"""
import numpy as np
import matplotlib.pyplot as plt
import math

h = 0.01
ta = 20
ti = 40
tf = 200
s= 1  #step size
diff = 10


def dzdx(t):
    global h, ta
    return h*(t-ta)

def RK4(ti,z0):
    
    z = np.array([z0])
    t = np.array([ti])
    global s
    
    for i in range(0, diff//s):
        
        kz1 = dzdx(t[i])
        kt1 = z[i]
        
        kz2 = dzdx(t[i] + kt1*s/2)
        kt2 = z[i] + kz1*s/2
        
        kz3 = dzdx(t[i] + kt2*s/2)
        kt3 = z[i] + kz2*s/2
        
        kz4 = dzdx(t[i] + kt3*s)
        kt4 = z[i] + kz3*s
        
        z = np.append(z, [z[i] + s*(kz1 + 2*kz2 + 2*kz3 + kz4)/6])
        t = np.append(t, [t[i] + s*(kt1 + 2*kt2 + 2*kt3 + kt4)/6])
                
    return t
        

def shooting():
    
    global ti,tf
    
    z1 = 10
    z2 = 20
    
    zi = z1 + ((z2-z1)*(tf - RK4(ti,z1)[diff//s]))/(RK4(ti,z2)[diff//s]-RK4(ti,z1)[diff//s])

    ans = RK4(ti, zi)
    
    return ans

def thomas(e,f,g,r):
    
    n = -1+diff//s
    x = np.zeros(-1+diff//s)
    
    for i in range(1,n):
        e[i] = (e[i]/f[i-1])
        f[i] = f[i] - e[i]*g[i-1]
    
        
    for i in range(1, n):
        r[i] = r[i] - e[i]*r[i-1]
        
    x[-2+diff//s] = r[n-1]/f[n-1]
    for i in range(n-2, -1, -1):
        x[i] = (r[i] - g[i]*x[i+1])/f[i]
            
    return x

def finite_difference():
    
    global h, ta, ti, tf
    
    e = np.full(-1+diff//s,-1.0)
    g= np.full(-1+diff//s, -1.0)
    f = np.full(-1+diff//s, 2+h*s*s)
    r = np.zeros(-1+diff//s)
    
    for i in range(0,-1+diff//s):
        
        r[i] = h*s*s*ta
        
        if(i == 0): r[i] = ti + h*s*s*ta
        if(i == -2+diff//s): r[i] = tf + h*s*s*ta
        
    a = np.array
    a = thomas(e,f,g,r)
    a = np.insert(a, 0, ti)
    a = np.append(a, [tf])
    
    return a


print("SHOOTING METHOD")
print(shooting())
print("FINITE DIFFERENCE METHOD")
print(finite_difference())

def plot():
    
    x = np.arange(1+diff//s)
    print(x)
    plt.plot(x, shooting().transpose(), label = 'Shooting Method', linewidth = 5)
    plt.plot(x, finite_difference(), label = 'Finite Difference')
    
    plt.xlabel('Lenght along the rod')
    plt.ylabel('Temperature')
    plt.legend(loc = 1)
    plt.title('Temperature vs Length curve')
    
    plt.show()
    
plot()


    
    




    
    
    
    
    
    
    
    
    
    
    
    
        
        
        
        
    
    
    
    