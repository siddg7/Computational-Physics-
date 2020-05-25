# -*- coding: utf-8 -*-
"""
Created on Sat May  9 20:45:30 2020

@author: Siddharth

Numerical Integration
    
    Specifically a multiple integration problem Chapra and Canale 21.9
"""

import numpy as np

def func(x,y):
    return 2*x*y + 2*x - x*x - 2*y*y + 72

def trapezoid(xi,xf,yi,yf,n):
    
    hx = (xf-xi)/n
    hy = (yf-yi)/n
    
    f = np.zeros([n+1,n+1])
    s = np.zeros([n+1,1])
    
    for i in range(0,n+1):
        for j in range(0,n+1):
            
            f[i,j] = func(hx*i,hy*j)
            
    m1 = 0
    
    for i in range(0,n+1):
        m1 = 0
        for j in range(1,n): m1 += f[j,i]
        s[i] = (xf-xi)*(f[0,i] + 2*m1 + f[n,i])/(2*n)
        
    m1 = 0
    
    for j in range(1,n): m1 += s[j]
    
    return (yf-yi)*(s[0]+2*m1+s[n])/(2*n)
        

def simpson3(xi,xf,yi,yf,n):
    
    if(n%2 != 0): return 0
    hx = (xf-xi)/n
    hy = (yf-yi)/n
    
    f = np.zeros([n+1,n+1])
    s = np.zeros([n+1,1])
    
    for i in range(0,n+1):
        for j in range(0,n+1):
            
            f[i,j] = func(hx*i,hy*j)
            
    m1 = 0
    m2 = 0
    
            
    for i in range(0,n+1):
        m1 = 0
        m2 = 0
        for j in range(1,n,2): m1 += f[j,i]
        for j in range(2,n-1,2): m2 += f[j,i]
        
        s[i] = (xf-xi)*(f[0,i] + 4*m1 + 2*m2 + f[n,i])/(3*n)
        
    m1 = 0
    m2 = 0
        
    for j in range(1,n,2): m1 += s[j]
    for j in range(2,n-1,2): m2 += s[j]
        
    return (yf-yi)*(s[0] + 4*m1 + 2*m2 + s[n])/(3*n)

print('Output of Simple Trapezoid Rule')
print(trapezoid(0,8,0,6,1)/48)
print('Output of Multiple Trapezoid Rule with n = 5')
print(trapezoid(0,8,0,6,5)/48)
print('Output of Multiple Trapezoid Rule with n = 10')
print(trapezoid(0,8,0,6,10)/48)
print('Output of Multiple Trapezoid Rule with n = 100')
print(trapezoid(0,8,0,6,100)/48)

print('\n')

print('Output of simple Simpsons 1/3 Rule')
print(simpson3(0,8,0,6,2)/48)
print('Output of Multiple Simpsons 1/3 Rule with n = 10')
print(simpson3(0,8,0,6,10)/48)
        
    
        
    
        
        
        
        
        
            
    
        
        
    
    
    
    
    
    
    
    
