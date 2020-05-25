# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 23:02:40 2020

Non Linear Regression

@author: Siddharth
"""

import numpy as np
import math
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

############# Linear Algebra Functions for Finding the inverse ################

def backward_substitution(A,B):
    
    x = np.zeros(B.shape)
    n = A.shape[0]
    
    x[n-1] = B[n-1]/A[n-1][n-1]
    for i in range(n-2, -1, -1):
        
        x[i] = (B[i] - np.sum(np.multiply(A[i],x)))/A[i][i]
        
    return x

def forward_substitution(A,B):
    
    x = np.zeros(B.shape)
    n = A.shape[0]

           
    x[0] = B[0]/A[0][0]

    for i in range(1,n):
        
        x[i] = (B[i] - np.sum(np.multiply(A[i],x)))/A[i][i]
        
    return x

def LUdecompose(A):
    
    L = np.identity(A.shape[0])
    U = np.zeros(A.shape)
    U = A
    n = A.shape[0]
    
    for i in range(0,n-1):
        
        for j in range(i+1,n):
            
            L[j][i] = (U[j][i]/U[i][i])
            U[j] = U[j] - (U[j][i]/U[i][i])*U[i]
            
    return np.array([L,U])

def solvewithLU(A,B):
    
    x = np.array(B.shape)
    d = np.array(B.shape)
    
    L = np.array(A.shape)
    U = np.array(A.shape)
    L = LUdecompose(A)[0]
    U = LUdecompose(A)[1]
    
    d = forward_substitution(L,B)
    x = backward_substitution(U,d)
    
    return x


def inverse(A):
    e1 = np.array(([1.0,0.0]))
    e2 = np.array(([0.0,1.0]))
    return np.array([solvewithLU(A,e1),solvewithLU(A,e2)])

###############################################################################

n = 5 # Number of data points

def f(x,a0,a1):
    
    return a0*(1-math.exp(-a1*x))

def dfda0(x,a1):

    return 1 - math.exp(-a1*x)

def dfda1(x,a0,a1):
    
    return a0*x*math.exp(-a1*x)

def gauss_newton(x,y,a0,a1):
    #x and y are numpy arrays(data points), a0 and a1 are the intital values
    
    a = np.zeros((2,1))
    dela = np.zeros((2,1))
    error = np.array([2,1])
        
    Z = np.zeros((n,2))
    D = np.zeros((n,1))
    
    a[0] = a0
    a[1] = a1
    
    while(error[0] > 0.0001 and error[1] > 0.0001):
        
        for i in range(0,n):
            
            Z[i][0] = dfda0(x[i],a[1])
            Z[i][1] = dfda1(x[i],a[0],a[1])
            
            D[i] = y[i] - f(x[i],a[0], a[1])
            
        dela = inverse(np.dot(Z.transpose(),Z))@(Z.transpose()@D)
        
        a = a + dela
        
        error[0] = abs((dela[0]/a[0]))*100
        error[1] = abs((dela[1]/a[1]))*100
        
    return a

x = np.array([0.25,0.75,1.25,1.75,2.25])
y = np.array([0.28,0.57,0.68,0.74,0.79])
a0 = 1
a1 = 1

print(gauss_newton(x,y,a0,a1))

####### Using Scipy to fit the given points and check the solution ##########

def test(x, a0, a1):
    one = np.ones(x.shape)
    return a0*(one-np.exp(-a1*x))

param, param_cov = curve_fit(test,x,y)
print(param)

        
def plot(x,y):
    
    plt.scatter(x,y)
    x_values = np.arange(0.25,2.5,0.1)
    plt.plot(x_values,test(x_values,gauss_newton(x,y,a0,a1)[0],gauss_newton(x,y,a0,a1)[1]), label = 'Gauss_Newton', linewidth = 5)
    
    plt.plot(x_values,test(x_values,param[0],param[1]), label = 'Scipy Curve Fit', color = 'red', linewidth = 2)
    plt.title('Non-Linear Regression')
    plt.xlabel('X Axis')
    plt.ylabel('Y Axis')
    plt.legend()
    plt.show()
    
plot(x,y)
    
    
    
    
        
    
        
    
        
        
        
        
    
    
    
    


    
    