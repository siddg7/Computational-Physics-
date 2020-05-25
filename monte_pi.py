# -*- coding: utf-8 -*-
"""
Created on Fri May 22 00:10:58 2020

@author: Siddharth
"""
import numpy as np
import random 


def linear_congruential(x0,a,c,m):
    x = np.zeros([m,1])
    x[0] = x0
    
    for i in range(0,m-1):
        #print(i)
        x[i+1] = (a*x[i] + c)%m
        
    r = np.array([])
    
    for i in range(1,m):
        r = np.append(r,[x[i]/m])
        
    return np.unique(r)

    
def plot(r1,r2):
    
    x = np.array([])
    y = np.array([])
    
    for i in r1:
        for j in r2:
            
            x = np.append(x,[i])
            y = np.append(y,[j])
            
    plt.scatter(x,y)

def calc_pi(r1,r2):
    x = r1
    y = r2

    n0 = 0
    n = 0
    
    for i in x:
        for j in y:
            n += 1
            
            if(((i-0.5)**2 + (j-0.5)**2) < 0.25):
                n0 += 1
    print(n0)
    print(n)
                
    return (4*n0/n)



            
    
            
            
    
                
                
            
            
            
    
    

        