# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 09:43:22 2020

Polynomial regression

@author: Siddharth
"""

import numpy as np
import matplotlib as plt

def cholesky(A):
    
    l = np.zeroes(A.shape)
    
    for k in range(0, l.shape[0]):
        
        for i in range(0, l.shape[1]):
            
            if (k==i):
                
                l[k][k] = (A[k][k] - np.sum(np.power(l,2), axis = 1))[k]
                continue
            
            l[k][i] = A[k][i] - np.sum(np.multiply(l[i],l[k]))
            
    return l


            

def pol_regression(m,n, x, y):

    if n < (m+1):
        
        print("Not possible reduce m")
        return 0
    
    A = np.zeros((m+1, m+1))
    B = np.zeros((m+1, 1))
    coeff = np.zeros((m+1, 1))
    
    for i in range(0,m+1):
        
        for j in range(0, m+1):
            
            A[i][j] = np.sum(np.power(x,i+j))
            
        B[i] = np.sum(np.multiply((np.power(x,i)),y))
        
    cholesky(A,B,coeff)
            
            