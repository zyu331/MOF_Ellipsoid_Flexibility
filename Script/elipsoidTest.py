#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 22:07:05 2020

@author: zyu
"""

import numpy as np


def ellipseGen(arrayIn,factor): 
    
    arrayT=[[arrayIn[0],arrayIn[3],arrayIn[4]],[arrayIn[3],arrayIn[1],arrayIn[5]],[arrayIn[4],arrayIn[5],arrayIn[2]]]
    array=np.linalg.inv(arrayT)
    flag=False
    
    while flag==False:
        x=float(2*np.sqrt(arrayIn[0])*(np.random.rand(1,1)-0.5))
        y=float(2*np.sqrt(arrayIn[0])*(np.random.rand(1,1)-0.5))
        z=float(2*np.sqrt(arrayIn[0])*(np.random.rand(1,1)-0.5))
    
        v=np.array([x,y,z])
        if np.matmul(np.matmul(v.T,array),v) < factor:
            flag=True
            cc=np.matmul(np.matmul(v.T,array),v)

    return v

a=np.zeros((100,3))

for i in range(0,100):
    a[i,:]=ellipseGen(np.array([0.05**2,0.02**2,0.01**2,0,0,0]),1)
    
    
    np.std(a)

cov=np.array(([0.05**2,0,0],[0,0.02**2,0],[0,0,0.01**2]))

a_cov=np.zeros((100,3))

for i in range(0,100):
    a_cov[i,:]=np.random.multivariate_normal([0,0,0],cov)