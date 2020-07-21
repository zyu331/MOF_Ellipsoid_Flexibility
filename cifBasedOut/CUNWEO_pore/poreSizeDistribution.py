#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 19:57:54 2020

@author: zyu
"""

import pandas as pd
import matplotlib.pyplot as plt
import os.path

dirname= os.path.dirname(os.path.realpath(__file__)).split('/')[-1]

rigid=pd.read_csv('Rigid/0.psd_histo',header=9,sep=' ')
rigid=rigid[rigid['Bin']<10]
flexCount=pd.read_csv('elipsoid/0.psd_histo',header=9,sep=' ')
flexCount=flexCount[flexCount['Bin']<10]['Count'].values

for i in range(1,10):
    data=pd.read_csv('elipsoid/'+str(i)+'.psd_histo',header=9,sep=' ')
    
    data=data[data['Bin']<10]
    flexCount=flexCount+data['Count'].values


flexCount2=pd.read_csv('Gaussian/'+'0.psd_histo',header=9,sep=' ')
flexCount2=flexCount2[flexCount2['Bin']<10]['Count'].values

for i in range(1,10):
    data=pd.read_csv('Gaussian/'+str(i)+'.psd_histo',header=9,sep=' ')
    
    data=data[data['Bin']<10]
    flexCount2=flexCount2+data['Count'].values


flexCount3=pd.read_csv('MD/'+'0.psd_histo',header=9,sep=' ')
flexCount3=flexCount3[flexCount3['Bin']<10]['Count'].values

for i in range(1,10):
    data=pd.read_csv('MD/'+str(i)+'.psd_histo',header=9,sep=' ')
    
    data=data[data['Bin']<10]
    flexCount3=flexCount3+data['Count'].values

fig,ax=plt.subplots()
plt.autoscale(enable=True, axis='x',)

ax.plot(rigid['Bin'].values,rigid['Count'].values,label='Rigid_MD')
ax.plot(rigid['Bin'].values,flexCount/10,label='eli')
ax.plot(rigid['Bin'].values,flexCount2/10,label='Gau')
ax.plot(rigid['Bin'].values,flexCount3/10,label='MD')

plt.xlabel('pore Size (A)')
plt.ylabel('count')
plt.legend()
plt.title(dirname)
plt.savefig(dirname,dpi=500)


