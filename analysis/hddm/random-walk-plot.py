#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  5 10:00:41 2022

@author: saeedeh
"""

import numpy as np
import matplotlib.pyplot as plt
def Randomwalk1D(n): #n here is the no. of steps that we require
   x = 0
   y = 0
   xposition = [0] #starting from origin (0,0)
   yposition = [0] 
   for i in range (1,n+1):
       step = np.random.uniform(0,1)
       if step < 0.50: # if step is less than 0.5 we move up    
           x += 1
           y += 1.1  #moving up in u direction
       if step > 0.50: # if step is greater than 0.5 we move down  
           x += 1
           y += -1 #moving down in y direction
 
       xposition.append(x)
       yposition.append(y)
   return [xposition,yposition]
Randwalk = Randomwalk1D(1000) #creating an object for the Randomwalk1D class and passing value of n as 100
plt.plot(Randwalk[0],Randwalk[1],'k-', label = "Randwalk1D") # 'r-' makes the color of the path red
plt.title("1-D Random Walks")
plt.show()