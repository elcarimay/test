# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    
    n = 50
    x = np.linspace(0,1,n)
    f = 1 
    y = np.sin(2*np.pi*f*x)
    
    plt.plot(x,y)