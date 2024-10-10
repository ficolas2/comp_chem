#!/usr/bin/env python 

import numpy as np
import matplotlib

# matplotlib.use('module://matplotlib-backend-kitty')
# matplotlib.use('TkAgg')
import matplotlib . pyplot as plt
from scipy import stats

#            benzene, naftalene, antracene, fenantrene, tetracene, pirene
x = np.array([2,      1.236,     0.8284,    1.2104,     0.59,      0.89])
y = np.array([50000,  36360,     26700,     34190,      21230,     29900])

sizey = 5
sizex = sizey * 1.618
fig, ax = plt.subplots(figsize=(sizex, sizey) )
ax.plot(x, y, linestyle = "None" , marker= "o" , markersize =8)
res = stats.linregress(x , y) # Regresion lineal
xx = np.array([min(x), max(x)])
beta = "$\\beta$ = %.1f $cm^{-1}$ " % res.slope
ax.plot (xx, res.intercept + res.slope * xx, label = beta, color = "tab:blue", lw =2)
print(res.intercept)

ax.set_xlabel("Theoretical $\\Delta E_{HOMO-LUMO}$ ($\\beta$)" , fontsize =18)
ax.set_ylabel("Experimental  $\\Delta E_{HOMO-LUMO}$ ($cm^{-1}$)" , fontsize =18)
ax.legend( fontsize =14 , frameon = False )
plt.savefig("report/beta2_plt.png")
# plt.show()
