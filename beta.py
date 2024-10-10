#!/usr/bin/env python 

import numpy as np
import matplotlib
# matplotlib.use('module://matplotlib-backend-kitty')
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from scipy import stats

# benzene, naftalene, antracene, fenantrene tetracene pirene
x = np.array([2, 3.6832, 5.3137, 5.4483]) # aqui las energias de deslocalizacion calculadas
y = np.array([37.0, 75.0, 105.0, 110.0])
sizey = 5
sizex = sizey * 1.618
fig, ax = plt.subplots(figsize=(sizex, sizey) )
ax.plot(x, y, linestyle = "None" , marker= "o" , markersize =8)
res = stats.linregress(x , y) # Regresion lineal
xx = np.array([min(x), max(x)])
beta = "$\\beta$ = %.1f kcal/mol " % res.slope
print(res.intercept)
ax.plot (xx, res.intercept + res.slope * xx, label = beta, color = "tab:blue", lw =2)

ax.set_xlabel("Theoretical DE ($\\beta$)" , fontsize =18)
ax.set_ylabel("Experimental DE (kcal/mol)" , fontsize =18)
ax.legend( fontsize =14 , frameon = False )
plt.savefig("report/beta1_plt.png")
# plt.show()
