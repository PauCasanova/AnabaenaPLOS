import pandas
import numpy as np
import matplotlib.pyplot as plt
import math
from matplotlib.ticker import MaxNLocator


plt.rcParams["figure.figsize"] = [9, 6]
plt.rcParams["figure.autolayout"] = True
plt.style.use('seaborn-dark-palette')
plt.rcParams['font.size'] = 15
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['font.weight'] ='bold'
plt.rcParams['xtick.major.size'] = 10
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['xtick.minor.size'] = 5
plt.rcParams['xtick.minor.width'] = 1
plt.rcParams['ytick.major.size'] = 10
plt.rcParams['ytick.major.width'] = 2
plt.rcParams['ytick.minor.size'] = 5
plt.rcParams['ytick.minor.width'] = 1

F=1.2578E-4
T=81.132
Q=162.03
K=7.36

x = np.linspace(0, 150, 100)
a=110
i=15.8
g=0
f=(F*x**2*(1+a/T))/(1+(F*x**2*(1+a/T))+(i/K)**2+(g/Q)**2)
plt.plot(x, f, color='orange', linewidth=3, label='First round (2h)')
a=20
i=40
g=120
f=(F*x**2*(1+a/T))/(1+(F*x**2*(1+a/T))+(i/K)**2+(g/Q)**2)
plt.plot(x, f, color='green', linewidth=3, label='Near to an heterocyst (50h)')
a=356
i=12
g=80
f=(F*x**2*(1+a/T))/(1+(F*x**2*(1+a/T))+(i/K)**2+(g/Q)**2)
plt.plot(x, f, color='red', linewidth=3, label='Far from an heterocyst (50h)')

plt.yticks(np.arange(0,0.9,0.2))
plt.ylim((0,1))
plt.xlim((0,150))
plt.ylabel(r'$g(R_j,A_j,I_j,G_i)$', fontweight='bold')
plt.xlabel('HetR concentration (nM)', fontweight='bold')
plt.title('wild type HetR dependent catalysis', y=1.0, pad=-22, fontweight='bold')
plt.legend(frameon=False,bbox_to_anchor=(0.55, 0.88),title='Conditions:')._legend_box.align = "left"
plt.tight_layout()
plt.savefig('Hill.png', transparent=True)
plt.close()
