import numpy as np
import matplotlib.pyplot as plt
import os


fig,ax = plt.subplots(1, 1, figsize=(8, 6))

k = 2
dn = os.path.dirname(os.path.realpath(__file__))

KE_array = np.loadtxt(dn+"/KEs_array.txt")
dt = 0.0125
n=8000
time = np.arange(0,n*dt,dt)

freq = 1.
nf = 6
for i in range(nf):
    ax.plot(time/2/np.pi/freq, KE_array[:,i], label='{:d}'.format(i+1))
ax.legend(loc='upper right', fontsize=12,ncol=nf)
ax.set_yscale('log')
fig.savefig("test.png", dpi=300, bbox_inches='tight')
#fig.savefig('{:s}/folder_{:d}/KE_plot.png'.format(dn,k), dpi=300, bbox_inches='tight')

