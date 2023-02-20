import numpy as np
import matplotlib.pyplot as plt

r, q, x, y, phi = np.loadtxt('out.txt', unpack=True)
with open('out.txt') as fl:
    comment = fl.readline()
words = comment.split()
nr = int(words[3])
ntheta = int(words[6])
x = x.reshape(nr,ntheta)
y = y.reshape(nr,ntheta)
phi = phi.reshape(nr,ntheta)

x = np.concatenate((x, x[:, 0, np.newaxis]), axis=1)
y = np.concatenate((y, y[:, 0, np.newaxis]), axis=1)
phi = np.concatenate((phi, phi[:, 0, np.newaxis]), axis=1)

fig, ax = plt.subplots(1,1)
clevels = np.linspace( phi.min(), phi.max(), 101)
im = ax.contourf( x, y, phi, clevels, cmap='seismic' )
for c in im.collections:
    c.set_edgecolor('face')
plt.colorbar( im )

n1 = nr//8
ax.plot( x[:,::n1], y[:,::n1], '-', color='lightgrey', lw=0.5 )
ax.plot( x[:,nr-1], y[:,nr-1], '-', color='lightgrey', lw=0.5 )
ax.plot( x.transpose()[:,::n1], y.transpose()[:,::n1], '-', color='lightgrey', lw=0.5 )
# style
ax.set_xlabel( r'$x$' )
ax.set_ylabel( r'$y$' )
#ax.set_title( r'Numerical solution: $\phi(x1,x2)$' )
ax.set_aspect( 'equal' )
fig.tight_layout()

plt.show()

