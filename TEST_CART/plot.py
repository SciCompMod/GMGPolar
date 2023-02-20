import numpy as np
import matplotlib.pyplot as plt
import rho_triangularity_logical as mod
import article_setup

r, q, x, y, f = np.loadtxt('sol_out.txt', unpack=True)
nq = 1
while r[nq]==r[nq-1]:
    nq += 1
nr = len(q)//nq
r = r.reshape((nr,nq))
q = q.reshape((nr,nq))
x = x.reshape((nr,nq))
y = y.reshape((nr,nq))
f = f.reshape((nr,nq))

x = np.concatenate([x, x[:, :1]], axis=1)
y = np.concatenate([y, y[:, :1]], axis=1)
r = np.concatenate([r, r[:, :1]], axis=1)
q = np.concatenate([q, q[:, :1]], axis=1)
f = np.concatenate([f, f[:, :1]], axis=1)

exact_func = np.vectorize(mod.phi_exact)
exact = exact_func(r, q, 0.3, 1.4)

err = f - exact

def plot_fig(val, title):
    fig, ax = plt.subplots(1,1)

    crop_val = max(-val.min(), val.max())
    clevels = np.linspace(-crop_val, crop_val, 101)
    im = ax.contourf( x[:,:], y[:,:], val[:,:], clevels, cmap='seismic' )
    for c in im.collections:
        c.set_edgecolor('face')
    ax.set_xlabel('x')
    ax.set_ylabel('y', rotation=0)
    cbar = plt.colorbar( im )
    cbar.formatter.set_powerlimits((0,0))
    n1 = 8
    nr = 8
    ax.plot( x[:,::nr], y[:,::nr], '-', color='lightgrey', lw=0.5 )
    ax.plot( x.transpose()[:,::nr], y.transpose()[:,::nr], '-', color='lightgrey', lw=0.5 )
    ax.set_title(title)
    plt.tight_layout()

plot_fig(err, title='Error')
plt.savefig('gmgpolar_czarny_cart_error.png')

#plot_fig(f, title = 'Calculated solution')
#plot_fig(exact, title = 'Exact solution')

plt.show()
