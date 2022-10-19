import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint



a,b,c,d = 1,3,1,5
I = 0.0

theta = 0.14
omega = 0.112
gamma = 1.5
epsilon = 0.01;
xmin = -0.5
xmax = 1.5
ymin = -0.15
ymax = 0.3;

_gamma = 3/(1 -  theta + theta**2)

def f(x,y):
    return y - a*x**3 + b*x**2 + I

def g(x,y):
    return c - d*x**2 - y

def dx_dt(u, t):
    x = u[0]
    y = u[1]
    return [f(x,y), g(x,y)]


fig, ax=plt.subplots(1)
    
#vector fields
X, Y = np.mgrid[xmin:xmax:20j,  ymin:ymax:20j ]
dx, dy =  f(X, Y), g(X, Y)
ax.quiver(X, Y, dx, dy, color = 'b', alpha=0.5)

# Trajectories in forward time.
ts = np.linspace(0, 100, 1000)
ic = [ (0.5, 0.09), (-0.5, 0.3), (0.3, 1.5), (0.5, 0.2), (-0.4, 0.08), (-0.25, -0.1) ]
for x0_y0 in ic :
    xs=odeint(dx_dt, x0_y0, ts)
    ax.plot(xs[:,0], xs[:,1], "r-")


# Label the axes and set fontsizes.
ax.set_xlabel('u', fontsize=15)
ax.set_ylabel('v', fontsize=15)
#ax.set_tick_params(labelsize=15)
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax);
ax.grid(lw=0.3)

# Plot the nullclines.
x=np.arange(xmin, xmax, 0.01)
ax.plot(x, x/gamma, 'g--', x, -x * (x - theta) * (x - 1) + omega, 'g--', alpha=0.5)

plt.show()


#-------- Plot  u=u(t) , w = w(t) --------
fig, ax=plt.subplots(2, 1)

ts = np.linspace(0, 500, 10000)
x0_y0 = (0.5, 0.09)
xs=odeint(dx_dt, x0_y0, ts)
ax[0].plot(ts, xs[:,0], "r-")
ax[1].plot(ts, xs[:,1], "r-")


ax[1].set_xlabel('time', fontsize=15)
ax[0].set_ylabel('u', fontsize=15)
ax[1].set_ylabel('w', fontsize=15)
ax[0].grid(lw=0.4)
ax[1].grid(lw=0.4)

plt.show()
     
####--------- Get the roots numerically
from scipy.optimize import root

def  fun_root (u) :
    '''return the objetive function along with the value of the jacobian''' 
    x = u[0]
    y = u[1]
    f = y - a*x**3 + b*x**2 + I
    g = c - d*x**2 - y
    eqn = np.array ( [f, g] )
    jacobiano = np.array([ [-x *(x-theta) - x * (x-1) - (x-theta) * (x-1), -1], \
                    [epsilon, - epsilon * gamma]] )
    return eqn, jacobiano

# Using the "symbolic" critical points as a seed 
c_points = [(0, 0)]
#c_points =[(0,0), (0.2, 0.02), (0.8, 0.08)]    
for cp in c_points : 
    sol = root(fun_root, cp, method='hybr', jac = True)
    c_point = sol.x
    print(f'Punto de equilibrio: {c_point}' ) 
    jacob = fun_root(c_point)[1] 
    # Get the eigenvalues
    eigen = np.linalg.eig (jacob)
    print(eigen)
 
    