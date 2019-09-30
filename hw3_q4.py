from matplotlib import pyplot as plt
import numpy as np
from sympy import *
from sympy.solvers import solve
from sympy import Symbol


E = Symbol('E')
TRUE_ANOMALY = 0.616371
N = 1000
R_MARS = 3390 #[km]
P = 8.6322 #period [hours]
mars_state = np.zeros([N, 2]) #list of position states of Mars' surface in form (x, y) [km]
sat_state = np.zeros([N, 2]) #list of position states of satellite in form (x, y) [km]

sat0 = np.array([3890., 0.]) #initial position of satellite; when t=0
sat_state[0] = sat0 #set first row in satellite state to initial state

for i in range(N):
    mars_state[i, 0] = R_MARS*np.cos(2*np.pi*i/N)
    mars_state[i, 1] = R_MARS*np.sin(2*np.pi*i/N)

    time = P * i / N
    Me = 2 * np.pi * time / P
    eccentric_anomaly = solve(E - TRUE_ANOMALY*sin(E) - Me, E)


plt.figure(figsize=[7,7])
plt.plot(mars_state[:,0], mars_state[:,1])
plt.xlim(-4000,4000)
plt.ylim(-4000,4000)
#plt.show()
