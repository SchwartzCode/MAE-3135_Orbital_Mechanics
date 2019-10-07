from matplotlib import pyplot as plt
import numpy as np
#from sympy import *
from sympy.solvers import solve
from sympy import Symbol

SPEC_ANGULAR_MOMENTUM = 1.637e10
TRUE_ANOMALY = 0.616371
MU = 4.26213e13
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

    total_anomaly = 2*np.pi * i / N
    r_sat = SPEC_ANGULAR_MOMENTUM**2 / ((MU) * (1 + TRUE_ANOMALY * np.cos(total_anomaly)))
    sat_state[i,0] = r_sat*np.cos(total_anomaly) / 1000 #divide by 1000 to convert from [m] to [km]
    sat_state[i,1] = r_sat*np.sin(total_anomaly) / 1000


plt.figure(figsize=[7,7])
plt.plot(mars_state[:,0], mars_state[:,1])
plt.plot(sat_state[:,0], sat_state[:,1])
plt.xlim(-20000,5000)
plt.ylim(-12500,12500)
plt.title("Satellite Orbit Around Mars")
plt.xlabel("X [km]")
plt.ylabel("Y [km]")
plt.show()
