import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D

mu = 398600 #[km^3/s^2]
I_vec = np.array([1.0, 0.0, 0.0])
J_vec = np.array([0.0, 1.0, 0.0])
K_vec = np.array([0.0, 0.0, 1.0])

def mag(a):
    return np.linalg.norm(a)

def state_to_orbital_elements(pos, vel):
  h = np.cross(pos,vel)
  h_mag = mag(h)
  i = np.arccos(h[2]/h_mag)
  e_vec = np.cross(vel, h) / mu - pos / mag(pos)
  e = mag(e_vec)

  N = np.cross(K_vec, h)

  small_omega = np.arctan2(np.dot(h, np.cross(N, e_vec)), h_mag*np.dot(e_vec, N))
  big_omega = np.arctan2( np.dot(N, J_vec), np.dot(N, I_vec) )
  true_anomaly = np.arctan2(np.dot(h, np.cross(e_vec, pos)), h_mag*np.dot(e_vec, pos) )


  specific_energy = 0.5 * mag(vel)**2 - (mu / mag(pos))
  a = -mu / (2*specific_energy)


  print("ORBITAL ELEMENTS FROM STATE:")
  print("Pos = ", pos,   "Vel = ", vel)
  print("=============================")
  print("a \t=", a)
  print("e \t=", e)
  print("i \t=", i)
  print("Omega \t=", big_omega)
  print("w \t=", small_omega)
  print("nu \t=", true_anomaly) #true anomaly
  print("h \t=", h, "\n") #specific angular momentum


  return a, e, i, big_omega, small_omega, true_anomaly

def orbital_elements_to_state(a, e, i, big_omega, small_omega, true_anom):
    N = np.array([np.cos(big_omega), np.sin(big_omega), 0.0])
    N = N / mag(N)


    h_norm = np.array([np.sin(i)*np.sin(big_omega), -np.sin(i)*np.cos(big_omega), np.cos(i)])
    N_t = np.cross(h_norm, N)
    r_mag = a*(1 - e**2) / (1 + e*np.cos(true_anom))

    radial_vec = np.cos(true_anom + small_omega)*N + np.sin(true_anom + small_omega) * N_t
    perp_vec = -np.sin(true_anom + small_omega)*N + np.cos(true_anom + small_omega) * N_t

    h_mag = np.sqrt(mu*a*(1 - e**2))

    r = r_mag*radial_vec
    v = mu/h_mag * ( (1+e*np.cos(true_anom))*perp_vec + e*np.sin(true_anom) * radial_vec )

    print("STATE VECTORS FROM ORBITAL ELEMENTS:")
    print("a = ", a, "  e = ", e, "  i = ", i)
    print("small omega = ", small_omega, "  big omega = ", big_omega)
    print("true anomaly = ", true_anom)
    print("====================================")
    print("posiiton = ", r)
    print("velocity = ", v)
    print("h \t = ", h_norm*h_mag, "\n")

def plot_3d_orbit():
    EARTH_RAD = 6378
    fig = plt.figure()
    ax = Axes3D(fig)

    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    x = EARTH_RAD*np.cos(u)*np.sin(v)
    y = EARTH_RAD*np.sin(u)*np.sin(v)
    z = EARTH_RAD*np.cos(v)
    ax.plot_wireframe(x, y, z, color="brown", label='Earth')


    plt.legend()

    plt.show()


position = np.array([-2465, 6040, 3413])    #[km]
velocity = np.array([1.727, 3.893, -5.883]) #[km/sec]
a, b, c, d, e, f = state_to_orbital_elements(position, velocity)
orbital_elements_to_state(a, b, c, d, e, f)

plot_3d_orbit()
