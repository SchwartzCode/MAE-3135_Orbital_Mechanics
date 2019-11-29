import numpy as np, math as m
np.set_printoptions(suppress=True, formatter={'float_kind':'{:0.3f}'.format})
arr = np.array
mag = np.linalg.norm
def pr(*args):
    if True:
        print(str(arr(args)).replace(' ', ', ').replace('[', '').replace(']', ''))


muU = 1.327e11
muE = 3.986e5
RE = 1.496e8
rE = 6378
muJ = 1.267e8
RJ = 7.786e8
rJ = 71490
muS = 3.793e7
RS = 1.433e9
rS = 60270

tt = 494.5 * 24 * 3600
e = .8431
PE = 2 * m.pi / m.sqrt(muU) * RE ** 1.5
PJ = 2 * m.pi / m.sqrt(muU) * RJ ** 1.5
nE = 2 * m.pi / PE
nJ = 2 * m.pi / PJ
a = RE / (1 - e)
T = 2 * m.pi / m.sqrt(muU) * a ** 1.5
Me = tt * 2 * m.pi / T


def E0(e, Me):
    if Me > np.pi:
        E = Me - e / 2
    else:
        E = Me + e / 2

    precision = 1e-5
    ratio = 1
    while abs(ratio) > precision:
        ratio = (E - e * np.sin(E) - Me) / (1 - e * np.cos(E))
        E = E - ratio
    return E

E = E0(e, Me)
nu = m.atan(m.tan(E / 2) * m.sqrt((1 + e) / (1 - e))) * 2
phi = nu - nJ * tt
#pr(m.degrees(phi))

rp = rE + 1000
ve = m.sqrt(muU / RE)
vo = m.sqrt(muU * (2 / RE - 1 / a))
vinf = vo - ve
vc = m.sqrt(muE / rp)
#print(m.sqrt(2 + (vinf / vc) **2))
dv = vc * (m.sqrt(2 + (vinf / vc) ** 2) - 1)
#pr(dv)

#part d

a2 = 1.079e9
e1 = 0.8431
e2 = 0.8379

aj = RE / (1 - e1)
#print("a_j:", aj)
vp = m.sqrt(muU * (2 / RE - 1 / aj))
h1 = RE * vp

nuJ = tt * m.sqrt(muU / RJ ** 3) + phi


vscr = muU / h1 * e1 * m.sin(nuJ)
vscp = muU / h1 * (1 + e1 * m.cos(nuJ))
vsc1 = mag(arr([vscr, vscp]))

vj = m.sqrt(muU / RJ)
vinf = m.sqrt(vscr ** 2 + (vscp - vj) **2)

h2 = m.sqrt(muU * a2 * (1 - e2 ** 2))

vsc2 = m.sqrt(2 * muU / RJ - muU / a2)
vscp2 = h2 / RJ
vscr2 = m.sqrt(vsc2 ** 2 - vscp2 ** 2)

phi1 = m.atan2(vscr, vscp - vj)
phi2 = m.acos((vscp2 - m.sqrt(muU / RJ)) / vinf)

delt = phi1 - phi2
beta = (m.pi - delt) / 2
eh = 1 / m.cos(beta)
rph = muJ * (eh - 1) / vinf ** 2

print(rph)

D = rph * m.sqrt(1 + 2*muJ / (rph * vinf ** 2))
pr(D)
