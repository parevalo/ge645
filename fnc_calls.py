import numpy as np
import matplotlib.pyplot as plt
from RTE_fnc import *

#  Variable declarations

ng = 8
thetaprime = np.radians(150.0)
phiprime = np.radians(0.0)
muprime = np.cos(thetaprime)
distnames = ('Planophile', 'Erectophile', 'Plagiphile', 'Extremophile', 'Uniform', 'Spherical')
dist = 0
nl = 100
epsilon = 0.00001
Ftot = 1.0
fdir = 0.7
Io = Ftot * (fdir / (abs(muprime)))
Id = Ftot * (1 - fdir) / np.pi
LAI = 5.0
dl = LAI/nl
rhold = 0.475
tauld = 0.45
R_s = 0.2

#  Get crosssections
xg, wg, gdir_out, gdif_out, gmdir_out, gmdif_out = xsections(ng, 1.0, dist, muprime, phiprime, rhold, tauld)

#  Calculate collided and uncollided up and down
Io_ucd = io_uncol_down(dl, nl, gdir_out, Io, muprime)
Fo_ucdsoil, Io_ucu = io_uncol_up(dl, nl, gdir_out, gdif_out, Io, muprime, ng, xg, R_s)
Id_ucd = id_uncol_down(dl, nl, gdif_out, Id, ng, xg)
Fd_ucdsoil, Id_ucu = id_uncol_up(dl, nl, gdif_out, Id, ng, xg, wg, R_s)
Q = fcs(nl, ng, wg, gmdir_out, gmdif_out, Io_ucd, Io_ucu, Id_ucd, Id_ucu)

#  Iterate over multiple collision source
S = np.zeros((nl, ng, ng))
s_bot = np.zeros((nl, ng, ng)) # Added to create Figure 1
s_top = np.zeros((nl, ng, ng)) # Added to create Figure 1
ic = np.zeros((nl+1, ng, ng))
for ims in range(100):
    print ims
    #S[:, :, :] += Q[:, :, :] # Vectorized version works fine too
    for i in range(ng):
        for j in range(ng):
            for k in range(nl):
                S[k, j, i] += Q[k, j, i]


    ic1 = sweep_down(nl, ng, xg, wg, gdif_out, dl, S, R_s, ic)
    cv, ic = sweep_up(nl, ng, xg, gdif_out, dl, S, ic1, epsilon)

    if cv:
        break

    S = multicoll_s(nl, ng, wg, gmdif_out, ic)
    s_bot[ims, :, :] = S[99, :, :]
    s_top[ims, :, :] = S[0, :, :]

#  Print theta_v, phi_v and RF values

RF = np.zeros((ng, ng/2)) # Reflectance factor, BRF if fdir=1, else HDRF)
for i in range(ng/2, ng):
    theta_v = np.degrees(xg[i])
    for j in range(ng):
        phi_v = np.degrees(xg[j] * np.pi + np.pi)
        RF = ic[0, j, i] * np.pi
        print theta_v, phi_v, RF

energy_bal(nl, ng, xg, wg, dl, R_s, rhold, tauld, gdir_out, gdif_out, Fo_ucdsoil,
           Fd_ucdsoil, Io_ucd, Id_ucd, Io_ucu, Id_ucu, ic)


# Put in plots script or something
def sconver(s_arr):
    convs = np.zeros(30)
    for i in range(30):
        mnsum = 0.0
        for j in range(8):
            mn = np.mean(s_arr[i,j,:])
            mnsum += mn

        convs[i] = mnsum
    return convs

cs_bot = sconver(s_bot)
cs_top = sconver(s_top)

plt.plot(cs_bot)
plt.plot(cs_top)
plt.ylim(0, 1)