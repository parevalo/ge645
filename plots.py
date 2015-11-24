import numpy as np
from RTE_fnc import *

#  PROVISIONAL CODE TO MAKE THE PLOTS

#  Variable declarations
ng = 8
thetaprime = np.radians(150.0)
phiprime = np.radians(0.0)
muprime = np.cos(thetaprime)
distnames = ('Planophile', 'Erectophile', 'Plagiphile', 'Extremophile', 'Uniform', 'Spherical')
dist = 0
nl = 100
epsilon = 0.0001
Ftot = 1.0
fdir = 0.7
Io = Ftot * (fdir / (abs(muprime)))
Id = Ftot * (1 - fdir) / np.pi
LAI = 3.0
dl = LAI/nl
rhold = 0.45
tauld = 0.45
R_s = 0.3

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
ic = np.zeros((nl+1, ng, ng))
for ims in range(100):
    print ims
    for i in range(ng):
        for j in range(ng):
            for k in range(nl):
                S[k, j, i] += Q[k, j, i]

    ic1 = sweep_down(nl, ng, xg, wg, gdif_out, dl, S, R_s, ic)
    cv, ic = sweep_up(nl, ng, xg, gdif_out, dl, S, ic1, epsilon)

    if cv:
        break

    S = multicoll_s(nl, ng, wg, gmdif_out, ic)

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


# PLOTS, convert into function or something
# Figure 3
# for i in range(6):
#     plt.plot(np.arange(0, 90, 8), pdfs[i, :])
#
# plt.xlabel('Leaf inclination in degrees')
# plt.ylabel('Leaf normal inclination distribution function')
#
# # Figure 8
# plt.figure()
# gfo = np.zeros((6, 90))
# for i in range(6):
#     for j in range(0, 90):
#         mp = np.cos(np.radians(j))
#         gfo[i, j] = g_dir(n, o, w, 1.0, pdfs[i, :], mp, phiprime)
#
#     plt.plot(np.arange(0,90), gfo[i, :])
#
# plt.xlabel('Projection Polar Angle in Degrees')
# plt.ylabel('Geometry function')
# plt.legend(distnames)
#
# # Figure 9
# gmfo = np.zeros((19, 6))
# t = np.arange(0.0, 0.6, 0.1)
# p = np.arange(1.0, 0.4, -0.1)
# bb = np.linspace(0, 180, 19)
# cosbb = np.cos(np.radians(bb))
#
# plt.figure()
# for i in range(6):
#     for j in range(0, 190, 10):
#         muprime = np.cos(np.radians(j))
#         gmfo[j/10, i] = gamma_d(n, o, w, 1.0, pdfs[5, :], p[i], t[i], muprime, phiprime, 1, 0.0)
#
#     plt.plot(cosbb, gmfo[:, i])
#
# plt.ylim(0, 0.4)
# plt.grid()
# plt.xlabel('Cosine of Scattering Angle')
# plt.ylabel('Area scattering phase function')
# lab = t/(t+p)
# plt.legend(lab, title="$\\tau \div \omega$ =")
#
# # Figure 11
# def sph2cart(r, theta, phi):
#     """
#     Function to convert between spherical to cartesian coordinates
#     :param r: Length of the vector
#     :param theta: Angle in xy plane, in degrees
#     :param phi: Angle with respect to Z axis, in degrees
#     :return: Position of the point in cartesian coordinates
#     """
#     xyz = np.zeros(3)
#     xyz[0] = r * np.sin(np.radians(phi)) * np.cos(np.radians(theta))
#     xyz[1] = r * np.sin(np.radians(phi)) * np.sin(np.radians(theta))
#     xyz[2] = r * np.cos(np.radians(phi))
#
#     return xyz
#
# thetaprime = np.radians(170.0)
# phiprime = np.radians(0.0)
# muprime = np.cos(thetaprime)
#
# plt.figure()
# omegaprime = sph2cart(1, 0, 170)
# mu_size = 180 / 5 + 1
# phi_size = 360 / 5 + 1
# totsize = mu_size * phi_size
# gmd = np.zeros((mu_size, phi_size))
# ords = np.zeros((mu_size, phi_size))
# for i in range(0, 185, 5):
#     mu = np.cos(np.radians(i))
#     for j in range(0, 365, 5):
#         phi = np.radians(j)
#         gmd[i/5, j/5] = gamma_d(n, o, w, 1.0, pdfs[0, :], rhold, tauld,  muprime, phiprime, mu, phi)
#         omega = sph2cart(1.0, j, i)
#         ords[i/5, j/5] = np.dot(omega, omegaprime)
#
# P = (4 * gmd)/gdir_out
# plt.scatter(ords, P)
#
# plt.ylim(0, 2)
# plt.xlim(-1, 1)
# plt.xticks(np.linspace(1.0, -1.0, 11))
# plt.yticks(np.linspace(0, 2, 11))
# plt.grid()
# plt.xlabel('Cosine of Scattering Angle')
# plt.ylabel('Normalized scattering phase function')
#
