import numpy as np
import collections
#import matplotlib.pyplot as plt


def gauss_quad(ng):
    '''
    Function for obtaining gauss quadrature of order ng

    inputs:
    ng : quadrature order (integer). Must be even number between 4 and 12.

    outputs:
    xg: ordinates (real)
    wg: weights (real)
    '''

    xg = np.zeros(ng)
    wg = np.zeros(ng)

    xx = (-0.861136312,-0.339981044,-0.9324695,-0.6612094,-0.2386192,-0.960289856,-0.796666477,-0.525532410,
         -0.183434642, -0.973906529,-0.865063367,-0.679409568,-0.433395394, -0.148874339, -0.981560634,
          -0.904117256,-0.769902674,-0.587317954, -0.367831499,-0.125233409)

    ww = (0.347854845, 0.652145155,0.1713245, 0.3607616, 0.4679139, 0.101228536, 0.222381034, 0.313706646,
         0.362683783, 0.066671344,0.149451349,0.219086363,0.269266719, 0.295524225,
         0.047175336,0.106939326,0.160078329,0.203167427, 0.233492537,0.249147046)

    shift = (0, 0, 2, 5, 9, 14)

    ng2 = (ng/2)

    for i in range(ng2):
        xg[i] = xx[i + shift[ng2-1]]
        wg[i] = ww[i + shift[ng2-1]]

    for i in range(ng2, ng):
        xg[i] = -xg[ng - 1 - i]
        wg[i] = wg[ng - 1  - i]

    return xg, wg


def check_quad(ng, xg, wg):
    '''
    Function to check if the quadrature weights sum to 2
    Checks if integral from 0 to 1 dx x is equal to 0.5
    '''

    ws = np.sum(wg)
    s = sum(wg[ng/2:ng] * xg[ng/2:ng])

    str = 'The weights sum to {0} and the integral equals {1}'.format(ws, s)

    return str


def convfactors(upperlim, lowerlim):
    """
    Return the conversion factors calculated from the upper and lower limits
    :param upperlim: Value of the desired upper limit
    :param lowerlim: Values of the desired lower limit
    :return: Two conversion factors
    """

    cf1 = (upperlim - lowerlim) / 2.0
    cf2 = (upperlim + lowerlim) / 2.0

    return cf1, cf2


def example_integral(ng, xg, wg):
    '''
    This function shows how to do an integral with lower and upper bounds A and B, respectively
    '''

    conv1, conv2 = convfactors(1.0, -1.0)

    sum = 0.00
    for i in range(ng):
        neword = conv1 * xg[i] + conv2
        sum = sum + abs(neword * wg[i])

    sum = sum * conv1

    return 'Integral has a value of {}'.format(sum)


def leaf_normal_pdf(ng, xg, wg):
    """
    Evaluates planophile, erectophile, plagiophile, extremophile, uniform and spherical pdf
    :param ng: Number of points in quadrature
    :param xg: Ordinates
    :param wg: Weights
    :return: Multidim. array with all the pdfs. Verifies only ONE of the PDFS (for now)
    """

    gL = np.zeros((6,ng))

    conv1, conv2 = convfactors(np.pi / 2.0, 0.0)
    sum = 0.00
    f = 2.0 / np.pi

    for i in range(ng):
        neword = conv1 * xg[i] + conv2
        cos2 = np.cos(2.0 * neword)
        cos4 = np.cos(4.0 * neword)
        gL[0, i] = f * (1.0 + cos2) # Planophile
        gL[1, i] = f * (1.0 - cos2) # Erectophile
        gL[2, i] = f * (1.0 - cos4) # Plagiophile
        gL[3, i] = f * (1.0 + cos4) # Extremophile
        gL[4, i] = f # Uniform
        gL[5, i] = np.sin(neword) # Spherical
        sum += (gL[0,i] * wg[i])

    sum *= conv1
    print 'The integral sums up to {}'.format(sum)

    return gL


def g_dir(ng, xg, wg, hL, gL, mu_p, phi_p):
    """
    Calculates the G function given a direction of photon travel OMEGA_prime and the leaf normal pdf
    :param ng: Number of points in quadrature
    :param xg: Ordinates
    :param wg: Weights
    :param hL: Leaf normal azimuth pdf(s)
    :param gL: Leaf inclination pdf(s)
    :param mu_p: mu prime
    :param phi_p: phi prime
    :return: Multidimensional array with the values of the G function for each leaf orientation angle and type
    """

    sin_tp = np.sqrt(1.0 - mu_p**2)

    # Define limits of integration and the conversion factors for integration over thetaL
    conv1_tl, conv2_tl = convfactors(np.pi/2.0, 0.0)

    # Limits and conv. factors for integration over phiL
    conv1_pl, conv2_pl = convfactors(np.pi*2.0, 0.0)

    # Integrate over theta_L
    sum_tl = 0.0
    for i in range(ng):
        neword_tl = conv1_tl * xg[i] + conv2_tl
        mu_tl = np.cos(neword_tl)
        sin_tl = np.sin(neword_tl)

        # Integrate over phi_L
        sum_pl = 0.00
        for j in range(ng):
            neword_pl = conv1_pl * xg[j] + conv2_pl
            dotproduct = abs(mu_tl * mu_p + sin_tl * sin_tp * np.cos(neword_pl - phi_p))
            sum_pl += wg[j] * hL / (2.0 * np.pi) * dotproduct  # hL is being taken as a single number here

        # Finish the phi_L integral
        sum_pl *= conv1_pl
        sum_tl += wg[i] * gL[i] * sum_pl

    # Finish the theta_L integral (the sum is Gdir)
    sum_tl *= conv1_tl

    return sum_tl


def g_dif(ng, xg, wg, hL, gL):
    """
    Evaluates the G function in all quadrature directions and checks for normalization
    :param ng: Number of points in quadrature
    :param xg: Ordinates
    :param wg: Weights
    :param hL: Leaf normal azimuth pdf(s)
    :param gL: Leaf inclination pdf(s)
    :return: 2D array with the values of G over different theta and phi angles (photon travel azimuth and zenith angles)
    """

    # Conversion factors to have ordinates simulate phiprime
    conv1_pp, conv2_pp = convfactors(np.pi*2.0, 0.0)

    # Get the gdif matrix direction by direction
    gf_out = np.zeros((ng, ng))
    for i in range(ng):
        mp = xg[i]
        for j in range(ng):
            pp = conv1_pp * xg[j] + conv2_pp
            gf_out[j, i] = g_dir(ng, xg, wg, hL, gL, mp, pp)

    # Check for normalization
    # (1/2PI) int_0^2PI dphi^prime int_0^1 dmu^prime G(OMEGA^prime) = 0.5
    sum_tp = 0.0
    for i in range(int(ng * 0.5), ng):
        sum_pp = 0.0
        for j in range(ng):
            sum_pp += wg[j] * gf_out[j, i]
        sum_pp *= conv1_pp
        sum_tp += wg[i] * sum_pp

    sum_tp /= (2.0 * np.pi)

    print "Gfunction check equals {}".format(sum_tp)

    return gf_out


def gamma_dir(ng, xg, wg, hL, gL, rho_ld, tau_ld, mp, pp):
    """
    Evaluates the GAMMA_d function (muprime, phiprime -> mu, phi) where (mu, phi) are the quadrature directions
     and checks for normalization
    :param ng: Number of points in quadrature
    :param xg: Ordinates
    :param wg: Weights
    :param hL: Leaf normal azimuth pdf(s)
    :param gL: Leaf inclination pdf(s)
    :param rho_ld: Hemispherical reflectance coef.
    :param tau_ld: Hemispherical transmittance coef.
    :param mp: Mu prime
    :param pp: Phi prime
    :return: 2D array gmdir_out EXPAND
    """

    #Conversion factors to simulate phi
    gmdir_out = np.zeros((ng, ng))
    conv1_p, conv2_p = convfactors(2.0 * np.pi, 0.0)

    # Get the gamma_d dir matrix direction by direction
    for i in range(ng):
        mu = xg[i]
        for j in range(ng):
            phi = conv1_p * xg[j] + conv2_p
            gmdir_out[j, i] = gamma_d(ng, xg, wg, hL, gL, rho_ld, tau_ld, mp, pp, mu, phi)

    return gmdir_out


def gamma_dif(ng, xg, wg, hL, gL, rho_ld, tau_ld, gdif):
    """
    Evaluates the gamma dif function for scattering from all quadrature directions to all exit quadrature directions
    :gdif: Output of the function g_dif
    :return: 4D array with the output of gamma dir for each of the
    """
    conv1_pp, conv2_pp = convfactors(2.0 * np.pi, 0.0)
    gmdif_out = np.zeros((ng, ng, ng, ng))
    # Get gamma_d_dif matrix direction by direction
    for i in range(ng):
        muprime = xg[i]
        for j in range(ng):
            phiprime = conv1_pp * xg[j] + conv2_pp
            gmdif_out[j, i, :, :] = gamma_dir(ng, xg, wg, hL, gL, rho_ld, tau_ld, muprime, phiprime)
            check_gammadir(ng, wg, rho_ld, tau_ld, gmdif_out[j, i, :, :], gdif[j, i])

    return gmdif_out


def gamma_d(ng, xg, wg, hL, gL, rho_ld, tau_ld, mp, pp, mu, phi):
    """
    Calculates the G function given a direction of the photon travel OMEGA PRIME and the leaf normal PDF
    :return: Value of the diffuse leaf scattering phase function? (Chapter 2, eq. 9)
    """

    sin_t = np.sqrt(1.0 - mu**2)
    sin_tp = np.sqrt(1.0 - mp**2)

    # Limits of integration over thetaL
    conv1_tl, conv2_tl = convfactors(np.pi/2.0, 0.0)
    # Limits of integration over phiL
    conv1_pl, conv2_pl = convfactors(np.pi*2.0, 0.0)

    # Integrate over thetaL
    sum_tl = 0.0
    for i in range(ng):
        neword_tl = conv1_tl * xg[i] + conv2_tl
        mu_tl = np.cos(neword_tl)
        sin_tl = np.sin(neword_tl)

        # Integrate over phiL
        sum_pl = 0.0
        for j in range(ng):
            neword_pl = conv1_pl * xg[j] + conv2_pl
            dotproduct1 = (mu_tl * mp + sin_tl * sin_tp * np.cos(neword_pl-pp))
            dotproduct2 = (mu_tl * mu + sin_tl * sin_t * np.cos(neword_pl-phi))

            if dotproduct1 * dotproduct2 <= 0.0:
                sum_pl += rho_ld * wg[j] * hL / (2.0 * np.pi) * abs(dotproduct1 * dotproduct2)
            else:
                sum_pl += tau_ld * wg[j] * hL / (2.0 * np.pi) * abs(dotproduct1 * dotproduct2)

        # Finish phiL integral
        sum_pl *= conv1_pl
        sum_tl += wg[i] * gL[i] * sum_pl

    # Finish thetaL integral
    sum_tl *= conv1_tl

    return sum_tl


def check_gammadir(ng, wg, rho_ld, tau_ld, gammadir, gdir):
    """
    Check the Gamma dir function for normalization. Original code requires xg but doesn't use
    :param ng:
    :param wg:
    :param rho_ld:
    :param tau_ld:
    :param gammadir:
    :param gdir:
    :return:
    """
    conv1_pp, conv2_pp = convfactors(2.0 * np.pi, 0.0)

    # Check for normalization
    sum_tp = 0.0
    for i in range(ng):
        sum_pp = 0.0
        for j in range(ng):
            sum_pp += wg[j] * gammadir[j, i]
        sum_pp *= conv1_pp
        sum_tp += wg[i] * sum_pp

    sum_tp /= np.pi
    sum_tp /= (gdir * (rho_ld + tau_ld))

    #print "Gamma_d check equals {}".format(sum_tp)


def xsections(ng, hL, dist, mp, pp, rho_ld, tau_ld):
    """
    Calculate all the crosssections and return arrays with their values
    :param ng: fill!
    :param hL: 
    :param gL: 
    :param dist: 
    :param mp: 
    :param pp: 
    :param rho_ld: 
    :param tau_ld: 
    :return: xg, wg, gdir_out, gdif_out, gmdir_out, gmdif_out
    """
    # Obtain Gauss quadrature
    xg, wg = gauss_quad(ng)

    # Check the quadrature
    check_quad(ng, xg, wg)

    # Get ALL pdfs of leaf normal orientation
    pdfs = leaf_normal_pdf(ng, xg, wg)

    # Get G functions
    gdir_out = g_dir(ng, xg, wg, 1.0, pdfs[dist, :], mp, pp)
    gdif_out = g_dif(ng, xg, wg, 1.0, pdfs[dist, :])

    # Get Gamma functions
    gmdir_out = gamma_dir(ng, xg, wg, 1.0, pdfs[dist, :], rho_ld, tau_ld, mp, pp)
    check_gammadir(ng, wg, rho_ld, tau_ld, gmdir_out, gdir_out)
    gmdif_out = gamma_dif(ng, xg, wg, 1.0, pdfs[dist, :], rho_ld, tau_ld, gdif_out)

    return xg, wg, gdir_out, gdif_out, gmdir_out, gmdif_out


def io_uncol_down(deltaL, nlayers, gdir, io, mp):
    """
    This function evaluates the downward uncollided direct solar intensity layer by layer"
    :param deltaL: Thickness of spatial cells (LAI/nlayers)
    :param nlayers: Number of layers in the canopy
    :param gdir: Output of the Gdir function
    :param io: intensity of the direct beam
    :param mp: cosine of theta_o
    :return: io_ucd: Downward uncollided direct solar radiation
    """

    L1 = 0.0
    io_ucd = np.zeros(nlayers)
    for k in range(nlayers):
        L2 = (k) * deltaL - (0.5 * deltaL) # eric has i+1
        prob = np.exp(-(1 / abs(mp)) * gdir * (L2 - L1))
        io_ucd[k] = io * prob

    return io_ucd


def io_uncol_up(deltaL, nlayers, gdir, gdif, io, mp, ng, xg, r_s):
    """
    Evaluates the upward uncollided direct solar radiation
    :param deltaL: Thickness of spatial cells (LAI/nlayers)
    :param nlayers: Number of layers in the canopy
    :param gdir: Output of the Gdir function
    :param gdif: Array with output of the gdif function
    :param io: intensity of the direct beam
    :param mp: cosine of theta_o
    :param ng: gauss quadrature order
    :param xg: gauss ordinates
    :param r_s: soil hemispherical reflectance (assumed Lambertian)
    :return: io_ucu: Upward uncollided direct solar radiation
             fo_ucdsoil: Downward uncollided flux density of direct solar radiation
             incident on the ground below the canopy.
    """

    # Uncollided direct solar radiation incident on the ground
    L1 = 0.0
    L2 = nlayers * deltaL
    fo_ucdsoil = abs(mp) * io * np.exp(-(1 / abs(mp)) * gdir * (L2-L1))

    # Upward uncollided intensity from reflection by soil of io_ucusoil
    io_ucusoil = (r_s / np.pi) * fo_ucdsoil

    #  Evaluate upward uncollided direct solar intensity layer by layer in all upw. dir.
    io_ucu = np.zeros((nlayers, ng, ng))
    for i in range(ng/2, ng):
        for j in range(ng):
            L2 = nlayers * deltaL
            for k in range(nlayers-1, -1, -1):
                L1 = (k) * deltaL - (0.5 * deltaL)  # Eric has k+1
                prob = np.exp(-(1 / np.abs(xg[i])) * gdif[j, i] * (L2 - L1))
                io_ucu[k, j, i] = io_ucusoil * prob

    return fo_ucdsoil, io_ucu


def id_uncol_down(deltaL, nlayers, gdif, i_d, ng, xg):
    """
    Evaluates the downward uncollided diffuse sky radiation
    :param deltaL: Thickness of spatial cells (LAI/nlayers)
    :param nlayers: Number of layers in the canopy
    :param gdif: Array with output of the gdif function
    :param i_d: intensity of the diffuse beam
    :param ng: gauss quadrature order
    :param xg: gauss ordinates
    :return: id_ucd: Upward uncollided diffuse sky radiation

    """
    #  Evaluate downward uncollided diffuse sky intensity layer by layer in all dwnw. dir.
    id_ucd = np.zeros((nlayers, ng, ng))
    for i in range(ng/2):
        for j in range(ng):
            L1 = 0.0
            for k in range(nlayers):
                L2 = (k) * deltaL - (0.5 * deltaL)  # Eric has k+1
                prob = np.exp(-(1 / np.abs(xg[i])) * gdif[j, i] * (L2 - L1))
                id_ucd[k, j, i] = i_d * prob

    return id_ucd


def id_uncol_up(deltaL, nlayers, gdif, i_d, ng, xg, wg, r_s):
    """
    Evaluates the upward uncollided diffuse sky radiation
    :param deltaL: Thickness of spatial cells (LAI/nlayers)
    :param nlayers: Number of layers in the canopy
    :param gdif: Array with output of the gdif function
    :param i_d: intensity of the diffuse beam
    :param ng: gauss quadrature order
    :param xg: gauss ordinates
    :param wg: gauss weights
    :param r_s: soil hemispherical reflectance (assumed Lambertian)
    :return: id_ucd: Upward uncollided diffuse sky radiation
             fd_ucdsoil: Downward uncollided flux density of diffuse sky radiation
             incident on the ground below the canopy.
    """

    #  Evaluate diffuse sky flux density incident on the soil below canopy
    L1 = 0.0
    L2 = nlayers * deltaL
    sum1 = 0.0
    conv = convfactors(2.0 * np.pi, 0.0)[0]   # We only need the first factor
    for i in range(ng/2):
        sum2 = 0.0
        for j in range(ng):
            prob = np.exp(-(1 / np.abs(xg[i])) * gdif[j, i] * (L2-L1))
            id_ucdsoil = i_d * prob
            sum2 += wg[j] * id_ucdsoil
        sum1 += wg[i] * np.abs(xg[i] * sum2 * conv)

    fd_ucdsoil = sum1

    #  Evaluate upward uncollided intensity due to reflection from soil
    id_ucusoil = fd_ucdsoil * (r_s / np.pi)

    #  Evaluate downward uncollided diffuse sky intensit layer by layer in all downward directions
    id_ucu = np.zeros((nlayers, ng, ng))
    for i in range(ng/2, ng):
        for j in range(ng):
            L2 = nlayers * deltaL
            for k in range(nlayers-1, -1, -1):
                L1 = (k) * deltaL - (0.5 * deltaL)  # Eric has k+1
                prob = np.exp(-(1 / np.abs(xg[i])) * gdif[j, i] * (L2 - L1))
                id_ucu[k, j, i] = id_ucusoil * prob

    return fd_ucdsoil, id_ucu


def fcs(nlayers, ng, wg, gmdir, gmdif, io_ucd, io_ucu, id_ucd, id_ucu):
    """
    Evaluates the first collision source Q
    :param nlayers:
    :param ng:
    :param wg:
    :param gmdir:
    :param gmdif:
    :param io_ucd:
    :param io_ucu:
    :param id_ucd:
    :param id_ucu:
    :return: Q: First collision source
    """

    #  Evaluate first collision source due to Io_ucd
    q = np.zeros((nlayers, ng, ng))
    for i in range(ng):
        for j in range(ng):
            for k in range(nlayers):
                q[k, j, i] = (1.0 / np.pi) * gmdir[j, i] * io_ucd[k]

    #  Evaluate first collision source due to Io_ucu, Id_ucd, Id_ucu

    conv1 = convfactors(1.0, -1.0)[0]  # We only need the first one
    conv2 = convfactors(2.0 * np.pi, 0.0)[0]
    for i in range(ng):
        for j in range(ng):
            for k in range(nlayers):
                sum1 = 0.0
                for n in range(ng):
                    sum2 = 0.0
                    for m in range(ng):
                        sum2 += wg[m] * (1.0 / np.pi) * gmdif[m, n, j, i] * \
                         (io_ucu[k, m, n] + id_ucd[k, m, n] + id_ucu[k, m, n])

                    sum1 += wg[n] * sum2 * conv2

                q[k, j, i] += sum1 * conv1

    return q


def sweep_down(nlayers, ng, xg, wg, gdif, deltaL, jj, r_s, ic):
    """
    This function sweeps downward in the phase-space mesh, handles the bottom boundary condition and
    evaluates the upward Ic at the ground
    :param nlayers:
    :param ng:
    :param xg:
    :param wg:
    :param gdif:
    :param deltaL:
    :param jj:
    :param r_s:
    :return: Value of collided intensity Ic on the ground
    """

    # Sweep downwards
    for i in range(ng/2):
        for j in range(ng):
            fij = (xg[i] / deltaL) - (0.5 * gdif[j, i])
            aij = ((0.5 * gdif[j, i]) + (xg[i] / deltaL)) / fij
            bij = 1.0 / fij
            for k in range(nlayers):
                ic[k+1, j, i] = aij * ic[k, j, i] - bij * jj[k, j, i]

    # Evaluate flux density incident on the ground
    sum1 = 0.0
    conv = convfactors(2.0 * np.pi, 0.0)[0]
    for i in range(ng/2):
        sum2 = 0.0
        for j in range(ng):
            sum2 += wg[j] * ic[nlayers+1, j, i]
        sum1 += wg[i] * np.abs(xg[i]) * sum2 * conv

    fc_soil = sum1

    #  Evaluate Ic upward at the ground
    for i in range(ng/2, ng):
        for j in range(ng):
            ic[nlayers+1, j, i] = (fc_soil * r_s) / np.pi  # Same as Eric

    return ic


def sweep_up(nlayers, ng, xg, gdif, deltaL, jj, ic, eps):
    """
    This routine sweeps upward in the phase-space mesh and checks for convergence
    :param nlayers:
    :param ng:
    :param xg:
    :param gdif:
    :param deltaL:
    :param jj:
    :param ic: Collided intensity field
    :param eps: Epsilon, convergence criterion for iteration of the scattering integral
    :return:
    """

    #  Save earlier iterate of ic[0, j, i] for convergence checking
    ic_old = np.zeros(ic.shape[1:3])
    for i in range(ng/2, ng):
        for j in range(ng):
            ic_old[j, i] = ic[0, j, i]

    #  Sweep upwards

    for i in range(ng/2, ng):
        for j in range(ng):
            fij = ((xg[i] / deltaL) + (0.5 * gdif[j, i]))
            cij = ((xg[i] / deltaL) - (0.5 * gdif[j, i])) / fij
            dij = 1.0 / fij
            for k in range(nlayers, -1, -1): #switch back to nlay-1
                ic[k, j, i] = cij * ic[k+1, j, i] + dij * jj[k, j, i]
                #  Eric has k instead of k+1 in the right side

    #  Check for convergence
    converg = True
    ct = 0
    for i in range(ng/2, ng):
        for j in range(ng):
            # ct += 1
            tv = abs(ic[0, j, i] - ic_old[j, i]) #changed back to 0
            # WHY indexed on 1?
            if tv > eps:
                converg = False
                #print 'Iter {0}: tv = {1}'.format(ct, tv)
                break
            else:
                converg = True
                break

    return converg, ic


def multicoll_s(nlayers, ng, wg, gmdif, ic):
    """
    Evaluates the multiple-collision source S
    :param nlayers:
    :param ng:
    :param wg:
    :param gmdif:
    :param ic:
    :return: S: Multiple collision source
    """

    #  Evaluate multiple collision source due to Ic
    conv1 = convfactors(1.0, -1.0)[0]
    conv2 = convfactors(2.0 * np.pi, 0.0)[0]
    s = np.zeros((nlayers+1, ng, ng))  # Eric has nlayers+1

    for i in range(ng):
        for j in range(ng):
            for k in range(nlayers):
                sum1 = 0.0
                for n in range(ng):
                    sum2 = 0.0
                    for m in range(ng):
                        ic_cellcenter = 0.5 * (ic[k, m, n] + ic[k+1, m, n])
                        sum2 += wg[m] * (1.0 / np.pi) * gmdif[m, n, j, i] * ic_cellcenter

                    sum1 += wg[n] * sum2 * conv2

                s[k, j, i] = sum1 * conv1

    return s


def energy_bal(nlayers, ng, xg, wg, deltaL, r_s, rhold, tauld, gdir, gdif, fo_ucdsoil,
               fd_ucdsoil, io_ucd, id_ucd, io_ucu, id_ucu, ic):
    """
    This function calculates the energy balance
    :param nlayers:
    :param ng:
    :param xg:
    :param wg:
    :param deltaL:
    :param r_s:
    :param rhold:
    :param tauld:
    :param gdir:
    :param gdif:
    :param fo_ucdsoil:
    :param fd_ucdsoil:
    :param io_ucd:
    :param id_ucd:
    :param io_ucu:
    :param id_ucu:
    :param ic:
    :return: Values of reflectance, transmittance and absorptance
    """

    conv = convfactors(2.0 * np.pi, 0.0)[0]

    #  Evaluate collided hemispherical transmittance
    ht_uc = fo_ucdsoil + fd_ucdsoil

    #  Evaluate collided hemispherical transmittance
    sum1 = 0.0
    for i in range(ng/2):
        sum2 = 0.0
        for j in range(ng):
            sum2 += wg[j] * ic[nlayers, j, i]
            #print sum2
            #CHECK THIS NLAYERS AND THOSE ON THE OTHER IC'S

        sum1 += wg[i] * abs(xg[i]) * sum2 * conv
    ht_c = sum1

    #  Evaluate uncollided hemispherical reflectance

    L1 = 0.0
    L2 = nlayers * deltaL
    sum1 = 0.0
    for i in range(ng/2, ng):
        sum2 = 0.00
        for j in range(ng):
            prob = np.exp(-(1 / np.abs(xg[i])) * gdif[j, i] * (L2 - L1))
            i_ucu = ((fo_ucdsoil + fd_ucdsoil) * r_s / np.pi) * prob
            sum2 += wg[j] * i_ucu

        sum1 += wg[i] * abs(xg[i]) * sum2 * conv

    hr_uc = sum1

    #  Evaluate collided hemispherical reflectance
    sum1 = 0.0
    for i in range(ng/2, ng):
        sum2 = 0.0
        for j in range(ng):
            sum2 += wg[j] * ic[0, j, i]

        sum1 += wg[i] * np.abs(xg[i]) * sum2 * conv
    hr_c = sum1

    #  Evaluate canopy absorption from Io_ucd
    sum1 = 0.0
    for k in range(nlayers):
        sum1 += io_ucd[k] * gdir

    abo_ucd = sum1 * (1.0 - (rhold + tauld)) * deltaL

    #  Evaluate canopy absorption from Io_ucu
    sum1 = 0.0
    for k in range(nlayers):
        for i in range(ng/2, ng):
            sum2 = 0.0
            for j in range(ng):
                sum2 += wg[j] * io_ucu[k, j, i] * gdif[j, i]

            sum1 += wg[i] * sum2 * conv

    abo_ucu = sum1 * (1.0 - (rhold + tauld)) * deltaL

    #  Evaluate canopy absorption from Id_ucd
    sum1 = 0.0
    for k in range(nlayers):
        for i in range(ng/2):
            sum2 = 0.0
            for j in range(ng):
                sum2 += wg[j] * id_ucd[k, j, i] * gdif[j, i]

            sum1 += wg[i] * sum2 * conv

    abd_ucd = sum1 * (1.0 - (rhold + tauld)) * deltaL

    #  Evaluate canopy absorption from Id_ucu
    sum1 = 0.0
    for k in range(nlayers):
        for i in range(ng/2, ng):
            sum2 = 0.0
            for j in range(ng):
                sum2 += wg[j] * id_ucu[k, j, i] * gdif[j, i]

            sum1 += wg[i] * sum2 * conv

    abd_ucu = sum1 * (1.0 - (rhold + tauld)) * deltaL

    #  Evaluate canopy absoption from Ic
    sum1 = 0.0
    for k in range(nlayers):
        for i in range(ng):
            sum2 = 0.0
            for j in range(ng):
                sum2 += wg[j] * ic[k, j, i] * gdif[j, i]

            sum1 += wg[i] * sum2 * conv

    ab_c = sum1 * (1.0 - (rhold + tauld)) * deltaL

    ab_uc = abo_ucd + abo_ucu + abd_ucd + abd_ucu

    # Calculating the totals
    hr = hr_uc + hr_c
    ht = ht_uc + ht_c
    ab = ab_uc + ab_c
    usa = (1 - r_s) * ht_uc
    csa = (1 - r_s) * ht_c
    sa = (1 - r_s) * ht
    ebal = hr + ab + (1 - r_s) * ht

    # return dict with variables
    vn = ('Uncollided Hemispherical Reflectance (hc_uc)',
            'Collided Hemispherical Reflectance (hr_c)',
            'Hemispherical Reflectance (hr)',
            'Uncollided Hemispherical Transmittance (ht_uc)',
            'Collided Hemispherical Transmittance (ht_c)',
            'Hemispherical Transmittance (ht)',
            'Uncollided Canopy Absorptance (ab_uc)',
            'Collided Canopy Absorptance (ab_c)',
            'Canopy Absorptance (ab)',
            'Uncollided Soil Absorptance (usa)',
            'Collided Soil Absorptance (csa)',
            'Soil Absorptance (sa)',
            'Energy Balance (=1.0)')

    var_collect = [(vn[0], hr_uc), (vn[1], hr_c), (vn[2], hr),
                   (vn[3], ht_uc), (vn[4], ht_c), (vn[5], ht),
                   (vn[6], ab_uc), (vn[7], ab_c), (vn[8], ab),
                   (vn[9], usa), (vn[10], csa), (vn[11], sa), (vn[12], ebal)]

    for v in var_collect:
        print v
