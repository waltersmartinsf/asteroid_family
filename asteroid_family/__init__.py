# Asteroid Family Package
"""
Asteroid Family Package

Author: Walter S. Martins Filho
email: walter at on.br
___

This package is usefull to simulate a synthethic asteroid family.

Functions:

homogeneus_family
differentiated_family
mass_distribution
mean_vej
mean_vej_distribution
velocity_field
catastrophic_energy
yarkovsky_dadt
yakovsky_change_unit
mag_absoluta

"""

import numpy as np
import time
import astropy.constants as const
from scipy.stats import maxwell
import os
from sympy import *

#Useful constants
pi = np.pi #valor de pi= 3.141592653589793
DEGRAD = 180.0/pi #change rad to degrees
GRAV = 6.67259e-8 #gravitational constant in  cgs
UA = 1.49597870e+11 #1 UA to m
AU = 1.49597870e13 # 1cm in AU
YEAR = 31557600. #1 year in seconds
km=100000 #1km in cm
AUyear = 2.10945021e-6 # 1cm/s in AU/year

#import data
#cts_velocidade = read_csv('./data/Resultados_Ajuste_Velocidade.csv', sep=' ')
#cts_massa = read_csv('./data/Resultados_Ajuste_Massa.csv',sep=' ')
#cts_energia = read_csv('./data/cts_energy.csv', sep=' ')

k3, k2, k1 = 0, 0.377843111, 1.2270257964
T3, T2, T1 = -0.376002915, 0.034284001, 1.4030194798
log_v1, log_v2 = 1.6118268805, 1.0858658107

#Data from cumulative distribution of fragments
b1, b2, b3 =  0.7856965123, 0.5894131708, 1.6077446609
W1, W2, W3 = -0.8236094358, -0.2323084658, -2.7077789041
log_m1, log_m2 = -2.9149258064, -2.4309082675

#Data from velocity-mass relationship
C1, C2 = 12.6685862091, 0.0278160911

B1, B2, B3 = 10**W1, 10**W2, 10**W3
r1 = (1-b1)/k1
r2 = (1-b2)/k2
m1 = 10**(log_m1)
m2 = 10**(log_m2)
v1 = 10**(log_v1)
v2 = 10**(log_v2)

# print b1, b2, b3
# print W1, W2, W3
# print k1, k2, k3
# print T1, T2, T3
# print log_m1, log_m2
# print log_v2, log_v1
# print C1, C2
# print r1, r2

#Test function:
def test_function():
    print("Baka-chan!")

#Functions for everything: differentiated and homogeneus families of asteroids
def mean_vej(norm_massa):
    """
    Calculate the mean velocity of ejection
    ___

    massa: normalized mass of a particle, double float number

    Example:

    In [1]: asteroid_family.mean_vej(0.5)
    Out[1]: 12.186130110528911

    """
    if norm_massa < m1:
        vej = C1*(norm_massa**(-r1))
    if norm_massa > m1 and norm_massa < m2:
        vej = C2*(norm_massa**(-r2))
    if norm_massa > m2:
        vej = v2
    return vej

def mass_distribution(maximum):
    """
    Create a normalized mass distribuiton, given a maximum and a minimun value of normalized mass.
    Return the normalized mass distribution for a asteroid family.
    ___

    maximum: the probability normalized major mass for the synthethic distribution of mass

    Example:

    In [1]: asteroid_family.mass_distribution(0.4)
    Out[1]:
    array([  3.03109848e-05,   3.03109848e-05,   3.03109848e-05, ...,
         1.37594174e-02,   1.39187440e-02,   1.40799154e-02])
    """
    log_massa = np.arange(-4.5, maximum, 0.001)
    log_N = np.zeros(len(log_massa))

    for i in range(len(log_massa)):
        if log_massa[i] <= log_m1:
            log_N[i] = (-b1) * log_massa[i] + float(W1)
        if log_massa[i] > log_m1 and log_massa[i] < log_m2:
            log_N[i] = (-b2) * log_massa[i] + float(W2)
        if log_massa[i] > log_m2:
            log_N[i] = (-b3) * log_massa[i] + float(W3)

    N = 10**(log_N)
    #print N[N >= 1], len(N[N >= 1])
    Prob = N[N >= 1]/sum(N[N >= 1])
    log_massa = log_massa[N >= 1]

    quant = 100
    massa_acumulada = 0
    while massa_acumulada < 1.0:
        sorteio = np.random.choice(log_massa, size=quant, p=Prob)
        massa = 10**(sorteio)
        massa_acumulada = sum(massa)
        quant = quant + 10

    #Normalizando
    massa = massa/massa_acumulada
    massa = np.sort(massa)

    return massa

def mean_vej_distribution(massa):
    """
    Calculate a mean ejection velocity field for a given normalized data distribution of mass
    ___

    massa: normalized mass distribuiton, array-like object

    Example:

    """
    VejM = np.zeros(len(massa))
    for i in range(len(VejM)):
        VejM[i] = mean_vej(massa[i])

    return VejM

def velocity_field(mass, show_time):
    """
    mass: normalized mass distribution
    show_time: string. If show_time=='YES', it's show the duration of the process.

    Return:

    Velocity field for the given normalized mass distribution, assmuming a Maxwellian Distribution for the velocity field for each mean velocity.
    ___

    Example:

    In [1]: asteroid_family.mass_distribution(0.3)

    In [2]: a = asteroid_family.mass_distribution(0.3)

    b = asteroid_family.mean_vej_distribution(a)

    In [3]: asteroid_family.velocity_field(b,'YES')
    Duration [seconds] =  0.412480831146
    Out[3]:
    array([  45.89043405,  118.93529017,   44.34338289, ...,    5.96927437,
          8.80950761,   17.1578176 ])
    """

    VejM = mean_vej_distribution(mass)
    min_Vej = 0
    tempo = time.time()

    ejecao = np.zeros(len(VejM))

    for i in range(len(VejM)):
        ejecao[i] = maxwell.rvs(scale=VejM[i]/np.sqrt(3.))
        if ejecao[i] < min_Vej:
            while ejecao[i] < min_Vej:
                ejecao[i] = maxwell.rvs(scale=VejM[i]/np.sqrt(3.))

    tempo = abs(tempo - time.time())
    if show_time == 'YES':
        print('\n Obtain velocity field: duration [seconds] = '+str(tempo)+'\n')

    return ejecao

def catastrophic_energy(rpb,rho,Vi):
    """
    Return the catastrofic energy in c.g.s.
    Based in Stewart and Leinhardt(2009)

    rpb: parental radius in km
    rho: mean density of the system in cgs
    Vi: impact velocity in km/s
    """
    # qs, qg, fi, mi: constants of material
    if rho < 4.0:
        qs, qg, fi, mi = 7.0e4, 1e-4, 8, 0.5
    else:
        qs, qg, fi, mi = 500.0, 1.0e-4, 6., 0.4
    Mcomb = (4/3.)*np.pi*((rpb*1e5)**3)*rho
    Vi = (Vi*1.e5)
    RC = ((3*Mcomb)/(4*pi))**(1/3.)
    QD = qs*(RC**(9*mi*(3-2*fi)))*(Vi**(2-3*mi))+qg*(RC**(3*mi))*(Vi**(2-3*mi))
    return QD

def min_Vej(E, Mtot, Mmin, v2, C1, C2, r1, r2, B1, b1, B2, b2, B3, b3):
    """
    Obtain the minimum velocity for given system. This function is used in other two functions: homogeneus_family and differentiated_family.

    E: kinetic energy.

    Mtot: total mass of the system in c.g.s. or normalized.

    Mmin: minimum mass of the system in c.g.s. or normalized.

    v2, C1, C2, r1, r2, B1, b1, B2, b2, B3, b3: constants in the model
    """
    parte1 = (((C1**(2)) * (Mtot**(2*r1)) * (b1*B1))/(2*(1-b1-2*r1)))*(((m1*Mtot)**(1-b1-2*r1))-(Mmin**(1-b1-2*r1)))
    parte2 = (((C2**(2)) * (Mtot**(2*r2)) * (b2*B2))/(2*(1-b2-2*r2)))*(((m2*Mtot)**(1-b2-2*r2))-((m1*Mtot)**(1-b2-2*r2)))
    parte3 = (((v2**2) * (b3*B3))/(2*(1-b3)))*((Mtot**(1-b3))-((m2*Mtot)**(1-b3)))
    Parcial = (parte1 + parte2 + parte3)
    return np.sqrt(E/Parcial)

def mag_absoluta(pv, D):
    '''
    pv: geometric albedo
    D: diameter in km

    Based in equation given by Parker et al.(2008).
    '''
    return 18.1 - 2.4 * np.log10(pv/0.1) - 5 * np.log10(D)

# Differentiated asteroid family
def differentiated_family(composition, rpb, rho_mantle, rho_core, Vi, fke, maximum, pv_core, pv_mantle, show_time):
    """
    composition: string. Given the initial meteoritic composition for estimate the core radius. Using data from Gaffey et al.(1993).
    Compostion can be use 'H' for condritic ordinary H composition, 'L' for ordinary L composition, 'LL' for ordinary LL composition,
    'CO' for condritic carbonaceus O composition, and 'CV' for condritic carbonaceus V composition.

    rpb: radius of the parental body in kilometers.

    rho_mantle: density of the mantle in c.g.s. units.

    rho_core: density of the core in c.g.s. units.

    Vi: impact velocity in km/s.

    fke: inelastic of the parental body.

    maximum: value of the maximum mass of the distrbution of the fragments.

    pv_core: geometric albedo for the fragments from the core of the parental body.

    pv_mantle: geometric albedo for the fragments from the mantle of the parental body.

    show_time: string. If show_time=='YES', it's show the duration of the process.
    ___

    Return:

    mass distribution [g] (array-like), velocity field (array-like), density in c.g.s. (array-like), radius in km (array-like), absolute magnitude (array-like)
    ___

    Example:

    In [1]: massa, vej, densidade, raio, H = asteroid_family.differentiated_family('H',130,3,7,5,0.01,0.5,0.2,0.2,'YES')

    Obtain velocity field: duration [seconds] =  0.377580881119

    """
    #Ratio between core radius and the parental radius given a initial compositon for the parental body.
    #Gaffey et al.(1993)
    H  = 0.51
    L  = 0.43
    LL = 0.39
    CV = (0.37+0.31)/2.
    #CO = (0.37+0.31)/2.

    #Parental body data
    rpb = rpb*1e5 #raio em cgs
    rho_nucleo = rho_core
    rho_manto = rho_mantle

    if composition == 'H':
        ratio = H
    if composition == 'L':
        ratio = L
    if composition == 'LL':
        ratio = LL
    if composition == 'CO' or composition == 'CV':
        ratio = CV

    rpb_nucleo = ratio*rpb

    #Total mass and volume
    volume_total = (4/3.)*np.pi*(rpb**3)
    volume_nucleo = (4/3.)*np.pi*(rpb_nucleo**3)
    massa_nucleo = rho_nucleo*volume_nucleo
    volume_manto = volume_total - volume_nucleo
    massa_manto = rho_manto*volume_manto
    Mtot  = massa_nucleo + massa_manto
    #Mmax = 0.5*Mtot
    #print 'Mtot[g] = ',Mtot
    #print 'Mmax[g] = ', Mmax
    Mmin = Mtot*(10**(-4.5))

    #Calculate the catastrophic energy for the system
    rho = (rho_manto+rho_nucleo)/2.
    QD = catastrophic_energy(rpb*1e-5, rho, Vi)

    #kinetic energy
    E = QD * Mtot * fke

    #Change the normalized system to the c.g.s. units
    B11, B22, B33 = float(B1*(Mtot**(b1))), float(B2*(Mtot**(b2))), float(B3*(Mtot**(b3)))
    #print 'B11, B22, B33 = ', B11, B22, B33

    Vmin = min_Vej(E, Mtot, Mmin, v2, C1, C2, r1, r2, B11, b1, B22, b2, B33, b3)

    #print 'Vmin [cm/s] = ', Vmin
    # print 'max V(m): ', (max(ejecao)*Vmin)*1e-5,' km/s'
    # print 'max <V(m)>: ', mean_vej(10**-4.5)*Vmin*1e-5,' km/s'
    # print 'V1(m1): ', (10**(log_v1))*Vmin*1e-5,' km/s'
    # print 'V2(m2): ', (10**(log_v2))*Vmin*1e-5,' km/s'

    C11 = C1 * Vmin * (Mtot**(r1))
    C22 = C2 * Vmin * (Mtot**(r2))
    #print 'C1, C2 = ', C11, C22

    #Create the mass distribution
    massa = mass_distribution(maximum)
    massa = massa*Mtot

    #Calculate the mean velocity
    ejecao = velocity_field(massa/sum(massa), show_time)
    Vej = ejecao*Vmin*AUyear

    #obtain the density of the fragments
    massa_acumulada = 0
    densidade = np.zeros(len(massa))
    i = 0
    for i in range(len(massa)):
        if i != max(range(len(massa))):
            if massa_acumulada < (massa_manto/Mtot):
                densidade[i] =rho_manto
                massa_acumulada = massa_acumulada + (massa[i]/Mtot)
                #print massa_acumulada
            else:
                densidade[i] = rho_nucleo
        else:
            densidade[i] = rho_nucleo

    #obtain the radius of the fragments
    raio = ((3./(4*np.pi))*(massa/densidade))**(1/3.)
    raio = raio*1e-5

    #obtain the absolute magnitude of the asteroid
    H = np.zeros(len(raio))
    for i in range(len(raio)):
        if densidade[i] == rho_nucleo:
            H[i] = mag_absoluta(pv_core,2*raio[i])
        if densidade[i] == rho_manto:
            H[i] = mag_absoluta(pv_mantle,2*raio[i])

    return massa, Vej, densidade, raio, H

def yarkovsky_dadt(D,obliq):
    """
    Calculate the Yarkovsky drift in semi-major axis in time [AU/Myr]. Based in Roig et al.(2008).
    ___
    D: diameter of the asteroid in km. (double-float or array-like object)
    obliq: obliquty of the asteroid in degrees (double-float or array-like object)
    ___

    Example:

    In [1]: mass, vej, density, radius, H = asteroid_family.differentiated_family('H',130,3,7,5,0.01,0.5,0.2,0.2,'YES')

    Obtain velocity field: duration [seconds] =  0.384955883026

    In [2]: import numpy as np

    In [3]: dadt = asteroid_family.yarkovsky_dadt(2*radius,np.zeros(len(raio)))

    In [4]: dadt
    Out[4]:
    array([  2.88357030e-05,   2.88357030e-05,   2.88357030e-05, ...,
                4.86320757e-06,   4.72341538e-06,   4.71617023e-06])

    In [5]: dadt = asteroid_family.yarkovsky_dadt(2*raio,np.random.uniform(0,np.pi,len(massa)))

    In [6]: print dadt
    [ -1.60343954e-05  -2.72928891e-05   1.11210017e-05 ...,  -3.34213765e-06
        -2.46130825e-06   4.51568450e-06]

    In [7]: import matplotlib.pyplot as plt
    In [8]: plt.hist(dadt)
    Out[8]:
    (array([ 169.,  282.,  301.,  362.,  320.,  299.,  326.,  310.,  264.,  177.]),
    array([ -2.88332378e-05,  -2.30789504e-05,  -1.73246630e-05,
            -1.15703756e-05,  -5.81608815e-06,  -6.18007481e-08,
            5.69248666e-06,   1.14467741e-05,   1.72010615e-05,
            2.29553489e-05,   2.87096363e-05]),
        <a list of 10 Patch objects>)
    In[9]: plt.show()

    """
    dadt = 2.5e-4*(1./D)*np.cos(obliq*(np.pi/180.))
    return dadt

def yakovsky_change_unit(dadt):
    """
    Change the drift of Yarkovsky effect from AU/Myr to AU/year. Essencial for some orbital integrators.
    """
    for i in range(len(dadt)):
        dadt [i] =dadt[i]/(1.e6)
    return dadt

#homogeneus family
def homogeneus_family(rpb, rho, Vi, fke, maximum, pv, show_time):
    """
    rpb: radius of the parental body in kilometers.

    rho: density of the fragments in c.g.s..

    Vi: impact velocity in km/s.

    fke: inelastic of the parental body.

    maximum: value of the maximum mass of the distrbution of the fragments.

    pv: geometric albedo for the fragments.

    show_time: string. If show_time=='YES', it's show the duration of the process.
    ___

    Return:

    mass distribution (array-like), velocity field (array-like), density (array-like), radius (array-like), absolute magnitude (array-like)
    ___

    Example:

    """
    #Parental body data
    rpb = rpb*1e5 #raio em cgs

    #Total mass and volume
    volume_total = (4/3.)*np.pi*(rpb**3)
    Mtot  = rho*volume_total
    Mmin = Mtot*(10**(-4.5))

    #Calculate the catastrophic energy for the system
    QD = catastrophic_energy(rpb*1e-5, rho, Vi)

    #kinetic energy
    E = QD * Mtot * fke

    #Change the normalized system to the c.g.s. units
    B11, B22, B33 = float(B1*(Mtot**(b1))), float(B2*(Mtot**(b2))), float(B3*(Mtot**(b3)))
    #print 'B11, B22, B33 = ', B11, B22, B33

    Vmin = min_Vej(E, Mtot, Mmin, v2, C1, C2, r1, r2, B11, b1, B22, b2, B33, b3)

    #print 'Vmin [cm/s] = ', Vmin
    # print 'max V(m): ', (max(ejecao)*Vmin)*1e-5,' km/s'
    # print 'max <V(m)>: ', mean_vej(10**-4.5)*Vmin*1e-5,' km/s'
    # print 'V1(m1): ', (10**(log_v1))*Vmin*1e-5,' km/s'
    # print 'V2(m2): ', (10**(log_v2))*Vmin*1e-5,' km/s'

    C11 = C1 * Vmin * (Mtot**(r1))
    C22 = C2 * Vmin * (Mtot**(r2))
    #print 'C1, C2 = ', C11, C22

    #Create the mass distribution
    massa = mass_distribution(maximum)
    massa = massa*Mtot

    #Calculate the mean velocity
    ejecao = velocity_field(massa/sum(massa), show_time)
    Vej = ejecao*Vmin*AUyear

    #obtain the density of the fragments
    densidade = np.zeros(len(massa))
    for i in range(len(massa)):
        densidade[i] = rho

    #obtain the radius of the fragments
    raio = ((3./(4*np.pi))*(massa/densidade))**(1/3.)
    raio = raio*1e-5

    #obtain the absolute magnitude of the asteroid
    H = np.zeros(len(raio))
    for i in range(len(raio)):
        H[i] = mag_absoluta(pv,2*raio[i])

    return massa, Vej, densidade, raio, H
###############################################################################

#Dynamic functions

def isotropic_velocity(VejM):
    """
    Assuming isotropic direction for the ejection field, this function obtain the radial, perpendicular and normal velocity of the fragment.
    ___
    VejM: velocity ejection of the fragment, double-float number
    ___

    Return: VR (array-like), VT (array-like), VW (array-like), V2 (array-like)

    VR: radial velocity (direction of the Sun)

    VT: tranversal velocity

    VW: normal velocity in the orbital plane

    V2: (VR^2 + VT^2 + VW^2)

    ___
    Example:

    """
    norm = np.sqrt(3)
    VejR=VejM*np.random.normal(0, 0.1)/norm
    VejT=VejM*np.random.normal(0, 0.1)/norm
    VejW=VejM*np.random.normal(0, 0.1)/norm
    Vej2=(VejR**2)+(VejT**2)+(VejW**2)
    return VejR, VejT, VejW, Vej2


def gauss_equations(Vej,a,e,i,period,show_time):
    """
    Calculate the Gauss equation for a given position of the parental body in proper space.
    ___

    Vej: velocity field of the fragments in AU/yr (array-like)

    a: semi-major axis of the parental body in astronomical units [AU]

    e: excentricity of the parental body

    i: inclination of the parental body in degrees

    period: orbital period of the parental body in years [yr]
    ___

    Return:
    VR: radial velocity of the fragments in AU/yr (array-like)

    VT: tranversal velocity of the fragments in AU/yr (array-like)

    VW: normal velocity of the fragments in AU/yr (array-like)

    A: semi-major axis in astronomical unit (array-like),

    E: excentricity (array-like),

    I: inclination I in degrees (array-like),

    dA: variations of semi-major axis in astronomical unit (array-like),

    dE: variations of excentricity (array-like)

    dI: variations of inclination in degrees (array-like).
    ___

    Example:

    #Creating a differentiated asteroid family in a = 2.26 AU, e = 0.09, and i = 6.28 degrees:

    massa, vej, rho, raio, mag = asteroid_family.differentiated_family('H',130,3,7,5,0.01,0.5,0.2,0.2,'YES')

    #Appling Gauss equations:
    VT, VR, VW, A, E, I, dA, dE, dI = asteroid_family.gauss_equations(vej, 2.263620, 0.096125, 6.287527, 3.41, 'YES')

    """
    #    f: true anomaly of the parental body
    #  wpf: true anomaly plus
    f = 95 #anomalia verdadeira (graus)
    wpf = 0 #relacao w+f .................................... Morbidelli et al.(1995)

    na = 2*np.pi*a/period #mean orbital velocity [AU/year]
    f = f/DEGRAD #Anomalia verdadeira: transformamos graus em radianos
    wpf = wpf/DEGRAD
    cosf = np.cos(f)
    sinf = np.sin(f)
    coswf = np.cos(wpf)
    eta1 = np.sqrt(1.0-(e**2))
    eta2 = 1.0+e*cosf

    tempo = time.time()
    A, E, I = [], [], []
    dA, dE, dI = [], [], []
    VR, VT, VW = [], [], []
    Vinf = 0
    contador = 0
    while contador < len(Vej):
        VejR, VejT, VejW, Vej2 = isotropic_velocity(Vej[contador])
        #print VejR, VejT, VejW
        VinfR = VejR
        VinfT = VejT
        VinfW = VejW
        #Calculando as variacoes em elementos orbitais_ eq.s de Gauss (Zappala et al., 1996)
        da = (a/na)*(2.0/eta1)*(eta2*VinfT+(e*sinf)*VinfR)
        de = ((e+2*cosf+e*(cosf)**2)/(eta2))*VinfT + sinf*VinfR
        de = (eta1/na)*de
        di = (eta1/na)*(coswf/eta2)*VinfW
        A.append(a+da)
        E.append(e+de)
        I.append(i+di*DEGRAD)
        dA.append(da)
        dE.append(de)
        dI.append(di*DEGRAD)
        VR.append(VinfR)
        VT.append(VinfT)
        VW.append(VinfW)
        #print 'Particula: ',contador+1
        contador = contador + 1


    tempo = time.time() - tempo
    if show_time == 'YES':
        print '\n Applied Gauss Equations: duration [seconds] = ', tempo,'\n'

    return VR, VT, VW, A, E, I, dA, dE, dI

def yarkovsky_drift(raio, a, obliq, P, A, rho, k, epsilon, PP):
    """
    Obtain the Yarkovsky drift as function of density, albedo and thermal inertia.
    Based on Walter S. Martins Filho's Bachelor Thesis.
    
    Thermal inertia is obtain based on the model of Delbo et al.(2007), Icarus, 190, 236-249.
    ___
    INPUT:
    radius: radius of the asteroid in km; numpy.ndarray object
    a: semimajor axis in astronomical unit; numpy.ndarray object
    obliq: obliquity in degrees; numpy.ndarray object
    A: bond albedo
    rho: density
    k: conductivity in MKS
    epsilon: 
    PP:
    P:
    
    OUTPUT:
    inercia: thermal inertia in MKS; numpy.ndarray object
    da: yarkovsky drift in AU/Myr; numpy.ndarray object
    """
    def drift(yarkovsky):
        F = yarkovsky
        diurno = ((-8.)*(1.-A))*phi*F*np.cos(obliq)/(9.*eta)
        sazonal = (4.*(1.-A))*phi*F*np.sin(obliq)**2/(9*eta)
        da = diurno + sazonal
        return da
    
    L = const.L_sun.value
    sigma = const.sigma_sb.value
    c = const.c.value
    obliq = np.radians(obliq)
    distancia = const.au.value*a
    eta = (2*np.pi)/(86400*PP)
    inercia = 300*(2*raio*1.e-3)**(-.48)
    phi = 3*L/(16*np.pi*(distancia**2)*rho*raio*c)
    temperature = (((1-A)*L)/(4*np.pi*(distancia**2)*epsilon*sigma))**(3/4.)
    omega = 2*np.pi/(3600*P)
    Theta = inercia*np.sqrt(omega)
    Theta = Theta/(epsilon*sigma)
    Theta = Theta/((((1-A)*L)/(4*np.pi*(distancia**2)*epsilon*sigma))**(3/4.))
    X = inercia*raio*np.sqrt(2*omega)/k
    
    ################# Sympy Routine ########################################
    #Obtainm the K-Functions
    x = Symbol('x')
    f1 = -(x-2)-exp(x)*((x-2)*cos(x)-x*sin(x))
    f2 = -x-exp(x)*(x*cos(x)+(x-2)*sin(x))
    f3 = x*(x+3)-exp(x)*(x*(x-3)*cos(x)-(3*(x-2))*sin(x))
    f4  = 3*(x+2)+exp(x)*((3*(x-2))*cos(x)+x*(x-3)*sin(x))
    j = (f1*f4+f2*f3)/(f1**2+f2**2)
    jj = (f4**2+f3**2)/(f1**2+f2**2)
    K1 = (1.+j)/x
    K2 = (1.+2.*j+j)/(x**2)
    K3 = (f1*f3-f2*f4)/((f1**2+f2**2)*x)
    #Yarkovsky's Force
    F = -K1*Theta/(1.+2.*K2*Theta+K3*Theta**2)
    yarkovsky = np.zeros(len(X))
    for i in range(len(X)):
        yarkovsky[i] = F[i].subs(x,X[i]).evalf() #obtain the yarkovsky force's value for the X-array
    ############# END of Sympy Routine #########################################

    da = drift(yarkovsky)
    return inercia, da