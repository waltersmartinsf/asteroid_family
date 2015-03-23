# Asteroid Family Package

#Importando bibliotecas que serao usadas e definindo janela grafica
from pylab import *
import matplotlib.pyplot as plt
import numpy as np
import time
from pandas import *
from astropy import constants as const
from scipy.stats import maxwell
import os
import seaborn
plt.rcParams['figure.figsize'] = (14.0,8.0) # change figure size

#Useful constants
pi = np.pi #valor de pi= 3.141592653589793
DEGRAD = 180.0/pi #change rad to degrees
GRAV = 6.67259e-8 #gravitational constant in  cgs
UA = 1.49597870e+11 #1 UA to m
AU = 1.49597870e13 # 1cm in AU
YEAR = 31557600. #1 year in seconds
km=100000 #1km in cm

#import data
cts_velocidade = read_csv('./data/Resultados_Ajuste_Velocidade.csv', sep=' ')
cts_massa = read_csv('./data/Resultados_Ajuste_Massa.csv',sep=' ')
cts_energia = read_csv('./data/cts_energy.csv', sep=' ')

ind, razaoV, chi_velocidade, k3, k2, k1, T3, T2, T1, log_v1, log_v2 = cts_velocidade.values[0]
ind, razaoM, chi_massa, b1, b2, b3, W1, W2, W3, log_m1, log_m2 =cts_massa.values[0]
ind, C1, C2 = cts_energia.values[0]

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
