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
