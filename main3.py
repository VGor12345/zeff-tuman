# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 19:26:13 2024

@author: home
main2 + main3 
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import string
import re


    #инициализация из файла
time1 = []
current = []
time2 = []
hf_sxr = []
time3 = []
ne5 = []
time4 = []
sxr0_50 = []

f = open (r'C:\Users\home\Documents\учебники\НИР\sxr_z_eff\разряды\24032708.txt')
with f:
    for _ in range(4):
        next(f)
    lines = f.readlines() #[line.rstrip() for line in f] 
for line in lines:
    line = line.split()
    if line[0] != '---------':
        time1.append(float(line[0]))
    if line[1] != '---------':
        current.append(float(line[1]))
    if line[2] != '---------':
        time2.append(float(line[2]))
    if line[3] != '---------':
        hf_sxr.append(float(line[3]))
    if line[4] != '---------':
        time3.append(float(line[4]))
    if line[5] != '---------':
        ne5.append(float(line[5]))
    if line[6] != '---------':
        time4.append(float(line[6]))
    if line[7] != '---------':
        sxr0_50.append(float(line[7]))

f.close();
    
    #интерполяция плотности time3 -> time2 
tstart = time3[0]
istart = 0 #15304
for i in range(len(time2)):
    if  time2[i] < tstart:
        continue
    elif time2[i] >= tstart:
        istart = i
        break
    break
tfin = time3[len(time3)-1]
ifin = 0 #96963
for i in range(len(time2)):
    if  time2[len(time2)-1-i] < tfin:
        continue
    elif time2[len(time2)-1-i] >= tfin:
        ifin = len(time2)-1-i-2900
        break
    break
timezeff = time2[istart:ifin]
hf_sxrzeff = hf_sxr[istart:ifin]
ne5zeff = np.interp(timezeff, time3, ne5)

"""
fig, ax = plt.subplots()
ax.set(xlabel='time (ms)', ylabel='n',
       title='n')
ax.grid()
ax.plot(timezeff,ne5zeff)
plt.show()
"""

    #преобразование сигнала hf_sxr 
from scipy.signal import savgol_filter
hf_sxr_smooth = savgol_filter(hf_sxr, 51, 3) # window size 51, polynomial order 3
hf_zero = sum(hf_sxr_smooth[2000:10000])/len(hf_sxr_smooth[2000:10000])
hf_sxr_smooth = hf_sxr_smooth - hf_zero #не забыть домножить потом до 1e-3 !!
hf_sxr_smooth_zeff = hf_sxr_smooth[istart:ifin]

"""
fig, ax = plt.subplots()
ax.set(xlabel='time (ms)', ylabel='signal',
       title='signals')
ax.grid()
#fig.savefig("test.png")#ax.plot(ri, Texpvec(lambdaMax, Ti), ri, Texpvec(lambda1, Ti)) #ax.plot(ri, Texpvec(lambda1, Ti))
#ax.plot(time2,hf_sxr,time2,hf_sxr_smooth,time2[2000:10000],hf_sxr_smooth[2000:10000])
ax.plot(time2,hf_sxr,time2,hf_sxr_smooth)
#ax.plot(time2,hf_sxr)
plt.show()    
"""

    #вычисление интеграла
import scipy.constants
import math
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import numpy as np

T0 = 300 #эВ
constant_formula = 1.9 * 10 ** -28 #константа перед интегралом

def Pp(T,n):
    # данные для апроксимации электронной температуры
    TeCentral = T
    TePeriphery = 1

    # данные для апроксимации электронной концентрации
    neCentral = n * 1.5 * pow(10, 13) #cm-3 3.* pow(10, 15)/50
    nePeriphery = pow(10, 9)

    gaunt_factor = 1.
    omega = 1.854e-3 #0.167 #телесный угол
    tgalpha = 2.437e-2 #tg половины угла раствора 
    radiusA = 25 #малый радиус такомака туман
    distanceFromDiagnosticsToTokomak = 26 #раастояние от диагностики до токамака
    
    #граници интегрирования по лямда
    lambda1 = 1
    lambdaMax = 6200 / TeCentral
    
    C = gaunt_factor * omega * scipy.pi * tgalpha ** 2 
    C_int = C / 12395
    
    ri = np.linspace(-1 * radiusA, radiusA, 1000)
    ni = nePeriphery + (neCentral - nePeriphery) * (1. - pow(ri / radiusA, 2)) ** 2  #neCentral + (neCentral - nePeriphery) * (1. - pow(ri / radiusA, 2)) ** 2 
    #ni = neCentral + (neCentral - nePeriphery) * (1. - pow(ri / radiusA, 2)) ** 2 
    #ni = (neCentral - nePeriphery) * ((1 - (ri/radiusA)**2)**2) + neCentral
    Ti = TePeriphery + (TeCentral - TePeriphery) * (1. - pow(ri / radiusA, 2)) ** 2  #TeCentral + (TeCentral - TePeriphery) * (1. - pow(ri / radiusA, 2)) ** 2 
    #Ti = TeCentral + (TeCentral - TePeriphery) * (1. - pow(ri / radiusA, 2)) ** 2 
    #Ti = (TeCentral - TePeriphery) * ((1 - (ri/radiusA)**2)**2) + TeCentral
    #plt.plot(ri,Ti)
    def Texp(lam, T_):
        return(math.exp(-12395. / (lam * T_)))
    
    Texpvec = np.vectorize(Texp, excluded=['lam'])
    
    tti = Texpvec(lambdaMax, Ti) - Texpvec(lambda1, Ti)
    yi = ni ** 2 * (ri + distanceFromDiagnosticsToTokomak + radiusA) ** 2 * Ti ** 0.5 * tti
    
    I = scipy.integrate.simps(yi, x=ri, dx=0.01)
    I = C_int * I
    return(I)


    #вычисление zeff в цикле     
zeff = []
a = 0.113
#print(Pp(T0, ne5zeff[10000]))
for i in range(len(timezeff)):
    zeff.append(a * hf_sxr_smooth_zeff[i] * 1e-3 / (1.9 * constant_formula * Pp(T0, ne5zeff[i])))
#i = 20000
#zeff_ = a * hf_sxr_smooth_zeff[i] * 1e-3 / (1.9 * constant_formula * Pp(T0, ne5zeff[i]))
#print(hf_sxr_smooth_zeff[i],ne5zeff[i],zeff_)
#plt.plot(timezeff,hf_sxr_smooth_zeff)
#plt.plot(timezeff, a * 1e-3 / (1.9 * constant_formula * Pp(T0, ne5zeff[i])))

for i in range(len(timezeff)):
    zeff.append(a * 1e-3 / (1.9 * constant_formula * Pp(T0, ne5zeff[i])))
plt.plot(timezeff, zeff)



fig, ax = plt.subplots(3)
ax[0].set(xlabel='time (ms)', ylabel='current (ka)',
       title='current')
ax[0].grid()
ax[0].plot(time1,current)
ax[1].set(xlabel='time (ms)', ylabel='density ($cm^{-2}$)')
ax[1].grid()
ax[1].plot(timezeff,ne5zeff)
ax[2].set(xlabel='time (ms)', ylabel='zeff',title='zeff')
ax[2].grid()
ax[2].plot(timezeff[10000:-20000],zeff[10000:-20000]) #
plt.show()
    
#plt.plot(timezeff[10000:-20000],zeff[10000:-20000])







