#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 22:42:26 2020
Updated on Mon March 2021, 08/25/2022

@author: Ensheng Weng
"""
import numpy as np
import matplotlib.pyplot as plt
from cycler import cycler

font = {'family' : 'serif',
        'color'  : 'darkred',
        'weight' : 'normal',
        'size'   : 14,
        }

linestyle_cycler = cycler('linestyle',['-','-.','--',':'])
ls_cycler = cycler('linestyle',
            [(0,()), # solid
             (0, (1, 10)), # loosely dotted
             (0, (1, 5)), # dotted
             (0, (1, 1)), # densely dotted
             (0, (5, 10)), # loosely dashed
             (0, (5, 5)), # dashed
             (0, (5, 1)), # densely dashed
             (0, (3, 10, 1, 10)), # loosely dashdotted
             (0, (3, 5, 1, 5)), # dashdotted
             (0, (3, 1, 1, 1)), # densely dashdotted
             (0, (3, 10, 1, 10, 1, 10)), # loosely dashdotdotted
             (0, (3, 5, 1, 5, 1, 5)), # dashdotdotted
             (0, (3, 1, 1, 1, 1, 1))] # densely dashdotdotted
                  )

#%% --- Global constants
nYears = 400 # Total years
PI = 3.1415926
# Environmental conditions
dPsi0  = 1.1 # mPa
TranspD = 0.5
initialR = 0.005 # m
TreeD0  = 0.25 # m
thetaW0 = 0.5

# Plant architecture and hydraulic parameters
alphaC = 150.0
alphaH = 40.0
thetaH = 0.5
thetaC = 1.5
rhoW   = 300 # kgC m-3
eta    = 0.575 # allocation to wood
WUE    = 1.2 # kgC/Ton H2O
Taper  = 0.7
kx0    = 9000 #Change this value doese not affect r_HU and r_SW
W0_0   = 1500.0 # m
r_DF   = 100 # 0.001 #1.0

#%%------- Main function ---------------
def PlantHydraulicModel():
    global do_increasingWTC0,do_constant_dPsi0, \
           W0,dDBH,TreeD,Acrown,Atrunk,Asw,Ktrunk,  \
           TreeHU,TreeW0,AgeSW,Transp,fareaM

    # ------ Local variables ----------
    TreeH  = np.zeros(nYears)
    dPsi   = np.zeros(nYears)
    Kx     = np.zeros(nYears) # xylem intrinsic conductivity, m per mPa per year
    Ax     = np.zeros(nYears) # new sapwood area
    Dsw    = np.zeros(nYears) # width of sapwood
    farea  = np.ones(nYears)
    Wyear  = np.zeros(nYears)
    Wtotal = np.zeros(nYears)
    HtotM  = np.zeros((nYears,nYears))
    KtrM  = np.zeros((nYears,nYears))

    # Initiate tree size
    TreeD[0] = 2* initialR
    TreeH[0] = alphaH * TreeD[0]**thetaH
    Atrunk[0] = PI * initialR**2
    Asw[0] = Atrunk[0]
    Ax[0]  = PI * initialR**2
    Acrown[0] = alphaC * TreeD[0]**thetaC

    # Initiate tree hydraulic capacity
    if do_increasingWTC0:
        W0[0] = W0_0 * (TreeD[0]/0.0025)**thetaW0
    else:
        W0[0] = W0_0 * (TreeD0/0.0025)**thetaH

    # Conductivity and tanso
    Kx [0]    = kx0 / TreeH[0]
    Ktrunk[0] = Kx[0] * Ax[0]

    #dPsi[0] = TranspD* Acrown[0]/Ktrunk[0]
    dPsi[0]   = dPsi0
    Wyear[0]  = Kx[0]* dPsi[0]
    Wtotal[0] = Wyear[0]
    farea[0]  = 1.0
    if do_constant_dPsi0:
        Transp[0] = Wyear[0]*Ax[0]/Acrown[0] #dPsi[0] * Ktrunk[0]/Acrown[0]
    else:
        Transp[0]= TranspD

    TreeHU[0] = Transp[0] * Acrown[0]
    TreeW0[0] = W0[0] * Ax[0]

    # Full matrix of age and rings
    HtotM[0,0]=Wtotal[0]
    fareaM[0,0]=farea[0]
    KtrM[0,0] = Ktrunk[0]

    # Model Run for nYears years
    farea[:] = 0.0
    for i in range(1,nYears): # yearly time step
        #New growth based on last year's transpiration
        dDBH[i] = Transp[i-1] * WUE * alphaC * eta /   \
                (0.25 * (2.0+thetaH)*PI*rhoW*Taper*alphaH)*  \
                (TreeD[i-1]**(thetaC-thetaH-1.0))

        # Initialize TreeHU and TreeW0
        TreeHU[i] = TreeHU[i-1]
        TreeW0[i] = TreeW0[i-1]

        #New sapwood growth and Tree structure update
        TreeD[i]  = TreeD[i-1] + dDBH[i]
        TreeH[i]  = alphaH * TreeD[i]**thetaH
        Acrown[i] = alphaC * TreeD[i]**thetaC
        Atrunk[i] = PI * 0.25 * TreeD[i]**2
        Ax[i]     = Atrunk[i] - Atrunk[i-1]
        Asw[i]    = Asw[i-1] + Ax[i]

        # Plant hydraulic properties
        farea[i]= 1.0 # New, 100% functional
        Kx[i] = kx0 /TreeH[i] # Constant conductivity, normalized by current height
        Ktrunk[i] = Ktrunk[i-1] + Kx[i]

        # W0 patterns
        if do_increasingWTC0:
            W0[i] = W0_0 * (TreeD[i]/0.0025)**thetaW0 # an increasing W0, a function of D
        else:
            W0[i] = W0_0 * (TreeD0/0.0025)**thetaH

        # Update hydraulic usage and trunk hydraulic status
        if do_constant_dPsi0:
            dPsi[i]   = dPsi0  # constant dPsi
        else:                  # constant Transp
            dPsi[i] = TranspD * Acrown[i]/Ktrunk[i] # increasing dPsi

        for j in range(i+1): # go through all wood rings
            r_lift = TreeH[j]/TreeH[i]
            if do_constant_dPsi0:
                Wyear[j] = Kx[j] * dPsi[i] * r_lift # Water transported in year i by ring j, over-used old tissues
            else:
                Wyear[j] = TranspD*Acrown[i]/Asw[i]
            # Update Wtotal and Transp
            Wtotal[j] = Wtotal[j] + Wyear[j]
            #Whole tree water flow
            Transp[i] = Transp[i] + Wyear[j] * Ax[j]*farea[j]
            #Update farea and Dsw
            #farea[j] = 1.0 - 1.0/(1.0 + np.exp(r_DF*(1.0-Wtotal[j]/W0[j])))
            farea[j] = 1. - np.exp(-r_DF * (1. - min(1.0,Wtotal[j]/W0[j])))
            Dsw[i]   = Dsw[i] + 0.5*dDBH[j]*farea[j]

            # Record each year's results
            fareaM[j,i] = farea[j]

        #Acrown-normalized transpiration rate (mm/yr)
        if do_constant_dPsi0:
            Transp[i] = Transp[i]/Acrown[i]
        else:
            Transp[i] = TranspD

        # Update TreeHU and TreeW0
        TreeHU[i] = TreeHU[i] + Transp[i] * Acrown[i]
        TreeW0[i] = TreeW0[i] + W0[i] * Ax[i]
        # Update whole tree hydraulic status
        Ktrunk[i] = 0.0
        Asw[i]    = 0.0
        for j in range(i+1): # go through all previous wood rings
            r_lift = TreeH[j]/TreeH[i]
            Ktrunk[i] = Ktrunk[i] + Kx[j]*Ax[j]*farea[j] #/TreeH[j]
            Asw[i] = Asw[i] + Ax[j]*farea[j]

    # Calculate sapwood age
    for i in range(nYears):
        nSW = 0
        totage = 0
        for j in range(i+1):
            if(fareaM[i,j]>0.5):
                nSW=nSW + 1
                totage = totage + (i-j)
        if(nSW>0):
           AgeSW[i]=totage/nSW
        else:
           AgeSW[i] = 0.0

#%%==============Main model run==================================
#----- Case run 1
#-----Variables passed to the main function -------------
do_increasingWTC0 = True # False #
do_constant_dPsi0 = True # False #
W0 = np.zeros(nYears)
dDBH = np.zeros(nYears)
TreeD = np.zeros(nYears)
Acrown = np.zeros(nYears)
Atrunk = np.zeros(nYears) # Cross sectional area
Asw = np.zeros(nYears) # area of sapwood
Ktrunk = np.zeros(nYears)
TreeHU = np.zeros(nYears) # Whole tree hydraulic usage
TreeW0 = np.zeros(nYears) # Whole tree WTC potential
AgeSW = np.zeros(nYears)
Transp = np.zeros(nYears) # transp/crown area
fareaM = np.zeros((nYears,nYears))
#fareaM[:,:] = np.nan

# -------- Variables for calculation
t = np.arange(nYears)

#------ Call main function
PlantHydraulicModel()

#------------- Keep results 1
W0_iW0     = W0.copy()
dDBH_iW0   = dDBH.copy()
TreeD_iW0  = TreeD.copy()
Acrown_iW0 = Acrown.copy()
Ktrunk_iW0 = Ktrunk.copy()
Atrunk_iW0 = Atrunk.copy()
Asw_iW0    = Asw.copy()
TreeHU_iW0 = TreeHU.copy()
TreeW0_iW0 = TreeW0.copy()
AgeSW_iW0  = AgeSW.copy()
Transp_iW0 = Transp.copy()
fareaM_iW0 = fareaM.copy()

#--New calculations
TrunkRadius_iW0 = 0.5 * TreeD
Hv_iW0   = 10000*Asw/(Acrown*3.5)
WTC_iW0  = (TreeW0-TreeHU)/(Acrown*Transp)
r_SW_iW0 = Asw/Atrunk
r_HU_iW0 = TreeHU/TreeW0

#%%---Case run 2
#--Variables passed to the main function
do_increasingWTC0 = False #
do_constant_dPsi0 = True # False #
W0 = np.zeros(nYears)
dDBH = np.zeros(nYears)
TreeD = np.zeros(nYears)
Acrown = np.zeros(nYears)
Atrunk = np.zeros(nYears) # Cross sectional area
Asw = np.zeros(nYears) # area of sapwood
Ktrunk = np.zeros(nYears)
TreeHU = np.zeros(nYears) # Whole tree hydraulic usage
TreeW0 = np.zeros(nYears) # Whole tree WTC potential
AgeSW = np.zeros(nYears)
Transp = np.zeros(nYears) # transp/crown area
fareaM = np.zeros((nYears,nYears))
#fareaM[:,:] = np.nan

PlantHydraulicModel()

#%% Keep results 2
W0_cW0     = W0.copy()
dDBH_cW0   = dDBH.copy()
TreeD_cW0  = TreeD.copy()
Acrown_cW0 = Acrown.copy()
Ktrunk_cW0 = Ktrunk.copy()
Atrunk_cW0 = Atrunk.copy()
Asw_cW0    = Asw.copy()
TreeHU_cW0 = TreeHU.copy()
TreeW0_cW0 = TreeW0.copy()
AgeSW_cW0  = AgeSW.copy()
Transp_cW0 = Transp.copy()
fareaM_cW0 = fareaM.copy()

#--New calculations
TrunkRadius_cW0 = 0.5 * TreeD
Hv_cW0   = 10000*Asw/(Acrown*3.5)
WTC_cW0  = (TreeW0-TreeHU)/(Acrown*Transp)
r_SW_cW0 = Asw/Atrunk
r_HU_cW0 = TreeHU/TreeW0

#%% ============ Plotting for two WTC0 schemes ==========
r=TrunkRadius_iW0
plt.figure(5)
plt.rc('axes', prop_cycle=(cycler('color', ['b', 'y', 'g', 'm', 'r'])+
    cycler('linestyle', ['-', '-.', '--', ':',(0, (3, 1, 1, 1, 1, 1))])))
plt.clf()

plt.subplot(3,3,1)
plt.plot(t, W0_iW0/1000, t, W0_cW0/1000)
plt.title('a', x=0.1, y=0.8,fontdict=font)
plt.ylabel('$WTC_{0}$(km)', fontdict=font)
#plt.xlabel('Tree age (year)', fontdict=font)

plt.subplot(3,3,2)
plt.plot(t,fareaM_iW0[20,:],t,fareaM_iW0[120,:]+0.01,t,fareaM_iW0[220,:]+0.02, t,fareaM_iW0[320,:]+0.03)
plt.legend( ('Year 20', 'Year 120', 'Year 220', 'Year 320'), \
              ncol=1, loc='lower right', \
              bbox_to_anchor=(0.8, 0.05) )
plt.title('b', x=0.1, y=0.8,fontdict=font)
plt.ylim((0,1.2))
plt.ylabel('$f_{function}$', fontdict=font)
#plt.xlabel('Tree age (year)', fontdict=font)

plt.subplot(3,3,3)
plt.plot(t,fareaM_cW0[20,:],t,fareaM_cW0[120,:]+0.01,t,fareaM_cW0[220,:]+0.02, t,fareaM[320,:]+0.03)
plt.legend( ('Year 20', 'Year 120', 'Year 220', 'Year 320'), \
              ncol=1, loc='lower right', \
              bbox_to_anchor=(0.8, 0.05) )
plt.title('c', x=0.1, y=0.8,fontdict=font)
plt.ylim((0,1.2))
plt.ylabel('$f_{function}$', fontdict=font)
#plt.xlabel('Tree age (year)', fontdict=font)

plt.subplot(3,3,4)
plt.plot(t, Ktrunk_iW0,t, Ktrunk_cW0)
plt.title('d', x=0.1, y=0.8,fontdict=font)
plt.ylabel('K$_{trunk}$\n(Mg $yr^{-1} MPa^{-1}$)', fontdict=font)
#plt.xlabel('Tree age (year)', fontdict=font)

plt.subplot(3,3,5)
plt.plot(t, r_SW_iW0,t, r_SW_cW0)
#plt.plot(t, x)
plt.ylim((0,1.2))
plt.title('e', x=0.1, y=0.8,fontdict=font)
#plt.legend( ('r=2.8'))
plt.tick_params(labelbottom='off')
plt.ylabel('$r_{SW/Trunk}$', fontdict=font)
#plt.xlabel('Tree age (year)', fontdict=font)

plt.subplot(3,3,6)
#plt.plot(t, WTC_iW0,t, WTC_cW0)
plt.plot(t, r_HU_iW0,t, r_HU_cW0)
#plt.plot(t, x)
#plt.ylim((0,1.2))
plt.title('f', x=0.1, y=0.8,fontdict=font)
#plt.legend( ('r=2.8'))
plt.tick_params(labelbottom='off')
plt.ylabel('$r_{HU/W0}$', fontdict=font)
#plt.xlabel('Tree age (year)', fontdict=font)

plt.subplot(3,3,7)
plt.plot(t, Transp_iW0,t, Transp_cW0)
plt.title('g', x=0.1, y=0.8,fontdict=font)
plt.ylim((0.0,1.5))
plt.ylabel('Transp(m $yr^{-1}$)', fontdict=font)
plt.xlabel('Tree age (year)', fontdict=font)

dDBH_iW0[0]=dDBH_iW0[1]
dDBH_cW0[0]=dDBH_cW0[1]
plt.subplot(3,3,8)
plt.plot(t, dDBH_iW0*1000.0, t, dDBH_cW0*1000.0)
plt.ylim((0.0,10.0))
plt.title('h', x=0.1, y=0.8,fontdict=font)
plt.ylabel('$dD$ (mm $yr^{-1}$)', fontdict=font)
plt.xlabel('Tree age (year)', fontdict=font)


plt.subplot(3,3,9)
plt.plot(t, Hv_iW0, t, Hv_cW0)
#plt.ylim((0.0,10.0))
plt.title('j', x=0.1, y=0.8,fontdict=font)
plt.ylabel('Hv ($cm^{2}$ $m^{-2}$)', fontdict=font)
plt.xlabel('Tree age (year)', fontdict=font)
