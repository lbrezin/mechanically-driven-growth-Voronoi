#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 09:16:41 2024

@author: louisbrezin
"""

import reader
import readFile
import math
import analysis
import numpy as np
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import collections
import statistics as stat
import tkinter as tk
from tkinter import filedialog as fd
import pickle
import seaborn as sns
from scipy import integrate
from scipy.stats import gaussian_kde
from numpy import ma
plt.rcParams["figure.autolayout"] = True
plt.rcParams['text.usetex'] = True
plt.style.use('seaborn-poster')
#plt.rcParams['figure.figsize'] = [10, 10]

#%% Read the files for pressure
root = tk.Tk()
root.withdraw()
filenamesPressure = fd.askopenfilenames();
NfileP = len(filenamesPressure);
DataP=[];
for i in range(NfileP):
    DataP.append(np.loadtxt(filenamesPressure[i]));
    
#%% Read the files for the parameters
rootP = tk.Tk()
rootP.withdraw()
filenamesParams = fd.askopenfilenames();
ParamsP=[];
for i in range(NfileP):
    ParamsP.append(reader.ReadParams(filenamesParams[i]));
#%% Extract nonzero elements
DataP=np.array(DataP);
#DataMask = DataP!=0;
ExtractDataP = ma.masked_equal(DataP,0);

#%% Defines the number of files that are in the solid state
solidcycle=50;

#Nfilecycle = NfileP
#A0 = np.zeros(NfileP);
#v0 = np.zeros(NfileP);

#%% Extract the homeostatic pressure as the last value of each run, and extract the paramater associated to it

P0 = np.zeros(NfileP);
HomeoPP = np.zeros(NfileP);

#%%

K0s = np.zeros(solidcycle);
K0l = np.zeros(NfileP-solidcycle);
HomeoPsK = np.zeros(solidcycle);
HomeoPlK = np.zeros(NfileP-solidcycle);
#%%
Gs = np.zeros(solidcycle);
Gl = np.zeros(NfileP-solidcycle);
HomeoPsG = np.zeros(solidcycle);
HomeoPlG = np.zeros(NfileP-solidcycle);
#%%
v0s = np.zeros(solidcycle);
v0l = np.zeros(NfileP-solidcycle);
HomeoPsV = np.zeros(solidcycle);
HomeoPlV = np.zeros(NfileP-solidcycle);
#%%
a0s = np.zeros(solidcycle);
a0l = np.zeros(NfileP-solidcycle);
HomeoPsA = np.zeros(solidcycle);
HomeoPlA = np.zeros(NfileP-solidcycle);
#%%
Drs = np.zeros(solidcycle);
Drl = np.zeros(NfileP-solidcycle);
HomeoPsD = np.zeros(solidcycle);
HomeoPlD = np.zeros(NfileP-solidcycle);

#%%
for i in range(NfileP):
    j=0;
    while(DataP[i,j,1]!=0 and j<np.shape(DataP)[1]-1):
        j = j+1;
    if(i<solidcycle):
        #K0s[i] = ParamsP[i][1][5];
        #HomeoPsK[i] = DataP[i,j-1,0];
        #v0s[i] = ParamsP[i][1][9];
        #HomeoPsV[i] = DataP[i,j-1,0];
        Drs[i] = ParamsP[i][1][10];
        HomeoPsD[i] = DataP[i,j-1,0];
        #Gs[i] = ParamsP[i][1][2];
        #HomeoPsG[i] = DataP[i,j-1,0];
        #a0s[i] = ParamsP[i][1][4];
        #HomeoPsA[i] = DataP[i,j-1,0];
    else:
        #K0l[i-solidcycle] = ParamsP[i][1][5];
        #HomeoPlK[i-solidcycle] = DataP[i,j-1,0];
        #v0l[i-solidcycle] = ParamsP[i][1][9];
        #HomeoPlV[i-solidcycle] = DataP[i,j-1,0];
        Drl[i-solidcycle] = ParamsP[i][1][10];
        HomeoPlD[i-solidcycle] = DataP[i,j-1,0];
        #Gl[i-solidcycle] = ParamsP[i][1][2];
        #HomeoPlG[i-solidcycle] = DataP[i,j-1,0];
        #a0l[i-solidcycle] = ParamsP[i][1][4];
        #HomeoPlA[i-solidcycle] = DataP[i,j-1,0];
    #P0[i] = ParamsP[i][1][0]; #For the change in preferred perimeter
    #HomeoPP[i] = DataP[i,j-1,0];
    #A0[i] = ParamsP[i][1][4]; #For the change in preferred area
    
    #v0[i] = ParamsP[i][1][9];
#%%Plot homeostatic pressure
#x = np.linspace(0.5, 1.5,11);# area
#x = np.linspace(3.2, 4.2,11);#perimeter
#x = np.linspace(0.6, 2.6,11); #perimeter elasticity
#x=[0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1,1.2,1.3,1.5];
#plt.scatter(A0,HomeoP)
#plt.scatter(P0,HomeoP,marker='o',color='blue')
#plt.axvline(x=3.81,color='blue',linewidth=0.75)
#plt.scatter(v0,HomeoP)
sizemarker = 40;
mliquid = 'v';
msolid = 'o';
fig,axs=plt.subplots(2,3, sharex=False,sharey=True)
axs[0,0].scatter(a0s,HomeoPsA,s=sizemarker,marker=msolid,label='solid',color='red',facecolor='none');
axs[0,0].scatter(a0l[:],HomeoPlA[:],s=sizemarker,marker=mliquid,label='liquid',color='red',facecolor='none');
axs[0,0].set_title('Preferred area');
axs[0,0].legend()
axs[0,0].set_ylabel(r'$P_h$', fontsize='xx-large')
axs[0,1].scatter(P0,HomeoPP,s=sizemarker,marker=msolid,color='blue',facecolor='none');
axs[0,1].axvline(x=3.81,color='blue',linewidth=0.75)
axs[0,1].set_title('Preferred perimeter');
axs[0,2].scatter(Gs,HomeoPsG,s=sizemarker,marker=msolid,label='solid',color='g',facecolor='none');
axs[0,2].scatter(Gl[:],HomeoPlG[:],s=sizemarker,marker=mliquid,label='liquid',color='g',facecolor='none');
axs[0,2].set_title('Perimeter elasticity');
#axs[0,2].legend()
axs[1,0].scatter(K0s,HomeoPsK,s=sizemarker,marker=msolid,label='solid',color='g',facecolor='none');
axs[1,0].scatter(K0l[:],HomeoPlK[:],s=sizemarker,marker=mliquid,label='liquid',color='g',facecolor='none');
axs[1,0].set_title('Area elasticity');
axs[1,0].set_ylabel(r'$P_h$', fontsize='xx-large')
#axs[1,0].legend()
axs[1,1].scatter(v0s,HomeoPsV,s=sizemarker,marker=msolid,label='solid',color='k',facecolor='none');
axs[1,1].scatter(v0l[:],HomeoPlV[:],s=sizemarker,marker=mliquid,label='liquid',color='k',facecolor='none');
axs[1,1].set_title('Self-advection');
#axs[1,1].legend()
axs[1,2].scatter(Drs,HomeoPsD,s=sizemarker,marker=msolid,label='solid',color='brown',facecolor='none');
axs[1,2].scatter(Drl[:],HomeoPlD[:],s=sizemarker,marker=mliquid,label='liquid',color='brown',facecolor='none');
axs[1,2].set_title('Rotational diffusion');
#axs[1,2].legend()
#axs[1,:].set_xlabel(r'$a_0$', fontsize='xx-large')



#plt.xlabel(r'$a_0$', fontsize='xx-large')
#plt.ylabel(r'$P_h$', fontsize='xx-large')
#plt.title('Rotational diffusion')
#plt.legend()

#%% Saving files
VariedParam = [a0s,a0l,P0,Gs,Gl,K0s,K0l,v0s,v0l,Drs,Drl];
HomeoRes = [HomeoPsA,HomeoPlA,HomeoPP,HomeoPsG,HomeoPlG,HomeoPsK,HomeoPlK,HomeoPsV,HomeoPlV,HomeoPsD,HomeoPlD];
#%% Dumping on a file
fileVariedParam = open('./Extracted_Data/VariedParamHomeoTable','w+b');
fileHomeoRes = open('./Extracted_Data/HomeoResTable','w+b');
#fileData = open('./Extracted_Data/Data_100exp_alpha03-2_p0_265-465','rb');
#fileParams = open('./Extracted_Data/Params_100exp_alpha03-2_p0_265-465','rb');
pickle.dump(VariedParam,fileVariedParam);
pickle.dump(HomeoRes,fileHomeoRes);