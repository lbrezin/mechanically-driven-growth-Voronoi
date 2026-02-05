#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A script to plot the figures of the Voronoi paper
Created on Wed Sep  4 16:50:41 2024

@author: louisbrezin
"""
#%% Load packages
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation as an
from matplotlib.ticker import MultipleLocator,FormatStrFormatter,MaxNLocator
import reader
import readFile
import tkinter as tk
from tkinter import filedialog as fd
import numpy as np
import analysis
from scipy.stats import gaussian_kde
import copy
import pickle
import time
import seaborn as sns
from mpltools import annotation
from matplotlib.gridspec import GridSpec
from scipy.stats import pearsonr
from scipy.stats import linregress
import matplotlib.colors as mcolors



    #%% Basic styling of plots
#plt.rcParams["figure.autolayout"] = True
_new_black = '#373737'
params = {
   'axes.labelsize': 9,
   'font.size': 9,
   'legend.fontsize': 9,
   'xtick.labelsize': 9,
   'ytick.labelsize': 9,
   'lines.markersize':2,
   'lines.linewidth': 0.75,
   'lines.markeredgewidth': 0.2,
   'markers.fillstyle' : 'none',
   'text.usetex': False,
   # Avoid black unless necessary
    'text.color': _new_black,
    'patch.edgecolor': _new_black,
    'patch.force_edgecolor': False, # Seaborn turns on edgecolors for histograms by default and I don't like it
    'hatch.color': _new_black,
    'axes.edgecolor': _new_black,
    # 'axes.titlecolor': _new_black # should fallback to text.color
    'axes.labelcolor': _new_black,
    'xtick.color': _new_black,
    'ytick.color': _new_black
   }
plt.style.use('default')
plt.rcParams.update(params)
#plt.style.use('seaborn-poster')
colorsolid = 'orange';
colorliq ='blue';
colorresident = 'grey';
colormutant= 'red';
#%% Load data fig1
file1a = open('./DataFigures/HomogeneousRef','rb');
Fig1a = pickle.load(file1a);
file1a.close()

file1aDivL = open('./DataFigures/HomogeneousRefLiquid-Division','rb')
Fig1aDivL = pickle.load(file1aDivL);
file1aDivL.close()

file1aDivS = open('./DataFigures/HomogeneousRefSolid-Division','rb')
Fig1aDivS = pickle.load(file1aDivS);
file1aDivS.close()
#%% Load data Fig 2
file2 = open('./DataFigures/MSD_20240902','rb');
file2b = open('./DataFigures/TimeMSD_20240902','rb');
MSDAllFiles = pickle.load(file2);
T = pickle.load(file2b)
file2.close()
file2b.close()


#%% Figure 1
"""
Plot of the state of the system in solid and liquid case (top line),
distribution of area, perimeter and pressure(bottom line)

"""
#fig1, axes1 = plt.subplots(nrows=2,ncols=3,figsize=(6,6),dpi=600)#,gridspec_kw={'wspace':0.5,'hspace':0.5},dpi=600)
fig1 = plt.figure(figsize=(6.5,3),dpi=600)
subfigs = fig1.subfigures(2, 1, hspace=0.2)
axes0 = subfigs[0].subplots(1,3)
reader.MakePlot(len(Fig1a[0])-1, Fig1a[0], axes0[1])
axes0[1].set_title('b)',fontfamily='serif',loc='left')

reader.MakePlot(len(Fig1a[1])-1, Fig1a[1], axes0[2])
axes0[2].set_title('c)',fontfamily='serif',loc='left')

axes0[0].set_xscale('log')
axes0[0].set_yscale('log')

legends= ['Liquid w/ division', 'Liquid w/o division', 'Solid w/ division', 'Solid w/o division']
colors = [colorliq, colorliq, colorsolid, colorsolid]
styles = ['dashed','solid','dashed','solid']
Nt=len(T);
start = Nt//2;
#ax2.plot(np.log(T[start:Nt]-T[start]),np.log(MSD[:,idc]));
for k in range(len(MSDAllFiles)):
    axes0[0].plot(0.1*(T[start:Nt]-T[start]),MSDAllFiles[k],label=legends[k],linestyle=styles[k],color=colors[k]);
axes0[0].plot(0.1*(T[start+60:Nt]-T[start]), np.power(0.1*(T[start+60:Nt]-T[start]),1),linewidth=1,color=_new_black)
axes0[0].set_xlabel(r'$t*r^0$')
axes0[0].set_ylabel(r'MSD (in $a^0$)')
axes0[0].get_xaxis().set_tick_params(which='minor', size=0)
axes0[0].get_xaxis().set_tick_params(which='minor', width=0) 
annotation.slope_marker((5, 7), 1, ax=axes0[0],invert=True,size_frac=0.15,pad_frac=0.05)
axes0[0].set_title('a)',fontfamily='serif',loc='left')

"""
Plot the distributions of area, perimeter and pressure
"""
d = .5  # proportion of vertical to horizontal extent of the slanted line
kwargs = dict(marker=[(-1, -d), (1, d)], markersize=7,
              linestyle="none", color=_new_black, mec=_new_black, mew=1, clip_on=False)
axes1 = subfigs[1].subplots(2,3,sharey=False,sharex=False)
subfigs[1].subplots_adjust(wspace=0.5,hspace=0.1)
subfigs[1].supylabel('Proba. density')
#% Plotting the distribution of perimeters for liquid and solid
AreaOTimeL, PeriOTimeL, Type = analysis.AvAreaAndPerim(Fig1a[1]);
AreaOTimeS, PeriOTimeS, Type = analysis.AvAreaAndPerim(Fig1a[0]);
AreaConcLlist, PeriConcLlist = analysis.AreaAndPerimConcatenated(Fig1aDivL[0])
AreaConcSlist, PeriConcSlist = analysis.AreaAndPerimConcatenated(Fig1aDivS[0])

AreaConcL = analysis.flatten(AreaConcLlist)
AreaConcS = analysis.flatten(AreaConcSlist)
PeriConcL = analysis.flatten(PeriConcLlist)
PeriConcS = analysis.flatten(PeriConcSlist)






#axes1[0].plot(x_perimL, kde_perimL(x_perimL)*bin_width_perimL ,label='Liquid-like',color='blue')
#axes1[0].plot(x_perimS, kde_perimS(x_perimS)*bin_width_perimS ,label='Solid-like',color='y')
sns.kdeplot(ax=axes1[0,0],data=PeriOTimeS,color=colorsolid)
sns.kdeplot(ax=axes1[0,0],data=PeriOTimeL,color=colorliq)
sns.kdeplot(ax=axes1[1,0],data=PeriOTimeS,color=colorsolid)
sns.kdeplot(ax=axes1[1,0],data=PeriOTimeL,color=colorliq)
#sns.kdeplot(ax=axes1[0,0],data=PeriConcS,color=colorsolid)
#sns.kdeplot(ax=axes1[0,0],data=PeriConcL,color=colorliq)
#sns.kdeplot(ax=axes1[1,0],data=PeriConcS,color=colorsolid)
#sns.kdeplot(ax=axes1[1,0],data=PeriConcL,color=colorliq)

axes1[1,0].set_xlabel('perimeter, $p$')
axes1[1,0].set(ylabel=None)
axes1[0,0].set(ylabel=None)
axes1[1,0].set_xticks([3.2,3.3, 3.7, 4,4.3])
axes1[1,0].set_xticklabels(['',r'${p_S^0}$',3.7,r'${p_L^0}$',4.3])
axes1[1,0].set_yticks([0,2])
axes1[0,0].set_yticks([112,114])
axes1[1,0].set_ylim(0,4)
axes1[0,0].set_ylim(110,114)
axes1[0,0].set_title('d)',fontfamily='serif',loc='left')
# hide the spines between ax and ax2
axes1[0,0].spines.bottom.set_visible(False)
axes1[1,0].spines.top.set_visible(False)
axes1[0,0].xaxis.tick_top()
axes1[0,0].tick_params(labeltop=False)  # don't put tick labels at the top
axes1[1,0].xaxis.tick_bottom()
axes1[0,0].plot([0, 1], [0, 0], transform=axes1[0,0].transAxes, **kwargs)
axes1[1,0].plot([0, 1], [1, 1], transform=axes1[1,0].transAxes, **kwargs)



#%Plotting the distribution of areas for liquid and solid

# sns.kdeplot(ax=axes1[0,1],data=AreaOTimeS,color=colorsolid)
# sns.kdeplot(ax=axes1[0,1],data=AreaOTimeL,color=colorliq)
# sns.kdeplot(ax=axes1[1,1],data=AreaOTimeS,color=colorsolid)
# sns.kdeplot(ax=axes1[1,1],data=AreaOTimeL,color=colorliq)

sns.kdeplot(ax=axes1[0,1],data=AreaConcS,color=colorsolid)
sns.kdeplot(ax=axes1[0,1],data=AreaConcL,color=colorliq)
sns.kdeplot(ax=axes1[1,1],data=AreaConcS,color=colorsolid)
sns.kdeplot(ax=axes1[1,1],data=AreaConcL,color=colorliq)

axes1[1,1].set_xlabel('area, $a$')
axes1[1,1].set(ylabel=None)
axes1[0,1].set(ylabel=None)
axes1[1,1].set_xticks([0.7,1, 1.3])
axes1[1,1].set_yticks([0,2])
axes1[0,1].set_yticks([80,85])
axes1[1,1].set_xticklabels([0.7,r'$a^0$', 1.3])
axes1[1,1].set_ylim(0,4)
axes1[0,1].set_ylim(75,85)
axes1[0,1].set_title('e)',fontfamily='serif',loc='left')
# hide the spines between ax and ax2
axes1[0,1].spines.bottom.set_visible(False)
axes1[1,1].spines.top.set_visible(False)
axes1[0,1].xaxis.tick_top()
axes1[0,1].tick_params(labeltop=False)  # don't put tick labels at the top
axes1[1,1].xaxis.tick_bottom()
axes1[0,1].plot([0, 1], [0, 0], transform=axes1[0,1].transAxes, **kwargs)
axes1[1,1].plot([0, 1], [1, 1], transform=axes1[1,1].transAxes, **kwargs)



#%Gather pressure computed in the simulation
PressureOTimeS, Type = analysis.AvPressure(Fig1a[0]);
PressureOTimeL, Type = analysis.AvPressure(Fig1a[1]);

sns.kdeplot(ax=axes1[0,2],data=PressureOTimeS,color=colorsolid)
sns.kdeplot(ax=axes1[0,2],data=PressureOTimeL,color=colorliq)
sns.kdeplot(ax=axes1[1,2],data=PressureOTimeS,color=colorsolid)
sns.kdeplot(ax=axes1[1,2],data=PressureOTimeL,color=colorliq)
#axes1[2].set_title('Pressure')
axes1[1,2].set_xlabel('pressure, $\mathcal{P}$')
axes1[1,2].set(ylabel=None)
axes1[0,2].set(ylabel=None)
axes1[1,2].set_ylim(0,4)
axes1[1,2].set_xlim(-2,0.5)
axes1[0,2].set_xlim(-2,0.5)
axes1[0,2].set_ylim(125,130)
axes1[1,2].set_yticks([0,2])
axes1[0,2].set_yticks([127,130])
axes1[0,2].set_title('f)',fontfamily='serif',loc='left')
# hide the spines between ax and ax2
axes1[0,2].spines.bottom.set_visible(False)
axes1[1,2].spines.top.set_visible(False)
axes1[0,2].xaxis.tick_top()
axes1[0,2].tick_params(labeltop=False)  # don't put tick labels at the top
axes1[1,2].xaxis.tick_bottom()
axes1[0,2].plot([0, 1], [0, 0], transform=axes1[0,2].transAxes, **kwargs)
axes1[1,2].plot([0, 1], [1, 1], transform=axes1[1,2].transAxes, **kwargs)

#fig1.savefig('Fig1-3x2.pdf',format='pdf',bbox_inches='tight')


#%% Figure 1 with all the quantities for growing homogeneous system
"""
Plot of the state of the system in solid and liquid case (top line),
distribution of area, perimeter and pressure(bottom line)

"""
#fig1, axes1 = plt.subplots(nrows=2,ncols=3,figsize=(6,6),dpi=600)#,gridspec_kw={'wspace':0.5,'hspace':0.5},dpi=600)
fig1 = plt.figure(figsize=(6.5,3),dpi=600)
subfigs = fig1.subfigures(2, 1, hspace=0.2)
axes0 = subfigs[0].subplots(1,3)
reader.MakePlot(len(Fig1a[0])-1, Fig1a[0], axes0[1])
axes0[1].set_title('b) Solid phase',fontfamily='serif',loc='left')

reader.MakePlot(len(Fig1a[1])-1, Fig1a[1], axes0[2])
axes0[2].set_title('c) Liquid phase',fontfamily='serif',loc='left')

axes0[0].set_xscale('log')
axes0[0].set_yscale('log')

legends= ['Liquid w/ division', 'Liquid w/o division', 'Solid w/ division', 'Solid w/o division']
colors = [colorliq, colorliq, colorsolid, colorsolid]
styles = ['dashed','solid','dashed','solid']
Nt=len(T);
start = Nt//2;
axes0[0].plot(0.1*(T[start:Nt]-T[start]),MSDAllFiles[0],label='Liquid',linestyle=styles[1],color=colors[0]);
axes0[0].plot(0.1*(T[start:Nt]-T[start]),MSDAllFiles[2],label='Solid',linestyle=styles[1],color=colors[2]);
axes0[0].plot(0.1*(T[start+60:Nt]-T[start]), np.power(0.1*(T[start+60:Nt]-T[start]),1),linewidth=1,color=_new_black)
axes0[0].set_xlabel(r'$r^0 t$')
axes0[0].legend(loc='best')
axes0[0].set_ylabel(r'MSD (in $a^0$)')
axes0[0].get_xaxis().set_tick_params(which='minor', size=0)
axes0[0].get_xaxis().set_tick_params(which='minor', width=0) 
annotation.slope_marker((5, 7), 1, ax=axes0[0],invert=True,size_frac=0.15,pad_frac=0.05)
axes0[0].set_title('a) Fluidization ',fontfamily='serif',loc='left')

"""
Plot the distributions of area, perimeter and pressure
"""
d = .5  # proportion of vertical to horizontal extent of the slanted line
kwargs = dict(marker=[(-1, -d), (1, d)], markersize=7,
              linestyle="none", color=_new_black, mec=_new_black, mew=1, clip_on=False)
axes1 = subfigs[1].subplots(1,3,sharey=False,sharex=False)
subfigs[1].subplots_adjust(wspace=0.5,hspace=0.1)
subfigs[1].supylabel('Proba. density')
AreaConcLlist, PeriConcLlist = analysis.AreaAndPerimConcatenated(Fig1aDivL[0])
AreaConcSlist, PeriConcSlist = analysis.AreaAndPerimConcatenated(Fig1aDivS[0])

AreaConcL = analysis.flatten(AreaConcLlist)
AreaConcS = analysis.flatten(AreaConcSlist)
PeriConcL = analysis.flatten(PeriConcLlist)
PeriConcS = analysis.flatten(PeriConcSlist)



sns.kdeplot(ax=axes1[0],data=PeriConcS,color=colorsolid)
sns.kdeplot(ax=axes1[0],data=PeriConcL,color=colorliq)


axes1[0].set_xlabel('perimeter, $p$')
axes1[0].set(ylabel=None)
axes1[0].set(ylabel=None)
axes1[0].set_xticks([2,3,4,5])
#axes1[0].set_xticklabels(['',r'${p_S^0}$',3.7,r'${p_L^0}$',4.3])
axes1[0].set_yticks([0,5,10,15,20])
axes1[0].axvline(x=4,color=colorliq,linestyle='dashed',linewidth=0.3);
axes1[0].axvline(x=3.3,color=colorsolid,linestyle='dashed',linewidth=0.3);
axes1[0].set_ylim(0,20)
axes1[0].set_xlim(2,5)
axes1[0].set_title('d)',fontfamily='serif',loc='left')
# hide the spines between ax and ax2
#axes1[0,0].spines.bottom.set_visible(False)
#axes1[1,0].spines.top.set_visible(False)
#axes1[0,0].xaxis.tick_top()
#axes1[0,0].tick_params(labeltop=False)  # don't put tick labels at the top
#axes1[1,0].xaxis.tick_bottom()
#axes1[0,0].plot([0, 1], [0, 0], transform=axes1[0,0].transAxes, **kwargs)
#axes1[1,0].plot([0, 1], [1, 1], transform=axes1[1,0].transAxes, **kwargs)



#%Plotting the distribution of areas for liquid and solid

# sns.kdeplot(ax=axes1[0,1],data=AreaOTimeS,color=colorsolid)
# sns.kdeplot(ax=axes1[0,1],data=AreaOTimeL,color=colorliq)
# sns.kdeplot(ax=axes1[1,1],data=AreaOTimeS,color=colorsolid)
# sns.kdeplot(ax=axes1[1,1],data=AreaOTimeL,color=colorliq)

sns.kdeplot(ax=axes1[1],data=AreaConcS,color=colorsolid)
sns.kdeplot(ax=axes1[1],data=AreaConcL,color=colorliq)


axes1[1].set_xlabel('area, $a$')
axes1[1].set(ylabel=None)
axes1[1].set(ylabel=None)
axes1[1].set_xticks([0.5,1, 1.5])
axes1[1].set_yticks([0,5,10])
#axes1[1].set_xticklabels([0.7,r'$a^0$', 1.3])
axes1[1].set_ylim(0,10)
axes1[1].set_xlim(0.5,1.5)
axes1[1].axvline(x=1,color=_new_black,linestyle='dashed',linewidth=0.3);
axes1[1].set_title('e)',fontfamily='serif',loc='left')
# hide the spines between ax and ax2
#axes1[0,1].spines.bottom.set_visible(False)
#axes1[1,1].spines.top.set_visible(False)
#axes1[0,1].xaxis.tick_top()
#axes1[0,1].tick_params(labeltop=False)  # don't put tick labels at the top
#axes1[1,1].xaxis.tick_bottom()
#axes1[0,1].plot([0, 1], [0, 0], transform=axes1[0,1].transAxes, **kwargs)
#axes1[1,1].plot([0, 1], [1, 1], transform=axes1[1,1].transAxes, **kwargs)



#%Gather pressure computed in the simulation
PressureConcLlist = analysis.PressureConcatenated(Fig1aDivL[0])
PressureConcSlist = analysis.PressureConcatenated(Fig1aDivS[0])

PressureConcL = analysis.flatten(PressureConcLlist)
PressureConcS = analysis.flatten(PressureConcSlist)

sns.kdeplot(ax=axes1[2],data=PressureConcS,color=colorsolid)
sns.kdeplot(ax=axes1[2],data=PressureConcL,color=colorliq)

#axes1[2].set_title('Pressure')
axes1[2].set_xlabel('pressure, $\mathcal{P}$')
axes1[2].set(ylabel=None)
axes1[2].set(ylabel=None)
axes1[2].set_ylim(0,10)
axes1[2].set_xlim(-2,0.5)
axes1[2].set_yticks([0,5,10])
axes1[2].set_title('f)',fontfamily='serif',loc='left')
# hide the spines between ax and ax2
# axes1[0,2].spines.bottom.set_visible(False)
# axes1[1,2].spines.top.set_visible(False)
# axes1[0,2].xaxis.tick_top()
# axes1[0,2].tick_params(labeltop=False)  # don't put tick labels at the top
# axes1[1,2].xaxis.tick_bottom()
# axes1[0,2].plot([0, 1], [0, 0], transform=axes1[0,2].transAxes, **kwargs)
# axes1[1,2].plot([0, 1], [1, 1], transform=axes1[1,2].transAxes, **kwargs)

# Also put variation of pressure as function of parameters
# axes2 = subfigs[2].subplots(1,3,sharey=False,sharex=False)
# subfigs[2].subplots_adjust(wspace=0.5,hspace=0.2)
# subfigs[2].supylabel(r'pressure, $\mathcal{P}$')

# filecycle=50;

# #Preferred perimeter
# categories = np.where(VariedParamSolid[0*filecycle:1*filecycle-1] < 3.81, 0,1);
# colormap =np.array([colorsolid,colorliq])
# axes2[0].scatter(VariedParamSolid[0*filecycle:1*filecycle-1],HomeoPressureSolid[0*filecycle:1*filecycle-1],color=colormap[categories]);
# axes2[0].set_xlabel(r'preferred perimeter, $p^0$')
# axes2[0].set_title('g)',fontfamily='serif',loc='left')

# #perimeter elasticity
# axes2[1].scatter(VariedParamSolid[1*filecycle:2*filecycle-1],HomeoPressureSolid[1*filecycle:2*filecycle-1],color=colorsolid);
# axes2[1].scatter(VariedParamLiquid[1*filecycle:2*filecycle-1],HomeoPressureLiquid[1*filecycle:2*filecycle-1],color=colorliq);
# axes2[1].set_xlabel(r'perimeter elasticity, $\Gamma$')
# axes2[1].set_title('h)',fontfamily='serif',loc='left')
# #axs[0,2].set_ylim(-2, 0.5)

# #area elasticity
# axes2[2].scatter(VariedParamSolid[2*filecycle:3*filecycle-1],HomeoPressureSolid[2*filecycle:3*filecycle-1],color=colorsolid);
# axes2[2].scatter(VariedParamLiquid[2*filecycle:3*filecycle-1],HomeoPressureLiquid[2*filecycle:3*filecycle-1],color=colorliq);
# axes2[2].set_xlabel(r'area elasticity, $K$')
# axes2[2].set_title('i)',fontfamily='serif',loc='left')

fig1.savefig('Fig1-AllDiv.pdf',format='pdf',bbox_inches='tight')
#%% Load data Fig 2
file2 = open('./DataFigures/MSD_20240902','rb');
file2b = open('./DataFigures/TimeMSD_20240902','rb');
MSDAllFiles = pickle.load(file2);
T = pickle.load(file2b)
file2.close()
file2b.close()

#%% Fig 2 
"""
Plot the MSD of solid and liquid populations, with and without division 
"""
fig2, axes2 = plt.subplots(nrows=1,ncols=1,figsize=(3.42,2),dpi=600)
axes2.set_xscale('log')
axes2.set_yscale('log')

legends= ['Liquid w/ division', 'Liquid w/o division', 'Solid w/ division', 'Solid w/o division']
colors = [colorliq, colorliq, colorsolid, colorsolid]
styles = ['dashed','solid','dashed','solid']
axes2.set_xscale('log')
axes2.set_yscale('log')
Nt=len(T);
start = Nt//2;
#ax2.plot(np.log(T[start:Nt]-T[start]),np.log(MSD[:,idc]));
for k in range(len(MSDAllFiles)):
    axes2.plot(0.1*(T[start:Nt]-T[start]),MSDAllFiles[k],label=legends[k],linestyle=styles[k],color=colors[k]);
axes2.plot(0.1*(T[start+60:Nt]-T[start]), np.power(0.1*(T[start+60:Nt]-T[start]),1),linewidth=1,color=_new_black)
axes2.set_xlabel(r'$t*r^0$')
axes2.set_ylabel(r'MSD (in $a^0$)')
axes2.get_xaxis().set_tick_params(which='minor', size=0)
axes2.get_xaxis().set_tick_params(which='minor', width=0) 
annotation.slope_marker((5, 7), 1, ax=axes2,invert=True,size_frac=0.15,pad_frac=0.05)
#axes2.tick_params()
#axes2.legend(fontsize=11)
fig2.savefig('Fig2.pdf',format='pdf',bbox_inches='tight')

#%% Load data Fig 3
#A0 data
filePressureA0 = open('./DataFigures/PressureA0','rb');
fileVariedParamA0 = open('./DataFigures/VariedParamA0','rb');
HomeoPressureA0 = pickle.load(filePressureA0);
VariedParamA0 = pickle.load(fileVariedParamA0);
filePressureA0.close();
fileVariedParamA0.close();

#Solid data
filePressureSolid = open('./DataFigures/PressureSolid','rb');
fileVariedParamSolid = open('./DataFigures/VariedParamSolid','rb');
HomeoPressureSolid = pickle.load(filePressureSolid);
VariedParamSolid = pickle.load(fileVariedParamSolid);
filePressureSolid.close();
fileVariedParamSolid.close();

#Liquid data
filePressureLiquid = open('./DataFigures/PressureLiquid','rb');
fileVariedParamLiquid = open('./DataFigures/VariedParamLiquid','rb');
HomeoPressureLiquid = pickle.load(filePressureLiquid);
VariedParamLiquid = pickle.load(fileVariedParamLiquid);
filePressureLiquid.close();
fileVariedParamLiquid.close();


#%% Plot Fig 3
"""
Plot a table of pressures as a function of various parameters, in the homogeneous system.
"""
coloring = ['red','blue','g','c','k','brown'];
labeling = ['Preferred area','Preferred perimeter','Perimeter elasticity','Area elasticity','Self-advection','Rotational diffusion'];
filecycle=50;
nline = 2;
ncol = 3;
fig,axs=plt.subplots(nline,ncol,figsize=(7,4),dpi=900,sharex=False,sharey=False)
fig.subplots_adjust(wspace=0.4,hspace=0.7)
#Nexp = int(Nfile/filecycle);
#Preferred area
axs[0,0].scatter(VariedParamA0[0:filecycle-1],HomeoPressureA0[0:filecycle-1],color=colorsolid);
axs[0,0].scatter(VariedParamA0[filecycle::],HomeoPressureA0[filecycle::],color=colorliq);
axs[0,0].set_xlabel(r'preferred area, $a^0$')
axs[0,0].set_ylabel(r'$\mathcal{P}$')
axs[0,0].set_title('a)',fontfamily='serif',loc='left')

#Preferred perimeter
axs[0,1].scatter(VariedParamSolid[0*filecycle:1*filecycle-1],HomeoPressureSolid[0*filecycle:1*filecycle-1],color=_new_black);
axs[0,1].set_xlabel(r'preferred perimeter, $p^0$')
axs[0,1].set_title('b)',fontfamily='serif',loc='left')

#perimeter elasticity
axs[0,2].scatter(VariedParamSolid[1*filecycle:2*filecycle-1],HomeoPressureSolid[1*filecycle:2*filecycle-1],color=colorsolid);
axs[0,2].scatter(VariedParamLiquid[1*filecycle:2*filecycle-1],HomeoPressureLiquid[1*filecycle:2*filecycle-1],color=colorliq);
axs[0,2].set_xlabel(r'perimeter elasticity, $\Gamma$')
axs[0,2].set_title('c)',fontfamily='serif',loc='left')
#axs[0,2].set_ylim(-2, 0.5)

#area elasticity
axs[1,0].scatter(VariedParamSolid[2*filecycle:3*filecycle-1],HomeoPressureSolid[2*filecycle:3*filecycle-1],color=colorsolid);
axs[1,0].scatter(VariedParamLiquid[2*filecycle:3*filecycle-1],HomeoPressureLiquid[2*filecycle:3*filecycle-1],color=colorliq);
axs[1,0].set_xlabel(r'area elasticity, $K$')
axs[1,0].set_ylabel(r'$\mathcal{P}$')
axs[1,0].set_title('d)',fontfamily='serif',loc='left')


#Self-advection
axs[1,1].scatter(VariedParamSolid[3*filecycle:4*filecycle-1],HomeoPressureSolid[3*filecycle:4*filecycle-1],color=colorsolid);
axs[1,1].scatter(VariedParamLiquid[3*filecycle:4*filecycle-1],HomeoPressureLiquid[3*filecycle:4*filecycle-1],color=colorliq);
axs[1,1].set_xlabel(r'self advection, $v^0$')
axs[1,1].set_title('e)',fontfamily='serif',loc='left')

#Rotational Diffusion
axs[1,2].scatter(VariedParamSolid[4*filecycle:5*filecycle-1],HomeoPressureSolid[4*filecycle:5*filecycle-1],color=colorsolid);
axs[1,2].scatter(VariedParamLiquid[4*filecycle:5*filecycle-1],HomeoPressureLiquid[4*filecycle:5*filecycle-1],color=colorliq);
axs[1,2].set_xlabel(r'rotational diffusion, $D_r$')
axs[1,2].set_title('f)',fontfamily='serif',loc='left')

fig.savefig('Fig3.pdf',format='pdf',bbox_inches='tight')

#%% Load data Fig 3 with division
fileHomeoPressure = open('./Extracted_Data/HomeoResTable','rb');
fileVariedParam = open('./Extracted_Data/VariedParamHomeoTable','rb');
HomeoPressure = pickle.load(fileHomeoPressure);
VariedParam = pickle.load(fileVariedParam);
fileHomeoPressure.close();
fileVariedParam.close();

#%% Fig 3 with division
coloring = ['red','blue','g','c','k','brown'];
labeling = ['Preferred area','Preferred perimeter','Perimeter elasticity','Area elasticity','Self-advection','Rotational diffusion'];
nline = 2;
ncol = 3;
fig,axs=plt.subplots(nline,ncol,figsize=(7,4),dpi=900,sharex=False,sharey=False)
fig.subplots_adjust(wspace=0.4,hspace=0.7)

#Preferred area
axs[0,0].scatter(VariedParam[0],HomeoPressure[0],color=colorsolid);
axs[0,0].scatter(VariedParam[1],HomeoPressure[1],color=colorliq);
axs[0,0].set_xlabel(r'preferred area, $a^0$')
axs[0,0].set_ylabel(r'$\mathcal{P}_h$')
axs[0,0].set_title('a)',fontfamily='serif',loc='left')

#Preferred perimeter
categories = np.where(VariedParam[2] < 3.81, 0,1);
colormap =np.array([colorsolid,colorliq])
axs[0,1].scatter(VariedParam[2],HomeoPressure[2],color=colormap[categories]);
axs[0,1].set_xlabel(r'preferred perimeter, $p^0$')
axs[0,1].set_title('b)',fontfamily='serif',loc='left')

#perimeter elasticity
axs[0,2].scatter(VariedParam[3],HomeoPressure[3],color=colorsolid);
axs[0,2].scatter(VariedParam[4],HomeoPressure[4],color=colorliq);
axs[0,2].set_xlabel(r'perimeter elasticity, $\Gamma$')
axs[0,2].set_title('c)',fontfamily='serif',loc='left')
#axs[0,2].set_ylim(-2, 0.5)

#area elasticity
axs[1,0].scatter(VariedParam[5],HomeoPressure[5],color=colorsolid);
axs[1,0].scatter(VariedParam[6],HomeoPressure[6],color=colorliq);
axs[1,0].set_xlabel(r'area elasticity, $K$')
axs[1,0].set_ylabel(r'$\mathcal{P}_h$')
axs[1,0].set_title('d)',fontfamily='serif',loc='left')


#Self-advection
axs[1,1].scatter(VariedParam[7],HomeoPressure[7],color=colorsolid);
axs[1,1].scatter(VariedParam[8],HomeoPressure[8],color=colorliq);
axs[1,1].set_xlabel(r'self advection, $v^0$')
axs[1,1].set_title('e)',fontfamily='serif',loc='left')

#Rotational Diffusion
axs[1,2].scatter(VariedParam[9],HomeoPressure[9],color=colorsolid);
axs[1,2].scatter(VariedParam[10],HomeoPressure[10],color=colorliq);
axs[1,2].set_xlabel(r'rotational diffusion, $D_r$')
axs[1,2].set_title('f)',fontfamily='serif',loc='left')

fig.savefig('Fig3-wDivision.pdf',format='pdf',bbox_inches='tight')
#%% Load data Fig 4 
#The data for the solid case is in the last of those simulation
fileFraction = open('./DataFigures/Fraction-InvasionSolid-newP','rb');
fileTime = open('./DataFigures/Time-InvasionSolid-newP','rb');
frac2Full = pickle.load(fileFraction);
T = pickle.load(fileTime);
fileFractionBR = open('./DataFigures/Fraction-InvasionSolid_BirthRate','rb');
fileTimeBR = open('./DataFigures/Time-InvasionSolid_BirthRate','rb');
frac2FullBR = pickle.load(fileFractionBR);
TBR = pickle.load(fileTimeBR);
#The data for the birth rate is in a separate file 
#The data for the liquid case is from another simulation
fileFractionLiquid = open('./DataFigures/Fraction-InvasionLiquid_NewP_20241113','rb');
fileTimeLiquid = open('./DataFigures/Time-InvasionLiquid_NewP_20241113','rb');
frac2FullLiquid = pickle.load(fileFractionLiquid);
TLiquid = pickle.load(fileTimeLiquid);

#Data for the opposite variables
fileFractionL_opposite = open('./DataFigures/Fraction-InvasionLiquid-PRM-opposite','rb');
fileFractionS_opposite = open('./DataFigures/Fraction-InvasionSolid-PRM-opposite','rb');

fileTimeL_opposite = open('./DataFigures/Time-InvasionLiquid-PRM-opposite','rb');

frac2FullL_opposite = pickle.load(fileFractionL_opposite);
frac2FullS_opposite = pickle.load(fileFractionS_opposite);


T_opposite = pickle.load(fileTimeL_opposite);


fileFractionL_opposite.close();
fileFractionS_opposite.close();
fileTimeL_opposite.close();


#%% Removes the zeros and 1 from fraction to not trigger bugs. Puts them at unpysically small values
frac2Full=np.where(frac2Full==0,0.00001,frac2Full)
frac2FullBR=np.where(frac2FullBR==0,0.00001,frac2FullBR)
frac2FullLiquid=np.where(frac2FullLiquid==0,0.00001,frac2FullLiquid)
fileFractionL_opposite = np.where(fileFractionL_opposite==0,0.00001,fileFractionL_opposite)
frac2Full=np.where(frac2Full==1,0.99999,frac2Full)
frac2FullBR=np.where(frac2FullBR==1,0.99999,frac2FullBR)
frac2FullLiquid=np.where(frac2FullLiquid==1,0.99999,frac2FullLiquid)
fileFractionS_opposite = np.where(fileFractionS_opposite==0,0.00001,fileFractionS_opposite)
#%% Plot Fig 4
coloring = ['grey','red','blue','g','c','m','k','brown'];
coloringbis = [colorliq,colorsolid];
labeling = [r'$b^0$',r'$a^0$',r'$p^0$',r'$\Gamma$',r'$K$',r'$\lambda$',r'$v^0$',r'$D_r$'];
state = ['Liquid', 'Solid'];
listyle = ['dotted','solid'];
label = ['a)','b)','c)','d)','e)','f)','g)','h)']
#labeling = ['Perimeter elasticity','Area elasticity','Pressure strength','Self-advection','Rotational diffusion'];
filecycle=10;
macrocycle=80;
nline = 2;
ncol = 4;
fig4,axs4=plt.subplots(nline,ncol, figsize=(7,3.5), dpi=900, sharex=True,sharey=True)
fig4.subplots_adjust(hspace=0.4)
tt=0;

state_legend = [
    plt.Line2D([0], [0], color=colorsolid, label='Solid state'),
    plt.Line2D([0], [0], color=colorliq, label='Liquid state')
]

# Legend for parameter change (linestyle)
change_legend = [
    plt.Line2D([0], [0], color= _new_black, linestyle='-', label='↑ param'),
    plt.Line2D([0], [0], color= _new_black, linestyle='--', label='↓ param')
]


for i in range(frac2Full.shape[1]):
    j = i%macrocycle;
    st = i//macrocycle;
    exp = j//filecycle;
    if j%filecycle==0:
        frac2FullAverageL = np.mean(frac2FullLiquid[:,st*macrocycle+filecycle*exp:st*macrocycle+filecycle*exp+filecycle-1],axis=1)
        #std = np.std(frac2Full[:,st*macrocycle+filecycle*exp:st*macrocycle+filecycle*exp+filecycle-1],axis=1)
        #y = np.divide(frac2Full[:,i],1-frac2Full[:,i]);
        #cond = y>10^(-3);
        if exp==1:
            axs4[exp//ncol,exp%ncol].plot(0.01*TLiquid,np.divide(frac2FullAverageL,1-frac2FullAverageL),color=colorliq,linestyle=(0, (5, 5)));
        else:
            axs4[exp//ncol,exp%ncol].plot(0.01*TLiquid,np.divide(frac2FullAverageL,1-frac2FullAverageL),color=colorliq);
        #Logitfrac2FullAverageL = np.mean(np.log(np.divide(frac2FullLiquid[:,st*macrocycle+filecycle*exp:st*macrocycle+filecycle*exp+filecycle-1],1-frac2FullLiquid[:,st*macrocycle+filecycle*exp:st*macrocycle+filecycle*exp+filecycle-1])),axis=1)
        
        #axs4[exp//ncol,exp%ncol].plot(0.01*TLiquid,Logitfrac2FullAverageL,color=colorliq);
        #Treat the specific case of the Birth rate in seperate file
        if exp==0:
            frac2FullAverageS = np.mean(frac2FullBR[:,st*macrocycle+filecycle*exp:st*macrocycle+filecycle*exp+filecycle-1],axis=1)
            axs4[exp//ncol,exp%ncol].plot(0.01*TBR[0:501],np.divide(frac2FullAverageS[0:501],1-frac2FullAverageS[0:501]),color=colorsolid);
            if st == 0:
                frac2FullAverageS_opposite = np.mean(frac2FullS_opposite[:,filecycle*exp:filecycle*exp+filecycle-1],axis=1)
                frac2FullAverageL_opposite = np.mean(frac2FullL_opposite[:,filecycle*exp:filecycle*exp+filecycle-1],axis=1)
                axs4[exp//ncol,exp%ncol].plot(0.01*T_opposite[0:501],np.divide(frac2FullAverageS_opposite[0:501],1-frac2FullAverageS_opposite[0:501]),color=colorsolid,linestyle=(0, (5, 5)));
                axs4[exp//ncol,exp%ncol].plot(0.01*T_opposite[0:501],np.divide(frac2FullAverageL_opposite[0:501],1-frac2FullAverageL_opposite[0:501]),color=colorliq,linestyle=(0, (5, 5)));
                #axs4[exp//ncol,exp%ncol].legend();
            #frac2FullAverageS = np.mean(np.log(np.divide(frac2FullBR[:,st*macrocycle+filecycle*exp:st*macrocycle+filecycle*exp+filecycle-1],1-frac2FullBR[:,st*macrocycle+filecycle*exp:st*macrocycle+filecycle*exp+filecycle-1])),axis=1)
            #axs4[exp//ncol,exp%ncol].plot(0.01*TBR[0:501],frac2FullAverageS[0:501],color=colorsolid);
        #Treat the case of the preferred area which has opposite variations compared to other parameters
        elif exp==1:
            frac2FullAverageS = np.mean(frac2Full[:,st*macrocycle+filecycle*exp:st*macrocycle+filecycle*exp+filecycle-1],axis=1)
            axs4[exp//ncol,exp%ncol].plot(0.01*T[0:102],np.divide(frac2FullAverageS[0:102],1-frac2FullAverageS[0:102]),color=colorsolid,linestyle=(0, (5, 5)));
            if st==0:
                frac2FullAverageS_opposite = np.mean(frac2FullS_opposite[:,filecycle*exp:filecycle*exp+filecycle-1],axis=1)
                frac2FullAverageL_opposite = np.mean(frac2FullL_opposite[:,filecycle*exp:filecycle*exp+filecycle-1],axis=1)
                axs4[exp//ncol,exp%ncol].plot(0.01*T_opposite[0:501],np.divide(frac2FullAverageS_opposite[0:501],1-frac2FullAverageS_opposite[0:501]),color=colorsolid);
                axs4[exp//ncol,exp%ncol].plot(0.01*T_opposite[0:501],np.divide(frac2FullAverageL_opposite[0:501],1-frac2FullAverageL_opposite[0:501]),color=colorliq);
        #All the other plots
        else:
            frac2FullAverageS = np.mean(frac2Full[:,st*macrocycle+filecycle*exp:st*macrocycle+filecycle*exp+filecycle-1],axis=1)
            axs4[exp//ncol,exp%ncol].plot(0.01*T[0:102],np.divide(frac2FullAverageS[0:102],1-frac2FullAverageS[0:102]),color=colorsolid);
            #frac2FullAverageS = np.mean(np.log(np.divide(frac2Full[:,st*macrocycle+filecycle*exp:st*macrocycle+filecycle*exp+filecycle-1],1-frac2Full[:,st*macrocycle+filecycle*exp:st*macrocycle+filecycle*exp+filecycle-1])),axis=1)
            #axs4[exp//ncol,exp%ncol].plot(0.01*T[0:102],frac2FullAverageS[0:102],color=colorsolid);
            # Treat the case of the opposite plots, where 
            if st == 0:
                frac2FullAverageS_opposite = np.mean(frac2FullS_opposite[:,filecycle*exp:filecycle*exp+filecycle-1],axis=1)
                frac2FullAverageL_opposite = np.mean(frac2FullL_opposite[:,filecycle*exp:filecycle*exp+filecycle-1],axis=1)
                axs4[exp//ncol,exp%ncol].plot(0.01*T_opposite[0:501],np.divide(frac2FullAverageS_opposite[0:501],1-frac2FullAverageS_opposite[0:501]),color=colorsolid,linestyle=(0, (5, 5)));
                axs4[exp//ncol,exp%ncol].plot(0.01*T_opposite[0:501],np.divide(frac2FullAverageL_opposite[0:501],1-frac2FullAverageL_opposite[0:501]),color=colorliq,linestyle=(0, (5, 5)));
        #axs4[exp//ncol,exp%ncol].fill_between(0.01*T,np.divide(frac2FullAverage,1-frac2FullAverage)-0.5*np.multiply(np.power(np.divide(np.ones(len(frac2FullAverage)),1-frac2FullAverage),2),std),np.divide(frac2FullAverage,1-frac2FullAverage)+0.5*np.multiply(np.power(np.divide(np.ones(len(frac2FullAverage)),1-frac2FullAverage),2),std),color=coloringbis[st],alpha=0.3);
        #axs[exp//ncol,exp%ncol].plot(100*T,y,marker=("+"),color=coloring[exp],alpha=0.5,linewidth=0.5);
        #axs[exp//3,exp%3].set_xscale('log');
        axs4[exp//ncol,exp%ncol].set_yscale('log');
        axs4[exp//ncol,exp%ncol].set_ylim(0.01,1000);
        axs4[exp//ncol,exp%ncol].set_yticks([0.1,10,1000]);
        #axs4[exp//ncol,exp%ncol].set_ylim(-10,10);
        #axs4[exp//ncol,exp%ncol].set_yticks([-10,-5,0,5,10]);
        axs4[exp//ncol,exp%ncol].minorticks_off();
        axs4[exp//ncol,exp%ncol].set_title(labeling[exp]);
        axs4[exp//ncol,exp%ncol].set_title(label[exp],fontfamily='serif',loc='left')
        #axs[exp//ncol,exp%ncol].legend();
        axs4[1,1].legend(handles=state_legend, fontsize='small', frameon=False)
        axs4[1,2].legend(handles=change_legend, fontsize='small', frameon=False)
        axs4[nline-1,exp%ncol].set_xlabel(r'$\tau_g$')
        axs4[exp//ncol,0].set_ylabel(r'$\frac{f}{1-f}$')

#fig4.savefig('Fig4-newP-opposite-20250922.pdf',format='pdf',bbox_inches='tight',pad_inches=0)
#fig4.savefig('Fig4-AverageLog.pdf',format='pdf',bbox_inches='tight')
#fig4.savefig('Fig4-Floor-newP.pdf',format='pdf',bbox_inches='tight')

#%% Load data Fig 5
fileMixLA0 = open('./DataFigures/MixtureLiquidA0-fig5','rb');
fileMixLP0 = open('./DataFigures/MixtureLiquidP0-fig5','rb');
fileMixSP0 = open('./DataFigures/MixtureSolidP0-fig5','rb');
fileMixSA0 = open('./DataFigures/MixtureSolidA0-fig5-redone-200','rb');
MixSA0= pickle.load(fileMixSA0);
MixSP0= pickle.load(fileMixSP0);
MixLA0= pickle.load(fileMixLA0);
MixLP0= pickle.load(fileMixLP0);
fileMixLA0.close();
fileMixLP0.close();
fileMixSP0.close();
fileMixSA0.close();

#%% Some function needed to plot
def PlotFig5(cellarea, cellperimeter,celltype, ax):
    perimeters = analysis.ValuePerType(cellperimeter, celltype);
    areas = analysis.ValuePerType(cellarea, celltype);
    SqPerimeter0 = np.power(perimeters[0],2);
    SqPerimeter1 = np.power(perimeters[1],2);
    ax.scatter(SqPerimeter0,areas[0],color=colorresident, label='Resident');
    ax.scatter(SqPerimeter1,areas[1],color=colormutant, label='Mutant');
    

#%% Plot Fig 5
AreaOTimeSA0, PeriOTimeSA0, TypeSA0 = analysis.AvAreaAndPerim(MixSA0[0]);
AreaOTimeSP0, PeriOTimeSP0, TypeSP0 = analysis.AvAreaAndPerim(MixSP0[0]);
AreaOTimeLA0, PeriOTimeLA0, TypeLA0 = analysis.AvAreaAndPerim(MixLA0[0]);
AreaOTimeLP0, PeriOTimeLP0, TypeLP0 = analysis.AvAreaAndPerim(MixLP0[0]);

fig5, axs5 = plt.subplots(2,3, figsize=(7,4.5), dpi=900, sharex='row',sharey='row')
fig5.subplots_adjust(hspace=0.4)
#Call the plots on each frame
PlotFig5(AreaOTimeSA0,PeriOTimeSA0,TypeSA0,axs5[0,1])
PlotFig5(AreaOTimeSP0,PeriOTimeSP0,TypeSP0,axs5[0,2])
PlotFig5(AreaOTimeLA0,PeriOTimeLA0,TypeLA0,axs5[1,1])
PlotFig5(AreaOTimeLP0,PeriOTimeLP0,TypeLP0,axs5[1,2])

#Solid
#Upper left
#axs5[0,0].scatter(np.power(PeriOTimeS,2),AreaOTimeS,color=colorsolid)
axs5[0,0].set_xlabel('squared perimeter')
axs5[0,0].set_ylabel(r"$a$")
axs5[0,0].set_title('a)',fontfamily='serif',loc='left')
axs5[0,0].set_ylabel('area')

#Upper middle
axs5[0,1].axhline(y=1,color=colorresident,linestyle='dashed');
axs5[0,1].axhline(y=0.9,color=colormutant,linestyle='dashed');
axs5[0,1].axvline(x=3.3**2,color=_new_black,linestyle='dashed');
axs5[0,1].set_yticks([0.8, 0.9, 1, 1.2]);
axs5[0,1].set_yticklabels([0.8, r'$a_m^0$',r'$a_r^0$', 1.2]);
axs5[0,1].set_xlabel('squared perimeter')
#axs5[0,1].set_ylabel('area')
axs5[0,1].set_title('b)',fontfamily='serif',loc='left')

#Upper right
axs5[0,2].set_yticks([0.8, 0.9, 1, 1.2]);
axs5[0,2].set_yticklabels([0.8, r'$a_m^0$',r'$a^0$', 1.2]);
axs5[0,2].set_xticks([10, 3.3**2, 3.63**2, 15,17])
axs5[0,2].set_xticklabels([10, r'${p_r^0}^2$', r'${p_m^0}^2$', 15, 17])
#axs5[0,2].plot(PeriOTimeSP0,0.565*PeriOTimeSP0,color='black')
axs5[0,2].axvline(x=3.3**2,color=colorresident,linestyle='dashed');
axs5[0,2].axvline(x=3.63**2,color=colormutant,linestyle='dashed');
axs5[0,2].axhline(y=1,color=_new_black,linestyle='dashed');
axs5[0,2].set_xlabel('squared perimeter')
#axs5[0,2].set_ylabel('area')
axs5[0,2].set_title('c)',fontfamily='serif',loc='left')


#Liquid

#Lower left
#axs5[1,0].scatter(np.power(PeriOTimeL,2),AreaOTimeL,color=colorliq)
axs5[1,0].set_xlabel('squared perimeter')
axs5[1,0].set_ylabel(r"$a$")
axs5[1,0].set_title('d)',fontfamily='serif',loc='left')
axs5[1,0].set_ylabel('area')



#Lower Middle
axs5[1,1].axhline(y=1,color=colorresident,linestyle='dashed');
axs5[1,1].axhline(y=0.9,color=colormutant,linestyle='dashed');
axs5[1,1].axvline(x=4**2,color=_new_black,linestyle='dashed');
axs5[1,1].set_yticks([0.8, 0.9, 1, 1.2]);
axs5[1,1].set_yticklabels([0.8, r'$a_m^0$',r'$a_r^0$', 1.2]);
axs5[1,1].set_xticks([4**2, 14,18,19])
axs5[1,1].set_xticklabels([r'${p_r^0}^2$',14,18,19])
axs5[1,1].set_xlabel('squared perimeter')
axs5[1,1].set_title('e)',fontfamily='serif',loc='left')


#Lower right
axs5[1,2].axvline(x=4**2,color=colorresident,linestyle='dashed');
axs5[1,2].axvline(x=4.4**2,color=colormutant,linestyle='dashed');
axs5[1,2].axhline(y=1,color=_new_black,linestyle='dashed');
axs5[1,2].set_xticks([4**2, 4.4**2, 14,18])
axs5[1,2].set_xticklabels([r'${p_r^0}^2$', r'${p_m^0}^2$', 14,18])
axs5[1,2].set_xlabel('squared perimeter')
axs5[1,2].set_title('f)',fontfamily='serif',loc='left')





#fig5.savefig('/home/louisbrezin/Documents/PostDoc/VertexVisual/Fig5-2x3-shareX.pdf',format='pdf',bbox_inches='tight')


#%% Load data Fig 5.bis
fileFracLA0 = open('./DataFigures/VariedFraction-Liquid-A0Mixture','rb');
fileFracLP0 = open('./DataFigures/VariedFraction-Liquid-P0Mixture','rb');
fileFracSP0 = open('./DataFigures/VariedFraction-Solid-P0Mixture','rb');
fileFracSA0 = open('./DataFigures/VariedFraction-Solid-A0Mixture','rb');
FracSA0= pickle.load(fileFracSA0);
FracSP0= pickle.load(fileFracSP0);
FracLA0= pickle.load(fileFracLA0);
FracLP0= pickle.load(fileFracLP0);
fileFracLA0.close();
fileFracLP0.close();
fileFracSP0.close();
fileFracSA0.close();

fileFracLA0_bis = open('./DataFigures/VariedFraction-Liquid-A0Mixture-bis','rb');
fileFracLP0_bis = open('./DataFigures/VariedFraction-Liquid-P0Mixture-bis','rb');
fileFracSP0_bis = open('./DataFigures/VariedFraction-Solid-P0Mixture-bis','rb');
fileFracSA0_bis = open('./DataFigures/VariedFraction-Solid-A0Mixture-bis','rb');
FracSA0_bis= pickle.load(fileFracSA0_bis);
FracSP0_bis= pickle.load(fileFracSP0_bis);
FracLA0_bis= pickle.load(fileFracLA0_bis);
FracLP0_bis= pickle.load(fileFracLP0_bis);
fileFracLA0_bis.close();
fileFracLP0_bis.close();
fileFracSP0_bis.close();
fileFracSA0_bis.close();

#%% Fig 5 bis -- Plot area vs perimeter for different mixture fractions
fig5b = plt.figure(figsize=(7,5),dpi=600)
subfigs = fig5b.subfigures(2, 1, hspace=0.1,height_ratios=[2, 1])
ax = subfigs[0].subplots(2,3,sharex=True, sharey=True)
subfigs[0].subplots_adjust(hspace=0.4)
axL = subfigs[1].subplots(1,2)
#fig5b,ax = plt.subplots(3,3,figsize=(7,5),sharex='row', sharey='row')
colors=[colorsolid,colorliq]
colors2 = [['goldenrod','darkgoldenrod'],['cornflowerblue','blue']]
titles = ['solid','liquid']
for index,items in enumerate([FracSP0, FracSP0_bis,FracLP0, FracLP0_bis]):
    row = index//2;
    Nfile = len(items)
    Nt = len(items[0])
    AvArea = [];
    AvPerimeter = []; 
    AvAreaErr = [];
    AvPerimeterErr = []; 
    Nsample = int(0.9*Nt) #Sets the number of sample as a fraction of the second half of the simulation
    #Nsample = 1
    seed = None
    rng = np.random.default_rng(seed)
    samples = rng.integers(Nt//2, Nt, size=Nsample)
    AvArea0 = np.zeros(Nfile); AvArea1 = np.zeros(Nfile);
    AvAreaError0 = np.zeros(Nfile); AvAreaError1 = np.zeros(Nfile);
    AvPeri0 = np.zeros(Nfile); AvPeri1 = np.zeros(Nfile);
    AvPeriError0 = np.zeros(Nfile); AvPeriError1 = np.zeros(Nfile);
    #DAreaDPeri0 = np.zeros(Nfile); DArea0 = np.zeros(Nfile); DPeri0 = np.zeros(Nfile);
    #DAreaDPeri1 = np.zeros(Nfile); DArea1 = np.zeros(Nfile); DPeri1 = np.zeros(Nfile);
    if index%2==0:
        fraction = np.linspace(0.1, 0.9,Nfile)
    elif index%2==1:
        fraction = np.linspace(0.05, 0.95,Nfile)
            
    for i in range(Nfile):
        AvAreaTemp = [];
        AvPeriTemp = [];
        AvAreaErrTemp = [];
        AvPeriErrTemp = [];
        for j in range(Nsample):
            cellarea = reader.CellArea(samples[j], items[i]);
            cellperimeter = reader.CellPeri(samples[j], items[i]);
            celltype = reader.CellType(samples[j], items[i]);
            AvAreaTemp.append(analysis.AverageValuePerType(cellarea,celltype)[:,0]);
            AvPeriTemp.append(analysis.AverageValuePerType(cellperimeter,celltype)[:,0]);
            AvAreaErrTemp.append(analysis.AverageValuePerType(cellarea,celltype)[:,1]);
            AvPeriErrTemp.append(analysis.AverageValuePerType(cellperimeter,celltype)[:,1]);
            
        AvAreaErr.append(np.sqrt(np.mean(np.power(AvAreaErrTemp,2),axis=0) / Nsample))
        AvPerimeterErr.append(np.sqrt(np.mean(np.power(AvPeriErrTemp,2),axis=0) / Nsample))
        AvArea.append(np.mean(AvAreaTemp,axis=0));
        AvPerimeter.append(np.mean(AvPeriTemp,axis=0));
        
        AvArea0[i] = AvArea[i][0];
        AvAreaError0[i] = AvAreaErr[i][0];
        AvArea1[i] = AvArea[i][1];
        AvAreaError1[i] = AvAreaErr[i][1];
        AvPeri0[i] = AvPerimeter[i][0];
        AvPeriError0[i] = AvPerimeterErr[i][0];
        AvPeri1[i] = AvPerimeter[i][1];
        AvPeriError1[i] = AvPerimeterErr[i][1];
        
            
    
    DArea0 = AvArea0 - (np.multiply(AvArea1,fraction) + np.multiply(AvArea0,1-fraction));
    DArea1 = AvArea1 - (np.multiply(AvArea1,fraction) + np.multiply(AvArea0,1-fraction));
    DPeri0 = AvPeri0 - (np.multiply(AvPeri1,fraction) + np.multiply(AvPeri0,1-fraction));
    DPeri1 = AvPeri1 - (np.multiply(AvPeri1,fraction) + np.multiply(AvPeri0,1-fraction));
    DAreaDPeri1 = np.divide(DArea1, DPeri1);
    DAreaDPeri0 = np.divide(DArea0, DPeri0);
    
    
    ax[row,0].scatter(fraction,DArea0,color=colors2[row][0],marker='o',label='resident')
    ax[row,0].scatter(fraction,DArea1,color=colors2[row][1],marker='^',label='mutant')
    ax[row,0].set_ylabel(r'$\Delta a $')
    
    
    
    ax[row,1].scatter(fraction,DPeri0,color=colors2[row][0],marker='o',label='resident')
    ax[row,1].scatter(fraction,DPeri1,color=colors2[row][1],marker='^',label='mutant')
    ax[row,1].set_ylabel(r'$\Delta p $')
    
    ax[row,2].scatter(fraction,DAreaDPeri0,color=colors2[row][0],marker='o',label='resident')
    ax[row,2].scatter(fraction,DAreaDPeri1,color=colors2[row][1],marker='^',label='mutant')
    if row ==1:
        ax[row,2].set_xlabel(r'fraction of mutant, $f$')
        ax[row,1].set_xlabel(r'fraction of mutant, $f$')
        ax[row,0].set_xlabel(r'fraction of mutant, $f$')
    ax[row,2].set_ylabel(r'$\Delta a / \Delta p $')
    if index%2 == 0:
        ax[row,2].legend()
    
#Plot the 50-50 mixture with the da/dp slope that was found.
DaDps = [0.565,0.037]
for index, items in enumerate([FracSP0[4],FracLP0[4]]):
    AreasFull, PerimetersFull, Celltypes = analysis.AvAreaAndPerim(items)
    Perimeters = analysis.ValuePerType(PerimetersFull, Celltypes);
    Areas = analysis.ValuePerType(AreasFull, Celltypes);
    p = np.linspace(min(PerimetersFull), max(PerimetersFull),50)
    axL[index].scatter(Perimeters[0],Areas[0],marker='o',color = colors2[index][0], label='resident');
    axL[index].scatter(Perimeters[1],Areas[1],marker='^', color = colors2[index][1], label='mutant');
    axL[index].plot(p,DaDps[index]*(p-np.mean(PerimetersFull))+1,color='black')
    axL[index].set_xlabel(r'perimeter, $p$')
    axL[index].legend()
axL[0].set_ylabel(r'area, $a$')

ax[0,0].set_title('a)',fontfamily='serif',loc='left')
ax[1,0].set_title('b)',fontfamily='serif',loc='left')
axL[0].set_title('c)',fontfamily='serif',loc='left')


fig5b.savefig('/home/louisbrezin/Documents/PostDoc/VertexVisual/Fig5-bis.pdf',format='pdf',bbox_inches='tight',pad_inches=0)


#%% Load data Fig 6
fileFractionr0 = open('./DataFigures/Fraction-homeostatisr0-fig6-20241008','rb');
fileFractionK = open('./DataFigures/Fraction-homeostatisK-fig6-20241008','rb');
fileFractionLambda = open('./DataFigures/Fraction-homeostatisLambda-fig6-20241008','rb');
fileFractionv0 = open('./DataFigures/Fraction-homeostatisv0-fig6','rb');
fileFractionDr = open('./DataFigures/Fraction-homeostatisDr-fig6','rb');
fileTime6 = open('./DataFigures/Time-fig6-20241008','rb');
fileTime6_v = open('./DataFigures/Time-fig6-v0-20241107','rb');
frac2Fullr0 = pickle.load(fileFractionr0);
frac2FullK = pickle.load(fileFractionK);
frac2FullLambda = pickle.load(fileFractionLambda);
frac2Fullv0 = pickle.load(fileFractionv0);
frac2FullDr = pickle.load(fileFractionDr);
T6 = pickle.load(fileTime6);
T6_v = pickle.load(fileTime6_v);

#%% Plot Fig 6
coloring = ['red','c','g','m','grey','k','c'];
labeling = ['0.1','1','10','100','4','Perim'];
frac2FullAverageVec = [];
fig6 = plt.figure(figsize=(7,4.5), dpi=600)
subfigs6 = fig6.subfigures(2, 1, hspace=0.1)
alphaplot = 0.3
#gs = GridSpec(1, 6, figure=subfigs6[0]);
#axs0 = subfigs6.add_subplot(gs[1:2])
#axs1 = subfigs6.add_subplot(gs[3:4])
axs = subfigs6[0].subplots(1,2,sharey=True,sharex=False)
axsL = subfigs6[1].subplots(1,3,sharey=True,sharex=True)
#fig6,axs = plt.subplots(1,3,figsize=(7,2), dpi=900,sharey=True,sharex=False);
filecycle=10;
linewidth6=1;
marker6='-';
#ax2 = ax.twinx();
#ax2 = ax.twiny();
for i in range(frac2Fullr0.shape[1]):
    exp = i//filecycle;
    #axs[0].plot(100*T6,np.divide(frac2Fullr0[:,i],1-frac2Fullr0[:,i]),color=coloring[int(i/filecycle)],alpha=0.5,linewidth=0.5)
    if i%filecycle==0:
        frac2FullAverage = np.mean(frac2Fullr0[:,filecycle*exp:filecycle*exp+filecycle-1],axis=1)
        std = np.std(frac2Fullr0[:,filecycle*exp:filecycle*exp+filecycle-1],axis=1)
        axs[0].plot(100*T6,np.divide(frac2FullAverage,1-frac2FullAverage),label=labeling[int(i/filecycle)],color=coloring[int(i/filecycle)],alpha=1,linewidth=linewidth6)
        #axs[0].fill_between(100*T6,np.divide(frac2FullAverage,1-frac2FullAverage)-0.5*std,np.divide(frac2FullAverage,1-frac2FullAverage)+0.5*std,color=coloring[int(i/filecycle)],alpha=alphaplot)
legend0=axs[0].legend(loc='upper right',prop={'size': 6},labelspacing=0.1);
legend0.set_title(r'$r_0$')
#axs[0].legend(prop={'size': 6})
axs[0].set_yscale('log');
axs[0].minorticks_off();
axs[0].set_yticks([0.1, 10,1000])
#axs[0].set_yticklabels([0.01, 1,100]);
axs[0].set_xlabel(r'timesteps')
axs[0].set_ylabel(r'$\frac{f}{1-f}$')
axs[0].set_ylim(0.1,1000)
axs[0].set_title('a)',fontfamily='serif',loc='left')


for i in range(frac2FullLambda.shape[1]):
    exp = i//filecycle;
    #axs[2].plot(0.01*T6,np.divide(frac2FullLambda[:,i],1-frac2FullLambda[:,i]),color=coloring[int(i/filecycle)],alpha=0.5,linewidth=0.5)
    if i%filecycle==0:
        frac2FullAverage = np.mean(frac2FullLambda[:,filecycle*exp:filecycle*exp+filecycle-1],axis=1)
        std = np.std(frac2FullLambda[:,filecycle*exp:filecycle*exp+filecycle-1],axis=1)
        axs[1].plot(0.01*T6,np.divide(frac2FullAverage,1-frac2FullAverage),label=labeling[int(i/filecycle)],color=coloring[int(i/filecycle)],alpha=1,linewidth=linewidth6)
        #axs[1].fill_between(0.01*T6,np.divide(frac2FullAverage,1-frac2FullAverage)-0.5*std,np.divide(frac2FullAverage,1-frac2FullAverage)+0.5*std,color=coloring[int(i/filecycle)],alpha=alphaplot)
legend2=axs[1].legend(loc='upper right',prop={'size': 6},labelspacing=0.1);
legend2.set_title(r'$\lambda$')
axs[1].set_yscale('log');
axs[1].minorticks_off();
axs[1].set_yticks([0.1, 10,1000])
axs[1].set_xlabel(r'$\tau_g$')
axs[1].set_title('b)',fontfamily='serif',loc='left')


for i in range(frac2FullK.shape[1]):
    exp = i//filecycle;
    #axs[1].plot(0.01*T6,np.divide(frac2FullK[:,i],1-frac2FullK[:,i]),color=coloring[int(i/filecycle)],alpha=0.5,linewidth=0.5)
    if i%filecycle==0:
        frac2FullAverage = np.mean(frac2FullK[:,filecycle*exp:filecycle*exp+filecycle-1],axis=1)
        axsL[0].plot(0.01*T6,np.divide(frac2FullAverage,1-frac2FullAverage),label=labeling[int(i/filecycle)],color=coloring[int(i/filecycle)],alpha=1,linewidth=linewidth6)
        std=np.std(frac2FullK[:,filecycle*exp:filecycle*exp+filecycle-1],axis=1)
        #axsL[0].fill_between(0.01*T6,np.divide(frac2FullAverage,1-frac2FullAverage)-0.5*std,np.divide(frac2FullAverage,1-frac2FullAverage)+0.5*std,color=coloring[int(i/filecycle)],alpha=alphaplot)
legend1=axsL[0].legend(loc='lower right',prop={'size': 6},labelspacing=0.1);
legend1.set_title(r'$K$')
axsL[0].set_yscale('log');
axsL[0].minorticks_off();
axsL[0].set_yticks([0.1, 10,1000])
axsL[0].set_ylabel(r'$\frac{f}{1-f}$')
axsL[0].set_xlabel(r'$\tau_g$')
axsL[0].set_title('c)',fontfamily='serif',loc='left')



labelv0 = [0,0.55,1.55]
for i in range(frac2Fullv0.shape[1]):
    exp = i//filecycle;
    #axs[0].plot(100*T6,np.divide(frac2Fullr0[:,i],1-frac2Fullr0[:,i]),color=coloring[int(i/filecycle)],alpha=0.5,linewidth=0.5)
    if i%filecycle==0:
        frac2FullAverage = np.mean(frac2Fullv0[:,filecycle*exp:filecycle*exp+filecycle-1],axis=1)
        axsL[1].plot(0.01*T6_v,np.divide(frac2FullAverage,1-frac2FullAverage),label=labelv0[int(i/filecycle)],color=coloring[int(i/filecycle)],alpha=1,linewidth=linewidth6)
        std=np.std(frac2Fullv0[:,filecycle*exp:filecycle*exp+filecycle-1],axis=1)
        #axsL[1].fill_between(0.01*T6_v,np.divide(frac2FullAverage,1-frac2FullAverage)-0.5*std,np.divide(frac2FullAverage,1-frac2FullAverage)+0.5*std,color=coloring[int(i/filecycle)],alpha=alphaplot)
legend3=axsL[1].legend(loc='lower right',prop={'size': 6},labelspacing=0.1);
legend3.set_title(r'$v_0$')
#axs[0].legend(prop={'size': 6})
axsL[1].set_yscale('log');
axsL[1].minorticks_off();
axsL[1].set_yticks([0.1, 10,1000])
#axs[0].set_yticklabels([0.01, 1,100]);
axsL[1].set_xlabel(r'$\tau_g$')

axsL[1].set_ylim(0.1,1000)
axsL[1].set_title('d)',fontfamily='serif',loc='left')



for i in range(frac2FullDr.shape[1]):
    exp = i//filecycle;
    #axs[1].plot(0.01*T6,np.divide(frac2FullK[:,i],1-frac2FullK[:,i]),color=coloring[int(i/filecycle)],alpha=0.5,linewidth=0.5)
    if i%filecycle==0:
        frac2FullAverage = np.mean(frac2FullDr[:,filecycle*exp:filecycle*exp+filecycle-1],axis=1)
        axsL[2].plot(0.01*T6_v,np.divide(frac2FullAverage,1-frac2FullAverage),label=labeling[int(i/filecycle)],color=coloring[int(i/filecycle)],alpha=1,linewidth=linewidth6)
        std=np.std(frac2FullDr[:,filecycle*exp:filecycle*exp+filecycle-1],axis=1)
        #axsL[2].fill_between(0.01*T6_v,np.divide(frac2FullAverage,1-frac2FullAverage)-0.5*std,np.divide(frac2FullAverage,1-frac2FullAverage)+0.5*std,color=coloring[int(i/filecycle)],alpha=alphaplot)
legend4=axsL[2].legend(loc='lower right',prop={'size': 6},labelspacing=0.1);
legend4.set_title(r'$D_r$')
axsL[2].set_yscale('log');
axsL[2].minorticks_off();
axsL[2].set_yticks([0.1, 10,1000])
axsL[2].set_xlabel(r'$\tau_g$')
axsL[2].set_title('e)',fontfamily='serif',loc='left')


#fig6.savefig('Fig6-AllTrajectories.pdf',format='pdf',bbox_inches='tight')
fig6.savefig('Fig6-5panels.pdf',format='pdf',bbox_inches='tight')

#%% Supplementary figure 4
# Load data SI Fig 4 
#The data for the solid case is in the last of those simulation
fileFractionSI = open('./DataFigures/Fraction-InvasionSolidLiquid-DRM','rb');
fileFractionL_opposite = open('./DataFigures/Fraction-InvasionLiquid-DRM-opposite','rb');
fileFractionS_opposite = open('./DataFigures/Fraction-InvasionSolid-DRM-opposite','rb');


fileTimeSI = open('./DataFigures/Time-InvasionSolidLiquid-DRM','rb');
fileTimeL_opposite = open('./DataFigures/Time-InvasionLiquid-DRM-opposite','rb');

frac2FullSI = pickle.load(fileFractionSI);
frac2FullL_opposite = pickle.load(fileFractionL_opposite);
frac2FullS_opposite = pickle.load(fileFractionS_opposite);


T_SI = pickle.load(fileTimeSI);
T_opposite = pickle.load(fileTimeL_opposite);

fileFractionSI.close();
fileTimeSI.close();
fileFractionL_opposite.close();
fileFractionS_opposite.close();
fileTimeL_opposite.close();



#%% Removes the zeros and 1 from fraction to not trigger bugs. Puts them at unpysically small values
frac2FullSI=np.where(frac2FullSI==0,0.00001,frac2FullSI)
frac2FullSI=np.where(frac2FullSI==1,0.99999,frac2FullSI)
#%% Plot Fig S4 - DRM invasion
coloring = ['grey','red','blue','g','c','m','k','brown'];
coloringbis = [colorliq,colorsolid];
labeling = [r'$b^0$',r'$a^0$',r'$p^0$',r'$\Gamma$',r'$K$',r'$\lambda$',r'$v^0$',r'$D_r$'];
state = ['Liquid', 'Solid'];
listyle = ['dotted','solid'];
label = ['a)','b)','c)','d)','e)','f)','g)','h)']
#labeling = ['Perimeter elasticity','Area elasticity','Pressure strength','Self-advection','Rotational diffusion'];
filecycle=10;
macrocycle=80;
nline = 2;
ncol = 4;
figS4,axsS4=plt.subplots(nline,ncol, figsize=(7,3.5), dpi=900, sharex=True,sharey=True)
figS4.subplots_adjust(hspace=0.4)
tt=0;

state_legend = [
    plt.Line2D([0], [0], color=colorsolid, label='Solid state'),
    plt.Line2D([0], [0], color=colorliq, label='Liquid state')
]

# Legend for parameter change (linestyle)
change_legend = [
    plt.Line2D([0], [0], color= _new_black, linestyle='-', label='↑ param'),
    plt.Line2D([0], [0], color= _new_black, linestyle='--', label='↓ param')
]

for i in range(frac2FullSI.shape[1]):
    j = i%macrocycle;
    st = i//macrocycle;
    exp = j//filecycle;
    if j%filecycle==0:
        if exp==1:
            frac2FullAverageSI = np.mean(frac2FullSI[:,st*macrocycle+filecycle*exp:st*macrocycle+filecycle*exp+filecycle-1],axis=1)
            axsS4[exp//ncol,exp%ncol].plot(0.01*T_SI,np.divide(frac2FullAverageSI,1-frac2FullAverageSI),color=coloringbis[st],linestyle=(0, (5, 5)));
            if st == 0:
                frac2FullAverageS_opposite = np.mean(frac2FullS_opposite[:,filecycle*exp:filecycle*exp+filecycle-1],axis=1)
                frac2FullAverageL_opposite = np.mean(frac2FullL_opposite[:,filecycle*exp:filecycle*exp+filecycle-1],axis=1)
                axsS4[exp//ncol,exp%ncol].plot(0.01*T_opposite,np.divide(frac2FullAverageS_opposite,1-frac2FullAverageS_opposite),color=colorsolid);
                axsS4[exp//ncol,exp%ncol].plot(0.01*T_opposite,np.divide(frac2FullAverageL_opposite,1-frac2FullAverageL_opposite),color=colorliq);
        else:
            frac2FullAverageSI = np.mean(frac2FullSI[:,st*macrocycle+filecycle*exp:st*macrocycle+filecycle*exp+filecycle-1],axis=1)
            axsS4[exp//ncol,exp%ncol].plot(0.01*T_SI,np.divide(frac2FullAverageSI,1-frac2FullAverageSI),color=coloringbis[st]);
            if st == 0:
                frac2FullAverageS_opposite = np.mean(frac2FullS_opposite[:,filecycle*exp:filecycle*exp+filecycle-1],axis=1)
                frac2FullAverageL_opposite = np.mean(frac2FullL_opposite[:,filecycle*exp:filecycle*exp+filecycle-1],axis=1)
                axsS4[exp//ncol,exp%ncol].plot(0.01*T_opposite,np.divide(frac2FullAverageS_opposite,1-frac2FullAverageS_opposite),color=colorsolid,linestyle=(0, (5, 5)));
                axsS4[exp//ncol,exp%ncol].plot(0.01*T_opposite,np.divide(frac2FullAverageL_opposite,1-frac2FullAverageL_opposite),color=colorliq,linestyle=(0, (5, 5)));
        #std = np.std(frac2Full[:,st*macrocycle+filecycle*exp:st*macrocycle+filecycle*exp+filecycle-1],axis=1)
        #y = np.divide(frac2Full[:,i],1-frac2Full[:,i]);
        #cond = y>10^(-3);
        #Logitfrac2FullAverageL = np.mean(np.log(np.divide(frac2FullLiquid[:,st*macrocycle+filecycle*exp:st*macrocycle+filecycle*exp+filecycle-1],1-frac2FullLiquid[:,st*macrocycle+filecycle*exp:st*macrocycle+filecycle*exp+filecycle-1])),axis=1)
        #axs4[exp//ncol,exp%ncol].plot(0.01*TLiquid,Logitfrac2FullAverageL,color=colorliq);
        axsS4[exp//ncol,exp%ncol].set_yscale('log');
        axsS4[exp//ncol,exp%ncol].set_ylim(0.01,1000);
        axsS4[exp//ncol,exp%ncol].set_yticks([0.1,10,1000]);
        #axs4[exp//ncol,exp%ncol].set_ylim(-10,10);
        #axs4[exp//ncol,exp%ncol].set_yticks([-10,-5,0,5,10]);
        axsS4[exp//ncol,exp%ncol].minorticks_off();
        axsS4[exp//ncol,exp%ncol].set_title(labeling[exp]);
        axsS4[exp//ncol,exp%ncol].set_title(label[exp],fontfamily='serif',loc='left')
        #axs[exp//ncol,exp%ncol].legend();
        axsS4[nline-1,exp%ncol].set_xlabel(r'$\tau_g$')
        axsS4[exp//ncol,0].set_ylabel(r'$\frac{f}{1-f}$')
        axsS4[1,1].legend(handles=state_legend, fontsize='small', frameon=False)
        axsS4[1,2].legend(handles=change_legend, fontsize='small', frameon=False)

#figS4.savefig('Fig-S4-opposite.pdf',format='pdf',bbox_inches='tight')
#fig4.savefig('Fig4-AverageLog.pdf',format='pdf',bbox_inches='tight')
#fig4.savefig('Fig4-Floor-newP.pdf',format='pdf',bbox_inches='tight')

#%% Supplementary fig 2
fileTime = open('./DataFigures/Time-HomeoPA0-figS2','rb');
T = pickle.load(fileTime);
fileTime.close();
fileData = open('./DataFigures/HomeoPA0-figS2','rb');
Data = pickle.load(fileData);
fileData.close();
fileParams = open('./DataFigures/HomeoPA0-figS2-Params','rb');
Params = pickle.load(fileParams);
fileParams.close();
Nfile = len(Data);


#%% Number of cells as function of homeostatic pressure
figS2, axesS2 = plt.subplots(nrows=1,ncols=1,figsize=(3.42,2),dpi=600)

AvOTimeNc = np.zeros(Nfile);
NIni = np.zeros(Nfile);
NRatio = np.zeros(Nfile);
Ph = np.zeros(Nfile);
RunAverage = np.zeros((Nfile,Nt));
for j in range(Nfile):
    NumCells = reader.CellNumOTime(Data[j])[1];
    RunAverage[j,:]=reader.RunningCellAverage(Data[j])[1];
    AvOTimeNc[j] = np.mean(NumCells[3:]);
    NRatio[j] = RunAverage[j,-1]/NumCells[0];
    Ph[j] = Params[j][1][7];
    
#%Fit and do a plot

coef = np.polyfit(Ph, NRatio, 1);
poly1d_fn = np.poly1d(coef);
Ph_homeo = (1/Params[0][1][4]-coef[1])/coef[0];
axesS2.scatter(Ph,NRatio);
axesS2.plot(Ph,poly1d_fn(Ph),'-k');
axesS2.set_xlabel(r'$P_h$')
axesS2.set_ylabel(r'$N/N_0$')
figS2.savefig('Fig-S2.pdf',format='pdf',bbox_inches='tight')


#%% Load data Supplementary Fig 6
fileFractionr0 = open('./DataFigures/Fraction-homeostatisr0-fig6-area-20241021','rb');
fileFractionK = open('./DataFigures/Fraction-Epistasis_DRM-K','rb');
fileFractionLambda = open('./DataFigures/Fraction-homeostatisLambda-fig6-area-20241021','rb');
fileFractionv0 = open('./DataFigures/Fraction-Epistasis_DRM-v0','rb');
fileFractionDr = open('./DataFigures/Fraction-Epistasis_DRM-Dr','rb');
fileFractionGamma = open('./DataFigures/Fraction-Epistasis_DRM-Gamma-a0-08','rb');
fileTime6 = open('./DataFigures/Time-Epistasis_DRM-Dr','rb');
fileTime6_v = open('./DataFigures/Time-fig6-area-20241021','rb');
frac2Fullr0 = pickle.load(fileFractionr0);
frac2FullK = pickle.load(fileFractionK);
frac2FullLambda = pickle.load(fileFractionLambda);
frac2Fullv0 = pickle.load(fileFractionv0);
frac2FullDr = pickle.load(fileFractionDr);
frac2FullGamma = pickle.load(fileFractionGamma);
T6 = pickle.load(fileTime6);
T6_area = pickle.load(fileTime6_v);

#%% Plot Supplementary Fig 6
coloring = ['red','c','g','m','grey','k','c'];
labeling = ['0.1','1','10','100','4','Perim'];
frac2FullAverageVec = [];
figS6 = plt.figure(figsize=(7,4.5), dpi=600)
subfigsS6 = figS6.subfigures(2, 1, hspace=0.1)
alphaplot = 0.3
#gs = GridSpec(1, 6, figure=subfigs6[0]);
#axs0 = subfigs6.add_subplot(gs[1:2])
#axs1 = subfigs6.add_subplot(gs[3:4])
axs = subfigsS6[0].subplots(1,3,sharey=True,sharex=False)
axsL = subfigsS6[1].subplots(1,3,sharey=True,sharex=True)
#fig6,axs = plt.subplots(1,3,figsize=(7,2), dpi=900,sharey=True,sharex=False);
filecycle=10;
linewidth6=1;
marker6='-';
#ax2 = ax.twinx();
#ax2 = ax.twiny();
for i in range(frac2Fullr0.shape[1]):
    exp = i//filecycle;
    #axs[0].plot(100*T6,np.divide(frac2Fullr0[:,i],1-frac2Fullr0[:,i]),color=coloring[int(i/filecycle)],alpha=0.5,linewidth=0.5)
    if i%filecycle==0:
        frac2FullAverage = np.mean(frac2Fullr0[:,filecycle*exp:filecycle*exp+filecycle-1],axis=1)
        std = np.std(frac2Fullr0[:,filecycle*exp:filecycle*exp+filecycle-1],axis=1)
        axs[0].plot(100*T6_area,np.divide(frac2FullAverage,1-frac2FullAverage),label=labeling[int(i/filecycle)],color=coloring[int(i/filecycle)],alpha=1,linewidth=linewidth6)
        #axs[0].fill_between(100*T6,np.divide(frac2FullAverage,1-frac2FullAverage)-0.5*std,np.divide(frac2FullAverage,1-frac2FullAverage)+0.5*std,color=coloring[int(i/filecycle)],alpha=alphaplot)
legend0=axs[0].legend(loc='upper right',prop={'size': 6},labelspacing=0.1);
legend0.set_title(r'$r_0$')
#axs[0].legend(prop={'size': 6})
axs[0].set_yscale('log');
axs[0].minorticks_off();
axs[0].set_yticks([0.1, 10,1000])
#axs[0].set_yticklabels([0.01, 1,100]);
axs[0].set_xlabel(r'timesteps')
axs[0].set_ylabel(r'$\frac{f}{1-f}$')
axs[0].set_ylim(0.02,1000)
axs[0].set_title('a)',fontfamily='serif',loc='left')


for i in range(frac2FullLambda.shape[1]):
    exp = i//filecycle;
    #axs[2].plot(0.01*T6,np.divide(frac2FullLambda[:,i],1-frac2FullLambda[:,i]),color=coloring[int(i/filecycle)],alpha=0.5,linewidth=0.5)
    if i%filecycle==0:
        frac2FullAverage = np.mean(frac2FullLambda[:,filecycle*exp:filecycle*exp+filecycle-1],axis=1)
        std = np.std(frac2FullLambda[:,filecycle*exp:filecycle*exp+filecycle-1],axis=1)
        axs[1].plot(0.01*T6_area,np.divide(frac2FullAverage,1-frac2FullAverage),label=labeling[int(i/filecycle)],color=coloring[int(i/filecycle)],alpha=1,linewidth=linewidth6)
        #axs[1].fill_between(0.01*T6,np.divide(frac2FullAverage,1-frac2FullAverage)-0.5*std,np.divide(frac2FullAverage,1-frac2FullAverage)+0.5*std,color=coloring[int(i/filecycle)],alpha=alphaplot)
legend2=axs[1].legend(loc='upper right',prop={'size': 6},labelspacing=0.1);
legend2.set_title(r'$\lambda$')
axs[1].set_yscale('log');
axs[1].minorticks_off();
axs[1].set_yticks([0.1, 10,1000])
axs[1].set_xlabel(r'$\tau_g$')
axs[1].set_title('b)',fontfamily='serif',loc='left')

for i in range(frac2FullGamma.shape[1]):
    exp = i//filecycle;
    #axs[2].plot(0.01*T6,np.divide(frac2FullLambda[:,i],1-frac2FullLambda[:,i]),color=coloring[int(i/filecycle)],alpha=0.5,linewidth=0.5)
    if i%filecycle==0:
        frac2FullAverage = np.mean(frac2FullGamma[:,filecycle*exp:filecycle*exp+filecycle-1],axis=1)
        std = np.std(frac2FullGamma[:,filecycle*exp:filecycle*exp+filecycle-1],axis=1)
        axs[2].plot(0.01*T6,np.divide(frac2FullAverage,1-frac2FullAverage),label=labeling[int(i/filecycle)],color=coloring[int(i/filecycle)],alpha=1,linewidth=linewidth6)
        #axs[1].fill_between(0.01*T6,np.divide(frac2FullAverage,1-frac2FullAverage)-0.5*std,np.divide(frac2FullAverage,1-frac2FullAverage)+0.5*std,color=coloring[int(i/filecycle)],alpha=alphaplot)
legend2=axs[2].legend(loc='upper right',prop={'size': 6},labelspacing=0.1);
legend2.set_title(r'$\Gamma$')
axs[2].set_yscale('log');
axs[2].minorticks_off();
axs[2].set_yticks([0.1, 10,1000])
axs[2].set_xlabel(r'$\tau_g$')
axs[2].set_title('c)',fontfamily='serif',loc='left')


for i in range(frac2FullK.shape[1]):
    exp = i//filecycle;
    #axs[1].plot(0.01*T6,np.divide(frac2FullK[:,i],1-frac2FullK[:,i]),color=coloring[int(i/filecycle)],alpha=0.5,linewidth=0.5)
    if i%filecycle==0:
        frac2FullAverage = np.mean(frac2FullK[:,filecycle*exp:filecycle*exp+filecycle-1],axis=1)
        axsL[0].plot(0.01*T6,np.divide(frac2FullAverage,1-frac2FullAverage),label=labeling[int(i/filecycle)],color=coloring[int(i/filecycle)],alpha=1,linewidth=linewidth6)
        std=np.std(frac2FullK[:,filecycle*exp:filecycle*exp+filecycle-1],axis=1)
        #axsL[0].fill_between(0.01*T6,np.divide(frac2FullAverage,1-frac2FullAverage)-0.5*std,np.divide(frac2FullAverage,1-frac2FullAverage)+0.5*std,color=coloring[int(i/filecycle)],alpha=alphaplot)
legend1=axsL[0].legend(loc='lower right',prop={'size': 6},labelspacing=0.1);
legend1.set_title(r'$K$')
axsL[0].set_yscale('log');
axsL[0].minorticks_off();
axsL[0].set_yticks([0.1, 10,1000])
axsL[0].set_ylabel(r'$\frac{f}{1-f}$')
axsL[0].set_xlabel(r'$\tau_g$')
axsL[0].set_title('d)',fontfamily='serif',loc='left')



labelv0 = [0,0.55,1.55]
for i in range(frac2Fullv0.shape[1]):
    exp = i//filecycle;
    #axs[0].plot(100*T6,np.divide(frac2Fullr0[:,i],1-frac2Fullr0[:,i]),color=coloring[int(i/filecycle)],alpha=0.5,linewidth=0.5)
    if i%filecycle==0:
        frac2FullAverage = np.mean(frac2Fullv0[:,filecycle*exp:filecycle*exp+filecycle-1],axis=1)
        axsL[1].plot(0.01*T6,np.divide(frac2FullAverage,1-frac2FullAverage),label=labelv0[int(i/filecycle)],color=coloring[int(i/filecycle)],alpha=1,linewidth=linewidth6)
        std=np.std(frac2Fullv0[:,filecycle*exp:filecycle*exp+filecycle-1],axis=1)
        #axsL[1].fill_between(0.01*T6_v,np.divide(frac2FullAverage,1-frac2FullAverage)-0.5*std,np.divide(frac2FullAverage,1-frac2FullAverage)+0.5*std,color=coloring[int(i/filecycle)],alpha=alphaplot)
legend3=axsL[1].legend(loc='lower right',prop={'size': 6},labelspacing=0.1);
legend3.set_title(r'$v_0$')
#axs[0].legend(prop={'size': 6})
axsL[1].set_yscale('log');
axsL[1].minorticks_off();
axsL[1].set_yticks([0.1, 10,1000])
#axs[0].set_yticklabels([0.01, 1,100]);
axsL[1].set_xlabel(r'$\tau_g$')

axsL[1].set_ylim(0.1,1000)
axsL[1].set_title('e)',fontfamily='serif',loc='left')



for i in range(frac2FullDr.shape[1]):
    exp = i//filecycle;
    #axs[1].plot(0.01*T6,np.divide(frac2FullK[:,i],1-frac2FullK[:,i]),color=coloring[int(i/filecycle)],alpha=0.5,linewidth=0.5)
    if i%filecycle==0:
        frac2FullAverage = np.mean(frac2FullDr[:,filecycle*exp:filecycle*exp+filecycle-1],axis=1)
        axsL[2].plot(0.01*T6,np.divide(frac2FullAverage,1-frac2FullAverage),label=labeling[int(i/filecycle)],color=coloring[int(i/filecycle)],alpha=1,linewidth=linewidth6)
        std=np.std(frac2FullDr[:,filecycle*exp:filecycle*exp+filecycle-1],axis=1)
        #axsL[2].fill_between(0.01*T6_v,np.divide(frac2FullAverage,1-frac2FullAverage)-0.5*std,np.divide(frac2FullAverage,1-frac2FullAverage)+0.5*std,color=coloring[int(i/filecycle)],alpha=alphaplot)
legend4=axsL[2].legend(loc='lower right',prop={'size': 6},labelspacing=0.1);
legend4.set_title(r'$D_r$')
axsL[2].set_yscale('log');
axsL[2].minorticks_off();
axsL[2].set_yticks([0.1, 10,1000])
axsL[2].set_xlabel(r'$\tau_g$')
axsL[2].set_title('f)',fontfamily='serif',loc='left')


#fig6.savefig('Fig6-AllTrajectories.pdf',format='pdf',bbox_inches='tight')
figS6.savefig('FigS6-5panels.pdf',format='pdf',bbox_inches='tight',pad_inches=0)

#%%Load data to plot 
fileDiffA = open('./DataFigures/AvDiff-MixLiq-DeltaA','rb');
DiffA = pickle.load(fileDiffA);
fileDiffA.close();
fileDiffP = open('./DataFigures/AvDiff-MixLiq-DeltaP','rb');
DiffP = pickle.load(fileDiffP);
fileDiffP.close();

#%%Plot Delta p as function of Delta a
#plt.scatter(AvArea0,AvPerimeter0)
fig, axes = plt.subplots(nrows=1,ncols=2,figsize=(6.5,2))
fig.subplots_adjust(wspace=0.4)
#axes[0,0].scatter(np.power(PeriOTimeS,2),AreaOTimeS,color=colorsolid)
#axes[0,0].set_xlabel(r"$p^2$")
#axes[0,0].set_ylabel(r"$a$")
#axes[0,0].set_title('a)',fontfamily='serif',loc='left')

#axes[0,1].scatter(np.power(PeriOTimeL,2),AreaOTimeL,color=colorliq)
#axes[0,1].set_xlabel(r"$p^2$")
#axes[0,1].set_ylabel(r"$a$")
#axes[0,1].set_title('b)',fontfamily='serif',loc='left')

axes[0].scatter(DiffA[0],DiffA[1],color=colorliq)
axes[0].set_xlabel(r"$\Delta a/a^0$")
axes[0].set_ylabel(r"$\Delta p/p^0$")
axes[0].ticklabel_format(style='sci',scilimits=(-2,2))
axes[0].set_title('a)',fontfamily='serif',loc='left')

axes[1].scatter(DiffP[1],DiffP[0],color=colorliq)
axes[1].set_ylabel(r"$\Delta a/a^0$")
axes[1].set_xlabel(r"$\Delta p/p^0$")
axes[1].ticklabel_format(style='sci',scilimits=(-2,2))
axes[1].set_title('b)',fontfamily='serif',loc='left')
fig.savefig("Homo-DeltaAP-mix-1x2.pdf",format='pdf',bbox_inches='tight')


#%% Supplementary figure that shows the realized perimeter as a function of preferred perimeter
fileRealizedPeri = open('./DataFigures/RealizedPerim-P0','rb');
fileVariedParam = open('./DataFigures/VariedParam-P0','rb');
RealizedPerim = pickle.load(fileRealizedPeri);
VariedParam = pickle.load(fileVariedParam);
fileRealizedPeri.close();
fileVariedParam.close();

#%%Plot the figure
fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(2,2))
categories = np.where(VariedParam < 3.81, 0,1);
mask = np.ma.masked_where(VariedParam < 3.81, VariedParam)
colormap =np.array([colorsolid,colorliq])
ax.scatter(VariedParam,RealizedPerim,color=colormap[categories])
ax.plot(mask,mask,color=_new_black)
ax.set_xlim(3,4.4)
ax.set_ylim(3.7,4.4)
ax.set_xlabel(r'$s^0$')
ax.set_ylabel(r'$s$')

fig.savefig("realizedPerim-P0",format='pdf',bbox_inches='tight',pad_inches=0)


#%% Figure about the new phase

# Load the data
fileNewPhase = open('./Extracted_Data/Data_NewPhase_P0','rb');
fileNewPhaseParams = open('./Extracted_Data/Params_NewPhase_P0','rb');
Data_NewPhase = pickle.load(fileNewPhase);
Params_NewPhase = pickle.load(fileNewPhaseParams);
fileNewPhase.close();
fileNewPhaseParams.close();

#%% Compute and plot data for the new phase
figNP, axesNP = plt.subplots(nrows=2,ncols=2,figsize=(6.5,4),sharex='col')
figNP.subplots_adjust(hspace=0.4)

HomeoPressure = np.zeros(len(Data_NewPhase));
VariedParam = np.zeros(len(Data_NewPhase));
slopes = np.zeros(len(Data_NewPhase));
for i in range(len(Data_NewPhase)):
    p02 = Params_NewPhase[i][1][1];
    VariedParam[i] = p02;
    #Average out last quarter of the simulation to compute pressure
    for j in range(int(3*Nt/4),Nt):
        cellpressure = reader.CellPressure(j, Data_NewPhase[i])
        HomeoPressure[i] = HomeoPressure[i] + np.average(cellpressure);
    HomeoPressure[i] = HomeoPressure[i]/(Nt-int(3*Nt/4));
    #Compute slope
    AreaOTime, PeriOTime, Type = analysis.AvAreaAndPerim(Data_NewPhase[i]);
    Peri2OTime = np.power(PeriOTime,2);
    resLinReg = linregress(Peri2OTime,AreaOTime);
    slopes[i] = resLinReg.slope
    
#Do the plots
indOP = 10;
indNP = 24;
reader.MakePlot(len(Data_NewPhase[indOP])-1,Data_NewPhase[indOP] , axesNP[0,1])
axesNP[0,1].set_title(r'b)    $p^0 = 4.76$',fontfamily='serif',loc='left')
reader.MakePlot(len(Data_NewPhase[indNP])-1,Data_NewPhase[indNP] , axesNP[1,1])
axesNP[1,1].set_title(r'd)    $p^0 = 5.61$',fontfamily='serif',loc='left')


axesNP[0,0].scatter(VariedParam,HomeoPressure)
axesNP[0,0].set_ylabel(r"Pressure, $\mathcal{P}$")
axesNP[0,0].set_title(r'a)',fontfamily='serif',loc='left')

#axesNP[0,0].set_xlabel(r"Preferred perimeter, $p^0$")

axesNP[1,0].scatter(VariedParam,slopes)
axesNP[1,0].set_ylabel(r"Area susceptibility, $\chi$")
axesNP[1,0].set_xlabel(r"Preferred perimeter, $p^0$")
axesNP[1,0].set_title(r'c)',fontfamily='serif',loc='left')



figNP.savefig("NewPhase.pdf",format='pdf',bbox_inches='tight')


#%% Load data supplementary figure on homogeneous populations
file1a = open('./DataFigures/HomogeneousRef','rb');
Fig1a = pickle.load(file1a);
file1a.close()

#%% Plot homogeneous area and perimeter with linear regression fit
Nfile = len(Fig1a)
slopes = np.zeros(Nfile);
fig = plt.figure(figsize=(6.5,2))
ax = fig.subplots(1,2)
fig.subplots_adjust(wspace=0.3)
colors = [colorsolid,colorliq]
for i in range(Nfile):
    AreaOTime, PeriOTime, Type = analysis.AvAreaAndPerim(Fig1a[i]);
    Peri2OTime = np.power(PeriOTime,2);
    p = np.linspace(min(PeriOTime),max(PeriOTime))
    p2 = np.linspace(min(Peri2OTime),max(Peri2OTime))
    pm = np.mean(PeriOTime)
    #print(pm)
    pm2 = np.mean(Peri2OTime)
    if i==0:
        resLinReg = linregress(Peri2OTime,AreaOTime);
        slopes[i] = resLinReg.slope
        ax[i].scatter(Peri2OTime,AreaOTime,color=colors[i])
        ax[i].set_xlabel(r'squared perimeter, $p^2$')
        ax[i].plot(p2,slopes[i]*(p2-pm2)+1,color=_new_black)
    elif i==1:
        resLinReg = linregress(PeriOTime,AreaOTime);
        slopes[i] = resLinReg.slope
        ax[i].scatter(PeriOTime,AreaOTime,color=colors[i])
        ax[i].set_xlabel(r'perimeter, $p$')
        ax[i].plot(p,slopes[i]*(p-pm)+1,color=_new_black)
    ax[i].set_ylabel(r'area, $a$')
    
ax[0].set_title(r'a)',fontfamily='serif',loc='left')
ax[1].set_title(r'b)',fontfamily='serif',loc='left')
#fig.savefig("Homogeneous-SM.pdf",format='pdf',bbox_inches='tight',pad_inches=0)

