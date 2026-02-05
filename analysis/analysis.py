#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Analysing the data

Created on Fri Mar 18 13:03:04 2022

@author: louis
"""
#%% Imports
import reader
import numpy as np
import multiprocessing
import collections
import statistics as stat
import matplotlib.pyplot as plt
import multiprocessing.pool
from IPython.core.debugger import Pdb;
from scipy.spatial import Voronoi, voronoi_plot_2d
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

#%%

def fraction(celltype):
    #Return an array clone with the unique cell types and an array frac with the associated
    #fraction of unique lineages in the population
    count_cell = collections.Counter(celltype);
    N = len(count_cell);
    Ncells = len(celltype);
    frac=np.zeros(N);
    clone=np.zeros(N);
    i=0
    for elem in count_cell.keys():
        clone[i] = elem;
        frac[i] = count_cell[elem]/Ncells;
        i = i+1;
    return clone,frac;

def cloneNumber(celltype):
    # Counts the number of different clones in the population
    count_cell = collections.Counter(celltype);
    number = len(count_cell);
    return number;

def frame_circle(cellpos,center,radius):
    #Extract the index of the cells of cellpos that are in the circle of center 'center' and of radius 'radius'
    Ncells = len(cellpos[:,0]);
    frame_data=[];
    for i in range(Ncells):
        r = np.power((cellpos[i,0]-center[0])**2+(cellpos[i,1]-center[1])**2,0.5);
        if r<=radius:
            frame_data.append(i);
    return frame_data;

def IndInBox(cellpos,center,x1,x2,y1,y2):
    #Extract the index of the cells of cellpos that are in the circle of center 'center' and of radius 'radius'
    Ncells = len(cellpos[:,0]);
    frame_data=[];
    for i in range(Ncells):
        x = cellpos[i,0]
        y = cellpos[i,1]
        if x > x1 and x < x2 and y>y1 and y<y2 :
            frame_data.append(i);
    return frame_data;

def IndInX(cellpos,xgrid):
    #Extract all the indices of cells according to a spatial xgrid
    Ngrid = len(xgrid);
    dx = (xgrid[0] - xgrid[-1])/Ngrid;
    cellindX = [[] for i in range(Ngrid)];
    for i in range(len(cellpos)):
        posx = cellpos[i,0];
        ind = int(posx//dx);
        cellindX[ind].append(i)
    return cellindX

def Pressure(cellarea,cellperi,celltype,K,G1,G2,A0,p01,p02):
    #Returns the pressure with the same list structure as the cells
    N = len(cellarea);
    P = np.zeros(N);
    for i in range(N):
        if celltype[i] == 0:
            P[i] = cellarea[i] - 1 + G1*cellperi[i]/2/cellarea[i]*(cellperi[i]-p01);
        elif celltype[i] == 1:
            P[i] = K*(cellarea[i] - A0) + G2*cellperi[i]/2/cellarea[i]*(cellperi[i]-p02);
    return P

def IndInY(cellpos,ygrid):
    #Extract all the indices of cells according to a spatial xgrid
    Ngrid = len(ygrid);
    dy = (ygrid[0] - ygrid[-1])/Ngrid;
    cellindY = [[] for i in range(Ngrid)];
    for i in range(len(cellpos)):
        posy = cellpos[i,1];
        ind = int(posy//dy);
        cellindY[ind].append(i)
    return cellindY

def AverageValuePerType(cellarea,celltype):
    #Computes the average area for cells of each type specified by celltype. Works only for 2 types 0 and 1 now
    AvArea = np.zeros((2,2));
    AreaType = ValuePerType(cellarea, celltype);
    AvArea[0,0] = np.average(AreaType[0]);
    AvArea[0,1] = np.std(AreaType[0]);
    AvArea[1,0] = np.average(AreaType[1]);
    AvArea[1,1] = np.std(AreaType[1]);
    return AvArea

def AverageAreaOTime(data):
    Ntime = len(data)
    AvArea = np.zeros((Ntime,2))
    for i in range(Ntime):
        cellarea = reader.CellArea(i, data)
        celltype = reader.CellType(i, data)
        AvAreaTemp = AverageValuePerType(cellarea, celltype)
        AvArea[i,0] = AvAreaTemp[0,0]
        AvArea[i,1] = AvAreaTemp[1,0]
    return AvArea

def AverageValue(val):
    AvVal = np.zeros(2);
    AvVal[0] = np.average(val);
    AvVal[1] = np.std(val);
    return AvVal;
    
    
    # N = len(cellarea);
    # Nc = cloneNumber(celltype);
    # AvArea = np.zeros(Nc);
    # j0=0; j1=0;
    # # count_cell = collections.Counter(celltype);
    # if (N != len(celltype)):
    #     print('cellarea and cell type must have the same size');
    # for i in range(N):
    #     if (celltype[i]==0):
    #         AvArea[0] = AvArea[0] + cellarea[i];
    #         j0 = j0+1;
    #     elif (celltype[i]==1):
    #         AvArea[1] = AvArea[1] + cellarea[i];
    #         j1 = j1 +1;
    # AvArea[0] = AvArea[0]/j0;
    # AvArea[1] = AvArea[1]/j1;
    # return AvArea

def ValuePerType(cellarea, celltype):
    #Returns a value per type, 0 and 1. Coded with cellarea but applies to everything
    N = len(cellarea);
    cellarea0=[]; cellarea1=[];
    for i in range(N):
        if (celltype[i] == 0):
            cellarea0.append(cellarea[i]);
        elif (celltype[i] == 1):
            cellarea1.append(cellarea[i]);
    return [cellarea0, cellarea1]

def PlotAreaVsPerim2(cellarea, cellperimeter, celltype,a1,a2,p1,p2):
    #For a given population of cellarea and celltype, plots area vs p^2, 
    #and compare it with the regular polygons
    fig=plt.figure();
    ax = fig.add_subplot(111);
    perimeters = ValuePerType(cellperimeter, celltype);
    areas = ValuePerType(cellarea, celltype);
    SqPerimeter0 = np.power(perimeters[0],2);
    SqPerimeter1 = np.power(perimeters[1],2);
    ax.scatter(SqPerimeter0,areas[0],color='blue', label='Resident',alpha=0.5);
    ax.scatter(SqPerimeter1,areas[1],color='red', label='Mutant',alpha=0.5);
    PerimeterAxis = np.linspace(min(min(SqPerimeter0),min(SqPerimeter1)),max(max(SqPerimeter0),max(SqPerimeter1)),100);
    AreaRegularPentagon = 0.0688191*PerimeterAxis;
    AreaRegularHexagon = 0.0721688*PerimeterAxis;
    ax.set_xlabel(r'$p^2$', fontsize='large');
    ax.set_ylabel(r'a', fontsize='large');
    ax.set_xticks([p1**2, p2**2, 14,16,18])
    ax.set_xticklabels([r'${p_r^0}^2$', r'${p_m^0}^2$', 14, 16,18],fontsize='large')
    ax.set_yticks([0.8, a1, a2, 1.2]);
    ax.set_yticklabels([0.8, r'$a_r^0$',r'$a_m^0$', 1.2],fontsize='large');
    ax.axhline(y=a1,color='blue',linestyle='dashed',alpha=0.5);
    ax.axhline(y=a2,color='red',linestyle='dashed',alpha=0.5);
    ax.axvline(x=p1**2,color='blue',linestyle='dashed',alpha=0.5);
    ax.axvline(x=p2**2,color='red',linestyle='dashed',alpha=0.5);
    #ax.plot(PerimeterAxis,AreaRegularPentagon,linestyle='dashed',color='grey',alpha=1,label=r'$a = \frac{1}{24 \tan{\pi/6}} p^2$');
    ax.legend(shadow=True,fontsize=16);
    ax.plot(PerimeterAxis,AreaRegularHexagon,color='blue');
    
def PlotNeighNum(cellneighnum):
    #Takes as input an array of the number of neighbors for each cells and plots a histogram
    labels, values = zip(*collections.Counter(cellneighnum).items())
    width=1
    plt.bar(labels,values,width)
    plt.xticks(labels,labels)
    plt.show()
    
def PlotAreaVsPerim(cellarea, cellperimeter, celltype,a1,a2,p1,p2):
    #For a given population of cellarea and celltype, plots area vs p^2, 
    #and compare it with the regular polygons
    fig=plt.figure();
    ax = fig.add_subplot(111);
    perimeters = ValuePerType(cellperimeter, celltype);
    areas = ValuePerType(cellarea, celltype);
    ax.scatter(perimeters[0],areas[0],color='blue', label='Resident');
    ax.scatter(perimeters[1],areas[1],color='red', label='Mutant');
    PerimeterAxis = np.linspace(min(min(perimeters[0]),min(perimeters[1])),max(max(perimeters[0]),max(perimeters[1])),100);
    AreaRegularPentagon = 0.0688191*PerimeterAxis;
    AreaRegularHexagon = 0.0721688*PerimeterAxis;
    ax.set_xlabel(r'$p^2$', fontsize=26);
    ax.set_ylabel(r'a', fontsize=26);
    ax.set_xticks([p1**2, p2**2, 14,16,18])
    ax.set_xticklabels([r'${p_r^0}^2$', r'${p_m^0}^2$', 14, 16,18])
    ax.set_yticks([0.8, a1, 1.2]);
    ax.set_yticklabels([0.8, r'$a^0$', 1.2]);
    ax.axhline(y=a1,color='blue',linestyle='dashed',alpha=0.5);
    ax.axhline(y=a2,color='red',linestyle='dashed',alpha=0.5);
    ax.axvline(x=p1**2,color='blue',linestyle='dashed',alpha=0.5);
    ax.axvline(x=p2**2,color='red',linestyle='dashed',alpha=0.5);
    #ax.plot(PerimeterAxis,AreaRegularPentagon,linestyle='dashed',color='grey',alpha=1,label=r'$a = \frac{1}{24 \tan{\pi/6}} p^2$');
    ax.legend(shadow=True,fontsize=16);
    #ax.plot(PerimeterAxis,AreaRegularHexagon,color='blue');
    
def AvAreaAndPerim(data):
    #Only works when there is no growth and a constant number of cells
    Ntime = len(data)
    Ncells = int(data[0][0][0]); #Assume that the number of cells is constant
    TotArea = np.zeros(Ncells)
    TotPeri = np.zeros(Ncells)
    for i in range(Ntime//2,Ntime): #Compute only for the second half of simulation
        cellarea = reader.CellArea(i, data)
        cellperimeter = reader.CellPeri(i, data)        
        TotArea[:] = TotArea[:] + cellarea
        TotPeri[:] = TotPeri[:] + cellperimeter
    celltype = reader.CellType(Ntime//2,data)
    AvArea = TotArea/(Ntime-Ntime//2)
    AvPeri = TotPeri/(Ntime-Ntime//2)
    return AvArea, AvPeri, celltype
    
def HasFixated(TimeFrac0):
    #Input is a 2d list TimeFraction[i][0,1], where i is the time point,
    #0 and 1 the type for which the fraction is given
    FIXATION = False;
    TypeFixated = -1; #Negative value if no type has fixated
    N = len(TimeFrac0);
    for i in range(N):
        if TimeFrac0[i] == 0:
            FIXATION = True;
            TypeFixated = 1;
            break
        elif TimeFrac0[i] == 1:
            FIXATION = True;
            TypeFixated = 0;
            break
    return FIXATION,TypeFixated
    

def HasPartiallyFixated(TimeFrac0,gate):
    #Input is a 2d list TimeFraction[i][0,1], where i is the time point,
    #0 and 1 the type for which the fraction is given
    FIXATION = False;
    TypeFixated = -1; #Negative value if no type has fixated
    N = len(TimeFrac0);
    for i in range(N):
        if TimeFrac0[i] < gate:
            FIXATION = True;
            TypeFixated = 1;
            break
        elif TimeFrac0[i] >= 1-gate:
            FIXATION = True;
            TypeFixated = 0;
            break
    return FIXATION,TypeFixated

def LinearFit(fracfull, T, timesteps):
    #Takes fraction of several simulations, the times of the simulation, and 
    #the number of identical simulations filecycle, 
    #and the final time to do the interpolation
    func = np.log(fracfull[0:timesteps]) - np.log(1-fracfull[0:timesteps]);
    coeff = np.polyfit(T[0:timesteps],func,1);
    return coeff