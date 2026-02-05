#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Analyse population dynamics data
Created on Wed Mar 16 14:09:06 2022

@author: louis
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
plt.rcParams["figure.autolayout"] = True
plt.rcParams['text.usetex'] = True
#plt.style.use('seaborn-poster')

#%% Read files
root = tk.Tk()
root.withdraw()
filenames = fd.askopenfilenames();
Nfile = len(filenames);
Data=[];
Params=[];
for i in range(Nfile):
    #Data.append(reader.ReadFile(filenames[i])[0]);
    #Params.append(reader.ReadFile(filenames[i])[1]);
    Data.append(readFile.ReadFileC(filenames[i])[0]);
    Params.append(readFile.ReadFileC(filenames[i])[1]);
    
#%% Pickle the data
fileData = open('./Data_NewPhase_P0','w+b');
fileParams = open('./Params_NewPhase_P0','w+b');
#fileData = open('./Extracted_Data/Data_100exp_alpha03-2_p0_265-465','rb');
#fileParams = open('./Extracted_Data/Params_100exp_alpha03-2_p0_265-465','rb');
pickle.dump(Data,fileData);
pickle.dump(Params,fileParams);
fileData.close()
fileParams.close()
#Data = pickle.load(fileData);
#Params = pickle.load(fileParams);

#%% Load pickled data
fileDataOther = open('./Extracted_Data/Invasion_SolLiq_pressure_20240626','rb');
fileParamsOther = open('./Extracted_Data/Params_SolLiq_pressure_20240626','rb');
DataOther= pickle.load(fileDataOther);
ParamsOther = pickle.load(fileParamsOther);
#Nfile=len(Data);

#%% Add old data to existing data
Data.extend(DataOther);
Params.extend(ParamsOther);
Nfile=len(Data);

#%% Choose wether to check values from one experiment, 0, or averaged values over multiple exp, !=0
switch = 1;


#%% Time in the simulation
if(Nfile == 1):
    filenum = 0;
    Nt = len(Data[filenum]);
    T=np.zeros(Nt);
    for i in range(Nt):
        T[i] = reader.SimTime(i, Data[filenum]);
        
        
    #%% Time loop to get fraction of cells in the entire population
    clone = [];
    frac = [];
    for i in range(Nt):
        #Compute the frequency
        celltype = reader.CellType(i, Data[filenum]);
        tempclone, tempfrac = analysis.fraction(celltype);
        clone.append(tempclone);
        frac.append(tempfrac)
    #%% Compute average area and perimeters, shape index, full pressure, and average birth rate
    filenum = 0;
    AvArea = np.zeros((Nt,2));
    AvPeri = np.zeros((Nt,2));
    AvShape = np.zeros((Nt,2));
    AvFullPressure = np.zeros((Nt,2));
    AvPressureArea = np.zeros((Nt,2));
    AvPressurePeri = np.zeros((Nt,2));
    AvArea2 = np.zeros((Nt,2));
    BirthRateDiff = np.zeros(Nt);
    K = Params[filenum][1][6];
    G1 = Params[filenum][1][2];
    G2 = Params[filenum][1][3];
    lamb = Params[filenum][1][8];
    divrate = Params[filenum][1][9];
    p01 = Params[filenum][1][0];
    p02 = Params[filenum][1][1];
    A0 = Params[filenum][1][4];
    for i in range(Nt):
        cellarea = reader.CellArea(i, Data[filenum]);
        cellperi = reader.CellPeri(i, Data[filenum]);
        celltype = reader.CellType(i, Data[filenum]);
        cellshape = np.divide(cellperi,np.power(cellarea,0.5));
        cellarea_1,cellarea_2 = analysis.ValuePerType(cellarea, celltype);
        cellperi_1,cellperi_2 = analysis.ValuePerType(cellperi, celltype);
        pressureperi1 = G1/2*np.divide(np.multiply(cellperi_1,cellperi_1-p01*np.ones(len(cellperi_1))),cellarea_1)
        pressurearea1 = K*(cellarea_1-np.ones(len(cellarea_1)))
        pressureFull_1 = pressurearea1 + pressureperi1
        pressureperi2 = G2/2*np.divide(np.multiply(cellperi_2,cellperi_2-p02*np.ones(len(cellperi_2))),cellarea_2)
        pressurearea2 = K*(cellarea_2-np.ones(len(cellarea_2)))
        pressureFull_2 = pressurearea2 + pressureperi2
        AvFullPressure[i,0] = np.average(pressureFull_1);
        AvPressureArea[i,0] = np.average(pressurearea1);
        AvPressurePeri[i,0] = np.average(pressureperi1);
        AvFullPressure[i,1] = np.average(pressureFull_2);
        AvPressureArea[i,1] = np.average(pressurearea2);
        AvPressurePeri[i,1] = np.average(pressureperi2);
        birthrate1_temp = divrate/0.01*np.exp(2*K*lamb*(cellarea_1-np.ones(len(cellarea_1))));
        birthrate2_temp = divrate/0.01*np.exp(2*K*lamb*(cellarea_2-A0*np.ones(len(cellarea_2))));
        BirthRateDiff[i] = np.average(birthrate2_temp)-np.average(birthrate1_temp);
        AvArea[i,:] = analysis.AverageValuePerType(cellarea, celltype)[:,0];
        AvPeri[i,:] = analysis.AverageValuePerType(cellperi, celltype)[:,0];
        AvShape[i,:] = analysis.AverageValuePerType(cellshape, celltype)[:,0];
        cell2_1 = np.power(cellarea_1-np.ones(len(cellarea_1)),2);
        cell2_2 = np.power(cellarea_2-Params[filenum][1][4]*np.ones(len(cellarea_2)),2);
        AvArea2[i,0] = np.average(cell2_1);
        AvArea2[i,1] = np.average(cell2_2);
        
    #%% Plotting the average full pressure over time.
    fig = plt.figure();
    ax = fig.add_subplot(1, 1, 1);
    ax.plot(T,-AvFullPressure[:,0],color='blue',label='Total pressure ',linewidth=0.75)
    ax.plot(T,-AvPressureArea[:,0],color='blue',linestyle= 'dotted', label='Area contribution')
    ax.plot(T,-AvPressurePeri[:,0],color='blue',linestyle= 'dashed',label='Perimeter contribution')
    #ax.plot(T,np.average(AvFullPressure[:,0])*np.ones(Nt),color='blue',label='Average resident ',marker='+')
    #ax.plot(T,-AvFullPressure[:,1],color='red',label='mutant',linewidth=0.75)
    #ax.plot(T,AvPressureArea[:,1],color='red',linestyle= 'dotted', label='mutant, area only ')
    #ax.plot(T,AvPressurePeri[:,1],color='red',linestyle= 'dashed',label='mutant, peri only ')
    #ax.plot(T,np.average(AvFullPressure[:,1])*np.ones(Nt),color='red',label='Average mutant ',marker='+')
    ax.set_title(r'Homogeneous population, random b-d process')
    ax.set_xlabel(r'Time (c.u)')
    ax.set_ylabel(r'Pressure (in units of $K A_0$)')
    ax.legend()
    fig.show()
    
    #%% Plotting the shape parameter over time
    fig = plt.figure();
    ax = fig.add_subplot(1, 1, 1);
    ax.plot(T,AvShape[:,0],color='blue',label=' ',linewidth=0.75)
    ax.set_title(r'Shape parameter')
    ax.set_xlabel(r'Time (c.u)')
    #ax.set_ylabel(r'')
    ax.legend()
    fig.show()
    
    #%% Testing things about spatial position
    Nx=21;
    K = Params[filenum][1][6];
    G1 = Params[filenum][1][2];
    G2 = Params[filenum][1][3];
    p01 = Params[filenum][1][0];
    p02 = Params[filenum][1][1];
    A0 = Params[filenum][1][4];
    #tstamp = 150;
    xgrid = np.linspace(0, Data[filenum][0][0][2],Nx);
    fractionX = np.zeros((Nt,Nx));
    PressureX = np.zeros((Nt,Nx));
    for j in range(Nt):
        cellpos = reader.CellPos(j, Data[filenum]);
        cellarea = reader.CellArea(j, Data[filenum]);
        cellperi = reader.CellPeri(j, Data[filenum]);
        celltype = reader.CellType(j, Data[filenum]);
        pressure = analysis.Pressure(cellarea, cellperi, celltype, K, G1, G2, A0, p01, p02)
        IndX = analysis.IndInX(cellpos, xgrid);
        for i in range(Nx):
            PressureX[j,i] = np.average(pressure[IndX[i]]);
            clone,frac = analysis.fraction(celltype[IndX[i]]);
            if len(clone) == 2:
                if clone[0] == 0:
                    fractionX[j,i] = frac[1] #Only keep the fraction of mutant
                elif clone[0] ==1:
                    fractionX[j,i] = frac[0]
            elif len(clone) == 1:
                if clone[0] == 0:
                    fractionX[j,i] = 0;
                elif clone[0] ==1:
                    fractionX[j,i] = 1;
    #%%Plotting the fraction as a function of space
    #for i in range(int(Nt/100)):
    #    plt.plot(xgrid,PressureX[i,:])
    plt.contourf(xgrid,T,fractionX)
    plt.xlabel('Position')
    plt.ylabel('time')
    plt.colorbar(label='frequency')
    
        

    #%% Non-linear (exponential) and Linear fitness advantage 
    #Moving average for area
    window = 1;
    w1 = np.zeros(Nt);
    w2 = np.zeros(Nt);
    s_timepoint = np.zeros(Nt);
    s_expT = np.zeros(Nt)
    s_simpson = np.zeros(Nt);
    BirthRateDiffSimp = np.zeros(Nt);
    dt = T[Nt-1]-T[Nt-2];
    for j in range(Nt-1):
        s_timepoint[j] = 2*K*lamb*divrate/0.01*(AvArea[j,1]-AvArea[j,0]+1-Params[filenum][1][4])*dt
        s_simpson[j] = 1/6*(s_timepoint[j] + s_timepoint[j+1] + 4*(s_timepoint[j] + s_timepoint[j+1])/2)
        BirthRateDiffSimp[j] = dt/6*(BirthRateDiff[j] + BirthRateDiff[j+1] + 4*(BirthRateDiff[j] + BirthRateDiff[j+1])/2)
    s_timepoint[Nt-1] = 2*K*lamb*divrate/0.01*(AvArea[Nt-1,1]-AvArea[Nt-1,0]+1-Params[filenum][1][4])*dt
    s_simpson[Nt-1] = 1/6*(s_timepoint[Nt-1] + s_timepoint[Nt-2] + 4*(s_timepoint[Nt-1] + s_timepoint[Nt-2])/2)
    BirthRateDiffSimp[Nt-1] = dt/6*(BirthRateDiff[Nt-1] + BirthRateDiff[Nt-2] + 4*(BirthRateDiff[Nt-1] + BirthRateDiff[Nt-2])/2)
    s_area = np.cumsum(s_simpson)
    s_birthrate = np.cumsum(BirthRateDiffSimp)
    #s_area[Nt-1] = np.nan;
    # Ncells = 500;
    # for j in range(Nt-window):
    #     w1[j] = np.average(AvArea[j:j+window,0]-1);
    #     w2[j] = np.average(AvArea[j:j+window,1])-Params[filenum][1][4]
    # for j in range(window):
    #     w1[Nt-1-j] = np.nan;
    #     w2[Nt-1-j] = np.nan;
    # s_area=0.05/1000*(-np.average(AvArea[:,0]-1)+np.average(AvArea[:,1])-Params[filenum][1][4])/(np.average(AvArea[:,0]-1)-Params[filenum][1][4]+np.average(AvArea[:,1]));
    # s_vect_area=100*2*K*lamb*divrate/2/Ncells*(-AvArea[:,0]+1+AvArea[:,1]-Params[filenum][1][4]);
    # s_vect_area2 = 100*2*np.power(K*lamb,2)*divrate/2/Ncells*(AvArea2[:,1] - AvArea2[:,0] )
    # #s_running_av = 2*K*lamb*divrate/2/Ncells*np.divide(w2-w1,w1+w2);
    # s_running_av = 100*2*K*lamb*divrate/2/Ncells*(w2-w1);
    # s_integral = integrate.cumulative_trapezoid(s_running_av,T,initial=np.log(0.1/0.9));
    #%% Plot
    fig = plt.figure();
    ax = fig.add_subplot(1, 1, 1);
    ax.plot(T,AvArea[:,0]-1,color='blue',label='Type 1 ')
    ax.plot(T,np.average(AvArea[:,0]-1)*np.ones(Nt),color='blue',label='AverageType 1 ',marker='+')
    ax.plot(T,AvArea[:,1]-Params[filenum][1][4],color='red',label='Type 2')
    ax.plot(T,np.average(AvArea[:,1]-Params[filenum][1][4])*np.ones(Nt),color='red',label='AverageType 2 ',marker='+')
    ax.legend()
    fig.show()
    print(-np.average(AvArea[:,0]-1)+np.average(AvArea[:,1]-Params[filenum][1][4]))
        
    #%% Plot the results.
    frac1 = np.zeros(Nt);
    frac2 = np.zeros(Nt);
    for i in range(Nt):
        if len(clone[i]) == 2:
            if clone[i][0] == 0:
                frac1[i] = frac[i][0];
                frac2[i] = frac[i][1];
            elif clone[i][0] ==1:
                frac1[i] = frac[i][1];
                frac2[i] = frac[i][0];
        elif len(clone[i]) == 1:
            if clone[i][0] == 0:
                frac1[i] = frac[i][0];
            elif clone[i][0] ==1:
                frac2[i] = frac[i][0];
    plt.scatter(T,np.log(np.divide(frac2,1-frac2)),color='red',label='Type 2')
    #plt.plot(T,s_area*T+np.log(0.1/0.9)*np.ones(Nt),color="black")
    #plt.scatter(T,s_area,color="blue")
    #plt.scatter(T,s_birthrate,color="blue")
    #plt.scatter(T,np.multiply(s_running_av,T)+np.log(0.1/0.9)*np.ones(Nt),color="black")
    #plt.scatter(T,np.multiply(s_vect_area + s_vect_area2,T)+np.log(0.1/0.9)*np.ones(Nt),color="black")
    #plt.scatter(T,frac2,color='red',label='Fluid-like')
    plt.legend()
    
    #%% Plot the number of cells
    NumCells = reader.CellNumOTime(Data[filenum])[1];
    plt.plot(T,NumCells)
    plt.title('Cell number')

    
    #%% Is there fixation ?
    fixe = analysis.HasFixated(frac1);
    
    #%% Time loop to get fraction of cells in a fixed circle
    
    #Define region of interest
    center = [10,10];
    radius = 10;
    
    IndOfInterest = [];
    clone = [];
    frac = [];
    for i in range(Nt):
        #Gets the indices of cells in circle
        position = reader.CellPos(i, Data[filenum])
        IndOfInterest.append(analysis.frame_circle(position, center, radius));
        #Compute the frequency
        celltype = reader.CellType(i, Data[filenum]);
        tempclone, tempfrac = analysis.fraction(celltype[IndOfInterest[i]]);
        clone.append(tempclone);
        frac.append(tempfrac)
    
    #%% Plot the results
    frac1 = np.zeros(Nt);
    frac2 = np.zeros(Nt);
    for i in range(Nt):
        frac1[i] = frac[i][0];
        frac2[i] = frac[i][1];
    plt.scatter(T,frac1,color='blue',label='Solid-like')
    plt.scatter(T,frac2,color='red',label='Fluid-like')
    plt.legend()
    
    #%% Plot area vs perimeter^2 with lines of regular pentagons
    timeframe = int((Nt-1)/8);
    filenum = 55;
    analysis.PlotAreaVsPerim2(reader.CellArea(timeframe, Data[filenum]), reader.CellPeri(timeframe, Data[filenum]), reader.CellType(timeframe, Data[filenum]),1,Params[filenum][1][4],Params[filenum][1][0],Params[filenum][1][1])
    plt.title(r'Realized area vs perimeter in a mixture')
    
    #%% Plot area vs perimeter^2 averaged over time
    for i in range(Nfile):
        filenum = i;
        AreaOTime, PeriOTime, Type = analysis.AvAreaAndPerim(Data[filenum]);
        analysis.PlotAreaVsPerim2(AreaOTime, PeriOTime, Type,1,Params[filenum][1][4],Params[filenum][1][0],Params[filenum][1][1])
        #plt.title(r'Solid-Solid mixture, $a_m^0=0.9$')
        
    #%%Gather pressure computed in the simulation
    cellPressureOTime = [];
    for i in range(Nt):
        cellPressureOTime.append(reader.CellPressure(i, Data[filenum]));
    cellPressureArray = np.array(cellPressureOTime);
    AvcellPressure = np.mean(cellPressureArray,axis=0);
    #Creates histogram
    bin_values_Comp, bin_edges_Comp = np.histogram(AvcellPressure)
    bin_width_Comp = bin_edges_Comp[1] - bin_edges_Comp[0];
    kde_Comp = gaussian_kde(AvcellPressure);
    x_Comp = np.linspace(bin_edges_Comp[0], bin_edges_Comp[-1], 200);
    
    #%% Computing distribution of pressures and plotting
    K = Params[filenum][1][6];
    A0 = Params[filenum][1][4];
    G1 = Params[filenum][1][2];
    p01 = Params[filenum][1][0];    
    
    #PressureHydroOTime = -K*(AreaOTime-A0*np.ones(len(AreaOTime)))
    PressureOTime0 = (-K*(AreaOTime-A0*np.ones(len(AreaOTime))) - G1/2*np.divide(np.multiply(PeriOTime,PeriOTime-p01*np.ones(len(PeriOTime))),AreaOTime)); #Density of pressure
    #sns.kdeplot(data=PressureOTime,x='Pressure')
    bin_values0, bin_edges0 = np.histogram(PressureOTime0)
    bin_width0 = bin_edges0[1] - bin_edges0[0];
    kde0 = gaussian_kde(PressureOTime0);
    x0 = np.linspace(bin_edges0[0], bin_edges0[-1], 200)
    
    PressureOTime1 = (-K*(AreaOTime-A0*np.ones(len(AreaOTime))) - G1/2*np.divide(np.multiply(PeriOTime,PeriOTime-p01*np.ones(len(PeriOTime))),AreaOTime)); #Density of pressure
    #sns.kdeplot(data=PressureOTime,x='Pressure')
    bin_values1, bin_edges1 = np.histogram(PressureOTime1)
    bin_width1 = bin_edges1[1] - bin_edges1[0];
    kde1 = gaussian_kde(PressureOTime1);
    x1 = np.linspace(bin_edges1[0], bin_edges1[-1], 200)
    
    #%% Scatter plot pressure with two different functions at some specific times
    K = Params[filenum][1][6];
    A0 = Params[filenum][1][4];
    G1 = Params[filenum][1][2];
    p01 = Params[filenum][1][0];
    filenum=0;
    timepoint = Nt-1;
    PressureInfered = (-K*(reader.CellArea(timepoint, Data[filenum])-A0) - G1/2*np.divide(np.multiply(reader.CellPeri(timepoint, Data[filenum]),reader.CellPeri(timepoint, Data[filenum])-p01),AreaOTime)); #Density of pressure
    fig = plt.figure();
    ax = fig.add_subplot(1, 1, 1);
    ax.scatter(np.arange(0,len(PressureInfered)),cellPressureArray[timepoint,:]-PressureInfered);

    
    #%%Plotting pressure distribution
    
    x = np.concatenate((x0, x1));
    
    fig = plt.figure();
    ax = fig.add_subplot(1, 1, 1);
    #ax.bar(x=bin_edges[:-1], height=bin_values / len(PressureOTime), width=bin_width)
    ax.plot(x0, kde0(x0)*bin_width0 ,linestyle='dotted',label='Liquid-like',color='black')
    ax.plot(x1, kde1(x1)*bin_width1 ,linestyle='solid',label='Solid-like',color='black')
    ax.plot(x_Comp,kde_Comp(x_Comp)*bin_width_Comp,color='red',label='Computed')
    #ax2.plot(x1, kde1(x1) * bin_width1,linestyle='solid')
    ax.set_ylabel('Probability')
    ax.set_xlabel('Pressure')
    ax.legend()
    #ax2.set_ylabel('Probability')
    #ax2.set_xlabel('Pressure')
   
    #sns.kdeplot(x=PressureOTime,common_norm=True)
    #plt.hist(PressureOTime,stacked=True,density=True,histtype='bar')
    #plt.xlabel('Pressure')

    #%% Get fitness from average over time of each cells
    K = Params[filenum][1][6];
    #lamb = Params[filenum][1][8];
    lamb =1;
    #divrate = Params[filenum][1][9];
    divrate = 0.000125;
    A0 = Params[filenum][1][4];
    AverageAreaTotal = analysis.AverageValuePerType(AreaOTime, Type)[:,0];
    s_static = lamb*divrate/0.01*(AverageAreaTotal[1]-A0-(AverageAreaTotal[0]-1));
    

    #%%Bar plot of number of neighbors
    filenum = 0;
    cellneighbors = reader.CellNeighborNum(Nt-1, Data[filenum])
    analysis.PlotNeighNum(cellneighbors)
    
    #%% Check the total average area doesn't change and corresponds to number changes
    AvTotalArea = np.zeros(Nt);
    filenum = 1;
    Nnc = reader.CellNumOTime(Data[filenum])[1];
    for i in range(Nt):
        AvTotalArea[i] = analysis.AverageValue(reader.CellArea(i, Data[filenum]))[0];
    
    #%% Plot total average area
    fig = plt.figure();
    ax = fig.add_subplot(1, 1, 1);
    ax.plot(T,AvTotalArea,color='blue',label='Average area')
    ax.plot(T,Nnc/Nnc[0],color='red',label='Nc/N_0')
    ax.legend()
    fig.show()
#%%    
#Multiple files
if(Nfile > 1):
    Nt = len(Data[0]);
    T=np.zeros(Nt);
    filecycle=50; # Sets the number of identical experiments
    for i in range(Nt):
        #We assume that all the simulation are done with the same time points
        T[i] = reader.SimTime(i, Data[0]);
        
#%% Pickle the time
fileTime = open('./Time-InvasionLiquid-PRM-opposite','w+b');
pickle.dump(T,fileTime);
fileTime.close()

#%% Number of cells as function of homeostatic pressure
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
#%%Fit and do a plot
coef = np.polyfit(Ph, NRatio, 1);
poly1d_fn = np.poly1d(coef);
Ph_homeo = (1/Params[0][1][4]-coef[1])/coef[0];
print(Ph_homeo)
plt.plot(Ph,NRatio,'yo',Ph,poly1d_fn(Ph),'-k');
plt.xlabel(r'$P_h$')
plt.ylabel(r'$N/N_0$')



#%% Gather total number of cells for each file
coloring = ['red','blue','g','c','m','grey','k','brown'];
labeling = ['Preferred area','Preferred perimeter','Perimeter elasticity','Area elasticity','Pressure strength','Birth rate','Self-advection','Rotational diffusion'];
filecycle=50;
N0 = reader.CellNumOTime(Data[0])[1][0];
Nf = np.zeros(Nfile//filecycle);
AvNumCells = np.zeros(Nt);
for j in range(Nfile):
    count = 0;
    NumCells = reader.CellNumOTime(Data[j])[1];
    if j%filecycle==filecycle-1:
        #plt.plot(T[0:len(NumCells)],NumCells,alpha=0.5,label=j);#,color=coloring[int(j/filecycle)]labeling[int(j/filecycle)]);
        Nf[j//filecycle-1] = N0/(Params[j][1][4]);
        #plt.plot(T,Nf[j//filecycle-1]*np.ones(len(T)),color='black',label='Theory')
        AvNumCells[:] = AvNumCells[:] + NumCells[:]
        AvNumCells[:] = AvNumCells[:]/filecycle
        plt.plot(100*T,AvNumCells,alpha=1);#,color=coloring[int(j/filecycle)],label='Sim Average');
        AvNumCells[:] = np.zeros(Nt);
    else:
        #plt.plot(T[0:len(NumCells)],NumCells,alpha=0.5,color=coloring[int(j/filecycle)]);#,color=coloring[int(j/filecycle)]);
        AvNumCells[:] = AvNumCells[:] + NumCells[:]
    #plt.plot(T,500/(1-0.1*j)*np.ones(Nt))
    #plt.title(r'')
    plt.xlabel(r't',fontsize='large');
    plt.ylabel(r'N',fontsize='large');
    plt.legend()
    
    #%% Computes homeostatic pressure of homogeneous tissue as a function of parameters
    #It is setup such that it is a homogeneous tissue of mutant population
    HomeoPressure = np.zeros(Nfile);
    AvPeri = np.zeros(Nfile)
    VariedParam = np.zeros(Nfile);
    filecycle=50;
    for i in range(Nfile):
        K = Params[i][1][6];
        A0 = Params[i][1][4];
        G2 = Params[i][1][3];
        p02 = Params[i][1][1];
        v0 = Params[i][1][11];
        Dr = Params[i][1][12];
        sw = i//filecycle;
        if sw == 0:
            #VariedParam[i] = A0;
            VariedParam[i] = p02;
            #VariedParam[i] = reader.CellNumOTime(Data[i])[1][Nt-1];
        elif sw ==1:
            VariedParam[i] = A0;
            #VariedParam[i] = G2
            #VariedParam[i] = v0+i%filecycle/10;
        elif sw ==2:
            VariedParam[i] = K;
        elif sw ==3:
            VariedParam[i] = v0;
        elif sw ==4:
            VariedParam[i] = Dr;
        elif sw ==5:
            VariedParam[i] = K;
        #Average out last quarter of the simulation
        for j in range(int(2*Nt/4),Nt):
            cellarea = reader.CellArea(j, Data[i]);
            cellperi = reader.CellPeri(j, Data[i]);
            pressureperi2 = -G2/2/K*np.divide(np.multiply(cellperi,cellperi-p02*np.ones(len(cellperi))),cellarea);
            pressurearea2 = -(cellarea-A0*np.ones(len(cellarea)));
            pressureFull_2 = pressurearea2 + pressureperi2;
            cellpressure = reader.CellPressure(j, Data[i]);
            #HomeoPressure[i] = HomeoPressure[i] + np.average(pressureFull_2);
            HomeoPressure[i] = HomeoPressure[i] + np.average(cellpressure);
            AvPeri[i] = AvPeri[i] + np.average(cellperi);
        HomeoPressure[i] = HomeoPressure[i]/(Nt-int(2*Nt/4));
        AvPeri[i] =  AvPeri[i]/(Nt-int(2*Nt/4));
        
    #%%Pickle the pressure
    fileVariedParam = open('./VariedParam-newphase-P0','w+b');
    filePressure = open('./Pressure-newphase-P0','w+b');
    pickle.dump(VariedParam,fileVariedParam);
    pickle.dump(HomeoPressure,filePressure);
    fileVariedParam.close()
    filePressure.close()

    #%% Pickle the realized perimeter
    fileVariedParam = open('./VariedParam-P0','w+b');
    filePeri = open('./RealizedPerim-P0','w+b');
    pickle.dump(VariedParam,fileVariedParam);
    pickle.dump(AvPeri,filePeri);
    fileVariedParam.close()
    filePeri.close()
    
    #%% Plotting homeostatic pressure in table of plots
    coloring = ['red','blue','g','c','k','brown'];
    labeling = ['Preferred area','Preferred perimeter','Perimeter elasticity','Area elasticity','Self-advection','Rotational diffusion'];
    filecycle=50;
    nline = 2;
    ncol = 3;
    fig,axs=plt.subplots(nline,ncol, sharex=False,sharey=False)
    Nexp = int(Nfile/filecycle);
    for i in range(Nexp):
        #axs.scatter(VariedParam[i*filecycle:(i+1)*filecycle],HomeoPressure[i*filecycle:(i+1)*filecycle],color='black');
        #axs.set_title('Cell number');
        #axs.set_xlabel(r'$N_c$')
        #axs.set_ylabel(r'Homeostatic Pressure')
        axs[i//ncol,i%ncol].scatter(VariedParam[i*filecycle:(i+1)*filecycle],HomeoPressure[i*filecycle:(i+1)*filecycle],color=coloring[i]);
        axs[i//ncol,i%ncol].set_title(labeling[i]);
        axs[nline-1,i%ncol].set_xlabel(r'$t$')
        axs[i//ncol,0].set_ylabel(r'Homeostatic Pressure')
        
    #%% Plotting homeostatic pressure in single plots
    coloring = ['blue','c','k','g','red','brown'];
    labeling = ['Perimeter elasticity','Self-advection','Rotational diffusion','Preferred perimeter','Preferred area','Area elasticity'];
    filecycle=50;
    Nexp = int(Nfile/filecycle);
    for i in range(Nexp):
        plt.figure()
        plt.scatter(VariedParam[i*filecycle:(i+1)*filecycle],AvPeri[i*filecycle:(i+1)*filecycle],color=coloring[i]);
        #plt.scatter(VariedParam[i*filecycle:(i+1)*filecycle],HomeoPressure[i*filecycle:(i+1)*filecycle],color=coloring[i]);
        #plt.scatter(np.divide(16*np.ones(filecycle),np.power(VariedParam[i*filecycle:(i+1)*filecycle],2)),HomeoPressure[i*filecycle:(i+1)*filecycle],color=coloring[i]);
        plt.ylabel('Homeostatic pressure')
        plt.title(labeling[i]);
        #axs[nline-1,i%ncol].set_xlabel(r'$t$')
        #axs[i//ncol,0].set_ylabel(r'Homeostatic Pressure')
        
            
            
    #%% Time loop to get fraction of cells in the entire population
    frac1Full = np.zeros((Nt,Nfile));
    frac2Full = np.zeros((Nt,Nfile));
    for j in range(Nfile):
        for i in range(Nt):
            #Compute the frequency
            celltype = reader.CellType(i, Data[j]);
            tempclone, tempfrac = analysis.fraction(celltype);
            if len(tempclone)>1:
                if tempclone[0]==0:
                    frac1Full[i,j]=tempfrac[0];
                    frac2Full[i,j]=tempfrac[1];
                elif tempclone[0]==1:
                    frac1Full[i,j]=tempfrac[1];                
                    frac2Full[i,j]=tempfrac[0];
            elif len(tempclone)==1:
                if tempclone[0]==0:
                    frac1Full[i,j]=1;
                    frac2Full[i,j]=0;
                elif tempclone[0]==1:
                    frac1Full[i,j]=0;                
                    frac2Full[i,j]=1;
            #clone[j].append(tempclone);
            #frac[j].append(tempfrac)
    #cloneAverage = np.mean(clone,axis=1);
    frac1Average = np.mean(frac1Full,axis=1);
    frac2Average = np.mean(frac2Full,axis=1);    
    
    #%% Pickle the fraction 
    filefraction = open('./Fraction-InvasionLiquid-PRM-opposite','w+b')
    pickle.dump(frac2Full,filefraction)
    filefraction.close()
        
        
    #%% Computes the average fitness amongst the different files, according to a file cycle
    Nexp = Nfile//filecycle;
    AvArea = np.zeros((Nt,2));
    AvArea2 = np.zeros((Nt,2));
    BirthRateDiffAv = np.zeros(Nt);
    BirthRateDiffSimp = np.zeros(Nt);
    AvFullPressure = np.zeros((Nexp,Nt,2));
    AvPressureArea = np.zeros((Nexp,Nt,2));
    AvPressurePeri = np.zeros((Nexp,Nt,2));
    TimeFullPressure = np.zeros((Nt,2));
    TimePressureArea = np.zeros((Nt,2));
    for k in range(Nexp):
        for l in range(filecycle):
            K = Params[k*filecycle+l][1][6];
            lamb = Params[k*filecycle+l][1][8];
            A0 = Params[k*filecycle+l][1][4];
            G1 = Params[k*filecycle+l][1][2];
            G2 = Params[k*filecycle+l][1][3];
            p01 = Params[k*filecycle+l][1][0];
            p02 = Params[k*filecycle+l][1][1];
            divrate = Params[k*filecycle+l][1][9];
            for i in range(Nt):
                cellarea = reader.CellArea(i, Data[k*filecycle+l]);
                cellperi = reader.CellPeri(i, Data[k*filecycle+l]);
                celltype = reader.CellType(i, Data[k*filecycle+l]);
                cellarea_1,cellarea_2 = analysis.ValuePerType(cellarea, celltype);
                cellperi_1,cellperi_2 = analysis.ValuePerType(cellperi, celltype);
                pressureperi1 = G1/2*np.divide(np.multiply(cellperi_1,cellperi_1-p01*np.ones(len(cellperi_1))),cellarea_1)
                pressurearea1 = K*(cellarea_1-np.ones(len(cellarea_1)))
                pressureFull_1 = -pressurearea1 - pressureperi1
                pressureperi2 = G2/2*np.divide(np.multiply(cellperi_2,cellperi_2-p02*np.ones(len(cellperi_2))),cellarea_2)
                pressurearea2 = K*(cellarea_2-A0*np.ones(len(cellarea_2)))
                pressureFull_2 = -pressurearea2 - pressureperi2
                TimeFullPressure[i,0] = np.average(pressureFull_1);
                TimePressureArea[i,0] = np.average(pressurearea1);
                TimeFullPressure[i,1] = np.average(pressureFull_2);
                TimePressureArea[i,1] = np.average(pressurearea2);
            AvFullPressure[k,:,:] = AvFullPressure[k,:,:] + TimeFullPressure[:,:];
            AvPressureArea[k,:,:] = AvPressureArea[k,:,:] + TimePressureArea[:,:];
        AvFullPressure[k,:,:] = AvFullPressure[k,:,:]/filecycle
        AvPressureArea[k,:,:] = AvPressureArea[k,:,:]/filecycle
        
    #%% Plot single pressure
    plt.figure()
    plt.plot(100*T,AvFullPressure[0,:,0],color='blue',linewidth=0.5);
    plt.plot(100*T,AvFullPressure[0,:,1],color='red',linewidth=0.5);
    plt.plot(100*T,(Params[0][1][7] + Params[0][1][10])*np.ones(Nt),color='red',linestyle='dashed')
            
    #%% Plot the pressures in a table of plots
    coloring = ['grey','red','blue','g','c','m','k','brown'];
    labeling = ['Birth rate','Preferred area','Preferred perimeter','Perimeter elasticity','Area elasticity','Pressure strength','Self-advection','Rotational diffusion'];
    #labeling = ['Perimeter elasticity','Area elasticity','Pressure strength','Self-advection','Rotational diffusion'];
    filecycle=10;
    nline = 2;
    ncol = 4;
    fig,axs=plt.subplots(nline,ncol, sharex=True,sharey=True)
    tt=0;
    for i in range(Nexp):
        #axs[i//ncol,i%ncol].plot(100*T,AvPressureArea[i,:,0],color='blue',linestyle='dashed',linewidth=0.5);
        axs[i//ncol,i%ncol].plot(100*T,AvFullPressure[i,:,0],color='blue',linewidth=0.5);
        axs[i//ncol,i%ncol].plot(100*T,AvFullPressure[i,:,1],color='red',linewidth=0.5);
        #axs[i//ncol,i%ncol].plot(100*T,AvPressureArea[i,:,1],color='red',linestyle='dashed',linewidth=0.5);
        axs[i//ncol,i%ncol].set_title(labeling[i]);
        axs[nline-1,i%ncol].set_xlabel(r'$t$')
        axs[i//ncol,0].set_ylabel(r'$\frac{f}{1-f}$')
        #axs[exp//ncol,exp%ncol].set_ylim(10**(-1.5),10**1.5)
        #axs[exp//ncol,exp%ncol].autoscale(enable=False)
        #plt.scatter(T,np.log(np.divide(frac2Full[:,i],1-frac2Full[:,i])),label=labeling[int(i/filecycle)],marker=("+"),color=coloring[int(i/filecycle)])

#%%    Computes the average fitness amongst the different files
    dt=0.01;
    for k in range(Nfile):
        BirthRateDiff = np.zeros(Nt);
        for i in range(Nt):
            cellarea = reader.CellArea(i, Data[filenum]);
            celltype = reader.CellType(i, Data[filenum]);
            cellarea_1,cellarea_2 = analysis.ValuePerType(cellarea, celltype);
            birthrate1_temp = divrate/0.01*np.exp(lamb*(cellarea_1-np.ones(len(cellarea_1))));
            birthrate2_temp = divrate/0.01*np.exp(lamb*(cellarea_2-A0*np.ones(len(cellarea_2))));
            BirthRateDiff[i] = np.average(birthrate2_temp)-np.average(birthrate1_temp);
            AvArea[i,:] = analysis.AverageValuePerType(cellarea, celltype)[:,0];
            cell2_1 = np.power(cellarea_1-np.ones(len(cellarea_1)),2);
            cell2_2 = np.power(cellarea_2-Params[filenum][1][4]*np.ones(len(cellarea_2)),2);
            AvArea2[i,0] = np.average(cell2_1);
            AvArea2[i,1] = np.average(cell2_2);
        BirthRateDiffAv[:] = BirthRateDiffAv[:] + BirthRateDiff[:];
    for j in range(1,Nt-1):
        BirthRateDiffSimp[j] = 1/6*(BirthRateDiffAv[j] + BirthRateDiffAv[j+1] + 4*(BirthRateDiffAv[j] + BirthRateDiffAv[j+1])/2)
    s_birth = np.cumsum(BirthRateDiffSimp);
    
    
    #%% Area over time
    filenum = 0;
    AreaOTime = analysis.AverageAreaOTime(Data[filenum]);
    fitnessAdvOTime = (AreaOTime[:,1] - AreaOTime[:,0] + 1 - Params[filenum][1][4])
    
    #%% Computes the theoretical fitness advantage in solid case
    s_theory = np.zeros(Nfile)
    for filenum in range(Nfile):
        K = Params[filenum][1][6];
        Gamma = Params[filenum][1][2];
        #lamb = Params[filenum][1][8];
        lamb =1;
        c = (6*math.tan(math.pi/6))/2;
        divrate = Params[filenum][1][9];
        A0 = Params[filenum][1][4];
        s_theory[filenum] = lamb*divrate/0.01*(1 - A0)*(Gamma/K*c**2)/(1+Gamma/K*c**2);
    
    #%% Plot the results
    #plt.scatter(T,np.log(np.divide(frac2Average,1-frac2Average)),color='black',label='Liquid-like, Average')
    coloring = ['red','c','g','m','grey','k','c'];
    labeling = ['0.1','1','10','2','4','Perim'];
    frac2FullAverageVec = [];
    fig,ax = plt.subplots(1,1,sharey=True,sharex=True);
    filecycle=10;
    #ax2 = ax.twinx();
    #ax2 = ax.twiny();
    for i in range(Nfile):
        exp = i//filecycle;
        if i%filecycle==0:
            #ax.plot(T,np.log(np.divide(frac2Full[:,i],1-frac2Full[:,i])),marker=("+"),color=coloring[int(i/filecycle)],alpha=0.5,linewidth=0.5)
            frac2FullAverage = np.mean(frac2Full[:,filecycle*exp:filecycle*exp+filecycle-1],axis=1)
            ax.plot(T,np.divide(frac2FullAverage,1-frac2FullAverage),label=labeling[exp],marker=("+"),color=coloring[int(i/filecycle)],alpha=1,linewidth=2)
            frac2FullAverageVec.append(frac2FullAverage);
        #else:
        #   ax.plot(T,np.divide(frac2Full[:,i],1-frac2Full[:,i]),marker=("+"),color=coloring[int(i/filecycle)],alpha=0.5,linewidth=0.5)
    #
    #Linear fit of the data
    # coeff = analysis.LinearFit(np.array(frac2FullAverageVec[1]), T, (Nt-1)//3);
    # b = np.exp(coeff[1]);
    # ax2.plot(T,coeff[0]*T+b,color='red')
    #ax.set_title(r'$\Gamma/K$')
    ax.legend(loc='upper left')
    ax.set_yscale('log');
    #ax2.legend(loc='upper right')
    #ax.tick_params(axis='both',labelsize=24);
    #ax2.set_xticks([]);
    #plt.title('Elasticity ratio')
    ax.set_xlabel(r'$t$',fontsize='large')
    ax.set_ylabel(r'$\ln \frac{f}{1-f}$',fontsize='large')
    

    #%% Linear fit of the data
    Nindep = Nfile//filecycle;
    coeffs = [];
    coeff_th = [];
    timesteps = (Nt-1)//1;
    for i in range(Nindep):
        coeffs.append(analysis.LinearFit(np.array(frac2FullAverageVec[i]), T, timesteps));
        coeff_th.append(Params[i*filecycle][1][8]*Params[i*filecycle][1][9]*Params[i*filecycle][1][10]/0.01)

    #%% Plot the results in a table of plots
    coloring = ['grey','red','blue','g','c','m','k','brown'];
    labeling = ['Birth rate','Preferred area','Preferred perimeter','Perimeter elasticity','Area elasticity','Pressure strength','Self-advection','Rotational diffusion'];
    state = ['Liquid', 'Solid'];
    listyle = ['dotted','solid'];
    #labeling = ['Perimeter elasticity','Area elasticity','Pressure strength','Self-advection','Rotational diffusion'];
    filecycle=10;
    macrocycle=80;
    nline = 2;
    ncol = 4;
    fig,axs=plt.subplots(nline,ncol, sharex=True,sharey=True)
    tt=0;
    for i in range(Nfile):
        j = i%macrocycle;
        st = i//macrocycle;
        exp = j//filecycle;
        if j%filecycle==0:
            frac2FullAverage = np.mean(frac2Full[:,st*macrocycle+filecycle*exp:st*macrocycle+filecycle*exp+filecycle-1],axis=1)
            y = np.divide(frac2Full[:,i],1-frac2Full[:,i]);
            cond = y>10^(-3);
            axs[exp//ncol,exp%ncol].plot(100*T,np.divide(frac2FullAverage,1-frac2FullAverage),marker=("+"),color=coloring[exp],alpha=1,linewidth=2,label=state[st],linestyle=listyle[st]);
            #axs[exp//ncol,exp%ncol].plot(100*T,y,marker=("+"),color=coloring[exp],alpha=0.5,linewidth=0.5);
            #axs[exp//3,exp%3].set_xscale('log');
            axs[exp//ncol,exp%ncol].set_yscale('log');
            axs[exp//ncol,exp%ncol].set_ylim(0.01,10);
            axs[exp//ncol,exp%ncol].set_yticks([0.1,10,10]);
            axs[exp//ncol,exp%ncol].set_title(labeling[exp]);
            #axs[exp//ncol,exp%ncol].legend();
            axs[nline-1,exp%ncol].set_xlabel(r'$t$')
            axs[exp//ncol,0].set_ylabel(r'$\frac{f}{1-f}$')
            #axs[exp//ncol,exp%ncol].set_ylim(10**(-1.5),10**1.5)
            #axs[exp//ncol,exp%ncol].autoscale(enable=False)
            #plt.scatter(T,np.log(np.divide(frac2Full[:,i],1-frac2Full[:,i])),label=labeling[int(i/filecycle)],marker=("+"),color=coloring[int(i/filecycle)])
        #else:
        #    axs[exp//ncol,exp%ncol].plot(100*T,np.divide(frac2Full[:,i],1-frac2Full[:,i]),marker=("+"),color=coloring[exp],alpha=0.5,linewidth=0.5,linestyle=listyle[st]);
            #plt.scatter(T,np.log(np.divide(frac2Full[:,i],1-frac2Full[:,i])),marker=("+"),color=coloring[exp])
    #plt.plot(T,2*Params[0][1][8]*Params[0][1][6]*Params[0][1][9]*(0.13302446585454544)*100*T+np.log(0.1/0.9)*np.ones(Nt),color="red")
    #plt.legend(loc='upper left')
    #plt.title('Liquid resident')
    
    #%% Storing fixation diagrams with mechanical parameters
    fixated0 = [];
    fixated1 = [];
    nonfixated = [];
    x0=[]; y0=[];x1=[]; y1=[];xn=[]; yn=[];
    
    for j in range(Nfile):
        temp = analysis.HasPartiallyFixated(frac1Full[:,j],0.05);
        #temp = analysis.HasFixated(frac1Full[:,j]);
        if temp[0] == True:
            if temp[1] == 0:
                #Perimeter and area
                fixated0.append([Params[j][1][1]-Params[j][1][0], Params[j][1][4]]);
                
                #Area and birth rate
                #fixated0.append([Params[j][1][4], Params[j][1][10]]);
                
                #Perimeter and birth rate
                #fixated0.append([Params[j][1][1]-Params[j][1][0], Params[j][1][10]]);
                
                #fixated0.append([Params[j][1][3], Params[j][1][6]]);
            elif temp[1] == 1:
                #Perimeter and area
                fixated1.append([Params[j][1][1]-Params[j][1][0], Params[j][1][4]]);
                
                #Area and birth rate
                #fixated1.append([Params[j][1][4], Params[j][1][10]]);
                
                #Perimeter and birth rate
                #fixated1.append([Params[j][1][1]-Params[j][1][0], Params[j][1][10]]);
                
                #fixated1.append([Params[j][1][3], Params[j][1][6]]);
        else:
            #Area and birth rate
            #nonfixated.append([Params[j][1][4], Params[j][1][10]]);
            
            #Perimeter and area
            nonfixated.append([Params[j][1][1]-Params[j][1][0], Params[j][1][4]]);
            
            #Perimeter and birth rate
            #nonfixated.append([Params[j][1][1]-Params[j][1][0], Params[j][1][10]]);
            
            #nonfixated.append([Params[j][1][3], Params[j][1][6]]);
    #Plotting
    if len(fixated0) >0 : x0,y0 = zip(*fixated0);
    if len(fixated1) >0 :x1,y1 = zip(*fixated1);
    if len(nonfixated) >0 :xn, yn = zip(*nonfixated);
    #%% Plotting
    fig = plt.figure();
    ax = fig.add_subplot(1, 1, 1);
    if len(x0) > 0 : ax.scatter(x0,y0,color='red',label='Extinction of mutant')
    if len(x1) > 0 : ax.scatter(x1,y1,color='blue',label='Fixation of mutant')
    if len(xn) > 0 : ax.scatter(xn,yn,color='grey',label='Coexistence')
    #ax.set_xlabel(r'$\gamma^2/\gamma^1$', fontsize=16);
    #ax.set_ylabel(r'$K^2/K^1$', fontsize=16);
    
    #Perimeter and area
    ax.set_ylabel(r'$a_m^0$', fontsize=20);
    ax.set_xlabel(r'$p_m^0-p_r^0$', fontsize=20);
    #And ploting the theoretical line
    perims = np.linspace(-0.5,0.5,100);
    areas = np.linspace(0.2,2,100);
    ax.plot(perims,np.ones(len(perims)) -1.86*perims,color='black',label='Theory')
    
    #Area and Birth Rate
    #ax.set_xlabel(r'$a_m^0$', fontsize=20);
    #ax.set_ylabel(r'$\Pi_r^h$', fontsize=20);
    
    
    #Perimeter and Birth Rate
    #ax.set_xlabel(r'$p_m^0-p_r^0$', fontsize=20);
    #ax.set_ylabel(r'$\Pi_r^h$', fontsize=20);
    
    ax.legend(shadow=True,fontsize=12)
    
    #%% Time loop to get fraction of cells in a fixed circle
    
    #Define region of interest
    center = [10,10];
    radius = 10;
    
    IndOfInterest = [];
    clone = [];
    frac = [];
    for i in range(Nt):
        #Gets the indices of cells in circle
        position = reader.CellPos(i, Data)
        IndOfInterest.append(analysis.frame_circle(position, center, radius));
        #Compute the frequency
        celltype = reader.CellType(i, Data);
        tempclone, tempfrac = analysis.fraction(celltype[IndOfInterest[i]]);
        clone.append(tempclone);
        frac.append(tempfrac)
    
    #%% Plot the results
    frac1 = np.zeros(Nt);
    frac2 = np.zeros(Nt);
    for i in range(Nt):
        frac1[i] = frac[i][0];
        frac2[i] = frac[i][1];
    plt.scatter(T,frac1,color='blue',label='Solid-like')
    plt.scatter(T,frac2,color='red',label='Fluid-like')
    plt.legend()
    
    #%% Compute average area and perimeters
    AvArea = [];
    AvPerimeter = []; 
    for i in range(Nfile):
        cellarea = reader.CellArea(Nt-1, Data[i]);
        cellperimeter = reader.CellPeri(Nt-1, Data[i]);
        celltype = reader.CellType(Nt-1, Data[i]);
        AvArea.append(analysis.AverageValuePerType(cellarea,celltype));
        AvPerimeter.append(analysis.AverageValuePerType(cellperimeter,celltype));
        
    Area0 = []; 
    Area1 = [];
    fArea0 = []; 
    fArea1 = [];
    Peri0 = [];
    Peri1 = [];
    AvAreaDiff = np.zeros(Nfile);
    AvPerimDiff = np.zeros(Nfile);

    f = np.linspace(0, 1, Nfile);
    for i in range(Nfile):
        preferredarea0 = Params[i][1][4]*np.ones(Nfile);
        #preferredarea0 = np.ones(Nfile);
        preferredarea1 = Params[i][1][4]*np.ones(Nfile);
        preferredperim0 = Params[i][1][0]*np.ones(Nfile);
        preferredperim1 = Params[i][1][1]*np.ones(Nfile);
        cellarea = reader.CellArea(Nt-1, Data[i]);
        cellperimeter = reader.CellPeri(Nt-1, Data[i]);
        celltype = reader.CellType(Nt-1, Data[i]);
        AvAreaDiff[i] = AvArea[i][1,0] - AvArea[i][0,0];
        AvPerimDiff[i] = AvPerimeter[i][1,0] - AvPerimeter[i][0,0];
        Area0.append(analysis.ValuePerType(cellarea,celltype)[0]-preferredarea0[i]);
        Area1.append(analysis.ValuePerType(cellarea,celltype)[1]-preferredarea1[i]);
        fArea0.append(np.multiply(f[i],analysis.ValuePerType(cellarea,celltype)[0])-f[i]);
        fArea1.append(np.multiply((1-f[i]),analysis.ValuePerType(cellarea,celltype)[1])-1+f[i]);
        Peri0.append(analysis.ValuePerType(cellperimeter,celltype)[0]-preferredperim0[i]);
        Peri1.append(analysis.ValuePerType(cellperimeter,celltype)[1]-preferredperim1[i]);
        
    #%%Plot Delta p as function of Delta a
    #plt.scatter(AvArea0,AvPerimeter0)
    plt.scatter(AvAreaDiff,AvPerimDiff)
    
        
    #%% Plot
    fig = plt.figure();
    ax = fig.add_subplot(1, 1, 1);
    ax.boxplot(Area0,showfliers=False, notch=True,widths=0.2,patch_artist=True,boxprops=dict(facecolor='blue', color='blue'),whiskerprops=dict(color='blue'))
    ax.boxplot(Area1,showfliers=False, notch=True,widths=0.2,patch_artist=True,boxprops=dict(facecolor='red', color='red'),whiskerprops=dict(color='red'))
    ax.set_title(r'$a-a_0^i$',fontsize='x-large')
    fig2 = plt.figure();
    ax2 = fig2.add_subplot(1, 1, 1);
    ax2.boxplot(Peri0,showfliers=False, notch=True,widths=0.2,patch_artist=True,boxprops=dict(facecolor='blue', color='blue'),whiskerprops=dict(color='blue'))
    ax2.boxplot(Peri1,showfliers=False, notch=True,widths=0.2,patch_artist=True,boxprops=dict(facecolor='red', color='red'),whiskerprops=dict(color='red'))
    ax2.set_title(r'$p-p_0^i$',fontsize='x-large')
    
    #%% Plot
    fig = plt.figure();
    ax = fig.add_subplot(1, 1, 1);
    f = np.linspace(0, 1, 11);
    p01 = 3.63;
    p02 = 3.9;
    AvArea0 = np.zeros(Nfile); AvArea1 = np.zeros(Nfile);
    AvAreaError0 = np.zeros(Nfile); AvAreaError1 = np.zeros(Nfile);
    AvPeri0 = np.zeros(Nfile); AvPeri1 = np.zeros(Nfile);
    AvPeriError0 = np.zeros(Nfile); AvPeriError1 = np.zeros(Nfile);
    for i in range(Nfile):
        AvArea0[i] = AvArea[i][0,0];
        AvAreaError0[i] = AvArea[i][0,1];
        AvArea1[i] = AvArea[i][1,0];
        AvAreaError1[i] = AvArea[i][1,1];
        AvPeri0[i] = AvPerimeter[i][0,0];
        AvPeriError0[i] = AvPerimeter[i][0,1];
        AvPeri1[i] = AvPerimeter[i][1,0];
        AvPeriError1[i] = AvPerimeter[i][1,1];
    ax.errorbar(f,AvArea0,yerr=AvAreaError0,fmt='o',color='blue',label='Type 1, p0=%.2f' %p01)
    ax.errorbar(f,AvArea1,yerr=AvAreaError1,fmt='o',color='red',label='Type 2, p0=%.2f' %p02)
    ax.plot(f,np.ones(len(f)),color='black',linestyle="--",linewidth=1)
    ## Comparison with theory
    # Check conservation of area
    conser = np.multiply(f,AvArea1) + np.multiply(1-f,AvArea0);
    #ax.plot(f,conser,color='green',linestyle="-",linewidth=2)
    # Check perimeter/area relationship
    #fracf = np.divide(f,1-f);
    #periTerm = np.multiply(AvPeri0,AvPeri0-p01)+np.multiply(AvPeri1,AvPeri1-p02)
    #areaTh = 1+np.power(np.multiply(-fracf,0.5*periTerm),0.5);
    #ax.plot(f,areaTh, color='black')
    ax.set_xlabel('Type 2 fraction');
    ax.set_ylabel('Average Area per last frame')
    ax.legend()
    fig.show()
    fig2 = plt.figure();
    ax2 = fig2.add_subplot(1, 1, 1);
    ax2.errorbar(f,AvPeri0,yerr=AvPeriError0,fmt='o',color='blue',label='Type 1, p0=%.2f' %p01)
    ax2.errorbar(f,AvPeri1,yerr=AvPeriError1,fmt='o',color='red',label='Type 2, p0=%.2f' %p02)
    ax2.plot(f,preferredperim0,color='blue',linestyle="--",linewidth=1);
    ax2.plot(f,preferredperim1,color='red',linestyle="--",linewidth=1);
    ax2.set_xlabel('Type 2 fraction');
    ax2.set_ylabel('Average Perimeter per last frame')
    ax2.legend()
    fig2.show()

