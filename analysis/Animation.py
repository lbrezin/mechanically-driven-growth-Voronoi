#!/usr/bin/env python
# coding: utf-8

import reader
import readFile
import matplotlib.animation as an
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import filedialog as fd
import numpy as np
import analysis
import copy
import pickle
import time
#plt.style.use('seaborn-poster')
plt.rcParams.update({'font.size': 40})
plt.rcParams["figure.autolayout"] = True
plt.rcParams['text.usetex'] = True
#plt.style.use('seaborn-poster')
plt.rcParams['figure.figsize'] = [10, 10]


#%% Reading files
root = tk.Tk()
root.withdraw()
filenames = fd.askopenfilenames();
Nfile = len(filenames);
Data=[];
Params=[];
start = time.time();
for i in range(Nfile):
    Data.append(reader.ReadFile(filenames[i])[0]);
    Params.append(reader.ReadFile(filenames[i])[1]);
end = time.time();
print(end-start);

#%% Comparing performances

root = tk.Tk()
root.withdraw()
filenames = fd.askopenfilenames();
Nfile = len(filenames);
Data=[];
Params=[];
start = time.time();
for i in range(Nfile):
    Data.append(readFile.ReadFileC(filenames[i])[0]);
    Params.append(readFile.ReadFileC(filenames[i])[1]);
end = time.time();
print(end-start);

#%% Pickle the data
fileData = open('./Extracted_Data/Data_mixture_frequency_20230406','w+b');
fileParams = open('./Extracted_Data/Params_mixture_frequency_20230406','w+b');
pickle.dump(Data,fileData);
pickle.dump(Params,fileParams);

#%% Load pickled data
fileDataOther = open('./Extracted_Data/Data_mixture_frequency_20230406','rb');
fileParamsOther = open('./Extracted_Data/Params_mixture_frequency_20230406','rb');
Data = pickle.load(fileDataOther);
Params = pickle.load(fileParamsOther);

#%% Selecting with file we are working with
filenum = 0;

#%% Animation by first plotting all the frames and sorting them
figs = [];
Nframe = len(Data[filenum]);
fig = plt.figure();
for i in range(2,int(Nframe/2)):
    tempax = fig.add_subplot(111)
    reader.MakeColorPlot(i, Data[filenum], tempax);
    figs.append([tempax])

#%% Create animation from the artists object and saves it
#figani = plt.figure()
ani = an.ArtistAnimation(fig=fig, artists=figs, interval=30)
writermp4 = an.FFMpegWriter(fps=20);
ani.save("tests.mp4",writer=writermp4)
#plt.show()

#%% Plotting animation of file filenum
step = 1; #Only show every step frame
fig2 = plt.figure();
anim = an.FuncAnimation(fig2, reader.MakeColorPlot, frames=len(Data[filenum][::step]), fargs=(Data[filenum][::step],fig2), interval=5)

#%% Saving animation
writergif = an.PillowWriter(fps=5);
writermp4 = an.FFMpegWriter(fps=10);
anim.save("invasion_solid_area09_example.mp4",writer=writermp4)

# %% Plots of individual frames 
#frameInd = len(Data[filenum])-1;
filenum=47
frameInd = int(8*len(Data[filenum])/8-1)
fig3,ax3 = plt.subplots(1,1,sharex=True);
#reader.MakeColorPlot(1,Data[filenum],ax3[0])
#fig4,ax4 = plt.subplots();
reader.MakeColorPlot(frameInd,Data[filenum],ax3)

#%% Plot with some cells colored
tocolor = [];
timeframe = len(Data[filenum])-1;
fig6 = plt.figure();
reader.ColorOneCell(timeframe,Data[filenum],fig6,tocolor)
fig6.show()
# tocolor = [66];
# timeframe = 905;
# fig7 = plt.figure();
# reader.ColorOneCell(timeframe,Data[filenum],fig7,tocolor)
# fig7.show()

#%% Plot of pressure, in a zoom
fig5 = plt.figure();
reader.MakePressurePlot(len(Data[filenum])-1,Data[filenum],fig5,1,0,1,Params[filenum][1][4]);
fig5.show()

#%% Plot number of cells over time
T, Ncells = reader.CellNumOTime(Data[0]);
plt.plot(T,Ncells);

#%% Pickle dump the time
fileTMSD = open('./TimeMSD_20240902','w+b');
pickle.dump(T,fileTMSD);


#%% Plot trajectory
NumTraj = 20;
start=10;
center = [12, 12];
radius = 6;
fractraj=4; #Fraction of total time we are looking at
T,Ncells =  reader.CellNumOTime(Data[filenum]);
indices = np.random.randint(0,Ncells[0],NumTraj);
indices_circle = analysis.frame_circle(reader.CellPos(len(Data)-1, Data[filenum]), center,radius)
Positions = reader.CellPosInd(Data[filenum],indices_circle,start);
#Get the types of cells in the circle
celltypes = reader.CellType(len(Data)-1, Data[filenum]);
celltypes_select = celltypes[indices_circle];
#Plot with color by type, switch=0, or color by cell, switch != 0
sw=1; 
fig = plt.figure();
ax = fig.add_subplot(1, 1, 1);
for j in range(len(indices_circle)):
    if sw==0:
        if celltypes_select[j]==0:
            ax.plot(Positions[0:int(len(Positions)/fractraj),j,0],Positions[0:int(len(Positions)/fractraj),j,1],color='blue');
        elif celltypes_select[j]==1:
            ax.plot(Positions[0:int(len(Positions)/fractraj),j,0],Positions[0:int(len(Positions)/fractraj),j,1],color='red');
    else:
        ax.plot(Positions[0:int(len(Positions)/fractraj),j,0],Positions[0:int(len(Positions)/fractraj),j,1]);
    
#%% Compute global diffusion constant as function of fraction
Diff = np.zeros(len(Data));
Diff0 = np.zeros(len(Data));
Diff1 = np.zeros(len(Data));
D0 = 0.1**2/2; # Diffusion of isolated individual cell
MSDf =[];
MSD0f = [];
MSD1f = [];
VariedParam = [];
for i in range(len(Data)):
    indices_circle = analysis.frame_circle(reader.CellPos(len(Data[i])-1, Data[i]), center,radius)
    Positions = reader.CellPosInd(Data[i],indices_circle,start);
    celltypes = reader.CellType(len(Data[i])-1, Data[i]);
    celltypes_select = celltypes[indices_circle];
    Positions0=[];
    Positions1=[];
    for l in range(len(Positions[0])):
        if celltypes_select[l]==0:
            Positions0.append(Positions[:,l,:]);
        elif celltypes_select[l]==1:
            Positions1.append(Positions[:,l,:]);
    Positions_ref= Positions.copy();
    Positions0_ref= copy.deepcopy(Positions0);
    Positions1_ref= copy.deepcopy(Positions1);
    MSD = np.zeros((int(len(Positions)/1),len(Positions[0])));
    for j in range(int(len(Positions)/1)):
        Positions_ref[j,:,:] = Positions_ref[j,:,:]-Positions[0,:,:];
        MSD[j,:] = np.power(Positions_ref[j,:,0],2) + np.power(Positions_ref[j,:,1],2);
    # Now compute the same thing but per type
    if len(Positions0)>0:
        MSD0 = np.zeros((len(Positions0[0]),len(Positions0)));
        for k in range(len(Positions0)):
            for j in range(len(Positions0[0])):
                Positions0_ref[k][j,:] = Positions0_ref[k][j,:]-Positions0[k][0,:];
                MSD0[j,k] = np.power(Positions0_ref[k][j,0],2) + np.power(Positions0_ref[k][j,1],2);
    if len(Positions1)>0:
        MSD1 = np.zeros((len(Positions1[0]),len(Positions1)));
        for k in range(len(Positions1)):
            for j in range(len(Positions1[0])):
                Positions1_ref[k][j,:] = Positions1_ref[k][j,:]-Positions1[k][0,:];
                MSD1[j,k] = np.power(Positions1_ref[k][j,0],2) + np.power(Positions1_ref[k][j,1],2);
    MSD_Av = np.mean(MSD,axis=1);
    if len(Positions0)>0:
        MSD0_Av = np.mean(MSD0,axis=1);
        D=0;
        for j in range(1,len(MSD0_Av)):
            D = D + MSD0_Av[j]/(j*100);
        Diff0[i] = 1/4*D/len(MSD0_Av)/D0;
    elif len(Positions0)==0:
        Diff0[i] = np.nan;
    if len(Positions1)>0:
        MSD1_Av = np.mean(MSD1,axis=1);
        D=0;
        for j in range(1,len(MSD1_Av)):
            D = D + MSD1_Av[j]/(j*100);
        Diff1[i] = 1/4*D/len(MSD1_Av)/D0;
    elif len(Positions1)==0:
        Diff1[i] = np.nan;
    D=0;
    for j in range(1,len(MSD_Av)):
        D = D + MSD_Av[j]/(j*100);
    Diff[i] = 1/4*D/len(MSD_Av)/D0;
    MSDf.append(MSD_Av);
    VariedParam.append(Params[i][1][1]);
    #MSD0f.append(MSD0_Av);
    #MSD1f.append(MSD_Av);

#%% Plotting as function of p0
plt.scatter(VariedParam,Diff);#plt.scatter(f, Diff0, color='blue',label='Solid')
plt.xlabel(r'$p^0$', fontsize='x-large');
plt.ylabel(r'D (a.u)', fontsize='x-large');
#%% Ploting
f = np.linspace(0, 1,50);
plt.scatter(f,Diff,color='black',label='Global');#plt.scatter(f, Diff0, color='blue',label='Solid')
#plt.scatter(f, Diff1, color='red',label='Liquid')
#plt.scatter(f, Diff0, color='blue',label='Solid')
plt.xlabel(r'f', fontsize='x-large');
plt.ylabel(r'D (a.u)', fontsize='x-large');
plt.legend()
#%%MSD
fig = plt.figure();
ax = fig.add_subplot(1, 1, 1);
for j in range(0,len(MSDf),2):
    fr = 0.05*j;
    ax.plot(np.log(100*T[start+1::]),np.log(MSDf[j][1::]),label='%.1f'%fr)
#ax.set_xscale('log')
#ax.set_yscale('log')
#ax.plot(100*T[start+1:start+20],MSDf[len(MSDf)-1][1]+100*T[start+1:start+20]-100*T[start+1],linestyle='dashed')
ax.plot(np.linspace(np.log(100*T[start+1]),np.log(100*T[start+1])+1,50),2*(np.linspace(np.log(100*T[start+1]),np.log(100*T[start+1])+1,50)-np.log(100*T[start+1]))+np.log(MSDf[len(MSDf)-1][1])+1.5,linestyle='dashed')
ax.plot(np.linspace(np.log(100*T[-1])-2,np.log(100*T[-1]),50),(np.linspace(np.log(100*T[-1])-2,np.log(100*T[-1]),50)-np.log(100*T[-1]))+np.log(MSDf[len(MSDf)-1][-1])+1,linestyle='dashed')
plt.legend()
#plt.plot(MSD0f[20])


#%% Computes trajectories, MSD and such when there is divisions, with tracker cells that do not divide.
MSDAllFiles = [];
for k in range(Nfile):
    filenum=k;
    Nt = len(Data[filenum]);
    #Initialize size of the lists
    celltypes = reader.CellType(0, Data[filenum]);
    ind0 = np.where(celltypes == 0)[0];
    start = Nt//2;
    Positions = np.zeros((Nt-start,len(ind0),2));
    for i in range(start,Nt):
        celltypes = reader.CellType(i, Data[filenum]);
        ind0 = np.where(celltypes == 0)[0];
        Positions[i-start,:,:] = reader.CellPosTimeInd(Data[filenum], ind0, i);
    
    #% Unwrapping the coordinates with PBC
    Displacements = np.zeros((Nt-start-1,len(ind0),2));
    TempPosArrival = np.zeros(2);
    UnwrapTraj = np.zeros((Nt-start,len(ind0),2));
    #Get size of the box
    Lx = Data[filenum][0][0][2]-Data[filenum][0][0][3];
    Ly = Data[filenum][0][0][5]-Data[filenum][0][0][4];
    for j in range(len(ind0)):
        #Initial position of the cell
        UnwrapTraj[0,j,:] = Positions[0,j,:];
        for i in range(1,Nt-start):
            #First add the box length if needed in x and y
            if(abs(Positions[i,j,0] - Positions[i-1,j,0])>Lx/2):
                if(Positions[i-1,j,0]>Lx/2):
                    TempPosArrival[0] = Positions[i,j,0] +Lx
                else:
                    TempPosArrival[0] = Positions[i,j,0] - Lx
            else:
                TempPosArrival[0] = Positions[i,j,0]
            if(abs(Positions[i,j,1] - Positions[i-1,j,1])>Ly/2):
                if(Positions[i-1,j,1]>Ly/2):
                    TempPosArrival[1] = Positions[i,j,1] +Ly
                else:
                    TempPosArrival[1] = Positions[i,j,1] - Ly
            else:
                TempPosArrival[1] = Positions[i,j,1]
            #Then compute the displacement vectors 
            Displacements[i-1,j,:] = TempPosArrival[:]- Positions[i-1,j,:];
            UnwrapTraj[i,j,:] = UnwrapTraj[i-1,j,:] + Displacements[i-1,j,:];
        
            
    
    
    #% Computing MSD for all cell over time
    Positions_ref = copy.deepcopy(UnwrapTraj);
    MSD = np.zeros((UnwrapTraj.shape[0],UnwrapTraj.shape[1]));
    for i in range(start,Nt):
        Positions_ref[i-start,:,:] = Positions_ref[i-start,:,:]-UnwrapTraj[0,:,:];
        MSD[i-start,:] = np.power(Positions_ref[i-start,:,0],2) + np.power(Positions_ref[i-start,:,1],2);
    
    
    #%Average MSD and append to the list
    MSDAv = np.mean(MSD,axis=1);
    MSDAllFiles.append(MSDAv);
    
#%% Save this data that took a while to do
fileMSD = open('./MSD_20240902','w+b');
pickle.dump(MSDAllFiles,fileMSD);
fileMSD.close()
#%%Plotting all the MSD from all files
fig3 = plt.figure();
ax3 = fig3.add_subplot(1, 1, 1);
alpha = 1;
legends= ['Liquid w/ division', 'Liquid w/o division', 'Solid w/ division', 'Solid w/o division']
colors = ['blue', 'blue', 'black', 'black']
styles = ['dashed','solid','dashed','solid']
ax3.set_xscale('log')
ax3.set_yscale('log')
#ax2.plot(np.log(T[start:Nt]-T[start]),np.log(MSD[:,idc]));
for k in range(Nfile):
    ax3.plot(T[start:Nt]-T[start],MSDAllFiles[k],label=legends[k],linestyle=styles[k],color=colors[k]);
ax3.set_xlabel('t (a.u.)',fontsize='large',fontweight='light')
ax3.set_ylabel('MSD (a.u.)',fontsize='large')
ax3.tick_params(labelsize=30)
ax3.legend()
#plt.savefig('MSDAll.pdf',format='pdf',bbox_inches='tight')

#%% Plotting MSD for a specific cell
fig2 = plt.figure();
ax2 = fig2.add_subplot(1, 1, 1);
alpha = 1;
idc=4;
ax2.set_xscale('log')
ax2.set_yscale('log')
#ax2.plot(np.log(T[start:Nt]-T[start]),np.log(MSD[:,idc]));
D = np.polyfit(T[start:Nt]-T[start], MSD[:,idc], 1)[0];
ax2.plot(T[start:Nt]-T[start],MSD[:,idc]);
#ax2.plot(T[start:Nt]-T[start],D*(T[start:Nt]-T[start]));
#ax2.plot(np.log(T[start:Nt]-T[start]),alpha*np.log(T[start:Nt]-T[start])+(np.log(MSD[1,idc]))*np.ones(Nt-start));


#%% Plotting the trajectories
fig = plt.figure();
ax = fig.add_subplot(1, 1, 1);
idc = 7
for j in range(len(ind0)):
    ax.plot(UnwrapTraj[:,j,0],UnwrapTraj[:,j,1])
#ax.plot(Positions[:,idc,0],Positions[:,idc,1])
#%% Extracting the diffusion coefficient and power law
coeff = np.polyfit(np.log(T[start+1:Nt]-T[start]), np.log(MSDAv[1::]), 1)
plt.plot(np.log(T[start+1:Nt]-T[start]),np.log(MSDAv[1::]),label='MSD');
plt.xlabel(r'$\log(t)$')
plt.ylabel(r'$\log(MSD)$')
plt.plot(np.log(T[start+1:Nt]-T[start]),coeff[0]*np.log(T[start+1:Nt]-T[start+1])+(coeff[1])*np.ones(Nt-start-1),label=r'$\alpha =$ %.2f' %coeff[0]);
plt.legend()
