#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 08:29:11 2023

@author: louisbrezin
"""

import numpy as np
import cython

def ReadFileC(str file):
    ## File must be a string of where the file to read is
    f = open(file, "r");
    cdef int i=0;
    cdef list params =[];
    cdef list data=[];
    cdef list Coord=[];
    while (True):
        line = f.readline();
        #Terminat loop at end of file
        if not line:
            break
        line = line.split();
        #Check if we enconter the line that separates time points
        if line[0]=='N,':
            #Skip the line with N, time, Box and extract the data
            data.append([]);
    i=0;
    f = open(file, "r");
    line = f.readline();
    line = line.split();
    if line[0] == "p0_0":
        params.append(line);
        line = f.readline();
        params.append([float(x) for x in line.split()]);
        line=f.readline();
        line = line.split();
    while (True):
        #Check if we enconter the line that separates time points
        if line[0]=='N,':
            #Skip the line with N, time, Box and extract the data
            #data[i].append([]);
            line = f.readline();
            data[i].append([float(x) for x in line.split()]);
        while line[0]!='N,':
            line = f.readline();
            line = line.split();
            #print(line);
            if not line:
                break
            if line[0]!='N,':
                Coord.append([float(x) for x in line]);
                #Pdb().set_trace();
        data[i].append(Coord);
        Coord=[];
        i = i+1;
        #Terminate loop at end of file
        if not line:
            break
    f.close()
    return data, params