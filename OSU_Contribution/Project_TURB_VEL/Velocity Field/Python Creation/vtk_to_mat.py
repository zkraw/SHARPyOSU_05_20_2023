#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 13:05:40 2022

@author: Rohit K S S Vuppala
         Graduate Student, 
         Mechanical and Aerospace Engineering,
         Oklahoma State University.

@email: rvuppal@okstate.edu

"""

import numpy as np
import scipy.io
from vtk import vtkStructuredPointsReader
from vtk.util import numpy_support as VN
import os 
from tqdm import tqdm

reader = vtkStructuredPointsReader()
reader.SetFileName("../VTK_resampled/boone-pickens_120.vtk")
reader.ReadAllVectorsOn()
reader.ReadAllScalarsOn()
reader.Update()

data = reader.GetOutput()


dim = data.GetDimensions()
vec = list(dim)
vec = [i-1 for i in dim]
vec.append(3)
#%%
vel = VN.vtk_to_numpy(data.GetPointData().GetArray('U'))
u   = vel[:,0] 
v   = vel[:,1]
w   = vel[:,2]

x = np.zeros(data.GetNumberOfPoints())
y = np.zeros(data.GetNumberOfPoints())
z = np.zeros(data.GetNumberOfPoints())

for i in range(data.GetNumberOfPoints()):
        x[i],y[i],z[i] = data.GetPoint(i)
    
nx = x[0:vec[0]+1].shape[0]
ny = y[0:(vec[0]+1)*(vec[1]+1):vec[0]+1].shape[0]
nz = z[0:(vec[0]+1)*(vec[1]+1)*(vec[2]+1):(vec[0]+1)*(vec[1]+1)].shape[0]


x_3d = np.reshape(x,(nx,ny,nz),order='F')
y_3d = np.reshape(y,(nx,ny,nz),order='F')
z_3d = np.reshape(z,(nx,ny,nz),order='F')

#%%
nstart = 0
nend   = 401

for n in tqdm(range(nstart,nend)):
    #
    path  = "../VTK_resampled/"
    fname = path+"boone-pickens_{0:03d}.vtk".format(n)
    reader = vtkStructuredPointsReader()
    reader.SetFileName(fname)
    reader.ReadAllVectorsOn()
    reader.ReadAllScalarsOn()
    reader.Update()

    data = reader.GetOutput()
    
    vel = VN.vtk_to_numpy(data.GetPointData().GetArray('U'))
    u   = vel[:,0] 
    v   = vel[:,1]
    w   = vel[:,2]

    u_3d = np.reshape(u,(nx,ny,nz),order='F')
    v_3d = np.reshape(v,(nx,ny,nz),order='F')
    w_3d = np.reshape(w,(nx,ny,nz),order='F')



    scipy.io.savemat('wind_{0:03d}.mat'.format(n),mdict={'u':u,'v':v,'w':w,'x':x,'y':y,'z':z})











