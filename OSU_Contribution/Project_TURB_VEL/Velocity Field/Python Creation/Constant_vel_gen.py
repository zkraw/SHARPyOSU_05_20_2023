import h5py
import matplotlib.pyplot as plt
import numpy as np
import os
import fnmatch
import scipy
import re

file_write = 'const_8_x_longer' # name of folder with data in same directory
set = {}

# Define the shape of the grid
grid_shape = (400, 40, 120)
set['grid_shape'] = grid_shape

# Lengths of grid (m)
x_dis = 400 
y_dis = 40 
z_dis = 120 

# in code dimensions are set 
dx =  x_dis / (grid_shape[0] - 1)
dy =  y_dis / (grid_shape[1] - 1)
dz =  z_dis / (grid_shape[2] - 1)

set['d_L'] = [dx, dy, dz]

# origin of grid
set['origin'] = [-x_dis-10, -y_dis/2, -z_dis/2]

# making a total of 3 files that specify the wind velocity filed at time = 0, 1, 2
time_step = 1 
for i in range(time_step):

    ux_data = np.ones(grid_shape)*8
    uy_data = np.zeros(grid_shape) 
    uz_data = np.zeros(grid_shape)

    h5_dirc = os.path.join(os.path.dirname(os.path.abspath(__file__)), file_write + '_h5_write')
    if not os.path.exists(h5_dirc):
        os.makedirs(h5_dirc)
    # Write ux.h5
    with h5py.File(os.path.join(h5_dirc, 'ux_{:04d}.h5'.format(i)), 'w') as f:
        f.create_dataset('data', data=ux_data, dtype= np.float64)

    # Write uy.h5
    with h5py.File(os.path.join(h5_dirc, 'uy_{:04d}.h5'.format(i)), 'w') as f:
        f.create_dataset('data', data=uy_data, dtype= np.float64)

    # Write uz.h5
    with h5py.File(os.path.join(h5_dirc, 'uz_{:04d}.h5'.format(i)), 'w') as f:
        f.create_dataset('data', data=uz_data, dtype= np.float64)

# Writing XDMF file for passing infomration to SHARPy
file = open('{}.xdmf'.format('RUN_' + file_write), 'w')
file.write('<?xml version="1.0" ?>\n')
file.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n')
file.write('<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="3.0">\n')
file.write('<Domain>\n')
# time information
file.write('<Time TimeType="HyperSlab">\n')
file.write('\t<DataItem Format="XML" NumberType="Float" Dimensions="3">\n')
file.write('\t0.0 1.0 0.0\n')
file.write('\t</DataItem>\n')
file.write('</Time>\n')
# Topology
key = set['grid_shape']
file.write('<Topology TopologyType="3DCORECTMesh" Dimensions="{} {} {}"/>\n'.format(key[0], key[1], key[2]))
# Geometry 
file.write('<Geometry GeometryType="ORIGIN_DXDYDZ">\n')
file.write('\t<DataItem Name="Origin" Dimensions="3" NumberType="Float" Precision="4" Format="XML">\n')
key = set['origin']
file.write('\t{} {} {}\n'.format(key[0], key[1], key[2]))
file.write('\t</DataItem>\n')
file.write('\t<DataItem Name="Spacing" Dimensions="3" NumberType="Float" Precision="4" Format="XML">\n')
key = set['d_L']
file.write('\t{:.4f} {:.4f} {:.4f}\n'.format(key[0], key[1], key[2]))
file.write('\t</DataItem>\n')
file.write('</Geometry>\n')
file.write('<Grid CollectionType="Temporal" GridType="Collection" Name="Collection">\n')

for i in range(time_step):
    file.write('\t<Grid Name="Grid{}">\n'.format(i))
    file.write('\t<Attribute Center="Node" ElementCell="" ElementDegree="0" ElementFamily="" ItemType="" Name="ux" Type="Scalar">\n')
    key = set['grid_shape']
    file.write('\t<DataItem DataType="Float" Dimensions="{} {} {} 3" Format="HDF" Precision="8">ux_{:04d}.h5</DataItem>\n'.format(key[0], key[1], key[2], i))
    file.write('\t</Attribute>\n')

    file.write('\t<Attribute Center="Node" ElementCell="" ElementDegree="0" ElementFamily="" ItemType="" Name="uy" Type="Scalar">\n')
    key = set['grid_shape']
    file.write('\t<DataItem DataType="Float" Dimensions="{} {} {} 3" Format="HDF" Precision="8">uy_{:04d}.h5</DataItem>\n'.format(key[0], key[1], key[2], i))
    file.write('\t</Attribute>\n')

    file.write('\t<Attribute Center="Node" ElementCell="" ElementDegree="0" ElementFamily="" ItemType="" Name="uz" Type="Scalar">\n')
    key = set['grid_shape']
    file.write('\t<DataItem DataType="Float" Dimensions="{} {} {} 3" Format="HDF" Precision="8">uz_{:04d}.h5</DataItem>\n'.format(key[0], key[1], key[2], i))
    file.write('\t</Attribute>\n')
    file.write('\t</Grid>\n')

file.write('</Grid>\n')
file.write('</Domain>\n')
file.write('</Xdmf>\n')