import h5py
import matplotlib.pyplot as plt
import numpy as np
import os
import re
import numpy as np
from vtk import vtkStructuredPointsReader
from vtk.util import numpy_support as VN

file_read = 'VTK_resampled_fine'# name of folder with data in same directory
write_name = 'Reversed_Fine' # name of the folder that will be written to
# frozen = True # flag specifying if wanting to only run first timestep in SHARPy
write_h5 = False
set = {}

# Define the shape of the grid
grid_shape = (750, 350, 75)
set['grid_shape'] = grid_shape

# Lengths of grid (m)
x_dis = 750
y_dis = 350 
z_dis = 75 

# in code dimensions are set 
dx =  x_dis / (grid_shape[0] - 1)
dy =  y_dis / (grid_shape[1] - 1)
dz =  z_dis / (grid_shape[2] - 1)

set['d_L'] = [dx, dy, dz]

set['time'] = [0.0, 1.0] # start time and time step

def proof():
    grid = (5, 5, 3)
    test_vals = [1, 2, 3, 4, 5] * (grid[1] * grid[2])
    test_final = np.zeros(grid) 

    count = 0
    for z in range(grid[2]):
        for y in range(grid[1]):
            for x in range(grid[0]):
                test_final[x, y, z] = test_vals[count]
                count += 1 
    with h5py.File('just a test.h5', 'w') as f:
        f.create_dataset('data', data=test_final, dtype= np.float64)

# origin of grid
# set['origin'] = [-x_dis, -y_dis/2, -z_dis/2]
set['origin'] = [-x_dis, -y_dis/2, 0]
test = np.array([[[1, 2], [3, 4], [5, 6]],
                [[7, 8], [9, 10], [11, 12]]])
test_nu = np.flip(test, axis=1)

# loading all .mat file names to access 
if write_h5:
    # Specify the directory to search in
    directory = os.path.join(os.path.dirname(os.path.abspath(__file__)), file_read)
    # Set the desired index to start searching from
    start_index = '41'
    stop_index = '56'
    # Get a list of all files in the directory
    files = os.listdir(directory)
    # Find files with the extension .mat and starting from the specified index
    pattern = re.compile(r'^bp_8ms_fine_1mres_([0-9]+)\.vtk$')
    vtk_files = [f for f in files if pattern.match(f) and int(start_index) <= int(pattern.match(f).group(1)) <= int(stop_index)]

    i_index = len(vtk_files) 
    
    for i, val in enumerate(vtk_files):
        reader = vtkStructuredPointsReader()
        reader.SetFileName(os.path.join(directory, val))
        reader.ReadAllVectorsOn()
        reader.ReadAllScalarsOn()
        reader.Update()

        data = reader.GetOutput()

        dim = data.GetDimensions()
        vec = list(dim)
        vec = [i-1 for i in dim]
        vec.append(3)

        vel = VN.vtk_to_numpy(data.GetPointData().GetArray('U'))
        u   = vel[:,0] 
        v   = vel[:,1]
        w   = vel[:,2]

        x = np.zeros(data.GetNumberOfPoints())
        y = np.zeros(data.GetNumberOfPoints())
        z = np.zeros(data.GetNumberOfPoints())

        for ii in range(data.GetNumberOfPoints()):
                x[ii],y[ii],z[ii] = data.GetPoint(ii)
            
        nx = x[0:vec[0]+1].shape[0]
        ny = y[0:(vec[0]+1)*(vec[1]+1):vec[0]+1].shape[0]
        nz = z[0:(vec[0]+1)*(vec[1]+1)*(vec[2]+1):(vec[0]+1)*(vec[1]+1)].shape[0]

        ux_data = np.zeros([nx, ny, nz]) 
        uy_data = np.zeros([nx, ny, nz]) 
        uz_data = np.zeros([nx, ny, nz])

        count = 0
        for z in range(nz):
            for y in range(ny):
                for x in range(nx):
                    # ux_data[x, y, z] = data['x'][0][count] assuming this is zero for now
                    uy_data[x, y, z] = v[count]
                    uz_data[x, y, z] = w[count]
                    count += 1
        # Flip the second axis (index 1)
        uy_data = np.flip(uy_data, axis=0)
        uz_data = np.flip(uz_data, axis=0)
        # u_3d = np.reshape(u,(nx,ny,nz),order='C')
        # v_3d = np.reshape(v,(nx,ny,nz),order='C')
        # w_3d = np.reshape(w,(nx,ny,nz),order='C')

        h5_dirc = os.path.join(os.path.dirname(os.path.abspath(__file__)), write_name + '_h5_write')
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

def xdmf_create(filename, set, time_step):
    # Writing XDMF file for passing infomration to SHARPy
    file = open('{}.xdmf'.format('RUN_' + filename), 'w')
    file.write('<?xml version="1.0" ?>\n')
    file.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n')
    file.write('<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="3.0">\n')
    file.write('<Domain>\n')
    # time information
    file.write('<Time TimeType="HyperSlab">\n')
    file.write('\t<DataItem Format="XML" NumberType="Float" Dimensions="3">\n')
    file.write('\t{} {} 0.0\n'.format(set['time'][0], set['time'][1]))
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

i_index = 15
xdmf_create(write_name, set, i_index)