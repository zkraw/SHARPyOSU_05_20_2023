import h5py
import matplotlib.pyplot as plt
import numpy as np

name = '1_cos'
set = {}

u_inf = 10
gust_intensity = 0.2
gust_length = 1 * u_inf

# Define the shape of the grid
# all variables are defined in xyz order it is then rearranged to 
# match convention in SHARPy
grid_shape = (40, 80, 80)
set['grid_shape']  = (grid_shape[2], grid_shape[1], grid_shape[0])

# Lengths of grid (m)
x_dis = 10
y_dis = 40
z_dis = 40

# in code dimensions are set 
dx =  x_dis / (grid_shape[0] - 1)
dy =  y_dis / (grid_shape[1] - 1)
dz =  z_dis / (grid_shape[2] - 1)
set['d_L'] = [dz, dy, dx]

# x_values to calculate z velocities
x = np.arange(0, -x_dis-dx, -dx).tolist()
x.reverse()
# to get x values in proper order as origin values will be furthest away

# origin of grid
set['origin'] = [-x_dis-5, -y_dis/2, -z_dis/2]

# gust offset is capture by grid offset in origin
def gust_vel (x, u_inf, gust_intensity, gust_length): 
    if x > 0 or x < -gust_length:
        return 0 
    vel = (1.0 - np.cos(2.0 * np.pi * (x) / gust_length)) * gust_intensity*(u_inf) * 0.5
    return vel

def one_through_ten (x): 
    return 1*-x

# Create the ux data (constant 15)
ux_data = np.zeros((40, 80, 80)) # actually x velocity 
uy_data = np.zeros((40, 80, 80)) # actually z velocity
uz_data = np.zeros((40, 80, 80)) # actually y velocity 

vel = []
for i, val in enumerate(x):
    vel.append(gust_vel(val, u_inf, gust_intensity, gust_length))
    uz_data[i, :, :] = vel[-1]
# count = 19
# for i in range(len(uy_data), 20, -1):
#     uy_data[:, :, i] = vel[count]
#     count += 1

# make first value the initial intensity
# uy_data[0][:][:] = gust_intensity

# Write ux.h5
with h5py.File('ux_{}.h5'.format(name), 'w') as f:
    f.create_dataset('data', data=ux_data, dtype= np.float64)

# Write uy.h5
with h5py.File('uy_{}.h5'.format(name), 'w') as f:
    f.create_dataset('data', data=uy_data, dtype= np.float64)

# Write uz.h5
with h5py.File('uz_{}.h5'.format(name), 'w') as f:
    f.create_dataset('data', data=uz_data, dtype= np.float64)

def xdmf_create(filename, set):
   file = open('{}.xdmf'.format(filename), 'w')
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
   file.write('\t {} {} {}\n'.format(key[0], key[1], key[2]))
   file.write('\t</DataItem>\n')
   file.write('\t<DataItem Name="Spacing" Dimensions="3" NumberType="Float" Precision="4" Format="XML">\n')
   key = set['d_L']
   file.write('\t {:.4f} {:.4f} {:.4f}\n'.format(key[0], key[1], key[2]))
   file.write('\t</DataItem>\n')
   file.write('</Geometry>\n')
   file.write('<Grid CollectionType="Temporal" GridType="Collection" Name="Collection">\n')
   file.write('\t<Grid Name="Grid">\n')

   file.write('\t<Attribute Center="Node" ElementCell="" ElementDegree="0" ElementFamily="" ItemType="" Name="ux" Type="Scalar>"\n')
   key = set['grid_shape']
   file.write('\t<DataItem DataType="Float" Dimensions="{} {} {} 3" Format="HDF" Precision="8">ux_{}.h5</DataItem>\n'.format(key[0], key[1], key[2], filename))
   file.write('\t</Attribute>\n')

   file.write('\t<Attribute Center="Node" ElementCell="" ElementDegree="0" ElementFamily="" ItemType="" Name="uy" Type="Scalar>"\n')
   key = set['grid_shape']
   file.write('\t<DataItem DataType="Float" Dimensions="{} {} {} 3" Format="HDF" Precision="8">uy_{}.h5</DataItem>\n'.format(key[0], key[1], key[2], filename))
   file.write('\t</Attribute>\n')

   file.write('\t<Attribute Center="Node" ElementCell="" ElementDegree="0" ElementFamily="" ItemType="" Name="uz" Type="Scalar>"\n')
   key = set['grid_shape']
   file.write('\t<DataItem DataType="Float" Dimensions="{} {} {} 3" Format="HDF" Precision="8">uz_{}.h5</DataItem>\n'.format(key[0], key[1], key[2], filename))
   file.write('\t</Attribute>\n')

   file.write('\t</Grid>\n')
   file.write('</Grid>\n')
   file.write('</Domain>\n')
   file.write('</Xdmf>\n')

xdmf_create(name, set)