import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

direct = __file__.split('\\')
del direct[-1]; del direct[-1]; del direct[-1]
direct =  '\\'.join(direct)
sys.path.append(direct)

# eigenvalues_trim = np.loadtxt(direct + '\\output\\horten_u_inf2800_M4N11Msf5\\stability\\eigenvalues.dat')
# eigenvect_i = np.loadtxt(direct + '\\output\\horten_u_inf2800_M4N11Msf5\\stability\\eigenvectors_i.dat')
# eigenvect_r = np.loadtxt(direct + '\\output\\horten_u_inf2800_M4N11Msf5\\stability\\eigenvectors_r.dat')

# forward = 'c:\\Users\\zkrawcz\\Desktop\\GitUploads\\SHARPy_Scaled_Skydweller\\OSU_Contribution\\output\\Stability Results\\Skydwell_Forward_500_Eig'

# eigen_real = pd.read_excel(forward + '\\eigenvectors_r.xlsx')
# eigen_vals = pd.read_excel(forward + '\\eigenvaluetable.xlsx')
# eigen_vals = eigen_vals.values[2:,]

# rigid_mode_locations = [8493, 8501] # visual inspection look for the 1's in the file. Use the row references in excel

# interest = eigen_real.values[rigid_mode_locations[0]-2:rigid_mode_locations[1]-1]


# sort = np.sort(interest[4][1:,])

# top_suspects = []
# for i in range(-2, -10, -1):
#     top_suspects.append(sort[i])

# index = []
# eigen_vals_res = []
# nat_freq = []
# for i in top_suspects:
#     for j in range(1, len(interest[4])):
#         if i == interest[4][j]:
#             index.append(j)
#             eigen_vals_res.append([eigen_vals[j-1][2], eigen_vals[j-1][3]])
#             nat_freq.append(eigen_vals[j-1][4])

# print(index)
# print(eigen_vals_res)
# print(nat_freq)



# C:\Users\zkrawcz\Desktop\GitUploads\SHARPy_Scaled_Skydweller\OSU_Contribution\output\Initial Velocity Set
eigenvalues_forward = np.loadtxt('c:\\Users\\zkrawcz\\Desktop\\GitUploads\\SHARPy_Scaled_Skydweller\\OSU_Contribution\\output\\Stability Results' + '\\Skydwell_forward_500_Eig\\eigenvalues.dat')
# eigenvalues_og = np.loadtxt('c:\\Users\\zkrawcz\\Desktop\\GitUploads\\SHARPy_Scaled_Skydweller\\OSU_Contribution\\output\\Stability Results' + '\\Skydwell_OG\\eigenvalues.dat')
eigenvalues_back = np.loadtxt('c:\\Users\\zkrawcz\\Desktop\\GitUploads\\SHARPy_Scaled_Skydweller\\OSU_Contribution\\output\\Stability Results' + '\\Skydwell_back_500_Eig\\eigenvalues.dat')
# eigenvect_r = np.loadtxt('c:\\Users\\zkrawcz\\Desktop\\GitUploads\\SHARPy_Scaled_Skydweller\\OSU_Contribution\\output\\Initial Velocity Set' + '\\stability\\eigenvectors_r.dat')

# fig = plt.figure()
# eigenvect_r = np.array(eigenvect_r)
# plt.imshow(eigenvect_r, interpolation= 'none')
# # fig.colorbar(eigenvect_r)
# plt.grid()
# plt.xlim(0, 100)
# plt.ylim(9230, 8230)
# plt.xlabel('Real Part, Re(lambda) [rad/s]')
# plt.ylabel('Imaginary Part, Im(lambda) [rad/s]')
# plt.show()


fig = plt.figure(1)
plt.scatter(eigenvalues_back[:, 0], eigenvalues_back[:, 1],
           marker='x',
           color='k')
plt.ylim(-10, 10)
plt.grid()
plt.title('Payload Back 9.5m (from 4.18m)')
plt.xlabel('Real Part, Re(lambda) [rad/s]')
plt.ylabel('Imaginary Part, Im(lambda) [rad/s]')


fig = plt.figure(2)
plt.scatter(eigenvalues_forward[:, 0], eigenvalues_forward[:, 1],
           marker='x',
           color='k')
plt.ylim(-10, 10)
plt.grid()
plt.title('Payload Forward 1m (from 4.18m)')
plt.xlabel('Real Part, $Re(\lambda)$ [rad/s]')
plt.ylabel('Imaginary Part, $Im(\lambda)$ [rad/s]')


fig = plt.figure(3)
plt.scatter(eigenvalues_back[:, 0], eigenvalues_back[:, 1],
           marker='x',
           color='blue')
plt.scatter(eigenvalues_forward[:, 0], eigenvalues_forward[:, 1],
           marker='x',
           color='red')
plt.xlim(-4, 1)
plt.ylim(-10, 10)
plt.grid()
plt.title('Combined')
plt.xlabel('Real Part, $Re(\lambda)$ [rad/s]')
plt.ylabel('Imaginary Part, $Im(\lambda)$ [rad/s]');
plt.show()