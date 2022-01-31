import sys, os
import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv

# 3D vector field plot based on the csv files specified below
dir = "/Users/hidekinakajima/Downloads/Radia-master/env/radia_python/"
head = "test_x"

entries = os.listdir(dir)
entries.sort()
list_file = []
for entry in entries:
    fname = os.path.basename(entry)
    #print(fname[0:len(head)])
    if fname[0:len(head)] == head and os.path.splitext(entry)[1] == '.csv':
        list_file.append(str(dir + os.sep + entry))
        #print(fname)
x, y, z, u, v, w, m = [], [], [], [], [], [], []

for file_path in list_file:
    f = os.path.basename(file_path)

    if f[f.find('_b')+2:f.find('.')] == 'x':
        f_x = file_path
        #print(f[len(head):f.find('_b')])
        """
        f_y = os.path.dirname(file_path) + os.sep + head + f[f.find('x')+1:f.find('y')] + 'y' + f[f.find('y')+1:f.find('_b')] + '_by.csv'
        f_z = os.path.dirname(file_path) + os.sep + head + f[f.find('x')+1:f.find('y')] + 'y' + f[f.find('y')+1:f.find('_b')] + '_bz.csv'
        """
        f_y = os.path.dirname(file_path) + os.sep + head + f[len(head):f.find('_b')] + '_by.csv'
        f_z = os.path.dirname(file_path) + os.sep + head + f[len(head):f.find('_b')] + '_bz.csv'

        df_x = np.loadtxt(f_x, delimiter=',', skiprows=1)
        df_y = np.loadtxt(f_y, delimiter=',', skiprows=1)
        df_z = np.loadtxt(f_z, delimiter=',', skiprows=1)
        dim = df_x[:,0]
        m_x = df_x[:,1]
        m_y = df_y[:,1]
        m_z = df_z[:,1]
        #print(f, len(y0), float(f[f.find('x')+1:f.find('z')]))
        for j in range(len(dim)):
            """
            x = np.append(x, float(f[f.find('x')+1:f.find('y')]))
            y = np.append(y, float(f[f.find('y')+1:f.find('_b')]))
            z = np.append(z, z0[j])
            """
            u = np.append(u, m_x[j])
            v = np.append(v, m_y[j])
            w = np.append(w, m_z[j])
            #m = np.append(m, np.sqrt(v**2 + u**2 + w**2))
            #print(float(f[f.find('x')+1:f.find('z')]), f[f.find('z')+1:f.find('_b')], f[f.find('_b')+2:f.find('.')], x0[0], y0[0])
    else:
        pass

#np.savetxt(head + '.csv', x, header='ME/eV,test', comments='', delimiter=",")
"""
# COMPUTE LENGTH OF VECTOR -> MAGNITUDE
c = np.sqrt(np.abs(v)**2 + np.abs(u)**2 + np.abs(w)**2)
c = (c.ravel() - c.min())/c.ptp()
# Repeat for each body line and two head lines
c = np.concatenate((c, np.repeat(c, 2)))
# Colormap
c = plt.cm.jet(c)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
#ax.quiver(x, y, z, u, v, w, color = c)
for x1,y1,z1,u1,v1,w1,m0,c0 in zip(x,y,z,u,v,w,m,c):
    ax.quiver(x1, y1, z1, u1, v1, w1, pivot = 'middle', length = m0*0.5, color = c0, normalize=True)
plt.show()
"""

x0,x1,dx,y0,y1,dy,z0,z1,dz = -100,101,5,-200,201,5,-15,16,5
ptsx, ptsy, ptsz = 41j, 81j, 7j
#x0,x1,dx,y0,y1,dy,z0,z1,dz = 0,11,5,0,11,5,0,11,5
#x = np.arange(x0, x1, dx)
#y = np.arange(y0, y1, dy)
#z = np.arange(z0, z1, dz)
#x, y, z = np.meshgrid(x, y, z)
x, y, z = np.mgrid[x0:x1:ptsx, y0:y1:ptsy, z0:z1:ptsz]
#x, y, z = np.mgrid[-50:50:21j, -200:200:21j, -15:15:7j]
#x, y, z = np.mgrid[0:10:3j, 0:10:3j, 0:10:3j]
#x, y, z = np.mgrid[np.min(x):np.max(x):len(x)j, np.min(y):np.max(y):len(y)j, np.min(z):np.max(z):len(z)j]
grid = pv.StructuredGrid(x, y, z)
B = np.column_stack((u.ravel(), v.ravel(), w.ravel()))
grid["ABC field magnitude"] = np.linalg.norm(B, axis=1)
grid["ABC field vectors"] = B
grid.set_active_vectors("ABC field vectors")
arrows_grid = grid.glyph(orient="ABC field vectors", factor=5.0)

p = pv.Plotter()
cmap = plt.cm.get_cmap("viridis")
#p.add_mesh(grid, cmap=cmap)
p.show_grid(color='black')
p.add_mesh(arrows_grid, cmap=cmap)
#p.add_volume(grid, cmap='viridis', mapper='gpu', show_scalar_bar=False)
p.show()
