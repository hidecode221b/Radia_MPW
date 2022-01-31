# Add paths to RadiaID and RadiaUtils
import sys, os
#sys.path.append('../PyRadia/RadiaID/') # RadiaInsertion Devices
#sys.path.append('../PyRadia/RadiaUtils/')# Magnetic materials, Radia export functions, etc.

# Import the radia_id module
import radia_id_odd_cham_rtg_mp_ep as rid
import pyvista as pv
import numpy as np

build = 'full'
# full, main, side, pole, side_pole, main_pole, main_side

# Main undulator parameters
period = 220 # Undulator period (mm)
gap = 20 # Udulator gap (mm)
#n_periods = 4 # Number of periods
n_poles = 9 # Number of poles

mat = 'ndfeb' # Should be 'cryo', 'ndfeb' or 'smco'

# Build a parameter object according to the material specification
if mat == 'smco':
    # For Sm2Co17 magnets
    params = rid.HybridWigParam(period, n_poles, gap, mag_mat='sm2co17', br=1.1, wig_build=build)
elif mat == 'ndfeb':
    # For room temperature NdFeB (Br may be a bit optimistic)
    params = rid.HybridWigParam(period, n_poles, gap, mag_mat='ndfeb', br=1.29, wig_build=build)
elif mat == 'cryo':
    # For PrFeB magnets at 80 K and Vanadium Permadur poles
    params = rid.HybridWigParam(period, n_poles, gap, mag_mat='ndfeb', br=1.6, pole_mat='fecov', wig_build=build)

# Build and solve
und = rid.HybridWiggler(params)
fileName = 'mpw' + str(n_poles) + 'cs_g' + str(gap) + '_' + build
"""
# Load model exported previously
fileName = 'mpw2cs_g20_full'
und = rid.Undulator()
und.load(fileName)
"""
und.print_wavelength(e=1.2, n=1, theta=0)

# Export the parameters and model
und.save(fileName) # save parameters and model
# Save VTK
und.exportGeometryToVTK(fileName)
"""
"""
# Plot geometry
if sys.platform == "win32":
    und.plot_geo('EdgeLines->True')
    #'EdgeLines->True|False','Faces->True|False','Axes->True|False'
else:
    # Plot VTK by PyVISTA
    grid = pv.read(fileName+'.vtk')
    grid.plot(cmap='viridis',show_scalar_bar=False,show_axes=False,show_edges=True,window_size = [1000, 580],lighting=True,component=1)

# Plot results (need to plot geometry above)
#und.plot_field(xyz_end=[0, 1000, 0], xyz_start=[0, -1000, 0], n=1000, b='bz', x_axis='d', plot_show=False, plot_title='field_y')
"""
# Save 3D vector field based on csv files
i = 0
x0,x1,dx,y0,y1,dy,z0,z1,dz = -100,101,5,-200,201,5,-15,16,5

for j in range(z0,z1,dz):
    for k in range(y0,y1,dy):
        x, y, z, d, ba = und.field(xyz_end=[x0, k, j], xyz_start=[x1, k, j], n=int((x1-x0-1)/dx+1), b='bx')
        np.savetxt("test_x"+str("{:04.0f}".format(i))+"_bx.csv", np.transpose([d,ba]), header='ME/eV,test', comments='', delimiter=",")
        x, y, z, d, ba = und.field(xyz_end=[x0, k, j], xyz_start=[x1, k, j], n=int((x1-x0-1)/dx+1), b='by')
        np.savetxt("test_x"+str("{:04.0f}".format(i))+"_by.csv", np.transpose([d,ba]), header='ME/eV,test', comments='', delimiter=",")
        x, y, z, d, ba = und.field(xyz_end=[x0, k, j], xyz_start=[x1, k, j], n=int((x1-x0-1)/dx+1), b='bz')
        np.savetxt("test_x"+str("{:04.0f}".format(i))+"_bz.csv", np.transpose([d,ba]), header='ME/eV,test', comments='', delimiter=",")
        i = i + 1
"""
"""
und.plot_field(xyz_end=[0, 1000, 0], xyz_start=[0, -1000, 0], n=1000, b='bz', x_axis='d', plot_show=True, plot_title='field_y')
und.plot_field([50, 0, 0], [-50, 0, 0], n=1000, b='bz', x_axis='d', plot_show=True, plot_title='field_x')
und.plot_traj(e=1.2, init_cond=[0, 0, 0, 0], y_range=None, n_points=1000, x_or_z='x', plot_show=True, plot_title='traj_y')
und.plot_field_int(xyz_end=[0, 2000, 0], xyz_start=[0, -2000, 0], dir_int=[1, 0, 0], n=1000, b='bz', x_axis='d', plot_show=True, plot_title='field_int')
und.plot_field_int(xyz_end=[400, 0, 0], xyz_start=[-400, 0, 0], dir_int=[0, 1, 0], n=1000, b='bz', x_axis='d', plot_show=True, plot_title='field_int')
und.print_wavelength(e=1.2, n=1, theta=0)

for j in range(-50,51,50):
    for k in range(-20,21,20):
        und.plot_field(xyz_end=[j, 1000, k], xyz_start=[j, -1000, k], n=1000, b='bx', x_axis='d', plot_show=False, plot_title='field_x')
        und.plot_field(xyz_end=[j, 1000, k], xyz_start=[j, -1000, k], n=1000, b='by', x_axis='d', plot_show=False, plot_title='field_y')
        und.plot_field(xyz_end=[j, 1000, k], xyz_start=[j, -1000, k], n=1000, b='bz', x_axis='d', plot_show=False, plot_title='field_z')
for j in range(-50,51,5):
    for k in range(-200,201,20):
        und.field(xyz_end=[j, k, -15], xyz_start=[j, k, 15], n=7, b='bx')
        und.field(xyz_end=[j, k, -15], xyz_start=[j, k, 15], n=7, b='by')
        und.field(xyz_end=[j, k, -15], xyz_start=[j, k, 15], n=7, b='bz')
for j in range(z0,z1,dz):
    for k in range(x0,x1,dx):
        x, y, z, d, ba = und.field(xyz_end=[k, y0, j], xyz_start=[k, y1, j], n=int((y1-y0-1)/dy+1), b='bx')
        np.savetxt("test_x"+str("{:04.0f}".format(i))+"_bx.csv", np.transpose([d,ba]), header='ME/eV,test', comments='', delimiter=",")
        x, y, z, d, ba = und.field(xyz_end=[k, y0, j], xyz_start=[k, y1, j], n=int((y1-y0-1)/dy+1), b='by')
        np.savetxt("test_x"+str("{:04.0f}".format(i))+"_by.csv", np.transpose([d,ba]), header='ME/eV,test', comments='', delimiter=",")
        x, y, z, d, ba = und.field(xyz_end=[k, y0, j], xyz_start=[k, y1, j], n=int((y1-y0-1)/dy+1), b='bz')
        np.savetxt("test_x"+str("{:04.0f}".format(i))+"_bz.csv", np.transpose([d,ba]), header='ME/eV,test', comments='', delimiter=",")
        i = i + 1
"""
