# Add paths to RadiaID and RadiaUtils
import sys
#sys.path.append('../PyRadia/RadiaID/') # RadiaInsertion Devices
#sys.path.append('../PyRadia/RadiaUtils/')# Magnetic materials, Radia export functions, etc.

# Import the radia_id module
import radia_id_odd_cham_rtg_mp_ep as rid
import pyvista as pv

build = 'full'
# full, main, side, pole, side_pole, main_pole, main_side
# full_half makes upper array

# Main undulator parameters
period = 220 # Undulator period (mm)
gap = 20 # Udulator gap (mm)
#n_periods = 4 # Number of periods
n_poles = 5 # Number of poles

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
fileName = 'mpw5cs_g20_full'
und = rid.Undulator()
und.load(fileName)
"""
und.print_wavelength(e=1.2, n=1, theta=0)

# Export the parameters and model
und.save(fileName) # save parameters and model
# Save VTK
und.exportGeometryToVTK(fileName)

# Plot geometry
if sys.platform == "win32":
    rad.ObjDrwOpenGL(self.obj, option)
    #rad.ObjDrwVTK(self.obj, option)
else:
    # Plot VTK by PyVISTA
    grid = pv.read(fileName+'.vtk')
    grid.plot(cmap='viridis',show_scalar_bar=False,show_axes=False,show_edges=True,window_size = [1000, 580],lighting=True,component=1)


# Plot results (need to plot geometry above)
#und.plot_field(xyz_end=[0, 1000, 0], xyz_start=[0, -1000, 0], n=1000, b='bz', x_axis='d', plot_show=True, plot_title='field_y')
# Plot 3D vector field by pyvista
#und.plot_vector_field(-100,101,5,-200,201,5,-15,16,5,3,plot_save=True, plot_title=fileName)
#und.plot_field([50, 0, 0], [-50, 0, 0], n=1000, b='bz', x_axis='d', plot_show=True, plot_title='field_x')
#und.plot_field_int(xyz_end=[50, 0, 0], xyz_start=[-50, 0, 0], dir_int=[0, 1, 0], n=1000, b='bz', x_axis='d', plot_show=True, plot_title='field_int')

"""
und.plot_field(xyz_end=[0, 1000, 0], xyz_start=[0, -1000, 0], n=1000, b='bz', x_axis='d', plot_show=True, plot_title='field_y')
und.plot_field([50, 0, 0], [-50, 0, 0], n=1000, b='bz', x_axis='d', plot_show=True, plot_title='field_x')
und.plot_traj(e=1.2, init_cond=[0, 0, 0, 0], y_range=None, n_points=1000, x_or_z='x', plot_show=True, plot_title='traj_y')
und.plot_field_int(xyz_end=[0, 2000, 0], xyz_start=[0, -2000, 0], dir_int=[1, 0, 0], n=1000, b='bz', x_axis='d', plot_show=True, plot_title='field_int')
und.plot_field_int(xyz_end=[200, 0, 0], xyz_start=[-200, 0, 0], dir_int=[0, 1, 0], n=1000, b='bz', x_axis='d', plot_show=True, plot_title='field_int')
und.print_wavelength(e=1.2, n=1, theta=0)
und.plot_vector_field(-100,101,5,-200,201,5,-15,16,5,3,plot_save=True, plot_title=fileName)

# setup force calc
poleName = 'mpw5cs_g20_pole_half'
pole = rid.Undulator()
pole.load(poleName)
grid = pv.read(poleName+'.vtk')
grid.plot(cmap='viridis',show_scalar_bar=False,show_axes=False,show_edges=True,window_size = [1000, 580],lighting=True,component=1)
# calc force
und.force(pole, und, normal_plane=None, point=None)
und.force2(pole, und)
"""
