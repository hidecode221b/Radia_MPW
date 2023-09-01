# Add paths to RadiaID and RadiaUtils
import sys
#sys.path.append('../PyRadia/RadiaID/') # RadiaInsertion Devices
#sys.path.append('../PyRadia/RadiaUtils/')# Magnetic materials, Radia export functions, etc.

# Import the radia_id module
import radia_id_mpw_hyb_vac79rd as rid
import pyvista as pv
import numpy as np

build = 'full'
# full, main, side, pole, side_pole, main_pole, main_side
# full_half makes upper array

# Main undulator parameters
period = 230 # Undulator period (mm)
gap = 15 # Udulator gap (mm)
#n_periods = 4 # Number of periods
n_poles = 5 # Number of poles

# Magnetization errors
sdr = 0 # standard deviation of amplitude errors (%/100)
sdt = 0 # standard deviation of angle errors (degrees)

mat = 'ndfeb' # Should be 'cryo', 'ndfeb' or 'smco'

# Build a parameter object according to the material specification
if mat == 'smco':
    # For Sm2Co17 magnets
    params = rid.HybridWigParam(period, n_poles, sdr, sdt, gap, mag_mat='sm2co17', br=1.1, wig_build=build)
elif mat == 'ndfeb':
    # For room temperature NdFeB (Br may be a bit optimistic)
    params = rid.HybridWigParam(period, n_poles, sdr, sdt, gap, mag_mat='ndfeb', br=1.29, wig_build=build)
elif mat == 'cryo':
    # For PrFeB magnets at 80 K and Vanadium Permadur poles
    params = rid.HybridWigParam(period, n_poles, sdr, sdt, gap, mag_mat='ndfeb', br=1.6, pole_mat='fecov', wig_build=build)

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
    und.plot_geo('EdgeLines->True')
    #'EdgeLines->True|False,Faces->True|False,Axes->True|False'
else:
    # Plot VTK by PyVISTA
    grid = pv.read(fileName+'.vtk')
    grid.plot(cmap='viridis',show_scalar_bar=False,show_axes=False,show_edges=True,window_size = [1000, 580],lighting=True,component=1)

# Plot results (need to plot geometry above)
#und.plot_field(xyz_end=[0, 1000, 0], xyz_start=[0, -1000, 0], n=1000, b='bz', x_axis='d', plot_show=True, plot_title='field_y')
"""
# Plot 3D vector field by pyvista
#und.plot_vector_field(-100,101,5,-200,201,5,-15,16,5,3,plot_save=True, plot_title=fileName)
#und.plot_field([50, 0, 0], [-50, 0, 0], n=1000, b='bz', x_axis='d', plot_show=True, plot_title='field_x')
#und.plot_field_int(xyz_end=[50, 0, 0], xyz_start=[-50, 0, 0], dir_int=[0, 1, 0], n=1000, b='bz', x_axis='d', plot_show=True, plot_title='field_int')
"""
# und.plot_traj(e=1.2, init_cond=[0, 0, 0, 0], y_range=None, n_points=1000, x_or_z='x', plot_show=True, plot_title='traj_y')
"""
und.plot_field(xyz_end=[0, 1000, 0], xyz_start=[0, -1000, 0], n=1000, b='bz', x_axis='d', plot_show=True, plot_title='field_y')
und.plot_field([50, 0, 0], [-50, 0, 0], n=1000, b='bz', x_axis='d', plot_show=True, plot_title='field_x')
und.plot_traj(e=1.2, init_cond=[0, 0, 0, 0], y_range=None, n_points=1000, x_or_z='x', plot_show=True, plot_title='traj_y')
und.plot_field_int(xyz_end=[0, 2000, 0], xyz_start=[0, -2000, 0], dir_int=[1, 0, 0], n=1000, b='bz', x_axis='d', plot_show=True, plot_title='field_int')
und.plot_field_int(xyz_end=[200, 0, 0], xyz_start=[-200, 0, 0], dir_int=[0, 1, 0], n=1000, b='bz', x_axis='d', plot_show=True, plot_title='field_int')
und.print_wavelength(e=1.2, n=1, theta=0)
und.plot_vector_field(-100,101,5,-200,201,5,-15,16,5,3,plot_save=True, plot_title=fileName)
"""
# setup force calc# calc force
#und.force(normal_plane=None, point=None)
#und.test_obj_export_vtk()
#und.force2(k=[1,1,4])
"""
for gap in [800,400,200,150,110,90,70,50,40,30,25,20,15,10,5,3]:
    params = rid.HybridWigParam(period, n_poles, sdr, sdt, gap, mag_mat='ndfeb', br=1.29, wig_build=build)
    und = rid.HybridWiggler(params)
    und.force2(k=[1,1,2])
"""

list_sd =[]
list_field_all =[]
for i in range(0, 11, 1):
    list_field =[]

    #for sdr in range(0, 11, 1):
    #    sdr = sdr/100
    #    if i==0:
    #        list_sd.append(sdr)

    for sdt in range(0, 21, 1):
        if i==0:
            list_sd.append(sdt)

        params = rid.HybridWigParam(period, n_poles, sdr, sdt, gap, mag_mat='ndfeb', br=1.29, wig_build=build)
        und = rid.HybridWiggler(params)
        #und.peak_field2(fileName)
        list_field.append(und.peak_field2(fileName))

    if i ==0:
        list_field_all.append(list_sd)
    list_field_all.append(list_field)
    print(i)
#print(list_field)
#list_sdr.append(list_field)
#print(list_sdr)
#print(list_field_all)
#np.savetxt(fileName+"_sdr001.csv",np.transpose(list_field_all),header='ME/eV,test',delimiter =",")
np.savetxt(fileName+"_sdt03.csv",np.transpose(list_field_all),header='ME/eV,test',delimiter =",")
