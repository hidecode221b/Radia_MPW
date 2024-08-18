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
    params = rid.HybridWigParam(period, n_poles, sdr, sdt, gap, ext_pole=[9.5, 0], mag_mat='ndfeb', br=1.20, wig_build=build)
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

# Export the geometric parameters
#und.save(fileName) 

# Save VTK (geometrical model)
#und.exportGeometryToVTK(fileName)
"""
# Plot geometry
if sys.platform == "win32":
    und.plot_geo('EdgeLines->True')
    #'EdgeLines->True|False,Faces->True|False,Axes->True|False'
else:
    # Plot VTK by PyVISTA (save vtk first)
    grid = pv.read(fileName+'.vtk')
    grid.plot(cmap='viridis',show_scalar_bar=False,show_axes=False,show_edges=True,window_size = [1000, 580],lighting=True,component=1)
"""
# Plot field (need to plot geometry above)
#und.plot_field(xyz_end=[0, 1000, 0], xyz_start=[0, -1000, 0], n=1000, b='bz', x_axis='d', plot_show=True, plot_title='field_y')
#und.plot_field([50, 0, 0], [-50, 0, 0], n=1000, b='bz', x_axis='d', plot_show=True, plot_title='field_x')

# Plot 3D vector field by pyvista
#und.plot_vector_field(-100,100,5,-400,400,5,-5,5,5,3,plot_save=True, plot_title=fileName)

# Plot kick maps (second field integral)
#und.plot_vector_field_int(-100,100,1,-600,600,1,-5,5,1,method='fld_int',plot_save=True,plot_title='Kickmap')

# Plot field integral
#und.plot_field_int(xyz_end=[0, 2000, 0], xyz_start=[0, -2000, 0], dir_int=[1, 0, 0], n=1000, b='bz', x_axis='d', plot_show=True, plot_title='field_int')
#und.plot_field_int(xyz_end=[200, 0, 0], xyz_start=[-200, 0, 0], dir_int=[0, 1, 0], n=1000, b='bz', x_axis='d', plot_show=True, plot_title='field_int')

# Plot beam trajectory
#und.plot_traj(e=1.2, init_cond=[0, 0, 0, 0], y_range=None, n_points=1000, x_or_z='x', plot_show=True, plot_title='traj_y')

# Magnetic force (unknown)
#und.force(normal_plane=None, point=None)
#und.test_obj_export_vtk()
#und.force2(k=[1,1,4])

"""
# Dependence of field and field integral upon gap and ext pole
list_save  = [] # to be saved
#list_gap = [5,10,15,20,25,30,35,40,45,50]
list_gap = [15] # undulator gap
list_epole = [8,8.5,9,9.5,10,10.5,11] # ending pole length
list_field = [] # peak field

for gap in list_gap:
    for epole in list_epole:
        params = rid.HybridWigParam(period, n_poles, sdr, sdt, gap, ext_pole=[epole, 0], mag_mat='ndfeb', br=1.2, wig_build=build)
        und = rid.HybridWiggler(params)
        #list_field.append(und.peak_field_tesla(fileName))
        list_field.append(und.field_int([0,1000,0], xyz_start=None, n=100, b='bz', method='fld_int')[4])

print(list_field)
#list_save.append(list_gap)
list_save.append(list_epole)
list_save.append(list_field)
np.savetxt(fileName+"fld_int_epole.csv",np.transpose(list_save),header='ME/eV,test',delimiter =",")
"""
"""
# Studies on standard errors
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
        #und.peak_field_tesla(fileName)
        list_field.append(und.peak_field_tesla(fileName))

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
"""