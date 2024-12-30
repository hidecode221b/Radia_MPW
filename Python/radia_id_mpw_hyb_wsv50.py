# ----------------------------------------------------------
# PyRadiaUndulators library
# radia_undulator.py
# Library of Radia models for undulators
#
# Gael Le Bec, ESRF, 2019
# ----------------------------------------------------------
# Modified for Hybrid wiggler and 3D vector field
# VTK preview by PyVISTA
#
# Hideki Nakajima, SLRI, 2022
#
# VAC/Kyma in-vacuum wiggler configration without side magnets
# ----------------------------------------------------------
# --- Imports
import radia as rad
import radia_util as rad_uti
import radia_mat_util as rad_mat
from matplotlib.pyplot import plot, show, xlabel, ylabel, figure, title
from itertools import accumulate
import pickle
import numpy as np
from datetime import datetime
import csv
import pyvista as pv

# Constants
h = 4.135667516e-15 # Plank constant (eV s)
c = 299792458 # Speed of light (m/s)

# -------------------------------------------
# Parameters for undulator Radia models
# -------------------------------------------
class UndParam():
    """
    Base class for undulator parameters
    """


# -------------------------------------------
# Parameters for hybrid wigglers
# -------------------------------------------
class HybridWigParam(UndParam):
    """
    Radia parameters for hybrid wigglers
    """
    def __init__(self, period, n_poles, gap, mag_width=70, mag_height=[50, 40, 0], mag_chamfer=[5, 0, 5], pole_length=8, pole_width=60, pole_height=40, pole_chamfer=[5, 0, 5], ext_pole=[8, 0], ext_mag=[8.5, 0, 5], mag_mat='ndfeb', br=1.29, sep_exp=[0, 0, 0], pole_mat='fecov', mag_area_max=200, pole_area_max=100, mag_long_sub=[2, 2], pole_long_sub=[1, 1], mag_color=[[0, 1, 1],[0, 0.5, 1]], pole_color=[1, 0, 1], wig_build='full'):
        """
        Parameters for hybrid wiggler -- SLRI BL1 type
        :param period: period (mm)
        :param n_periods: number of periods
        :param gap: gap (mm)
        :param mag_width=120: width of the magnet blocks (mm)
        :param mag_height=[130, 110, 80]: height of the magnet blocks [main, side, ext_side] (mm)
        :param mag_chamfer=[10, 9, 10]: magnet chamfer [x, y, z] (mm)
        :param pole_length=35: length of the poles (mm)
        :param pole_width=70: width of the poles (mm)
        :param pole_height=90: height of the poles (mm)
        :param pole_chamfer=[9, 9, 10]: pole chamfer [x, y, z] (mm)
        :param ext_pole=[20,0]: extremity pole params [pole_length (mm), distance to previous obj (mm)]
        :param ext_mag=[37.5,0,10]: extremity mag params [mag_length (mm), distance to extremity pole (mm),extremity chamfer (mm)]]
        :param mag_mat='ndfeb': magnet block material
        :param br=1.3: remanent field (T)
        :param pole_mat='xc6': pole material
        :param mag_area_max=200: max triangle area in magnet blocks
        :param pole_area_max=100: max triangle area in poles
        :param mag_long_sub=[4,1]: long. subdivision in magnets [4,1]
        :param pole_long_sub=[6,2]: long. subdivision in poles [6,2]
        :param mag_color=[[0,1,1],[0,0.5,1]]: color of magnet blocks [main, side]
        :param pole_color=[1,0,1]: color of poles
        :param build='full': building elements
        """

        self.period = period
        self.n_poles = n_poles
        self.gap=gap
        self.mag_width = mag_width
        self.mag_height = mag_height
        self.mag_chamfer = mag_chamfer
        self.pole_length = pole_length
        self.pole_width = pole_width
        self.pole_height = pole_height
        self.pole_chamfer = pole_chamfer
        self.ext_pole = ext_pole
        self.ext_mag = ext_mag
        self.mag_mat = mag_mat
        self.br = br
        self.sep = sep_exp
        self.pole_mat = pole_mat
        self.mag_area_max = mag_area_max
        self.pole_area_max = pole_area_max
        self.mag_long_sub = mag_long_sub
        self.pole_long_sub = pole_long_sub
        self.mag_color = mag_color
        self.pole_color = pole_color
        self.wig_build = wig_build


# -------------------------------------------
# Precision parameters for Radia solver
# -------------------------------------------
class PrecParams():
    def __init__(self, tolerance=1e-5, max_iter=10000):
        """
        Initialize precision parameters for the Radia solver
        :param tolerance: tolerance (default = 1e-5)
        :param max_iter:  maximum number of iterations (default=10000)
        """
        self.tolerance = tolerance
        self.max_iter = max_iter

# -------------------------------------------
# Radia undulator class
# -------------------------------------------
class Undulator():
    """
    Radia model undulator
    """
    def set_gap(self, gap):
        """
        Build an undulator with a new gap (for PPM and hybrid planar undulators)
        :param gap: new gap
        """
        # --- Hold the new gap
        self.radia_und_param.gap = gap

        # --- Build a new undulator
        self.obj = self.build_undulator(self.sym)

        # --- Solve
        self.solve()

    def solve(self,  print_switch=True):
        """
        Radia solver
        :param  print_switch=True: quiet mode if False
        """
        # --- Get tolerance parameters
        tol = self.radia_prec_param.tolerance
        max_iter = self.radia_prec_param.max_iter
        # --- Solve
        self.time_start = datetime.today()
        if print_switch:
            t = self.time_start
            print(t.day, '/', t.month, '/', t.year, ' at ', t.hour, ':', t.minute, ':', t.second,
                  ' Start to solve the magnetization problem...')
        self.result = rad.Solve(self.obj, tol, max_iter)
        self.time_end = datetime.today()
        if print_switch:
            t = self.time_end
            print(t.day, '/', t.month, '/', t.year, ' at ', t.hour, ':', t.minute, ':', t.second,
                  ' Magnetization problem solved.')
        # --- Peak field
        self.peak_field = abs(rad.Fld(self.obj, 'bz', [0, 0, 0]))
        # --- Deflection parameter
        self.k = 0.09337 * self.radia_und_param.period * self.peak_field

    def field(self, xyz_end, xyz_start=[0, 0, 0], n=100, b='bz'):
        """
        Compute the field along a straight line
        :param xyz_end: end point [x, y, z]
        :param xyz_start: starting point [x, y, z] (default = [0, 0, 0])
        :param n: number of points (default = 100)
        :param b: field component (default = 'bz')
        :return: x, y, z, d, bz: positions, distance to initial point and field
        """
        # --- Sampling
        x0, y0, z0 = xyz_start
        x1, y1, z1 = xyz_end
        # Steps
        dx = (x1 - x0) / (n - 1)
        dy = (y1 - y0) / (n - 1)
        dz = (z1 - z0) / (n - 1)
        # Positions
        x = [x0 + k * dx for k in range(n)]
        y = [y0 + k * dy for k in range(n)]
        z = [z0 + k * dz for k in range(n)]
        xyz = [[x[k], y[k], z[k]] for k in range(n)]
        # Distance to initial point
        d = [((x[k] - x[k - 1]) ** 2 + (y[k] - y[k - 1]) ** 2 + (z[k] - z[k - 1]) ** 2) ** 0.5 for k in range(1, n)]
        d = list(accumulate(d))
        d.insert(0, 0)

        # --- Field computation
        ba = rad.Fld(self.obj, b, xyz)
        #np.savetxt("test_x"+str("{:+03.0f}".format(x[0]))+"y"+str("{:+04.0f}".format(y[0]))+"_"+b+".csv", np.transpose([d,ba]), header='ME/eV,test', comments='', delimiter=",")
        # --- Return
        return x, y, z, d, ba

    def plot_field(self, xyz_end, xyz_start=[0, 0, 0], n=100, b='bz', x_axis='d', plot_show=True, plot_title=''):
        """"
        Compute and plot the field along a straight line
        :param xyz_end: end point [x, y, z]
        :param xyz_start: starting point [x, y, z] (default = [0, 0, 0])
        :param n: number of points (default = 100)
        :param b: field component (default = 'bz')
        :param xaxis: defines the x axis of the plot, x_axis = 'x', 'y', 'z' or 'd' (default = 'd', i.e. distance)
        :param show: show the plot if True
        :param plot_title: plot title
        :return: a matplotlib figure
        """

        # --- Compute the field
        x, y, z, d, ba = self.field(xyz_end=xyz_end, xyz_start=xyz_start, n=n, b=b)
        # --- Plot
        fig = figure()
        if x_axis == 'x':
            l = x
        elif x_axis == 'y':
            l = y
        elif x_axis == 'z':
            l = z
        else:
            l = d
        plot(l, ba)
        xlabel('Distance (mm)')
        ylabel('Field (T)')
        title(plot_title)
        #np.savetxt("test_x"+str("{:+02.0f}".format(x[0]))+"z"+str("{:+02.0f}".format(z[0]))+"_"+b+".csv", np.transpose([l,ba]), header='ME/eV,test', comments='', delimiter=",")
        if plot_show:
            show()
        return fig

    def plot_vector_field(self,x0,x1,dx,y0,y1,dy,z0,z1,dz,fac=3,plot_save=True,plot_title='test'):
        """"
        Compute and plot the 3D field specfied in (x0,y0,z0:x1,y1,z1)
        :param x0,x1,dx: start, end, step for x axis, and so on
        :param fac: factor of arrow head size
        :return pyvista plot
        """

        nx = int((x1-x0-1)/dx) + 1
        ny = int((y1-y0-1)/dy) + 1
        nz = int((z1-z0-1)/dz) + 1
        ptsx = complex(0,nx)
        ptsy = complex(0,ny)
        ptsz = complex(0,nz)
        bx,by,bz = [],[],[]

        for j in range(z0,z1,dz):
            for k in range(y0,y1,dy):
                x, y, z, d, ba = self.field(xyz_end=[x0, k, j], xyz_start=[x1, k, j], n=nx, b='bx')
                bx = np.append(bx, ba)
                x, y, z, d, ba = self.field(xyz_end=[x0, k, j], xyz_start=[x1, k, j], n=nx, b='by')
                by = np.append(by, ba)
                x, y, z, d, ba = self.field(xyz_end=[x0, k, j], xyz_start=[x1, k, j], n=nx, b='bz')
                bz = np.append(bz, ba)

        x, y, z = np.mgrid[x0:x1:ptsx, y0:y1:ptsy, z0:z1:ptsz]
        grid = pv.StructuredGrid(x, y, z)
        B = np.column_stack((bx.ravel(), by.ravel(), bz.ravel()))
        #np.savetxt(fileName + "_B.csv", B, header='Bx,By,Bz', comments='', delimiter=",")
        grid["ABC field magnitude"] = np.linalg.norm(B, axis=1)
        grid["ABC field vectors"] = B
        grid.set_active_vectors("ABC field vectors")
        arrows_grid = grid.glyph(orient="ABC field vectors", factor=fac)
        if plot_save:
            arrows_grid.save(plot_title + '_v3d.vtk')

        p = pv.Plotter()
        #p.add_mesh(grid, cmap=cmap)
        p.show_grid(color='black')
        p.add_mesh(arrows_grid, cmap="viridis")
        p.show()

    def traj(self, e, init_cond=[0, 0, 0, 0], y_range=None, n_points=100):
        """
        Compute the trajectory of an electron in the undulator
        :param e: energy of the electron (GeV)
        :param init_cond: initial coordinates (mm) and angles (rad) [x0, dx/dy0, z0, dz/dy0] (default: [0, 0, 0, 0])
        :param y_range: initial and final value of the longitudinal coordinate (mm)
                                    (default: +/- 3 periods before the limits of the object)
        :param n_points: number of points (default: 100)
        :return: the trajectory: n_points lists of [y, x, dxdy, z, dzdy]
        """

        # --- Limits
        if y_range is None:
            lim = rad.ObjGeoLim(self.obj)[3] + 3 * self.radia_und_param.period
            #lim = self.radia_und_param.period * self.radia_und_param.n_periods / 2 + 3 * self.radia_und_param.period
            y_range = [-lim, lim]

        # --- Compute the trajectory
        trj = rad.FldPtcTrj(self.obj, e, init_cond, y_range, n_points)

        # --- Transpose
        trj_np = np.array(trj)
        trj_tr = np.ndarray.tolist(trj_np.transpose())

        # --- Return the trajectory
        return trj_tr

    def plot_traj(self, e, init_cond=[0, 0, 0, 0], y_range=None, n_points=100, x_or_z='x', plot_show=True, plot_title=''):
        """
        Compute and plot the trajectory of a particle in the undulator
        :param e: energy of the electron (GeV)
        :param init_cond: initial coordinates (mm) and angles (rad) [x0, dx/dy0, z0, dz/dy0] (default: [0, 0, 0, 0])
        :param y_range: initial and final value of the longitudinal coordinate (mm)
                            (default: +/- 3 periods before the limits of the object)
        :param n_points: number of points (default: 100)
        :param x_or_z: select the trajectory component to plot
                        -- 'x': horizontal position (default)
                        -- 'z': vertical position
                        -- 'dxdy': horizontal angle
                        -- 'dzdy': vertical angle
        :param plot_show: show the plot if True
        :param plot_title: title of the plot
        :return: a matplotlib figure
        """

        n_poles = self.radia_und_param.n_poles
        gap = self.radia_und_param.gap
        sub_div = self.radia_und_param.pole_long_sub[0]
        # --- Compute the trajectory
        y, x, dxdy, z, dzdy = self.traj(e, init_cond, y_range, n_points)
        # --- Plot
        fig = figure()
        if x_or_z == 'x':
            trj = x
            ylab = 'Horiz. position (mm)'
        elif x_or_z == 'z':
            trj = z
            ylab = 'Vert. position (mm)'
        elif x_or_z == 'dxdy':
            trj = dxdy
            ylab = 'Horiz. angle (rad)'
        elif x_or_z == 'dzdy':
            trj = dzdy
            ylab = 'Vert. angle (rad)'
        else:
            trj = y
            ylab = 'Long. position (mm)'
        plot(y, trj)

        np.savetxt("traj_"+ x_or_z + "_g"+ str("{:02.0f}".format(gap)) + "_p" +  str("{:02.0f}".format(n_poles)) + "_n"+ str("{:01.0f}".format(sub_div)) + ".csv", np.transpose([y,trj]), header='ME/eV,test', comments='', delimiter=",")

        xlabel('Long. position (mm)')
        ylabel(ylab)
        title(plot_title)
        if plot_show:
            show()
        return fig

    def field_int(self, xyz_end, xyz_start=None, n=100, b='bz', method='fld_int'):
        """
        Compute the field integral along a straight line
        :param xyz_end: end point [x, y, z]
        :param xyz_start=None: starting point [x, y, z] (-xyz_end if None)
        :param n: number of points (default = 100)
        :param b: field component (default = 'bz')
        :param method: computation method:  -- 'fld_int': use rad.FldInt
                                            -- 'fld_int_fin': use rad.FldInt('fin')
                                            -- 'fld': numerical integration of field from rad.Fld
        :return: x, y, z, d, ib, bz: positions, distance to initial point, integrated field and field
        """
        if xyz_start is None:
            xyz_start = [-xyz for xyz in xyz_end]
        if method == 'fld_int':
            ib = rad.FldInt(self.obj, 'inf', b, xyz_start, xyz_end)
            x, y, z, d, b = None, None, None, None, None
        if method == 'fld_int_fin':
            ib = rad.FldInt(self.obj, 'fin', b, xyz_start, xyz_end)
            x, y, z, d, b = None, None, None, None, None
        if method == 'fld':
            # --- Compute the field
            x, y, z, d, b = self.field(xyz_end, xyz_start=xyz_start, n=n, b=b)
            # --- Integrate
            dxyz = ( (xyz_end[0] - xyz_start[0]) ** 2
                   + (xyz_end[1] - xyz_start[1]) ** 2
                   + (xyz_end[2] - xyz_start[2]) ** 2 ) ** 0.5 / (n - 1)
            ib = dxyz * np.cumsum(b)
        # --- return
        return x, y, z, d, ib, b

    def plot_field_int(self, xyz_end, xyz_start=[0, 0, 0], dir_int=[0, 1, 0], n=100, b='bz', x_axis='d', plot_show=True, plot_title=''):
        """"
        Compute and plot the field along a straight line
        :param xyz_end: end point [x, y, z]
        :param xyz_start: starting point [x, y, z] (default = [0, 0, 0])
        :param dir_int: direction of the integration (default = [0, 1, 0])
        :param n: number of points (default = 100)
        :param b: field component (default = 'bz')
        :param method: computation method:  -- 'fld_int': use rad.FldInt
                                            -- 'fld_int_fin': use rad.FldInt('fin')
                                            -- 'fld': numerical integration of field from rad.Fld
        :param xaxis: defines the x axis of the plot, x_axis = 'x', 'y', 'z' or 'd' (default = 'd', i.e. distance)
        :param show: show the plot if True
        :param plot_title: plot title
        :return: a matplotlib figure
        """

        # --- Init
        dxyz = [(xyz_end[k] - xyz_start[k]) / (n - 1) for k in range(3)]
        x, y, z, d = np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n)
        ib = np.zeros(n)
        # --- Compute the field
        for k in range(n):
            x[k], y[k], z[k] = xyz_start[0] + k * dxyz[0], xyz_start[1] + k * dxyz[1], xyz_start[2] + k * dxyz[2]
            xyz_end_k = [x[k] + 1000 * dir_int[0], y[k] + 1000 * dir_int[1], z[k] + 1000 * dir_int[2]]
            xyz_start_k = [x[k] - 1000 * dir_int[0], y[k] - 1000 * dir_int[1], z[k] - 1000 * dir_int[2]]
            xk, yk, zk, dk, ib[k], bk = self.field_int(xyz_end=xyz_end_k, xyz_start=xyz_start_k, b=b)

        for k in range(n-1):
            d[k+1] = d[k]  + ((x[k+1] - x[k]) ** 2 + (y[k+1] - y[k]) ** 2 + (z[k+1] - z[k]) ** 2) ** 0.5

        # --- Plot
        fig = figure()
        if x_axis == 'x':
            l = x
        elif x_axis == 'y':
            l = y
        elif x_axis == 'z':
            l = z
        else:
            l = d
        plot(l, ib)
        xlabel('Distance (mm)')
        ylabel('Field integral (T mm)')
        title(plot_title)
        if plot_show:
            show()
        return fig

    def phase_int(self, xyz_end, xyz_start=None, n=100, b='bz'):
        """
        Compute the phase integral of the device
        :param xyz_end: end point [x, y, z]
        :param xyz_start=None: starting point [x, y, z] (-xyz_end if None)
        :param n: number of points (default = 100)
        :param b: field component (default = 'bz')
        :return: x, y, z, d, phase_int: positions, distance to initial point, phase integral (T^2mm^3)
        """

        # --- Field integral
        x, y, z, d, ib, bz = self.field_int(xyz_end, xyz_start=xyz_start, n=n, b=b, method='fld')

        # --- b^2 integral
        if xyz_start is not None:
            dxyz = ((xyz_end[0] - xyz_start[0]) ** 2
                    + (xyz_end[1] - xyz_start[1]) ** 2
                    + (xyz_end[2] - xyz_start[2]) ** 2) ** 0.5 / (n - 1)
        else:
            dxyz = ((2 * xyz_end[0]) ** 2 + (2 * xyz_end[1]) ** 2 + (2 * xyz_end[2] ) ** 2) ** 0.5 / (n - 1)
        ib2 = dxyz * np.cumsum(ib * ib)

        # --- Return
        return x, y, z, d, ib2

    def plot_geo(self, option):
        """
        Plot the undulator geometry
        """
        rad.ObjDrwOpenGL(self.obj, option)
        #rad.ObjDrwVTK(self.obj, option)

    def plot_geo2(self, obj, option):
        """
        Plot the specific geometry
        """
        rad.ObjDrwOpenGL(obj, option)
        #rad.ObjDrwVTK(self.obj, option)

    def save(self, filename, path=''):
        """
        Save the undulator object to filename.rad and its parameters to filename.radp.
        Use the load() method for loading the undulator.

        The filename.rad can be opened with the Mathematica interface, using the RadUtiLoad[] function

        :param filename: file name without extension
        :param path: absolute path if specified (default = '', i.e. relative path)
        """

        # --- Save the undulator
        # Dump the undulator
        dmp = rad.UtiDmp(self.obj, 'bin')
        # Write to a file
        f = open(path + filename + '.rad', 'wb')
        f.write(dmp)
        f.close()

        # --- Save the parametres
        f = open(path + filename + '.radp', 'wb')
        pickle.dump(self.radia_und_param, f)
        pickle.dump(self.radia_prec_param, f)
        f.close()

    def export(self, filename, path=''):
        """
        Save the Radia model of the undulator to filename.rad
        :param filename: file name without extension
        :param path: absolute path if specified (default = '', i.e. relative path)
        :return: True
        """
        rad_uti.save(self.obj, filename, path=path)

    def load(self, filename, path=''):
        """"
        Load the undulator to filename.rad and its parameters to filename.radp.

        The filename.und can be also opened with the Mathematica interface, using the RadUtiLoad[] function

        :param filename: file name without extension
        :param path: absolute path if specified (default = '', i.e. relative path)
        :return: True if no error, else False
        """

        # --- Load the parametres
        try:
            f = open(path + filename + '.radp', 'rb')
            self.radia_und_param = pickle.load(f)
            self.radia_prec_param = pickle.load(f)
            f.close()
        except FileNotFoundError:
            print('Error: Para file not found')
            return

        # --- Load the undulator
        try:
            f = open(path + filename + '.rad', 'rb')
            und_bin = f.read()
            f.close()
            self.obj = rad.UtiDmpPrs(und_bin)
        except FileNotFoundError:
            print('Error: Model file not found')
            return

    def force(self, normal_plane=None, point=None):
        """
        Compute the forces from one part of the object to the other part
        The object is cut in two subparts defined by a plane and a point
        :param normal_plane=None: normal vector [nx, ny, nz] ([0, 0, 1] if None)
        :param point=None: [0, 0, 0] if None
        :return: fx, fy, fz : magnetic forces (N)
        """

        # --- Initialize
        if normal_plane is None:
            normal_plane = [0, 0, 1]
            inverse_plane = [0, 0, -1]
        if point is None:
            point = [0, 0, 0]
        # --- Duplicate and cut
        obj_dpl = rad.ObjDpl(self.obj, 'FreeSym->True')
        src = rad.ObjCutMag(obj_dpl, point, normal_plane)
        dis = rad.ObjCutMag(obj_dpl, point, inverse_plane)
        #obj_0, obj_1 = rad.ObjCutMag(obj_dpl, point, normal_plane)
        # --- Compute the forces
        # FORCES COMPUTATIONS NOT YET AVAILABLE WITH RADIA PYTHON !!!!
        fx, fy, fz = rad.FldEnrFrc(src[0], self.obj,'fx|fy|fz',[1, 1, 2])
        print('Magnetic forces (fx): ', fx, ' N')
        print('Magnetic forces (fy): ', fy, ' N')
        print('Magnetic forces (fz): ', fz, ' N')
        #return fx, fy, fz

    def force2(self, k):
        """
        Compute the forces from whole structure to half of pole array
        Get half of pole array from the global variable in class
        self.obj is full (whole) structure
        k specifies subdivision of destination object [1,1,2]
        because force above does not work
        :return: fx, fy, fz : magnetic forces (N)
        """

        # --- Compute the forces
        # FORCES COMPUTATIONS NOT YET AVAILABLE WITH RADIA PYTHON !!!!
        ft= rad.FldEnrFrc(und_frc, self.obj,'f',k)
        print('Magnetic forces (fx,fy,fz): ', ft, ' Newton')
        #return fx, fy, fz

    def test_obj_export_vtk(self):
        """
        export und_frc in vtk to check geometry
        """

        self.exportGeometryToVTK2(und_frc, "und_frc")
        self.plot_geo2(und_frc, 'EdgeLines->True')

    def wavelength(self, e=6, n=1, theta=0):
        """
        Return the radiation wavelength and energy
        :param e=6: energy of the electron beam (GeV)
        :param n=1: harmonic number
        :param theta: observer angle (rad)
        :return: wavelength (nm), energy (keV)
        """
        try:
            sep_y = self.radia_und_param.sep[1]
            period = self.radia_und_param.period + 6 * sep_y
            n_poles = self.radia_und_param.n_poles
            n_half_poles = int((n_poles - 1) / 2)
            p_len = self.radia_und_param.pole_length
            ly_mag = (period - 2 * p_len) / 4
            pole_ext_length = self.radia_und_param.ext_pole[0]
            mag_ext_length = self.radia_und_param.ext_mag[0]
            d_pole = self.radia_und_param.ext_pole[1]
            d_mag = self.radia_und_param.ext_mag[1]
            d_pole = sep_y
            d_mag = sep_y
            if n_poles % 2 == 0:
                total_length = period * n_poles / 2 + (ly_mag + pole_ext_length + d_pole + d_mag + mag_ext_length) * 2
            else:
                total_length = period * n_half_poles + (p_len / 2 + ly_mag * 2 + pole_ext_length + d_pole + d_mag + mag_ext_length) * 2
            # --- Peak field
            if n_poles % 2 == 0:
                peak_field = abs(rad.Fld(self.obj, 'bz', [0, -1 * (ly_mag + p_len / 2), 0]))
            else:
                peak_field = abs(rad.Fld(self.obj, 'bz', [0, 0, 0]))
            # --- Deflection parameter
            k = 0.09337 * period * peak_field

            gamma = 1957 * e
            lambda_r = 0.001 * period / (2 * n * gamma ** 2) * (1 + k ** 2 / 2 + theta ** 2) # (m)
            e_r = h * c / lambda_r # (eV)
            return 1e9 * lambda_r, 1e-3 * e_r, total_length, peak_field, k
        except AttributeError:
            print('Error: Missing attribute in wavelength()')
            return None, None, None, None, None

    def print_wavelength(self, e=6, n=1, theta=0):
        """
        Compute and print the wavelength and energy
        :param e=6: energy of the electron beam (GeV)
        :param n=1: harmonic number
        :param theta: observer angle (rad)
        :return: None
        """

        lambda_r, e_r, total_length, peak_field, k = self.wavelength(e=e, n=n, theta=theta)

        if lambda_r is not None:
            print('Wavelength: ', lambda_r, ' nm')
            print('Energy: ', e_r, 'keV')

        print('Peak field: ', peak_field, ' T')
        print('K: ', k)
        print('Total length: ', total_length, ' mm')

# https://github.com/ochubar/Radia/issues/17
    def chunks(self, lst, n):
        """Yield successive n-sized chunks from a list called 'lst'."""
        for i in range(0, len(lst), n):
            yield lst[i:i + n]

    def exportGeometryToVTK(self, fileName='radia_Geometry'):
        '''
        Writes the geometry of RADIA object "obj" to file fileName.vtk for use in Paraview. The format is VTK legacy, because it's simple. The file consists of polygons only (no cells).
        '''

        vtkData = rad.ObjDrwVTK(self.obj, 'Axes->False')

        lengths = vtkData['polygons']['lengths']
        nPoly = len(lengths)
        offsets = list(accumulate(lengths))
        offsets.insert(0, 0) # prepend list with a zero
        points = vtkData['polygons']['vertices']
        nPnts = int(len(points)/3)

        # format the points array to be floats rather than double
        points = [round(num, 8) for num in points]

        # define the connectivity list
        conn = list(range(nPnts))

        # define colours array
        colors = vtkData['polygons']['colors']

        # pre-process the output lists to have chunkLength items per line
        chunkLength = 9 # this writes 9 numbers per line (9 is the number used in Paraview if data is saved as the VTK Legacy format)
        offsets = list(self.chunks(offsets, chunkLength))
        points = list(self.chunks(points, chunkLength))
        conn = list(self.chunks(conn, chunkLength))
        colors = list(self.chunks(colors, chunkLength))

        # write the data to file
        with open(fileName + ".vtk", "w", newline="") as f:
            f.write('# vtk DataFile Version 5.1\n')
            f.write('vtk output\nASCII\nDATASET POLYDATA\n')
            f.write('POINTS ' + str(nPnts) + ' float\n')

            writer = csv.writer(f, delimiter=" ")
            writer.writerows(points)
            f.write('\n')
            f.write('POLYGONS ' + str(nPoly+1) + ' ' + str(nPnts) + '\n')
            f.write('OFFSETS vtktypeint64\n')
            writer.writerows(offsets)
            f.write('CONNECTIVITY vtktypeint64\n')
            writer.writerows(conn)
            f.write('\n')
            f.write('CELL_DATA ' + str(nPoly) + '\n')
            f.write('COLOR_SCALARS Radia_colours 3\n')
            writer.writerows(colors)

    def exportGeometryToVTK2(self, obj, fileName='radia_Geometry'):
        '''
        Writes the geometry of RADIA object "obj" to file fileName.vtk for use in Paraview. The format is VTK legacy, because it's simple. The file consists of polygons only (no cells).
        '''

        vtkData = rad.ObjDrwVTK(obj, 'Axes->False')

        lengths = vtkData['polygons']['lengths']
        nPoly = len(lengths)
        offsets = list(accumulate(lengths))
        offsets.insert(0, 0) # prepend list with a zero
        points = vtkData['polygons']['vertices']
        nPnts = int(len(points)/3)

        # format the points array to be floats rather than double
        points = [round(num, 8) for num in points]

        # define the connectivity list
        conn = list(range(nPnts))

        # define colours array
        colors = vtkData['polygons']['colors']

        # pre-process the output lists to have chunkLength items per line
        chunkLength = 9 # this writes 9 numbers per line (9 is the number used in Paraview if data is saved as the VTK Legacy format)
        offsets = list(self.chunks(offsets, chunkLength))
        points = list(self.chunks(points, chunkLength))
        conn = list(self.chunks(conn, chunkLength))
        colors = list(self.chunks(colors, chunkLength))

        # write the data to file
        with open(fileName + ".vtk", "w", newline="") as f:
            f.write('# vtk DataFile Version 5.1\n')
            f.write('vtk output\nASCII\nDATASET POLYDATA\n')
            f.write('POINTS ' + str(nPnts) + ' float\n')

            writer = csv.writer(f, delimiter=" ")
            writer.writerows(points)
            f.write('\n')
            f.write('POLYGONS ' + str(nPoly+1) + ' ' + str(nPnts) + '\n')
            f.write('OFFSETS vtktypeint64\n')
            writer.writerows(offsets)
            f.write('CONNECTIVITY vtktypeint64\n')
            writer.writerows(conn)
            f.write('\n')
            f.write('CELL_DATA ' + str(nPoly) + '\n')
            f.write('COLOR_SCALARS Radia_colours 3\n')
            writer.writerows(colors)

# -------------------------------------------
# Planar Hybrid wiggler class
# -------------------------------------------
class HybridWiggler(Undulator):
    def __init__(self, radia_und_param, radia_prec_param=None, sym=True, solve_switch=True, print_switch=False):
        # Hold parameters
        self.radia_und_param = radia_und_param

        # Precision parameters
        if radia_prec_param is None:
            self.radia_prec_param = PrecParams()
        else:
            self.radia_prec_param = radia_prec_param

        # Build an undulator
        self.obj = self.build_undulator(sym)

        # Solve
        if solve_switch:
            self.solve(print_switch)

    def build_block_main(self, mag_id=0, half_period_id=0):
        """
        Build a pm block for a hybrid wiggler
        :param mag_id: identifies the PM block in the period: 0: +H, 1: -H, 2: -H, 3: +H (default: 0)
        :param period_id: identifies the period in the undulator
        :return: the PM block
        """

        period = self.radia_und_param.period
        n_poles = self.radia_und_param.n_poles
        ly_pole = self.radia_und_param.pole_length
        # --- Dimensions
        lx_mag = self.radia_und_param.mag_width
        ly_mag = (period - 2 * ly_pole) / 4
        lz_mag = self.radia_und_param.mag_height[0]
        cx_mag = self.radia_und_param.mag_chamfer[0]
        cy_mag = self.radia_und_param.mag_chamfer[1]
        cz_mag = self.radia_und_param.mag_chamfer[2]
        sep_x = self.radia_und_param.sep[0]
        sep_y = self.radia_und_param.sep[1]
        sep_z = self.radia_und_param.sep[2]
        period = period + 6 * sep_y
        if n_poles % 2 == 0 and n_poles > 0:
            sh_ep = ly_mag + ly_pole / 2
        else:
            sh_ep = 0
        # --- Position
        if mag_id == 0:
            pos_y = ly_pole / 2 + ly_mag / 2 + half_period_id * period / 2 + sep_y - sh_ep
        elif mag_id == 1:
            if n_poles % 2 == 0:
                if half_period_id == 0:
                    pos_y = ly_pole / 2 + 1.5 * ly_mag + half_period_id * period / 2 + sep_y / 2 - sh_ep
                else:
                    pos_y = ly_pole / 2 + 1.5 * ly_mag + half_period_id * period / 2 + sep_y - sh_ep
            else:
                pos_y = ly_pole / 2 + 1.5 * ly_mag + half_period_id * period / 2 + sep_y * 2 - sh_ep
        elif mag_id == 2:
            if n_poles % 2 == 0:
                pos_y = 1.5 * ly_pole + 2.5 * ly_mag + half_period_id * period / 2 + sep_y * 2.5 - sh_ep
            else:
                pos_y = 1.5 * ly_pole + 2.5 * ly_mag + half_period_id * period / 2 + sep_y * 4 - sh_ep
        elif mag_id == 3:
            pos_y = 1.5 * ly_pole + 3.5 * ly_mag + half_period_id * period / 2 + sep_y * 5 - sh_ep
        # --- Material
        mat = rad_mat.set_pm_mat(self.radia_und_param.mag_mat, self.radia_und_param.br)
        if n_poles % 2 == 0 and n_poles > 0:
            if mag_id == 1:
                if half_period_id % 2 == 0:
                    axis = [0, 1, 0]
                else:
                    axis = [0, -1, 0]
            else:
                if half_period_id % 2 == 0:
                    axis = [0, -1, 0]
                else:
                    axis = [0, 1, 0]
        else:
            if half_period_id % 2 == 0:
                axis = [0, 1, 0]
            else:
                axis = [0, -1, 0]
        # --- Build the block (extrusion)
        if mag_id == 0 or mag_id == 2:
            m1=[[lx_mag/4-cx_mag/2,pos_y-cy_mag/2, sep_z],[lx_mag/2-cx_mag,ly_mag+cy_mag]]
        elif mag_id == 1 or mag_id == 3:
            m1=[[lx_mag/4-cx_mag/2,pos_y+cy_mag/2, sep_z],[lx_mag/2-cx_mag,ly_mag+cy_mag]]
        m2=[[lx_mag/4,pos_y,cz_mag + sep_z],[lx_mag/2,ly_mag]]
        m3=[[lx_mag/4,pos_y,lz_mag + sep_z],[lx_mag/2,ly_mag]]
        pm=rad.ObjMltExtRtg([m1,m2,m3],axis)
        # --- Longitudinal subdivision
        rad.ObjDivMag(pm, [self.radia_und_param.mag_long_sub[1], self.radia_und_param.mag_long_sub[0], self.radia_und_param.mag_long_sub[1]])
        # --- Chamfers
        # --- Set the material
        rad.MatApl(pm, mat)
        # --- Set the color
        color = self.radia_und_param.mag_color[0]
        #if axis[1] > 0:
            #color = [0, 1, 1]
        #else:
            #color = [0, 1, 0.8]
        rad.ObjDrwAtr(pm, color)
        # --- Return the block
        return pm

    def build_block_side(self, mag_id=0, half_period_id=0):
        """
        Build a side PM block for hybrid wiggler
        :param mag_id: block type : 0 (-X) or 1 (+X)
        :param period_id: period number
        :return: pm block
        """

        period = self.radia_und_param.period
        n_poles = self.radia_und_param.n_poles
        ly_pole = self.radia_und_param.pole_length
        lx_pole = self.radia_und_param.pole_width
        # --- Dimensions
        lx_mag = self.radia_und_param.mag_width
        ly_mag = (period - 2 * ly_pole) / 4
        lz_mag = self.radia_und_param.mag_height[1]
        cx_mag = self.radia_und_param.mag_chamfer[0]
        cy_mag = self.radia_und_param.mag_chamfer[1]
        cz_mag = self.radia_und_param.mag_chamfer[2]
        cx_pole = self.radia_und_param.pole_chamfer[0]
        cy_pole = self.radia_und_param.pole_chamfer[1]
        cz_pole = self.radia_und_param.pole_chamfer[2]
        sep_x = self.radia_und_param.sep[0]
        sep_y = self.radia_und_param.sep[1]
        sep_z = self.radia_und_param.sep[2]
        period = period + 6 * sep_y
        n_half_poles = int((n_poles - 1) / 2)
        if n_poles % 2 == 0 and n_poles > 0:
            sh_ep = ly_mag + ly_pole / 2
        else:
            sh_ep = 0
        # --- Position and axis
        #print(mag_id, str(half_period_id), str(half_period_id % 2))
        # odd n_poles
        if mag_id == 0:
            pos_y = ly_pole / 4 + half_period_id * period / 2
            if half_period_id % 2 == 0:
                axis = [1, 0, 0]
            else:
                axis = [-1, 0, 0]
        elif mag_id == 1:
            pos_y = (period / 2 - ly_pole / 4) + half_period_id * period / 2
            if half_period_id % 2 == 0:
                axis = [-1, 0, 0]
            else:
                axis = [1, 0, 0]
        # even n_poles
        elif mag_id == 2:
            pos_y = (period - sep_y * 6) / 2 + half_period_id * period / 2 + sep_y * 1.5 - sh_ep
            if half_period_id % 2 == 0:
                axis = [-1, 0, 0]
            else:
                axis = [1, 0, 0]
        # --- Material
        mat = rad_mat.set_pm_mat(self.radia_und_param.mag_mat, self.radia_und_param.br)
        # --- Build the block (extrusion)
        if mag_id == 0 or mag_id == 1:
            if mag_id == 0:
                s1=[[(lx_pole-cx_pole+(lx_mag-lx_pole)/2-cx_mag)/2+sep_x,pos_y-cy_pole/2,sep_z],[cx_pole+(lx_mag-lx_pole)/2-cx_mag,ly_pole/2-cy_pole]]
            else:
                s1=[[(lx_pole-cx_pole+(lx_mag-lx_pole)/2-cx_mag)/2+sep_x,pos_y+cy_pole/2,sep_z],[cx_pole+(lx_mag-lx_pole)/2-cx_mag,ly_pole/2-cy_pole]]
            s2=[[(lx_pole+(lx_mag-lx_pole)/2)/2 + sep_x,pos_y,cz_mag + sep_z],[(lx_mag-lx_pole)/2,ly_pole/2]]
            s3=[[(lx_pole+(lx_mag-lx_pole)/2)/2 + sep_x,pos_y,lz_mag + sep_z],[(lx_mag-lx_pole)/2,ly_pole/2]]
        else:
            s1=[[(lx_pole-cx_pole+(lx_mag-lx_pole)/2-cx_mag)/2+sep_x,pos_y,sep_z],[cx_pole+(lx_mag-lx_pole)/2-cx_mag,ly_pole-2*cy_pole]]
            s2=[[(lx_pole+(lx_mag-lx_pole)/2)/2+sep_x,pos_y,cz_mag+sep_z],[(lx_mag-lx_pole)/2,ly_pole]]
            s3=[[(lx_pole+(lx_mag-lx_pole)/2)/2+sep_x,pos_y,lz_mag+sep_z],[(lx_mag-lx_pole)/2,ly_pole]]
        pm=rad.ObjMltExtRtg([s1,s2,s3],axis)
        # --- Longitudinal subdivision
        if mag_id == 0 or mag_id == 1:
            rad.ObjDivMag(pm, [self.radia_und_param.mag_long_sub[0], self.radia_und_param.mag_long_sub[1], self.radia_und_param.mag_long_sub[1]])
        else:
            rad.ObjDivMag(pm, [self.radia_und_param.mag_long_sub[0], self.radia_und_param.mag_long_sub[1] * 2, self.radia_und_param.mag_long_sub[1]])
        # --- Set the material
        rad.MatApl(pm, mat)
        # --- Set the color
        color = self.radia_und_param.mag_color[1]
        #if axis[0] > 0:
            #color = [0, 0.5, 1]
        #else:
            #color = [0, 0.8, 1]
        rad.ObjDrwAtr(pm, color)
        # --- Return the block
        return pm

    def build_pole(self, pole_id=0, half_period_id=0):
        """
        Build a iron pole for a hybrid undulator
        :param pole_id: identifies the pole in the period: 0: first pole, 1: second pole (default: 0)
        :param period_id: identifies the period in the undulator
        :return: the pole
        """

        period = self.radia_und_param.period
        n_poles = self.radia_und_param.n_poles
        # --- Dimensions
        lx_pole = self.radia_und_param.pole_width
        ly_pole = self.radia_und_param.pole_length
        lz_pole = self.radia_und_param.pole_height
        cx_pole = self.radia_und_param.pole_chamfer[0]
        cy_pole = self.radia_und_param.pole_chamfer[1]
        cz_pole = self.radia_und_param.pole_chamfer[2]
        sep_x = self.radia_und_param.sep[0]
        sep_y = self.radia_und_param.sep[1]
        sep_z = self.radia_und_param.sep[2]
        ly_mag = (period - 2 * ly_pole) / 4
        period = period + 6 * sep_y
        n_half_poles = int((n_poles - 1) / 2)
        if n_poles % 2 == 0 and n_poles > 0:
            sh_ep = ly_mag + ly_pole / 2
        else:
            sh_ep = 0
        # --- Longitudinal position
        if pole_id == 0:
            pos_y = ly_pole / 4 + half_period_id * period / 2
        elif pole_id == 1:
            pos_y = (period / 2 - ly_pole / 4) + half_period_id * period / 2
        elif pole_id == 2:
            pos_y = (period - sep_y * 6) / 2 + half_period_id * period / 2 + sep_y * 1.5 - sh_ep
        # --- Material
        mat = rad_mat.set_soft_mat(self.radia_und_param.pole_mat)
        # --- Build the block (extrusion)
        # odd n_poles
        if pole_id == 0 or pole_id == 1:
            if pole_id == 0:
                p1=[[lx_pole/4-cx_pole/2,pos_y-cy_pole/2,sep_z],[lx_pole/2-cx_pole,ly_pole/2-cy_pole]]
            else:
                p1=[[lx_pole/4-cx_pole/2,pos_y+cy_pole/2,sep_z],[lx_pole/2-cx_pole,ly_pole/2-cy_pole]]
            p2=[[lx_pole/4,pos_y,cz_pole+sep_z],[lx_pole/2,ly_pole/2]]
            p3=[[lx_pole/4,pos_y,lz_pole+sep_z],[lx_pole/2,ly_pole/2]]
        # even n_poles
        elif pole_id == 2:
            p1=[[lx_pole/4-cx_pole/2,pos_y,sep_z],[lx_pole/2-cx_pole,ly_pole-2*cy_pole]]
            p2=[[lx_pole/4,pos_y,cz_pole+sep_z],[lx_pole/2,ly_pole]]
            p3=[[lx_pole/4,pos_y,lz_pole+sep_z],[lx_pole/2,ly_pole]]
        pole=rad.ObjMltExtRtg([p1,p2,p3])
        # --- Longitudinal subdivision
        if pole_id == 0 or pole_id == 1:
            rad.ObjDivMag(pole, [self.radia_und_param.pole_long_sub[1], self.radia_und_param.pole_long_sub[1], self.radia_und_param.pole_long_sub[0]])
        else:
            rad.ObjDivMag(pole, [self.radia_und_param.pole_long_sub[1], self.radia_und_param.pole_long_sub[1] * 2, self.radia_und_param.pole_long_sub[0]])
        # --- Set the material
        rad.MatApl(pole, mat)
        # --- Set the color
        color = self.radia_und_param.pole_color
        rad.ObjDrwAtr(pole, color)
        # --- Return the block
        return pole

    def build_block_side_ext(self, mag_id=0):
        """
        Build a extremity side PM block for hybrid wiggler
        :return: pm block
        """

        period = self.radia_und_param.period
        n_poles = self.radia_und_param.n_poles
        ly_pole = self.radia_und_param.pole_length
        lx_pole = self.radia_und_param.pole_width
        ly_epole = self.radia_und_param.ext_pole[0]
        d_pole = self.radia_und_param.ext_pole[1]
        # --- Dimensions
        lx_mag = self.radia_und_param.mag_width
        lz_mag = self.radia_und_param.mag_height[1]
        lz_side = self.radia_und_param.mag_height[2]
        cx_mag = self.radia_und_param.mag_chamfer[0]
        cy_mag = self.radia_und_param.mag_chamfer[1]
        cz_mag = self.radia_und_param.mag_chamfer[2]
        cx_pole = self.radia_und_param.pole_chamfer[0]
        cy_pole = self.radia_und_param.pole_chamfer[1]
        cz_pole = self.radia_und_param.pole_chamfer[2]
        sep_x = self.radia_und_param.sep[0]
        sep_y = self.radia_und_param.sep[1]
        sep_z = self.radia_und_param.sep[2]
        d_pole = sep_y
        ly_mag = (period - 2 * ly_pole) / 4
        period = period + 6 * sep_y
        n_half_poles = int((n_poles - 1) / 2)
        if n_poles % 2 == 0 and n_poles > 0:
            sh_ep = (ly_mag + ly_pole / 2) * 1
        else:
            sh_ep = 0
        # --- Position and axis
        if mag_id == 0:   # pole at the end compensation
            pos_y = ly_pole / 4 + (n_half_poles * period / 2)
            if n_half_poles % 2 == 1:
                axis = [-1, 0, 0]
            else:
                axis = [1, 0, 0]
        elif mag_id == 1:   # extremity pole
            if n_poles % 2 == 0:
                pos_y = ly_pole / 2 + ly_mag * 2 + ly_epole / 2 + (n_half_poles * period / 2) + d_pole + sep_y * 3.5 + sh_ep
            else:
                pos_y = ly_pole / 2 + ly_mag * 2 + ly_epole / 2 + (n_half_poles * period / 2) + d_pole + sep_y * 2 + sh_ep
            if n_poles % 2 == 0:
                if int(n_poles / 2) % 2 == 1:
                    axis = [1, 0, 0]
                else:
                    axis = [-1, 0, 0]
            else:
                if n_half_poles % 2 == 1:
                    axis = [1, 0, 0]
                else:
                    axis = [-1, 0, 0]
        # --- Material
        mat = rad_mat.set_pm_mat(self.radia_und_param.mag_mat, self.radia_und_param.br)
        # --- Build the block (extrusion)
        if mag_id == 0:
            s1=[[(lx_pole-cx_pole+(lx_mag-lx_pole)/2-cx_mag)/2 + sep_x,pos_y-cy_pole/2, sep_z],[cx_pole+(lx_mag-lx_pole)/2-cx_mag,ly_pole/2-cy_pole]]
            s2=[[(lx_pole+(lx_mag-lx_pole)/2)/2 + sep_x,pos_y,cz_mag + sep_z],[(lx_mag-lx_pole)/2,ly_pole/2]]
            s3=[[(lx_pole+(lx_mag-lx_pole)/2)/2 + sep_x,pos_y,lz_mag + sep_z],[(lx_mag-lx_pole)/2,ly_pole/2]]
        elif mag_id == 1:
            s1=[[(lx_pole-cx_pole+(lx_mag-lx_pole)/2-cx_mag)/2 + sep_x,pos_y, sep_z],[cx_pole+(lx_mag-lx_pole)/2-cx_mag,ly_epole-2*cy_pole]]
            s2=[[(lx_pole+(lx_mag-lx_pole)/2)/2 + sep_x,pos_y,cz_mag + sep_z],[(lx_mag-lx_pole)/2,ly_epole]]
            s3=[[(lx_pole+(lx_mag-lx_pole)/2)/2 + sep_x,pos_y,lz_side + sep_z],[(lx_mag-lx_pole)/2,ly_epole]]
        pm=rad.ObjMltExtRtg([s1,s2,s3],axis)
        # --- Longitudinal subdivision
        rad.ObjDivMag(pm, [self.radia_und_param.mag_long_sub[0], self.radia_und_param.mag_long_sub[1], self.radia_und_param.mag_long_sub[1]])
        # --- Set the material
        rad.MatApl(pm, mat)
        # --- Set the color
        color = self.radia_und_param.mag_color[1]
        #if axis[0] > 0:
            #color = [0, 0.5, 1]
        #else:
            #color = [0, 0.8, 1]
        rad.ObjDrwAtr(pm, color)
        # --- Return the block
        return pm

    def build_pole_ext(self, pole_id=0):
        """
        Build an extremity iron pole for a hybrid undulator
        :return: the pole
        """

        period = self.radia_und_param.period
        n_poles = self.radia_und_param.n_poles
        ly_epole = self.radia_und_param.ext_pole[0]
        d_pole = self.radia_und_param.ext_pole[1]
        # --- Dimensions
        lx_pole = self.radia_und_param.pole_width
        ly_pole = self.radia_und_param.pole_length
        lz_pole = self.radia_und_param.pole_height
        cx_pole = self.radia_und_param.pole_chamfer[0]
        cy_pole = self.radia_und_param.pole_chamfer[1]
        cz_pole = self.radia_und_param.pole_chamfer[2]
        sep_x = self.radia_und_param.sep[0]
        sep_y = self.radia_und_param.sep[1]
        sep_z = self.radia_und_param.sep[2]
        d_pole = sep_y
        ly_mag = (period - 2 * ly_pole) / 4
        period = period + 6 * sep_y
        n_half_poles = int((n_poles - 1) / 2)
        if n_poles % 2 == 0 and n_poles > 0:
            sh_ep = (ly_mag + ly_pole / 2) * 1
        else:
            sh_ep = 0
        # --- Longitudinal position
        if pole_id == 0:   # pole at the end compensation
            pos_y = ly_pole / 4 + (n_half_poles * period / 2)
        elif pole_id == 1:   # extremity pole
            if n_poles % 2 == 0:
                pos_y = ly_pole / 2 + ly_mag * 2 + ly_epole / 2 + (n_half_poles * period / 2) + d_pole + sep_y * 3.5 + sh_ep
            else:
                pos_y = ly_pole / 2 + ly_mag * 2 + ly_epole / 2 + (n_half_poles * period / 2) + d_pole + sep_y * 2 + sh_ep
        # --- Material
        mat = rad_mat.set_soft_mat(self.radia_und_param.pole_mat)
        # --- Build the block (extrusion)
        if pole_id == 0:
            p1=[[lx_pole/4-cx_pole/2,pos_y-cy_pole/2, sep_z],[lx_pole/2-cx_pole,ly_pole/2-cy_pole]]
            p2=[[lx_pole/4,pos_y,cz_pole + sep_z],[lx_pole/2,ly_pole/2]]
            p3=[[lx_pole/4,pos_y,lz_pole + sep_z],[lx_pole/2,ly_pole/2]]
        elif pole_id == 1:
            p1=[[lx_pole/4-cx_pole/2,pos_y, sep_z],[lx_pole/2-cx_pole,ly_epole-2*cy_pole]]
            p2=[[lx_pole/4,pos_y,cz_pole + sep_z],[lx_pole/2,ly_epole]]
            p3=[[lx_pole/4,pos_y,lz_pole + sep_z],[lx_pole/2,ly_epole]]
        pole=rad.ObjMltExtRtg([p1,p2,p3]);
        # --- Longitudinal subdivision
        rad.ObjDivMag(pole, [self.radia_und_param.pole_long_sub[1], int(self.radia_und_param.pole_long_sub[1]), self.radia_und_param.pole_long_sub[0]])
        # --- Set the material
        rad.MatApl(pole, mat)
        # --- Set the color
        color = self.radia_und_param.pole_color
        rad.ObjDrwAtr(pole, color)
        # --- Return the block
        return pole

    def build_block_ext(self):
        """
        Build a extremity PM block for hybrid wiggler
        :return: pm block
        """

        period = self.radia_und_param.period
        n_poles = self.radia_und_param.n_poles
        ly_pole = self.radia_und_param.pole_length
        ly_epole = self.radia_und_param.ext_pole[0]
        ly_emag = self.radia_und_param.ext_mag[0]
        d_pole = self.radia_und_param.ext_pole[1]
        d_mag = self.radia_und_param.ext_mag[1]
        # --- Dimensions
        lx_mag = self.radia_und_param.mag_width
        lz_mag = self.radia_und_param.mag_height[0]
        cx_mag = self.radia_und_param.mag_chamfer[0]
        cy_mag = self.radia_und_param.mag_chamfer[1]
        cz_mag = self.radia_und_param.mag_chamfer[2]
        sep_x = self.radia_und_param.sep[0]
        sep_y = self.radia_und_param.sep[1]
        sep_z = self.radia_und_param.sep[2]
        d_pole = sep_y
        d_mag = sep_y
        ly_mag = (period - 2 * ly_pole) / 4
        period = period + 6 * sep_y
        if n_poles % 2 == 0 and n_poles > 0:
            sh_ep = (ly_mag + ly_pole / 2) * 1
        else:
            sh_ep = 0
        # --- Position and axis
        n_half_poles = int((n_poles - 1) / 2)
        if n_poles % 2 == 0:
            pos_y = ly_pole / 2 + ly_mag * 2 + ly_epole  + (n_half_poles * period / 2) + d_pole + d_mag + ly_emag / 2 + sep_y * 3.5 + sh_ep
        else:
            pos_y = ly_pole / 2 + ly_mag * 2 + ly_epole  + (n_half_poles * period / 2) + d_pole + d_mag + ly_emag / 2 + sep_y * 2 + sh_ep
        if n_poles % 2 == 0:
            if int(n_poles / 2) % 2 == 1:
                axis = [0, 1, 0]
            else:
                axis = [0, -1, 0]
        else:
            if n_half_poles % 2 == 1:
                axis = [0, 1, 0]
            else:
                axis = [0, -1, 0]
        # --- Material
        mat = rad_mat.set_pm_mat(self.radia_und_param.mag_mat, self.radia_und_param.br)
        # --- Build the block (extrusion)
        m1=[[lx_mag/4-cx_mag/2,pos_y-cy_mag/2, sep_z],[lx_mag/2-cx_mag,ly_emag+cy_mag]]
        m2=[[lx_mag/4,pos_y,cz_mag + sep_z],[lx_mag/2,ly_emag]]
        m3=[[lx_mag/4,pos_y,lz_mag + sep_z],[lx_mag/2,ly_emag]]
        pm=rad.ObjMltExtRtg([m1,m2,m3],axis)
        # --- Longitudinal subdivision
        #rad.ObjDivMag(pm, [self.radia_und_param.mag_long_sub[1], int(self.radia_und_param.mag_long_sub[0]), self.radia_und_param.mag_long_sub[1]])
        rad.ObjDivMag(pm, [self.radia_und_param.mag_long_sub[1], self.radia_und_param.mag_long_sub[0], self.radia_und_param.mag_long_sub[1]])
        # --- Set the material
        rad.MatApl(pm, mat)
        # --- Set the color
        color = self.radia_und_param.mag_color[0]
        #if axis[1] > 0:
            #color = [0, 1, 1]
        #else:
            #color = [0, 1, 0.8]
        rad.ObjDrwAtr(pm, color)
        # --- Return the block
        return pm

    def build_undulator(self, sym=True):
        """
        Build an hybrid undulator
        :param sym: apply the symmetries (default = True). Build 1/8 undulator if False.
        :return: the undulator
        """

        # --- Hold the symmetry
        self.sym = sym
        # global quantities for force calc (half of pole array)
        global und_frc
        # --- Build an empty container
        und = rad.ObjCnt([])
        und_frc = rad.ObjCnt([])
        n_poles = self.radia_und_param.n_poles
        # --- Build all the periods
        if n_poles < 2:
            mag_0 = self.build_block_main(0, 0)
            #mag_side_0 = self.build_block_side(0, 0)
            pole_0 = self.build_pole(0, 0)
            rad.ObjAddToCnt(und, [mag_0, pole_0]) # full
            rad.ObjAddToCnt(und_frc, [pole_0]) # und_frc
            if n_poles == 1:
                # --- Undulator extremity
                mag_main_ext_1 = self.build_block_main(1, 0) # Last 2 standard magnet
                #mag_side_ext_1 = self.build_block_side_ext(1) # Extremity side magnet block
                pole_ext_1 = self.build_pole_ext(1) # Extremity pole
                mag_ext = self.build_block_ext() # Extremity magnet
                if self.radia_und_param.wig_build[:4] == 'full':
                    rad.ObjAddToCnt(und, [mag_main_ext_1, pole_ext_1, mag_ext]) # full
                    #rad.ObjAddToCnt(und_frc, [pole_ext_1]) # und_frc
                elif self.radia_und_param.wig_build[:9] == 'side_pole':
                    rad.ObjAddToCnt(und, [pole_ext_1]) # side & pole
                elif self.radia_und_param.wig_build[:9] == 'main_pole':
                    rad.ObjAddToCnt(und, [mag_main_ext_1, mag_ext, pole_ext_1]) # main & pole
                elif self.radia_und_param.wig_build[:9] == 'main_side':
                    rad.ObjAddToCnt(und, [mag_main_ext_1, mag_ext]) # main & side
                elif self.radia_und_param.wig_build[:4] == 'main':
                    rad.ObjAddToCnt(und, [mag_main_ext_1, mag_ext]) # main
                elif self.radia_und_param.wig_build[:4] == 'side':
                    rad.ObjAddToCnt(und, []) # side
                elif self.radia_und_param.wig_build[:4] == 'pole':
                    rad.ObjAddToCnt(und, [pole_ext_1]) # pole

        elif n_poles % 2 == 0:
            n_half_poles = int(n_poles / 2)
            for k in range(n_half_poles):
                #print('even', k)
                mag_0 = self.build_block_main(1, k)
                mag_1 = self.build_block_main(2, k)
                #mag_side_0 = self.build_block_side(2, k)
                pole_0 = self.build_pole(2, k)
                # Add to the container
                if self.radia_und_param.wig_build[:4] == 'full':
                    rad.ObjAddToCnt(und, [mag_0, mag_1, pole_0]) # full
                    #rad.ObjAddToCnt(und_frc, [pole_0]) # und_frc
                elif self.radia_und_param.wig_build[:9] == 'side_pole':
                    rad.ObjAddToCnt(und, [pole_0]) # side & pole
                elif self.radia_und_param.wig_build[:9] == 'main_pole':
                    rad.ObjAddToCnt(und, [mag_0, mag_1, pole_0]) # main & pole
                elif self.radia_und_param.wig_build[:9] == 'main_side':
                    rad.ObjAddToCnt(und, [mag_0, mag_1]) # main & side
                elif self.radia_und_param.wig_build[:4] == 'main':
                    rad.ObjAddToCnt(und, [mag_0, mag_1]) # main
                elif self.radia_und_param.wig_build[:4] == 'side':
                    rad.ObjAddToCnt(und, []) # side
                elif self.radia_und_param.wig_build[:4] == 'pole':
                    rad.ObjAddToCnt(und, [pole_0]) # pole

            # --- Undulator extremity
            mag_main_ext_1 = self.build_block_main(1, n_half_poles) # Last 2 standard magnet
            #mag_side_ext_1 = self.build_block_side_ext(1) # Extremity side magnet block
            pole_ext_1 = self.build_pole_ext(1) # Extremity pole
            mag_ext = self.build_block_ext() # Extremity magnet
            #mag_ext_cla = self.build_block_ext_clamp() # Extremity magnet clamp

            if self.radia_und_param.wig_build[:4] == 'full':
                #pass
                rad.ObjAddToCnt(und, [mag_main_ext_1, pole_ext_1, mag_ext]) # full
                #rad.ObjAddToCnt(und_frc, [pole_ext_1]) # und_frc
            elif self.radia_und_param.wig_build[:9] == 'side_pole':
                rad.ObjAddToCnt(und, [pole_ext_1]) # side & pole
            elif self.radia_und_param.wig_build[:9] == 'main_pole':
                rad.ObjAddToCnt(und, [mag_main_ext_1, mag_ext, pole_ext_1]) # main & pole
            elif self.radia_und_param.wig_build[:9] == 'main_side':
                rad.ObjAddToCnt(und, [mag_main_ext_1, mag_ext]) # main & side
            elif self.radia_und_param.wig_build[:4] == 'main':
                rad.ObjAddToCnt(und, [mag_main_ext_1, mag_ext]) # main
            elif self.radia_und_param.wig_build[:4] == 'side':
                rad.ObjAddToCnt(und, []) # side
            elif self.radia_und_param.wig_build[:4] == 'pole':
                rad.ObjAddToCnt(und, [pole_ext_1]) # pole

        else:
            n_half_poles = int((n_poles - 1) / 2)
            for k in range(n_half_poles):
                #print('odd',k)
                mag_0 = self.build_block_main(0, k)
                mag_1 = self.build_block_main(1, k)
                #mag_side_0 = self.build_block_side(0, k)
                #mag_side_1 = self.build_block_side(1, k)
                pole_0 = self.build_pole(0, k)
                pole_1 = self.build_pole(1, k)
                # Add to the container
                if self.radia_und_param.wig_build[:4] == 'full':
                    rad.ObjAddToCnt(und, [mag_0, mag_1, pole_0, pole_1]) # full
                    #rad.ObjAddToCnt(und_frc, [pole_0, pole_1]) # und_frc
                elif self.radia_und_param.wig_build[:9] == 'side_pole':
                    rad.ObjAddToCnt(und, [pole_0, pole_1]) # side & pole
                elif self.radia_und_param.wig_build[:9] == 'main_pole':
                    rad.ObjAddToCnt(und, [mag_0, mag_1, pole_0, pole_1]) # main & pole
                elif self.radia_und_param.wig_build[:9] == 'main_side':
                    rad.ObjAddToCnt(und, [mag_0, mag_1]) # main & side
                elif self.radia_und_param.wig_build[:4] == 'main':
                    rad.ObjAddToCnt(und, [mag_0, mag_1]) # main
                elif self.radia_und_param.wig_build[:4] == 'side':
                    rad.ObjAddToCnt(und, []) # side
                elif self.radia_und_param.wig_build[:4] == 'pole':
                    rad.ObjAddToCnt(und, [pole_0, pole_1]) # pole

            # --- Undulator extremity
            mag_main_ext_0 = self.build_block_main(0, n_half_poles) # Last 1 standard magnet
            mag_main_ext_1 = self.build_block_main(1, n_half_poles) # Last 2 standard magnet
            #mag_side_ext_0 = self.build_block_side_ext(0) # End side magnet block
            #mag_side_ext_1 = self.build_block_side_ext(1) # Extremity side magnet block
            pole_ext_0 = self.build_pole_ext(0) # End pole
            pole_ext_1 = self.build_pole_ext(1) # Extremity pole
            mag_ext = self.build_block_ext() # Extremity magnet
            #mag_ext_cla = self.build_block_ext_clamp() # Extremity magnet clamp

            if self.radia_und_param.wig_build[:4] == 'full':
                rad.ObjAddToCnt(und, [mag_main_ext_0, mag_main_ext_1, pole_ext_0, pole_ext_1, mag_ext]) # full
                #rad.ObjAddToCnt(und_frc, [pole_ext_0, pole_ext_1]) # und_frc
            elif self.radia_und_param.wig_build[:9] == 'side_pole':
                rad.ObjAddToCnt(und, [pole_ext_0, pole_ext_1]) # side & pole
            elif self.radia_und_param.wig_build[:9] == 'main_pole':
                rad.ObjAddToCnt(und, [mag_main_ext_0, mag_main_ext_1, mag_ext, pole_ext_0, pole_ext_1]) # main & pole
            elif self.radia_und_param.wig_build[:9] == 'main_side':
                rad.ObjAddToCnt(und, [mag_main_ext_0, mag_main_ext_1, mag_ext]) # main & side
            elif self.radia_und_param.wig_build[:4] == 'main':
                rad.ObjAddToCnt(und, [mag_main_ext_0, mag_main_ext_1, mag_ext]) # main
            elif self.radia_und_param.wig_build[:4] == 'side':
                rad.ObjAddToCnt(und, []) # side
            elif self.radia_und_param.wig_build[:4] == 'pole':
                rad.ObjAddToCnt(und, [pole_ext_0, pole_ext_1]) # pole

        # --- Move to the specified gap
        if n_poles >= 0:
            tr = rad.TrfTrsl([0, 0, self.radia_und_param.gap / 2])
            rad.TrfOrnt(und, tr)
            rad.TrfOrnt(und_frc, tr)

        # --- Symmetries
        if sym:
            rad.TrfZerPerp(und, [0, 0, 0], [1, 0, 0])
            rad.TrfZerPerp(und_frc, [0, 0, 0], [1, 0, 0])
            if n_poles > 0:
                if self.radia_und_param.wig_build[len(self.radia_und_param.wig_build)-4:] == 'half':
                    pass
                else:
                    rad.TrfZerPara(und, [0, 0, 0], [0, 0, 1])
            if n_poles % 2 == 0 and n_poles > 0:
                rad.TrfZerPara(und, [0, 0, 0], [0, 1, 0])
                rad.TrfZerPara(und_frc, [0, 0, 0], [0, 1, 0])
            else:
                rad.TrfZerPerp(und, [0, 0, 0], [0, 1, 0])
                rad.TrfZerPerp(und_frc, [0, 0, 0], [0, 1, 0])

        # --- Cut the half for force2
        und_frc = rad.ObjCutMag(und, [0,0,0], [0,0,1], 'Frame->Lab')[0]

        # --- Return the undulator
        return und

    def print_data(self, short=False):
        """
        Print undulator data
        :param short=False: print a short summary if True
        """
        sep_y = self.radia_und_param.sep[1]
        period = self.radia_und_param.period + 6 * sep_y
        n_half_poles = int((self.radia_und_param.n_poles - 1) / 2)
        p_len = self.radia_und_param.pole_length
        ly_mag = (period - 2 * p_len) / 4
        pole_ext_length = self.radia_und_param.ext_pole[0]
        mag_ext_length = self.radia_und_param.ext_mag[0]
        d_pole = self.radia_und_param.ext_pole[1]
        d_mag = self.radia_und_param.ext_mag[1]
        d_pole = sep_y
        d_mag = sep_y
        total_length = period * n_half_poles + (p_len / 2+ ly_mag * 2 + pole_ext_length + d_pole + d_mag + mag_ext_length) * 2
        print('Period: ', self.radia_und_param.period, ' mm')
        print('Number of poles: ', self.radia_und_param.n_poles)
        print('Peak field: ', self.peak_field, ' T')
        print('K: ', self.k)
        print('Total length: ', self.total_length, ' mm')
        if not short:
            print('PM Material: ', self.radia_und_param.mag_mat)
            print('Br: ', self.radia_und_param.br, ' T')
            print('PM block length: ', ly_mag , ' mm')
            print('Pole material: ', self.radia_und_param.pole_mat)
            print('Gap: ', self.radia_und_param.gap, ' mm')
