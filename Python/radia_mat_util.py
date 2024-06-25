# ----------------------------------------------------------
# PyRadiaUndulators library
# radia_mat_util.py
# Useful materials for Radia models
#
# Gael Le Bec, ESRF, 2019
# ----------------------------------------------------------

import radia as rad
from math import pi

# -------------------------------------------
# Materials
# -------------------------------------------
# XC6 / AISI 1006 low carbon steel
# xc6 = rad.MatSatIsoFrm([2118, 1.362], [63.06, 0.2605], [17.138, 0.4917])
# AFK1 FeCo Alloy from MetalImphy (Fe : 74.2, Co: 25%, Cr: 0.3%, Mn: 0.5%)
# fe_co = rad.MatSatIsoFrm([2001., 1.704], [38.56, 0.493], [1.24, 0.152])
# AFK502 Vanadium Permendur alloy from MetalImphy (Fe : 49%, Co: 49%, V: 2%) similar to Vacoflux50 from VacuumSchmelze
# fe_co_v = rad.MatSatIsoFrm([10485., 1.788], [241.5, 0.437], [7.43, 0.115])
# Armco pure iron
# armco = rad.MatSatIsoFrm([4236., 1.385], [126.2, 0.2649], [34., 0.5])

def set_soft_mat(mat_name):
    """"
    Returns a soft magnetic material according to mat_name
    :param mat_name: material specification
            -- 'xc6': XC6 / AISI 1006 low carbon steel
            -- 'feco': FeCo Alloy
            -- 'fecov': Vanadium Permendur alloy
            -- 'armco': Armco Pure Iron
            -- '1300-100a': Isovac 1300-100A electric steel from VostAlpine
            -- 'ni': Pure nickel
    :return the material
    """
    mu_0 = 4. * pi * 1e-07
    if mat_name == 'xc6':
        mat = rad.MatSatIsoFrm([2118, 1.362], [63.06, 0.2605], [17.138, 0.4917])
    elif mat_name == 'feco':
        mat = rad.MatSatIsoFrm([2001., 1.704], [38.56, 0.493], [1.24, 0.152])
    elif mat_name == 'fecov':
        mat = rad.MatSatIsoFrm([10485., 1.788], [241.5, 0.437], [7.43, 0.115])
    elif mat_name == 'armco':
        mat = rad.MatSatIsoFrm([4236., 1.385], [126.2, 0.2649], [34., 0.5])
    elif mat_name == '1300-100a':
        h = [48, 56, 65, 74, 82, 91, 100, 108, 117, 125, 134, 143, 152, 161, 171, 182, 193, 217, 230, 243, 255, 271,
             294, 318, 344, 401, 517, 698, 983, 1559, 2557, 3930, 7500, 10000, 25000, 40000, 100000]
        j = [0.100, 0.150, 0.200, 0.250, 0.300, 0.350, 0.400, 0.450, 0.500, 0.550, 0.600, 0.650, 0.700, 0.750, 0.800,
             0.850, 0.900, 1.000, 1.050, 1.100, 1.150, 1.200, 1.250, 1.300, 1.350, 1.400, 1.450, 1.500, 1.550, 1.600,
             1.650, 1.700, 1.799, 1.847, 2.000, 2.050, 2.060]
        # Note: Points (25000, 2000) and higher fields extrapolated
        mat = rad.MatSatIsoTab([[mu_0 * h[i], 1. * j[i]] for i in range(len(h))])
    elif mat_name == 'ni':
        h = [0, 4.66354, 8.69715, 15.4198, 26.5123, 45.3358, 69.5375, 113.907, 156.932, 250, 500, 10000] 
        j = [0, 0.17756, 0.28847, 0.38594, 0.44979, 0.49683, 0.53714, 0.56733, 0.58745, 0.59568, 0.60937, 0.62]
        # Note: Points (250, ) and higher fields extrapolated
        mat = rad.MatSatIsoTab([[mu_0 * h[i], 1. * j[i]] for i in range(len(h))])
    else:
        print('Unknown material set to xc6')
        mat = rad.MatSatIsoFrm([2118, 1.362], [63.06, 0.2605], [17.138, 0.4917])
    return mat

def set_pm_mat(mat_name, br=None):
    """"
    Returns a hard PM material
    :param mat_name: material type
            -- 'ndfeb': NdFeB material
            -- 'prfeb80K': Cryo PrFeB material
            -- 'sm2co17': Sm2Co17 material
            -- 'smco5': SmCo5 material
            -- 'ferrite': ferrite material
    :param br=None: remanent field
    :return the material
    """
    if mat_name == 'ndfeb':
        if br is None:
            br = 1.3
        mat = rad.MatStd('NdFeB', br)
    elif mat_name == 'prfeb80K':
        if br is None:
            br = 1.6
    elif mat_name == 'sm2co17':
        if br is None:
            br = 1.1
        mat = rad.MatStd('Sm2Co17', br)
    elif mat_name == 'smco5':
        if br is None:
            br = 0.85
        mat = rad.MatStd('SmCo5', br)
    elif mat_name == 'ferrite':
        if br is None:
            br = 0.35
        mat = rad.MatStd('Ferrite', br)
    else:
        mat = None

    return mat