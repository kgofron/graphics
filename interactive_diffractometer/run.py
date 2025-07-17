# in shell: source /epics/iocs/ioc-hkl/iochkl/bin/activate
# export GI_TYPELIB_PATH=/usr/local/lib/girepository-1.0
import numpy as np
import pandas as pd
import math
import gi
from gi.repository import GLib
gi.require_version('Hkl', '5.0')
from gi.repository import Hkl
from hkl2dfhkl import hkl2dfhkl
from dfhkl2dfhklaxes import dfhkl2dfhklaxes
from detectorposcalc import real2det
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

'''
hkl intensities plotted on a full cylindrical detector with 4-circle geometry
'''

hkl_path = '/epics/crystals/SiO2/EntryWithCollCode176.hkl'
min_intensity = 50

R = 10
zmin = -4
zmax = 4
cyl_center = (4,0)
ray_origin = (0,0,0)
if __name__=="__main__":
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    user = Hkl.UnitEnum.USER
    detector = Hkl.Detector.factory_new(Hkl.DetectorType(0))
    factory  = Hkl.factories()['E4CV']
    geometry = factory.create_new_geometry()
    geometry.wavelength_set(0.5, Hkl.UnitEnum.USER)
    sample = Hkl.Sample.new("toto")
    lattice_list = [5.41,5.41,5.41,90.,90.,90.]
    a,b,c,alpha,beta,gamma=lattice_list

    alpha = math.radians(alpha)
    beta  = math.radians(beta)
    gamma = math.radians(gamma)
    lattice = Hkl.Lattice.new(a,b,c,alpha,beta,gamma)
    sample.lattice_set(lattice)

    # go from hkl file output by cif2hkl to a dataframe of reflections/intensities
    df = hkl2dfhkl(hkl_path)
    # add columns for real axes motor positions to reflection df
    df2 = dfhkl2dfhklaxes(df, factory, geometry, detector, sample, user)

    xs, ys, zs, intensities = [], [], [], []

    for idx, refl in df2.iterrows():
        o = refl['omega']
        c = refl['chi']
        p = refl['phi']
        t = refl['tth']
        inten = refl['intensity']
        detxyz = real2det(o, c, p, t, R, cyl_center, ray_origin)
        if (detxyz is not None) and (inten>min_intensity):
            x, y, z = detxyz
            if (zmin<=z<=zmax):
                print(f'[{o}, {c}, {p}, {t}], {detxyz}')
                xs.append(x)
                ys.append(y)
                zs.append(z)
                intensities.append(inten)

    xs = np.array(xs)
    ys = np.array(ys)
    zs = np.array(zs)
    intensities = np.array(intensities)

    # plot all points, colored by intensity
    sc = ax.scatter(xs, ys, zs, c=intensities, cmap='viridis', marker='o')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    # colorbar
    cb = plt.colorbar(sc, ax=ax, shrink=0.5, aspect=10)
    cb.set_label('Intensity')

    plt.savefig('test.png')
