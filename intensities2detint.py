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

def intensities2detint(hkl_path, min_intensity, R, cyl_center, ray_origin, xmin,ymin,ymax,zmin,zmax):
    lst = []
    user = Hkl.UnitEnum.USER
    detector = Hkl.Detector.factory_new(Hkl.DetectorType(0))
    factory  = Hkl.factories()['E4CV']
    geometry = factory.create_new_geometry()
    geometry.wavelength_set(0.5, Hkl.UnitEnum.USER)
    sample = Hkl.Sample.new("toto") # sample. tab to check attributes

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
            x = float(detxyz[0])
            y = float(detxyz[1])
            z = float(detxyz[2])
            if (x>xmin) and (y<ymax) and (y>ymin) and (z<zmax) and (z>zmin):
                lst.append((x, y, z, inten))
    if lst is not None:
        return lst
    else:
        print("no points found")
        return None
