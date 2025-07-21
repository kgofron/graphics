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
import subprocess

def intensities2detint(cif_path, hkl_path, wavelength, min_intensity, R, cyl_center, ray_origin, zmin, zmax, gamma_axis, delta_axis):
    lst = []
    #generate hkl file with given cif file, wavelength
    cif2hkl_bin = '/usr/bin/cif2hkl'
    cmd = [cif2hkl_bin, '--mode', 'NUC', '--out', hkl_path, '--lambda', str(wavelength), '--xtal', cif_path]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = proc.communicate()

    # go from hkl file output by cif2hkl to a dataframe of reflections/intensities
    latt, df = hkl2dfhkl(hkl_path)
    print(latt)

    a = latt['a']
    b = latt['b']
    c = latt['c']
    alpha = latt['alpha']
    beta = latt['beta']
    gamma = latt['gamma']

    user = Hkl.UnitEnum.USER
    detector = Hkl.Detector.factory_new(Hkl.DetectorType(0))
    factory  = Hkl.factories()['E6C']
    geometry = factory.create_new_geometry()
    geometry.wavelength_set(wavelength, Hkl.UnitEnum.USER)
    sample = Hkl.Sample.new("toto") # sample. tab to check attributes

    alpha = math.radians(alpha)
    beta  = math.radians(beta)
    gamma = math.radians(gamma)
    lattice = Hkl.Lattice.new(a,b,c,alpha,beta,gamma)
    sample.lattice_set(lattice)

    # add columns for real axes motor positions to reflection df
    df2 = dfhkl2dfhklaxes(df, min_intensity, factory, geometry, detector, sample, user)
    theta, z, intensities = [], [], []
    #df2.to_csv('refls2.csv')
    for idx, refl in df2.iterrows():
        mu = refl['mu']
        omega = refl['omega']
        chi = refl['chi']
        phi = refl['phi']
        gamma = refl['gamma']
        delta = refl['delta']
        h = refl['h']
        k = refl['k']
        l = refl['l']
        inten = refl['intensity']
        detthetaz = real2det(gamma_axis, delta_axis, gamma, delta, R, \
            cyl_center, ray_origin)
        if (detthetaz is not None):
            theta = float(detthetaz[0])
            z = float(detthetaz[1])
            if (z<zmax) and (z>zmin):
                lst.append((theta, z, inten, h, k, l, mu, omega, chi, phi, gamma, delta))
    if lst is not None:
        return lst
    else:
        print("no points found")
        return None
