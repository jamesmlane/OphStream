# ----------------------------------------------------------------------------
#
# TITLE - io.py
# AUTHOR - James Lane
# PROJECT - Ophiuchus Stream Honours
# CONTENTS:
#	1. ReadSnapFile
#   2. SnaptoBinary
#   3.
#
# ----------------------------------------------------------------------------
#
# Docstrings and metadata:
'''Input/Output module for general GADGET-2 snapshots'''
__author__ = "James Lane"

#Imports
from snapdata import Snapdata
import numpy as np
import os
import csv
from astropy.io import fits

def ReadSnapFile(filename):
    '''
    ReadSnapFile:

    Reads a snapshot and return the data.

    Args:
        filename (string) - Filename

    Returns:
        data (array) - Snapshot data
    '''
    abs_path = os.path.abspath(filename)
    return Snapdata(abs_path)
#def

def SnapToBinary(files,fileout,velocities):
    '''
    SnapToBinary:

    Takes a series of snapshot filenames, reads them in and transforms them
    into numpy binary format.

    Args:
        filename (string array) - Array of filenames, sorted.
        fileout (string) - Output filename
        velocities (boolean) - Include velocities or not

    Returns:
        None
    '''
    ns = len(files)
    npt = ReadSnapFile(files[0]).npart
    x = np.empty((ns,npt))
    y = np.empty((ns,npt))
    z = np.empty((ns,npt))
    vx = np.empty((ns,npt))
    vy = np.empty((ns,npt))
    vz = np.empty((ns,npt))
    r = np.empty((ns,npt))
    p = np.empty((ns,npt))
    t = np.empty(ns)
    p_mass = 0
    count = 0
    for i in range(ns):
        data = ReadSnapFile(files[i])
        x[i] = data.x
        y[i] = data.y
        z[i] = data.z
        vx[i] = data.vx
        vy[i] = data.vy
        vz[i] = data.vz
        p[i] = data.p
        t[i] = data.time
        if p_mass == 0: p_mass = data.p_mass
        print('Done snap '+str(count))
        count += 1
    ###i
    if velocities == True: kinematics = np.array([x,y,z,vx,vy,vz,p])
    else: kinematics = np.array([x,y,z,p])

    if fileout[-4:] != '.npy': fileout = fileout+'.npy'
    np.save(fileout, np.array([kinematics, t, p_mass]))
#def

def ReadHarrisGCs(filename):
    '''
    ReadHarrisGCs:

    Read in the Harris Milky Way Globular Cluster catalog. Extract absolute
    magnitude [mags], 2D half-mass radius [pc], distance [kpc], velocity
    dispersion [km/s].

    Args:
        filename (array) - Filename of the Harris globular cluster catalog.

    Returns:
        output (array) - ( Mv [mags], Rh [pc], dist [kpc], sig [km/s] )
    '''

    #Get absolute V-band magnitude and half-light radius
    data = fits.getdata(filename,0)
    Mv = data['M_V,t']
    Rh_pc = data['r_h'] * 60 * data['R_Sun'] * 1000 / 206264.8
    sigma_v = data['sigma_v']
    sigma_err = data['sigma_v_err']
    R_gc = data['R_gc']
    name = data['ID']
    return Mv,Rh_pc,R_gc,sigma_v,sigma_err,name
#def

def ReadHuxorGCs(filename):
    '''
    ReadHuxorGCs:

    Read in the Huxor M31 Globular Cluster catalog. Extract absolute magnitude
    [mags], 2D half-mass radius [pc]

    Args:
        filename (array) - Filename of the Huxor globular cluster catalog.

    Returns:
        output (array) - ( Mv [mags], Rh [pc] )
    '''

    #Get absolute magnitudes and half-light radii
    data = np.genfromtxt(filename, dtype=[('Mv','f8'),('Rh','f8')],
                            usecols=[6,8], delimiter='\t', skip_header=2,
                            missing_values='-')
    return np.array([-data['Mv'],data['Rh']])
#def

def ReadPeacockGCs(filename):
    '''
    ReadPeacockGCs:

    Read in the Peacock M31 Globular Cluster catalog. Extract absolute magnitude
    [mags], 2D half-mass radius [pc]

    Args:
        filename (array) - Filename of the Peacock globular cluster catalog.

    Returns:
        output (array) - ( Mv [mags], Rh [pc] )
    '''

    gmag = np.array([])
    grcol = np.array([])
    rhalf = np.array([])
    with open(filename) as f:
        reader = csv.reader(f)
        i = 0
        for row in reader:
            if i < 18:
                i+=1
                continue
            ##fi
            #Only choose good globular clusters and those that have g, g-r, rh
            if int(row[3]) == 1 and row[5] != '' and row[7] != '' and row[12] != '':
                gmag = np.append(gmag,float(row[5]))
                grcol = np.append(grcol, float(row[7]))
                rhalf = np.append(rhalf, float(row[12]))
            ##fi
            i+=1
        ###i
    #wth
    app_vmag = gmag-0.59*grcol-0.01
    abs_vmag = app_vmag - 5*np.log10(780000/10)
    return np.array([abs_vmag,rhalf])
#def
