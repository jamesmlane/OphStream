# ----------------------------------------------------------------------------
#
# TITLE - units.py
# AUTHOR - James Lane
# PROJECT - Ophiuchus Stream Honours
# CONTENTS:
#   1. Phys2Code
#   2. Code2Phys
#	3. GetConversions
#
# ----------------------------------------------------------------------------
#
# Docstrings and metadata:
'''Unit handling module for general GADGET-2 snapshots'''
__author__ = "James Lane"

#Imports
import numpy as np
import pdb

def Phys2Code(coords,unitsys):
    '''
    Converts 3D position and velocity information from physical to code units.

    Args:
        coords (array) -    6 x N array of coordinates to convert
                            (3 position [kpc], 3 velocity [km/s])
        unitsys (string) - Simulation number based coordinate system

    Returns:
        coords_code (array) -   6 x N array of coordinates in code units
    '''
    _,code2kpc,code2kms,_ = GetConversions(unitsys)
    coords_code = np.array(coords)
    coords_code[:3] /= code2kpc
    coords_code[3:] /= code2kms
    return coords_code
#def

def Code2Phys(coords,unitsys):
    '''
    Converts 3D position and velocity information from code to physical units.

    Args:
        coords (array) -    6 x N array of coordinates to convert
                            (3 position [kpc], 3 velocity [km/s])
        unitsys (string) - Simulation number based coordinate system

    Returns:
        coords_phys (array) -   6 x N array of coordinates in physical units
    '''
    _,code2kpc,code2kms,_ = GetConversions(unitsys)
    coords_phys = np.array(coords)
    coords_phys[:3] *= code2kpc
    coords_phys[3:] *= code2kms
    return coords_phys
#def

def GetConversions(unitsys):
    '''
    Gets the factors to convert code to physical values.

    Args:
        unitsys (string) -  The string denoting the unit system, based on the
                            simulation number --> ie: 'GC##'

    Returns:
        conversion_factors (array) - The conversion factors code2Myr (time),
                                    code2kpc (Distance), code2kms (velocity),
                                    code2Msol (mass)
    '''
    # Deprecated
    # if unitsys == 'GC2':
    #     code2Myr = 0.1643
    #     code2kpc = 0.0023
    #     code2kms = 13.67
    #     code2Msol = 1E5
    # elif unitsys == 'GC3':
    #     code2Myr = 0.5195
    #     code2kpc = 0.0023
    #     code2kms = 4.32
    #     code2Msol = 1E4
    # elif unitsys == 'GC5':
    #     code2Myr = 0.8143
    #     code2kpc = 0.0023
    #     code2kms = 3.72
    #     code2Msol = 1E4
    # elif unitsys == 'GC6':
    #     code2Myr = 1.1052
    #     code2kpc = 0.0038
    #     code2kms = 3.36
    #     code2Msol = 1E4
    # elif unitsys == 'GC7':
    #     code2Myr = 1.4719
    #     code2kpc = 0.0046
    #     code2kms = 3.06
    #     code2Msol = 1E4
    # elif unitsys == 'GC8':
    # 	code2Myr = 1.6457
    # 	code2kpc = 0.0023
    # 	code2kms = 1.37
    # 	code2Msol = 1E3
    if unitsys == 'GC9':
        rh = 36
        code2Msol = 2E4
        code2pc = rh/1.3
        code2kms = np.sqrt( 0.004301 * code2Msol / code2pc )
        code2Myr = 0.978 * code2pc / code2kms
        code2kpc = code2pc / 1000
    elif unitsys == 'GC10':
        rh = 18
        code2Msol = 1E4
        code2pc = rh/1.3
        code2kms = np.sqrt( 0.004301 * code2Msol / code2pc )
        code2Myr = 0.978 * code2pc / code2kms
        code2kpc = code2pc / 1000
    elif unitsys == 'GC11':
        rh = 18
        code2Msol = 5E3
        code2pc = rh/1.3
        code2kms = np.sqrt( 0.004301 * code2Msol / code2pc )
        code2Myr = 0.978 * code2pc / code2kms
        code2kpc = code2pc / 1000
    elif unitsys == 'GC12':
        rh = 72
        code2Msol = 2E4
        code2pc = rh/1.3
        code2kms = np.sqrt( 0.004301 * code2Msol / code2pc )
        code2Myr = 0.978 * code2pc / code2kms
        code2kpc = code2pc / 1000
    elif unitsys == 'GC13':
        rh = 18
        code2Msol = 2.5E3
        code2pc = rh/1.3
        code2kms = np.sqrt( 0.004301 * code2Msol / code2pc )
        code2Myr = 0.978 * code2pc / code2kms
        code2kpc = code2pc / 1000
    elif unitsys == 'GC14':
        rh = 9
        code2Msol = 2.5E3
        code2pc = rh/1.3
        code2kms = np.sqrt( 0.004301 * code2Msol / code2pc )
        code2Myr = 0.978 * code2pc / code2kms
        code2kpc = code2pc / 1000
    elif unitsys == 'GC15':
        rh = 36
        code2Msol = 1E4
        code2pc = rh/1.3
        code2kms = np.sqrt( 0.004301 * code2Msol / code2pc )
        code2Myr = 0.978 * code2pc / code2kms
        code2kpc = code2pc / 1000
    elif unitsys == 'GC16':
        rh = 14
        code2Msol = 2.0E3
        code2pc = rh/1.3
        code2kms = np.sqrt( 0.004301 * code2Msol / code2pc )
        code2Myr = 0.978 * code2pc / code2kms
        code2kpc = code2pc / 1000
    elif unitsys == 'GC20':
        rh = 90 #3D
        code2Msol = 2.0E4
        code2pc = rh/1.3
        code2kms = np.sqrt( 0.004301 * code2Msol / code2pc )
        code2Myr = 0.978 * code2pc / code2kms
        code2kpc = code2pc / 1000
    elif unitsys == 'GC21':
        rh = 28 #3D
        code2Msol = 1.0E4
        code2pc = rh/1.3
        code2kms = np.sqrt( 0.004301 * code2Msol / code2pc )
        code2Myr = 0.978 * code2pc / code2kms
        code2kpc = code2pc / 1000
    elif unitsys == 'GC22':
        rh = 90 #3D
        code2Msol = 2.0E3
        code2pc = rh/1.3
        code2kms = np.sqrt( 0.004301 * code2Msol / code2pc )
        code2Myr = 0.978 * code2pc / code2kms
        code2kpc = code2pc / 1000
    else: return None
    return [code2Myr,code2kpc,code2kms,code2Msol]
