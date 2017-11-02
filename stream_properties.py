# ----------------------------------------------------------------------------
#
# TITLE - ophiuchus_properties
# AUTHOR - James Lane
# PROJECT - ophstream
# CONTENTS:
#	1.
#
# ----------------------------------------------------------------------------
#
# Docstrings and metadata:
'''Returns properties of the Ophiuchus Stream for easy access in Python'''
__author__ = "James Lane"

#Imports
import numpy as np

def GetGalCenXYZ():
    '''
    GetGalCenXYZ:

    Returns the galactocentric rectangular positions of the Ophiuchus Stream
    members as measured by Sesar et al. (2015). The coordinates are defined
    such that X is positive towards galactic center, Y is positive towards
    direction of galactic rotation, and z is positive towards GNP.

    Args: None

    Returns:
        x (array)
        y (array)
        z (array)

    '''
    x = np.array([  -0.8389508937,
                    -0.8046981415,
                    -1.183070456,
                    -1.093302254,
                    -1.272679717,
                    -1.599882823,
                    -0.8410581994,
                    -1.094970034,
                    -1.09207824,
                    -1.207962562,
                    -1.219576929,
                    -1.126949616,
                    -1.158278112,
                    -1.215369165])

    y = np.array([  0.5538807862,
                    0.5520129534,
                    0.5990711539,
                    0.588207452,
                    0.6124682588,
                    0.6489713382,
                    0.5568587692,
                    0.5913880388,
                    0.5885220972,
                    0.6053546714,
                    0.6063691437,
                    0.589048907,
                    0.595121864,
                    0.6049068491])

    z = np.array([  4.641615491,
                    4.663744668,
                    4.364544594,
                    4.442606605,
                    4.288764042,
                    4.03499888,
                    4.644662738,
                    4.44490282,
                    4.44056634,
                    4.330518409,
                    4.349458163,
                    4.422719619,
                    4.398678464,
                    4.342763314])

    return x, y, z
#def

def GetVRad():
    '''
    GetVRad:

    Returns the radial velocity of the Ophiuchus Stream members as measured by
    Sesar et al. (2015).

    Args: None

    Returns:
        vrad (array)

    '''
    vrad = np.array([   286.7,
                        285.3,
                        290,
                        291.3,
                        290.3,
                        289.8,
                        286,
                        286.7,
                        287.5,
                        288.8,
                        288,
                        289.4,
                        291.8,
                        286.4])
    return vrad
#def

def Getlb():
    '''
    Getlb:

    Returns the galactic longitude and latitude of Ophiuchus Stream members as
    measured by Sesar et al. (2015).

    Args: None

    Returns:
        l (array)
        b (array)
    '''
    l = np.array([  4.2542172,
                    4.2206215,
                    4.8211255,
                    4.6754368,
                    4.9909027,
                    5.5431722,
                    4.2782209,
                    4.7016991,
                    4.677135,
                    4.8884178,
                    4.9046188,
                    4.7039609,
                    4.7729883,
                    4.8899332])

    b = np.array([  31.852146,
                    31.858472,
                    31.45703,
                    31.59739,
                    31.324753,
                    30.955722,
                    31.875476,
                    31.615576,
                    31.581219,
                    31.343475,
                    31.496055,
                    31.600574,
                    31.569801,
                    31.442099])

    return l, b
#def

def GetFannedlb():
    '''
    GetFannedlb:

    Returns the galactic longitude and latitude of Ophiuchus Stream members
    identified as fanned candidates by Sesar et al. (2016).

    Args: None

    Returns:
        l (array)
        b (array)
    '''
    l = np.array([  5.70969653,
                    4.60461146,
                    6.44034091,
                    6.77319205,
                    7.91476459,
                    6.56873542])

    b = np.array([  30.71312853,
                    31.79888546,
                    30.48611263,
                    29.97629008,
                    27.04767995,
                    28.77314511])

    return l, b
#def

def GetFannedVRad():
    '''
    GetFannedVRad:

    Returns the radial velocity of the Ophiuchus Stream members identified as
    fanned candidates by Sesar et al. (2016).

    Args: None

    Returns:
        vrad (array)
    '''
    vrad = np.array([   289.9,
                        287.2,
                        271.4,
                        318.2,
                        236.6,
                        258.6])

    return vrad
#def
