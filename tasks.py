# ----------------------------------------------------------------------------
#
# TITLE - tasks.py
# AUTHOR - James Lane
# PROJECT - Ophiuchus Stream Honours
# CONTENTS:
#	1.CalcDensity3D
#
# ----------------------------------------------------------------------------
#
# Docstrings and metadata:
'''Analysis module for GADGET-2 globular cluster simulations'''
__author__ = "James Lane"


#Imports
import numpy as np
import scipy.linalg
from pykdtree.kdtree import KDTree

def CalcDensity2D(x,y,mass,nbranch=32):
    '''
    CalcDensity2D

    Calculates the 2D density of a series of particles.

    Args:
        x [array] - X particle positions.
        y [array] - Y particle positions.
        mass [float] - particle mass.
        nbranch [int] - number of branches to use for the density calculation.

    Returns:
        dens2d [array] - 2-dimensional density.
    '''
    z = np.zeros(len(x))
    kd_tree = KDTree( np.array([x,y,z]).T )
    dists, idx = kd_tree.query( np.array([x,y,z]).T, k=nbranch)
    return ( nbranch * mass ) / ( np.pi * np.power(dists[:,-1],2) )
#def

def CalcDensity3D(x,y,z,mass,nbranch=32):
    '''
    CalcDensity3D

    Calculates the 3D density of a series of particles.

    Args:
        x [array] - X particle positions.
        y [array] - Y particle positions.
        z [array] - Z particle positions.
        mass [float] - particle mass.
        nbranch [int] - number of branches to use for the density calculation.

    Returns:
        dens3d [array] - 3-dimensional density.
    '''
    kd_tree = KDTree( np.array([x,y,z]).T )
    dists, idx = kd_tree.query( np.array([x,y,z]).T, k=nbranch)
    return ( nbranch * mass ) / ( (4/3) * np.pi * np.power(dists[:,-1],3) )
#def

def CalcGalacticVRad(x,y,z,vx,vy,vz):
    '''
    CalcGalacticVRad

    Calculates heliocentric radial velocity using galactocentric rectangular
    coordinates.

    Args:
        x [array] - X particle positions.
        y [array] - Y particle positions.
        z [array] - Z particle positions.
        vx [array] - X particle velocities.
        vy [array] - Y particle velocities.
        vz [array] - Z particle velocities.

    Returns:
        vrad [array] - Radial velocities.
    '''

    # Create the velocity vector, but subtract out the rotational motion of the
    # sun and its peculiar velocities.
    # Peculiar velocity from Schonrich+ 2012, Rotational velocity from Bovy+ 2012
    # UVW frame is right handed, with U positive towards galactic center.
    helio_vx = 11.1
    helio_vy = 12.2 + 218.0
    helio_vz = 7.25

    vx -= helio_vx
    vy -= helio_vy
    vz -= helio_vz

    v_magnitude = np.sqrt(np.square(vx)+np.square(vy)+np.square(vz))

    # Create the radial vector from the sun to the stream.
    # Standard values, see Astropy documentation.
    x_radial = x - (-8.3)
    y_radial = y - (0.0)
    z_radial = z - (0.027)
    radial_magnitude = np.sqrt(np.square(x_radial)+np.square(y_radial)+np.square(z_radial))
    dist = radial_magnitude
    x_radial /= radial_magnitude
    y_radial /= radial_magnitude
    z_radial /= radial_magnitude
    x_tangent = vx / v_magnitude
    y_tangent = vy / v_magnitude
    z_tangent = vz / v_magnitude
    radial_dot_tangent = (  np.multiply(x_radial,x_tangent)+
                            np.multiply(y_radial,y_tangent)+
                            np.multiply(z_radial,z_tangent))
    vrad = v_magnitude * radial_dot_tangent
    return vrad
#def

def CalcTdynamical():
    '''
    Calculate dynamical properties of a Plummer sphere globular cluster

    Args:
        M [float] - Mass of the globular cluster.
        a [float] - Scale radius of the globular cluster.

    Returns:
        Vmax [float] - Maximum circular velocity.
        RVmax [float] - Radius of maximum circular velocity.
    '''




    return
#def


def FitGreatCircle(x, y, z):
    '''
    FitGreatCircle

    Fits a great circle to points in 3D space.

    Args:
        x [array] - X position

    ***INCOMPLETE***
    '''

    # Add (0, 0, 0) to the data, as the great circle should go through the origin
    x = np.append(x, 0)
    y = np.append(y, 0)
    z = np.append(z, 0)

    # Fit a linear plane through the data points
    A = np.c_[x, y, np.ones(x.shape[0])]
    C,_,_,_ = scipy.linalg.lstsq(A, z)

    # Calculate the great circle parameters
    z2 = C[0]**2 + C[1]**2

    theta0 = np.arcsin(z2/np.sqrt(z2 + z2**2))
    phi0 = np.arctan2(C[1], C[0])

    return C, theta0, phi0
#def

def GreatCircle(t, theta0, phi0):
    '''
    Calculates the point on a great circle defined my theta0 and phi0 in Cartesian coordinates.

    Sources:
        - http://demonstrations.wolfram.com/ParametricEquationOfACircleIn3D/

    Args:
        t: [float or 1D ndarray] phase angle of the point in the great circle
        theta0: [float] inclination of the great circle ('roll' in respect to phase=0, aka. red dot)
        phi0: [float] nodal angle of the great circle ('yaw' in respect to phase=0, aka. red dot)
    Return:
        [tuple or 2D ndarray] a tuple of (X, Y, Z) coordinates in 3D space
            (becomes a 2D ndarray if the input parameter t is also a ndarray)
    '''


    # Calculate individual cartesian components of the great circle points
    x = -np.cos(t)*np.sin(phi0) + np.sin(t)*np.cos(theta0)*np.cos(phi0)
    y =  np.cos(t)*np.cos(phi0) + np.sin(t)*np.cos(theta0)*np.sin(phi0)
    z =  np.sin(t)*np.sin(theta0)

    return (x, y, z)
#def

# def FitGreatCircleAlt(gl,gb,weights):
#
#     gl_rad = gl * (np.pi/180)
#     gb_rad = np.pi/2 - gb * (np.pi/180)
#
#     def func(v,c1,c2):
#         a = np.sqrt( np.square(1/c1) - 1 )
#         return np.arctan( np.power(np.sin(u + c2) * a , -1) )
#
#     popt = scipy.optimize.curve_fit(func,gl_rad,gb_rad,sigma=weights)
#
#     return popt
#
# def GreatCircleAlt()
