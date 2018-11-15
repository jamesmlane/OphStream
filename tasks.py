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
from astropy import table
# from ophstream.autogadget import OphstreamRun
from astropy import coordinates as coords
from astropy import units
from scipy.optimize import curve_fit,newton
from scipy.integrate import quad as quad_integrate
from snapdata import Snapdata
import ophstream.units
import ophstream.tasks
import sys
import os
import pdb

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

def StackCatalog(filelist,prog_names=None,prog_mass=None,prog_radius=None,prog_vdisp=None):
    '''
    Stack catalogs of observed stream information.

    Args:
        filelist (array) - The list of files to stack
        prog_names (str array) - N-d array of progenitor names
        prog_mass (flt array) - N-d array of progenitor masses in Msol
        prog_radius (flt arry) - N-d array of progenitor 3d r1/2 in pc
        prog_vdisp (flt array) - N-d array of progenitor 1d vdisp in km/s

    Returns:
        None

    Outputs:
        catalog (.txt) - Catalog of properties for all inputs
    '''

    n_files = len(filelist)
    if (len(prog_names) != n_files) or (len(prog_mass) != n_files) or (len(prog_radius) != n_files) or (len(prog_vdisp) != n_files):
        sys.exit('Not all input arrays of the same length as filelist')
    ##fi

    meas_data = np.empty((n_files,10))
    import pdb
    for i in range(n_files):
        meas_data[i,:] = np.genfromtxt(filelist[i])
    ###i
    fmt_string = ('%4s','%4.2f','%4.2f','%4.2f','%4.2f','%4.2f','%4.2f','%4.2f','%4.2f','%4.2f','%4.2f','%4.2f','%4.2f','%4.2f')
    pdb.set_trace()
    np.savetxt('catalog_out.txt', np.array([prog_names,prog_mass,prog_radius,prog_vdisp,meas_data.T[:]]).T, fmt=fmt_string)
    pdb.set_trace()
#def

def StackCatalogAP(filelist,prog_names=None,prog_mass=None,prog_radius=None,prog_vdisp=None):
    '''
    Stack catalogs of observed stream information with astropy.

    Args:
        filelist (array) - The list of files to stack
        prog_names (str array) - N-d array of progenitor names
        prog_mass (flt array) - N-d array of progenitor masses in Msol
        prog_radius (flt arry) - N-d array of progenitor 3d r1/2 in pc
        prog_vdisp (flt array) - N-d array of progenitor 1d vdisp in km/s

    Returns:
        None

    Outputs:
        catalog (.txt) - Catalog of properties for all inputs
    '''

    n_files = len(filelist)
    if (len(prog_names) != n_files) or (len(prog_mass) != n_files) or (len(prog_radius) != n_files) or (len(prog_vdisp) != n_files):
        sys.exit('Not all input arrays of the same length as filelist')
    ##fi

    meas_data = np.empty((n_files,10))



    import pdb
    for i in range(n_files):
        meas_data[i,:] = np.genfromtxt(filelist[i])
    ###i

    meas_names = (  'Name','Mass (Msol)','3D R1/2 (pc)','1D Vdisp (km/s)','Length (deg)','Length Error',
                    'Length N (Nmsto)','Length N Error','Width (arcmin)','Width Error','Width N (Nmsto)',
                    'Width N Error','Vdisp (km/s)','Vdisp Error')
    # meas_fmts = {   'Name':'%.2f',
    #                 'Mass':
    #                 '3D R1/2':
    #                 '1D Vdisp':
    #                 'Length':
    #                 'Length Error':
    #                 'Length N':
    #                 'Length N Error':'Width','Width Error','Width N',
    #                                 'Width N Error','Vdisp','Vdisp Error'}
    meas_table = Table( [ prog_names, prog_mass, prog_radius, prog_vdisp, meas_data.T ]
                        , names=meas_names)
    meas_table.write('Catalog_out.txt', format='ascii.fixed_width',
                    overwrite=True)
#def

def StackCatalogsAutoGadget(filelist,
                            output='./stream_params/autogadget_combined_catalog.FIT'):
    '''
    StackCatalogsAutoGadget:
        Combine multiple autogadget catalogs of stream properties into one
        stream property catalog.

    Args:
        filelist (str array) - Array of strings of filenames

    Returns:
        None

    Outputs:
        combined_catalog (astropy .FIT) - The combined catalog of stream
            properties. Will by default be named autogadget_combined_catalog.FIT

    '''
    # Empty list for new catalogs array
    catalogs = []
    for i,filename in enumerate(filelist):
        catalogs.append( table.Table.read(filename) )
    ####

    # Make the new catalog
    new_catalog = table.vstack(catalogs)
    new_catalog.write('./autogadget_combined_catalog.FIT', format='fits',
                        overwrite=True)

#def

def GetStreamLatLong(snap_file, conversions, output=False):
    '''
    GetStreamLatLong:
        Convert a snapshot into binned histograms of the data in stream-spline
        latitude and longitude. Stream equator is defined by a quadratic
        best-fit.

    Args:
        snap_file (hdf5) - The snapshot to convert
        conversions (arr) - Array of code2physical conversion factors of form:
            [code2kpc, code2kms, code2Msol, code2Myr]
        output (bool) - Output to file or return (2 variables) [False]

    Returns:
        None

    Output:
        longlat (2xn arr) - Array of stream-based longitudes and latitudes
            (respectively in that order)

    '''

    # First get the conversion factors
    if len(conversions) != 4:
        sys.exit('The conversions array is not length 4')
    ##fi
    code2kpc, code2kms, code2Msol, code2Myr = conversions[:]

    # Read in the specified file and extract particle data:
    abs_path = os.path.abspath(snap_file)
    data = Snapdata(abs_path)
    n_part = data.npart
    p_mass = data.p_mass * code2Msol
    x = data.x * code2kpc
    y = data.y * code2kpc
    z = data.z * code2kpc
    vx = data.vx * code2kms
    vy = data.vy * code2kms
    vz = data.vz * code2kms
    time = data.time * code2Myr

    # Convert positions to l,b
    galactocen_coords = coords.Galactocentric(x=x * units.kpc,
                                              y=y * units.kpc,
                                              z=z * units.kpc)
    galactic_coords = galactocen_coords.transform_to(coords.Galactic)
    gl = np.array(galactic_coords.l)
    gb = np.array(galactic_coords.b)

    # Calculate distance to sun in kpc, manually
    dist = np.sqrt(np.square(x+8.3)+np.square(y)+np.square(z-0.027))

    # Get surface density
    surfdens = CalcDensity2D(gl,gb,p_mass)
    surfdens_weights = np.power(surfdens,-3)

    ################################################################################

    # Fit the stream with a weighted quadratic:

    # The quadratic function
    def quad_func(x,a,b,c):
        return a*np.power(x,2)+b*x+c
    #def

    # Weighted quadratic
    where_quad_fit = np.where( (gl > 2) &
                               (gl < 8) &
                               (gb > 28) &
                               (gb < 33) )[0]
    wquad_popt,_ = curve_fit(quad_func,gl[where_quad_fit],gb[where_quad_fit],
                            sigma=surfdens_weights[where_quad_fit])
    wquada, wquadb, wquadc = wquad_popt
    wquadgl = np.arange(0,10,0.1)
    wquadgb = quad_func(wquadgl,wquada,wquadb,wquadc)

    ################################################################################

    # Use Newtons method to get the nearest neighbour point on the quadratic line to
    # each point in the stream, then get arc length along the stream.

    # Function to return dot product between quadratic point and stream point
    def stream_quad_dot(l,streamx,streamy,a,b,c):
        # Get the derivative at the quadratic point
        quad_b = a*(l**2) + b*l + c
        quad_slope = 2*a*l + b
        quad_yint = quad_b-l*quad_slope
        return l*(streamx-l) + (quad_b-quad_yint)*(streamy-quad_b)
    #def

    # Function to return arc length along the curve:
    def stream_quad_arc(l,a,b):
        return np.sqrt( 1 + np.square( 2*a*l + b ) )
    #def

    # Find the closest point on the quadratic to each stream point
    close_l = np.zeros(n_part)
    for i in range(n_part):
        close_l[i] = newton(stream_quad_dot, gl[i], tol=0.001, maxiter=100,
                            args=(gl[i],gb[i],wquada,wquadb,wquadc))
    ###i
    close_b = wquada*np.square(close_l) + wquadb*close_l + wquadc

    # Latitude in the new coordinate system is the distance from the curve:
    stream_lat = np.sqrt(np.square(close_l-gl)+np.square(close_b-gb))

    # Get the sign of the latitude by seeing if the new point is lower than the
    # curve.
    for i in range(n_part):
        if gb[i] < ( wquada*(gl[i])**2 + wquadb*gl[i] + wquadc ):
            stream_lat[i] *= -1
    ###i

    # Longitude in the new coordinate system is the arc length along the curve
    initial_long = np.median(gl)
    stream_long = np.zeros(n_part)
    arclen_args = (wquada,wquadb)
    for i in range(n_part):
        stream_long[i] = quad_integrate(stream_quad_arc, initial_long, gl[i],
                                        args=arclen_args )[0]
    ###i

    # Now output to file or return depending on user keywords.
    if output == True:
        outarr = np.array([ stream_long, stream_lat ])
        np.save(snap_file[:-5]+'_spline_latlon.npy', outarr)
    else:
        return np.array([ stream_long, stream_lat ])
    ##ie
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
