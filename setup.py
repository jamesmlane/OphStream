# ----------------------------------------------------------------------------
#
# TITLE - setup.py
# AUTHOR - James Lane
# PROJECT - Ophiuchus Stream Honours
# CONTENTS:
#	1. CheckRecentering
#   2. CheckSigma
#
# ----------------------------------------------------------------------------
#
# Docstrings and metadata:
'''Setup module for GADGET-2 globular cluster simulations'''
__author__ = "James Lane"


#Imports
import numpy as np
import ophstream.units
import subprocess
import os
import sys
import pdb
from snapdata import Snapdata

def MakeRecenter(ic, gc, apo):
    '''
    MakeRecenter:

    Create recentering files and prepares Gadget-2 for compiling. Has to be
    split up into multiple scipts for some reason... Because subprocess is
    dumb!!

    Args:
        ic (int) - Initial condition class to use.
        gc (string) - Unit system to use (ie: 'GC##')
        apo (int) - The apocenter number to use for these initial conditions

    Returns:
        None

    Output:
        file (.c) - new main.c file for compiling Gadget-2
    '''

    # Start by writing the initial conditions to the main.c file for Gadget-2
    main_in = open('parms/main_shell.c','r')
    main_out = open('main.c','w')

    # Get the correct class of initial conditions
    if ic == 2:
        _apk_temp = GetClass2ICs(gc)
    elif ic == 3:
        _apk_temp = GetClass3ICs(gc)
    elif ic == 4:
        _apk_temp = GetClass4ICs(gc)
    ##fi


#def

def RunRecenter(gc,run):
    '''
    RunRecenter:

    Run Gadget-2 then examine the output. Has to be split up into multiple
    scipts for some reason... Because subprocess is dumb!!

    Args:
        gc (string) - Unit system to use (ie: 'GC##')
        run (string) - 3-digit run number (ie: '###')

    Returns:
        None

    Output:
        file (hdf5) - Recentered plummer model initial conditions.
    '''
    # Bind raw_input to input for Python 2.
    try:
    	input = raw_input
    except NameError:
    	pass
    #try

    print('\n================')
    print('RUNNING GADGET-2')
    print('================\n')

    subprocess.call('source ./gadget_run',shell=True)

    CheckRecentering('out/recentered_999.hdf5')
    CheckSigma('out/recentered_999.hdf5')

    recenter_good = input('\nRecentering successfull? name '+gc+'_run_'+run+'? [y/n] ... ')
    if recenter_good == 'y':
        subprocess.call('cp out/recentered_999.hdf5 ../../../gadget/initial_conditions/plummer/IC4/'+gc+'_run_'+run+'.hdf5', shell=True)
    ##fi
    subprocess.call('rm -rf out/*', shell=True)
    print('Done...')
#def

def CheckRecentering(filename):
    '''
    Confirm that the recentering program placed the globular cluster at the
    correct position.

    Args:
        filename (string) - Name of the .hdf5 recentered file.
    '''

    path = os.path.abspath(filename)
    data = Snapdata(path)

    print('Recentered median positions (code):')
    print('x:  '+str(np.median(data.x)))
    print('y:  '+str(np.median(data.y)))
    print('z:  '+str(np.median(data.z)))
    print('vx: '+str(np.median(data.vx)))
    print('vy: '+str(np.median(data.vy)))
    print('vz: '+str(np.median(data.vz)))
#def


def CheckSigma(filename):
    '''
    Confirm that the dynamically correct globular cluster has been created.

    Args:
        filename (string) - Name of the .hdf5 recentered file.
    '''

    path = os.path.abspath(filename)
    data = Snapdata(path)

    print('Recentered velocity dispersions (code):')
    print('sigma vx: '+str(np.std(data.vx)))
    print('sigma vy: '+str(np.std(data.vy)))
    print('sigma vz: '+str(np.std(data.vz)))
#def

def GetRhalf(filename, code2kpc=None, gc=None):
    '''
    Get the half-mass radius in pc for a recentered globular cluster.

    Args:
        filename (string) - Name of the .hdf5 recentered file.
        gc (string) - Globular cluster model.

    Output:
        The measured half-mass radius in pc.
    '''

    if code2kpc == None and code2kpc == None:
        sys.exit('No parameters supplied to GetTimesIC4')
    elif code2kpc == None:
        code2kpc,_,_,_ = ophstream.units.GetConversions(gc)
    ##fi

    path = os.path.abspath(filename)
    data = Snapdata(path)
    r_center = np.sqrt( np.square( data.x - np.median(data.x) ) +
                        np.square( data.y - np.median(data.y) ) +
                        np.square( data.z - np.median(data.z) ) )
    r_sorted = np.sort(r_center) * code2kpc * 1000
    print('Rhalf: '+str(r_sorted[ int(len(r_sorted)/2) ])+' pc')
#def

def CheckSigmaRHalf(filename,gc):
    '''
    Check the, more accurate, velocity dispersion within the half-mass radius.

    Args:
        filename (string) - Name of the .hdf5 recentered file.
        gc (string) - Globular cluster model.
    '''

    path = os.path.abspath(filename)
    data = Snapdata(path)
    rhalf = 1.3 # In our coordinate system because 'a' is the scale
    r_center = np.sqrt( np.square( data.x - np.median(data.x) ) +
                        np.square( data.y - np.median(data.y) ) +
                        np.square( data.z - np.median(data.z) ) )
    inside_rhalf = np.where( r_center < rhalf )

    print('Recentered velocity dispersions within half-mass radius (code):')
    print('sigma vx: '+str(np.std(data.vx[inside_rhalf])))
    print('sigma vy: '+str(np.std(data.vy[inside_rhalf])))
    print('sigma vz: '+str(np.std(data.vz[inside_rhalf])))

def GenPointMassIC(filename,coords):
    '''
    Generate the initial conditions for a single point mass run of GADGET-2,
    Flips the sign of the velocities.

    Args:
        filename (string) - filename of the initial conditions
        coords (array) - array of position and velocity [x,y,z,vx,vy,vz]
                            in code units.
    '''
    coords[3:] *= -1
    fileout = open(filename,'w')
    fileout.write('1 1E-10\n') # Single particle, point mass
    for coord in coords:
        fileout.write(str(coord)+' ')
    ###i
    fileout.close()

def GenRecenterParms(ic,gc):
    '''
    GetRecenterParms:
    Generate a formatted recenter coordinate string block for main.c

    Args:
        ic (int) - The class of apocenter to return (2,3)
        gc (string) - The globular cluster unit system.
        return (boolean) - Return the parameters
    Returns:
        None
    Outputs:
        Formatted string blocks
    '''

    if ic == 2:
        _,apo2,apo3,apo4,_ = GetClass2ICs(gc)
    elif ic == 3:
        _,apo2,apo3,apo4 = GetClass3ICs(gc)
    elif ic == 4:
        _,apo2,apo3,apo4 = GetClass4ICs(gc)
    ##fi

    print('Apocenter 2:')
    print('dxi_off[0] = '+str( apo2[0] )+';' )
    print('dxi_off[1] = '+str( apo2[1] )+';' )
    print('dxi_off[2] = '+str( apo2[2] )+';' )
    print('dvi_off[0] = '+str( apo2[3] )+';' )
    print('dvi_off[1] = '+str( apo2[4] )+';' )
    print('dvi_off[2] = '+str( apo2[5] )+';' )
    print('\n')

    print('Apocenter 3:')
    print('dxi_off[0] = '+str( apo3[0] )+';' )
    print('dxi_off[1] = '+str( apo3[1] )+';' )
    print('dxi_off[2] = '+str( apo3[2] )+';' )
    print('dvi_off[0] = '+str( apo3[3] )+';' )
    print('dvi_off[1] = '+str( apo3[4] )+';' )
    print('dvi_off[2] = '+str( apo3[5] )+';' )
    print('\n')

    print('Apocenter 4:')
    print('dxi_off[0] = '+str( apo4[0] )+';' )
    print('dxi_off[1] = '+str( apo4[1] )+';' )
    print('dxi_off[2] = '+str( apo4[2] )+';' )
    print('dvi_off[0] = '+str( apo4[3] )+';' )
    print('dvi_off[1] = '+str( apo4[4] )+';' )
    print('dvi_off[2] = '+str( apo4[5] )+';' )
#def

def GenPotParms(gc):
    '''
    GetPotParms:
    Generate a formatted potential parameter string block for the Gadget-2
    parameter file.

    Args:
        gc (string) - The globular cluster unit system.
    Returns:
        None
    Outputs:
        Formatted string blocks
    '''

    _,code2kpc,_,code2Msol = ophstream.units.GetConversions(gc)

    # Masses in Msol
    Md = 5E10 / code2Msol
    Mb = 2.75E10 / code2Msol
    Mh = 9.5E12 / code2Msol

    # Scale distances in kpc
    Ad = 3.5 / code2kpc
    Bd = 0.525 / code2kpc
    Ab = 0.7 / code2kpc
    Ah = 140 / code2kpc

    print('Potential Parameters: '+gc)
    print( 'MdExt '+str(Md) )
    print( 'MbExt '+str(Mb) )
    print( 'MhExt '+str(Mh) )
    print( 'aExt '+str(Ad) )
    print( 'bExt '+str(Bd) )
    print( 'cExt '+str(Ab) )
    print( 'dExt '+str(Ah) )
#def

def GenPotParmsIC4(code2kpc=None,code2Msol=None,gc=None,output=False):
    '''
    GetPotParms:
    Generate a formatted potential parameter string block for the Gadget-2
    parameter file specifically for IC4.

    Args:
        gc (string) - The globular cluster unit system.
        output (Bool) - Output the values instead?
    Returns:
        None
    Outputs:
        Formatted string blocks
    '''

    if code2kpc == None and code2Msol == None and gc == None:
        sys.exit('No arguments supplied for GenPotParmsIC4')
    ##fi

    if gc == None and ( code2Msol == None or code2kpc == None ):
        sys.exit('Must supply arguments to both code2kpc and code2Msol if not supplying gc')
    ##fi

    if code2kpc == None and code2Msol == None:
        _,code2kpc,_,code2Msol = ophstream.units.GetConversions(gc)
    ##fi

    # Masses in Msol
    Md = 6.8E10 / code2Msol
    Mb = 0.5E10 / code2Msol
    Mhvir = 8.0E11 / code2Msol

    # Scale distances in kpc
    Ad = 3.0 / code2kpc
    Bd = 0.28 / code2kpc
    Rc = 1.9 / code2kpc
    Rs = 16 / code2kpc

    # Unchanging parameters
    conc = 15.3 # Fixed from Bovy+2015
    alpha = 1.8

    if output == False:
        print('Potential Parameters:')
        print( 'MdExt    '+str(Md) )
        print( 'MbExt    '+str(Mb) )
        print( 'MhvirExt '+str(Mhvir) )
        print( 'aExt     '+str(Ad) )
        print( 'bExt     '+str(Bd) )
        print( 'cExt     '+str(Rc) )
        print( 'dExt     '+str(Rs) )
        print( 'concExt  '+str(conc) )
        print( 'alphaExt '+str(alpha) )
    ##fi
    if output == True:
        return Md,Mb,Mhvir,Ad,Bd,Rc,Rs,conc,alpha
    ##fi
#def

def GetTimesIC4(apo,code2Myr=None,gc=None):
    '''
    GetTimesIC4:
    Generate class 4 simulation lengths.

    Args:
        gc (string) - The globular cluster unit system.
        apo (string) - The IC4 apocenter number.

    Returns
        t (float) - The simulation length in code units
    '''
    if code2Myr == None and code2Myr == None:
        sys.exit('No parameters supplied to GetTimesIC4')
    elif code2Myr == None:
        code2Myr,_,_,_ = ophstream.units.GetConversions(gc)
    ##fi

    ic4_apo_times = np.array([123.9691,360.9102,600.8504,836.7917])
    return ic4_apo_times[ int(apo-1) ] / code2Myr
#def

def GetClass2ICs(gc):
    '''
    GetApo2ICs:
    Generate class 2 apocenters in the correct units.

    Args:
        gc (string) - The globular cluster unit system.
    Returns:
        apos (list) - list of ic numpy arrays [x,y,z,vx,vy,vz], going forward
                        in time, in kpc and km/s
    '''
    _,code2kpc,code2kms,_ = ophstream.units.GetConversions(gc)

    apo1 = np.array([-0.5961,-0.9372,-11.6343,-97.1332,21.4329,4.5783])
    apo2 = np.array([-10.3206,2.9708,4.1791,46.4406,-3.3091,113.0465])
    apo3 = np.array([5.9037,-0.7927,10.0107,84.4855,-28.9282,-54.0116])
    apo4 = np.array([6.4276,-2.9608,-9.229,-81.2,21.2537,-65.909])
    apo5 = np.array([-9.6442,3.0358,-5.5803,-51.7081,27.0406,102.4739])

    apo1[:3] /= code2kpc
    apo1[3:] /= code2kms
    apo2[:3] /= code2kpc
    apo2[3:] /= code2kms
    apo3[:3] /= code2kpc
    apo3[3:] /= code2kms
    apo4[:3] /= code2kpc
    apo4[3:] /= code2kms
    apo5[:3] /= code2kpc
    apo5[3:] /= code2kms

    return [apo1,apo2,apo3,apo4,apo5]
#def

def GetClass3ICs(gc):
    '''
    GetApo2ICs:
    Generate class 2 apocenters in the correct units.

    Args:
        gc (string) - The globular cluster unit system.
    Returns:
        apos (list) - list of ic numpy arrays [x,y,z,vx,vy,vz], going forward
                        in time, in kpc and km/s
    '''
    _,code2kpc,code2kms,_ = ophstream.units.GetConversions(gc)

    apo1 = np.array([-2.0986,  -0.8438, -14.8627, -81.8248,  18.9293,	 11.0887, 108.9729])
    apo2 = np.array([-12.1537,  3.725,   7.8313,   51.0661, -6.7023,	 83.5431, 314.9216])
    apo3 = np.array([10.4652,  -2.0605,  10.5005,  61.8756, -22.5759,	-65.6819, 523.8696])
    apo4 = np.array([4.2073,   -2.6032, -14.1864, -78.2368,  22.5567, -27.5625, 728.8186])

    apo1[:3] /= code2kpc
    apo1[3:] /= code2kms
    apo2[:3] /= code2kpc
    apo2[3:] /= code2kms
    apo3[:3] /= code2kpc
    apo3[3:] /= code2kms
    apo4[:3] /= code2kpc
    apo4[3:] /= code2kms

    return [apo1,apo2,apo3,apo4]
#def

def GetClass4ICs(code2kpc,code2kms):
    '''
    GetApo2ICs:
    Generate class 2 apocenters in the correct units.

    Args:
        gc (string) - The globular cluster unit system.
    Returns:
        apos (list) - list of ic numpy arrays [x,y,z,vx,vy,vz], going forward
                        in time, in kpc and km/s
    '''

    apo1 = np.array([-2.7111,  -0.866,  -16.4423, -71.5236,	 17.2722,  11.7234])
    apo2 = np.array([-13.0558,  4.3714,	 9.1938,   49.1522, -8.1264,   74.042])
    apo3 = np.array([11.2067,  -2.4948,  11.993,   55.9778, -22.1673, -58.0949])
    apo4 = np.array([4.6442,   -3.1648, -15.6974, -68.3264,	 23.1416, -24.048])

    apo1[:3] /= code2kpc
    apo1[3:] /= code2kms
    apo2[:3] /= code2kpc
    apo2[3:] /= code2kms
    apo3[:3] /= code2kpc
    apo3[3:] /= code2kms
    apo4[:3] /= code2kpc
    apo4[3:] /= code2kms

    return [apo1,apo2,apo3,apo4]
#def
