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

    # Get the correct apocenter number
    apk = _apk_temp[apo-1]

    print('\n==============')
    print('WRITING MAIN.C')
    print('==============\n')

    for i,line in enumerate(main_in):
        # Add the offsets.
        if i == 90:
            main_out.write('  dxi_off[0] = '+str( apk[0] )+';\n')
            main_out.write('  dxi_off[1] = '+str( apk[1] )+';\n')
            main_out.write('  dxi_off[2] = '+str( apk[2] )+';\n')
            main_out.write('  dvi_off[0] = '+str( apk[3] )+';\n')
            main_out.write('  dvi_off[1] = '+str( apk[4] )+';\n')
            main_out.write('  dvi_off[2] = '+str( apk[5] )+';\n')
            print('Apocenter '+str(apo)+':')
            print('  dxi_off[0] = '+str( apk[0] )+';' )
            print('  dxi_off[1] = '+str( apk[1] )+';' )
            print('  dxi_off[2] = '+str( apk[2] )+';' )
            print('  dvi_off[0] = '+str( apk[3] )+';' )
            print('  dvi_off[1] = '+str( apk[4] )+';' )
            print('  dvi_off[2] = '+str( apk[5] )+';\n' )
        ##fi
        main_out.write(line)
    ###i

    # Move new main.c into place.
    subprocess.call('mv main.c source/',shell=True)
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
        subprocess.call('cp out/recentered_999.hdf5 ../Transfer/'+gc+'_run_'+run+'.hdf5', shell=True)
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

    import os
    from snapdata import Snapdata

    path = os.path.abspath(filename)
    data = Snapdata(path)

    print('Recentered median positions (code):')
    print('x: '+str(np.median(data.x)))
    print('y: '+str(np.median(data.y)))
    print('z: '+str(np.median(data.z)))
    print('vx: '+str(np.median(data.vx)))
    print('vy: '+str(np.median(data.vy)))
    print('vz: '+str(np.median(data.vz)))
#def


def CheckSigma(filename):
    '''
    Confirm that the recentering program placed the globular cluster at the
    correct position.

    Args:
        filename (string) - Name of the .hdf5 recentered file.
    '''

    import os
    from snapdata import Snapdata

    path = os.path.abspath(filename)
    data = Snapdata(path)

    print('Recentered velocity dispersions (code):')
    print('sigma vx: '+str(np.std(data.vx)))
    print('sigma vy: '+str(np.std(data.vy)))
    print('sigma vz: '+str(np.std(data.vz)))
#def

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

def GenPotParmsIC4(gc):
    '''
    GetPotParms:
    Generate a formatted potential parameter string block for the Gadget-2
    parameter file specifically for IC4.

    Args:
        gc (string) - The globular cluster unit system.
    Returns:
        None
    Outputs:
        Formatted string blocks
    '''

    _,code2kpc,_,code2Msol = ophstream.units.GetConversions(gc)

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
    conc = 15.3
    alpha = 1.8

    print('Potential Parameters: '+gc)
    print( 'MdExt    '+str(Md) )
    print( 'MbExt    '+str(Mb) )
    print( 'MhvirExt '+str(Mhvir) )
    print( 'aExt     '+str(Ad) )
    print( 'bExt     '+str(Bd) )
    print( 'cExt     '+str(Rc) )
    print( 'dExt     '+str(Rs) )
    print( 'concExt  '+str(conc) )
    print( 'alphaExt '+str(alpha) )
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

def GetClass4ICs(gc):
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
