# ----------------------------------------------------------------------------
#
# TITLE - autogadget
# AUTHOR - James Lane
# PROJECT - ophstream
# CONTENTS:
#	1.
#
# ----------------------------------------------------------------------------
#
# Docstrings and metadata:
'''
This code will run auto-gadget, a python shell around (most of) the running
of GADGET-2 for the purpose of simulating the Ophiuchus stream.
'''
__author__ = "James Lane"

### Imports

# Primary use imports
import numpy as np
import ophstream.setup as ophsetup
import ophstream.misc as ophmisc
import ophstream.tasks as ophtasks
import ophstream.stream_properties as ophstream_props
from astropy import coordinates as coords
from astropy import units
from astropy import table
from snapdata import Snapdata
from matplotlib import pyplot as plt
import matplotlib.cm
import matplotlib as mpl
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import binned_statistic_2d as bin2d
from scipy.stats import linregress
from scipy.optimize import curve_fit,newton
from scipy.integrate import quad as quad_integrate
from scipy.integrate import dblquad as quad_integrate2d
import scipy.interpolate


# Secondary use imports
import subprocess
import shutil
import os
import sys
import pdb
import time
import glob
import psutil

class OphstreamRun:

    '''
    OphstreamRun:

    A class that handles properties of a run of the Ophiuchus stream progenitor.

    Attributes:
        mass (float) - The mass of the progenitor [Msol]
        magnitude (float) - The V-band absolute magnitude of the progenitor. [Msol]
        radius_2d (float) - The projected half-mass radius of the progenitor. [pc]
        radius_3d (float) - The half-mass radius of the progenitor. [pc]
        radius_scale (float) - The Plummer scale radius of the progenitor. [pc]
        apocenter (int) - The apocenter number from which to start this simulations.
        ic (int) - The class of initial condition to use for this simulations.

    Args:
        mass (float) - The mass of the progenitor [solar masses]
        radius (float) - The half-mass (3d) of the progenitor [pc]
        apocenter (int) - The apocenter number at which the orbit stars
        magnitude (bool) - If True the mass argument is given in
                            magnitudes [False]
        hmr2d (bool) - If True the half-mass radius argument is the projected
                        half-mass radius of the system [False]
        ic (int) - The initial condition class, likely 4. Support for others
                        likely doesn't exist [4]
        MtoL (float) - The mass to light ratio. Assume 1.45 (Messier 5)

    Returns:
        None

    Raises:
        None

    '''

    def __init__(   self,
                    mass,
                    radius,
                    apocenter,
                    magnitude=False,
                    hmr2d=False,
                    ic=4,
                    MtoL=1.45,
                    Nparticles=1E4
                    ):

        # Determine the mass and magnitude, universally assume M/L = 1.45
        if magnitude == False:
            self.mass = mass
            self.magnitude = 4.77-(2.5*np.log10(mass/MtoL))
        else:
            self.mass = MtoL * np.power(10,(4.77-mass)/2.5)
            self.magnitude = mass
        ##fi

        # Determine the radius
        if hmr2d == False:
            self.radius_2d = (3.0/4.0) * radius
            self.radius_3d = radius
            self.radius_scale = radius / 1.3
        if hmr2d == True:
            self.radius_2d = radius
            self.radius_3d =  (4.0/3.0) * radius
            self.radius_scale = radius * (4.0/3.0) / 1.3
        ##fi

        # Set the scale factors.
        self.code2pc = self.radius_3d / 1.3
        self.code2kpc = self.radius_3d / 1.3 / 1000
        self.code2Msol = self.mass * 1
        self.code2kms = np.sqrt( 0.004301 * self.code2Msol / self.code2pc )
        self.code2Myr = 0.978 * self.code2pc / self.code2kms

        # Set the dynamic values
        self.vmax = 0.04069 * np.sqrt( self.code2Msol / self.code2pc )
        self.rmax = np.sqrt(2) * self.code2pc
        self.tdyn = 0.98 * self.rmax / self.vmax # maximum velocity / radius

        # Set supporting information
        self.apocenter = apocenter
        self.ic = ic
        self.matched_snap = None    # The snap that best-matches the real data.
                                    # Will be filled in by FindBestSnap

        # Declare the filename that will accompany this progenitor. Note has
        # no ending so will need to be added in application.
        self.filename = (   'M'+'{:.1E}'.format(self.mass)+'_Rh'+str(round(self.radius_3d,1))+'_apo'+str(self.apocenter) )

    #def

    def GetConversions(self):
        '''
        GetConversions:
            Convenience function that returns the factors to convert simulation
            units to physical units.

        Args:
            None

        Returns:
            factors (arr) - Array of conversion factors: code2kpc, code2kms,
                code2Msol, code2Myr

        Outputs:
            recentered
        '''
        return [self.code2kpc, self.code2kms, self.code2Msol, self.code2Myr]
    #def

    def Recenter(   self,
                    nprocs=5,
                    sim_dir='./recenter',
                    dest_dir='./initial_conditions',
                    storage_dir='./storage/recenter/'):
        '''
        Recenter:
            Generates the initial conditions of the simulation by running
            GADGET-RECENTER to place the progenitor at the proper distance.
            This method anticipates being in a directory with parms/
            containing the default parameter files.

        Args:
            nprocs (int) - Number of processors to use when recentering [5]
            sim_dir (str) - Location to perform recenter simulations relative to
                home_dir [./recenter]
            dest_dir (str) - Location to save the recentered file relative
                to home_dir [./initial_conditions]
            storage_dir (str) - The directory to store misc. parameter files
                [./sim_storage]

        Returns:
            None

        Outputs:
            recentered
        '''

        print('\n=================================')
        print('Recentering '+self.filename)
        print('=================================')

        # Append self.filename automatically to directories
        if dest_dir[-len(self.filename):] != self.filename:
            dest_dir += '/'+self.filename
        if storage_dir[-len(self.filename):] != self.filename:
            storage_dir += '/'+self.filename
        ##fi

        # Get the argument directory root paths
        sim_path = os.path.abspath(sim_dir)
        dest_path = os.path.abspath(dest_dir)
        storage_path = os.path.abspath(storage_dir)

        # Make sure directory structures exist
        if os.path.isdir(dest_path) == False:
            os.makedirs(dest_path)
        ##fi

        # Make sure directory structures exist
        if os.path.isdir(storage_path) == False:
            os.makedirs(storage_path)
        ##fi

        # Set home directory, descend into the recentering directory.
        home_path = os.getcwd()
        os.chdir(sim_path)

        # Get the initial conditions for the requested apocenter. Only support
        # class 4 initial conditions right now (2 and 3 defunct so don't
        # bother).
        print('\nWriting parameters...')
        apk = ophsetup.GetClass4ICs(self.code2kpc,self.code2kms)[self.apocenter-1]

        # First create the position-velocity initial condition parameter file.
        dxv_out = open('./parms/dxv.param', 'w')
        dxv_out.write( str( apk[0] )+' ' ) # Positions
        dxv_out.write( str( apk[1] )+' ' )
        dxv_out.write( str( apk[2] )+' ' )
        dxv_out.write( str( apk[3] )+' ' ) # Velocities
        dxv_out.write( str( apk[4] )+' ' )
        dxv_out.write( str( apk[5] ) )
        dxv_out.close()

        # Calculate the dynamical timescale in code units. This is ~2.28 in
        # code units. Make the evolution time 5 times the dynamical time.
        evolution_time = 5 * self.tdyn / self.code2Myr

        # # Create the file that runs the program:
        # gadget_run = open('./gadget_run','w')
        # gadget_run.write('#!/bin/bash\n') # generic bash
        # gadget_run.write('rm -f log.txt\n') # Remove prior log
        # gadget_run.write('mpirun -np '+str(nprocs)) # Run mpi with nprocs processors
        # gadget_run.write(' ./source/Gadget2-RECENTER') # Gadget-2 binary
        # gadget_run.write(' parms/recenter.param') # Gadget-2 parameters
        # gadget_run.write(' '+str(round(evolution_time,2))) # t_relax (virialize)
        # gadget_run.write(' '+str(round(evolution_time,2))) # t_evolve
        # gadget_run.write(' ./parms/dxv.param') # position-velocity parameters
        # gadget_run.write(' >> log.txt 2>&1') # Write to log.txt

        # Check to make sure we have an 'out' directory, otherwise make one.
        if 'out' not in os.listdir('./'):
            os.mkdir('out')
        ##

        # Write the command to run GADGET-2
        recenter_command = 'mpirun -np '+str(nprocs) # Run mpi with nprocs processors
        recenter_command += ' ./source/Gadget2-RECENTER' # Gadget-2 binary
        recenter_command += ' parms/recenter.param' # Gadget-2 parameters
        recenter_command += ' '+str(round(evolution_time,2)) # t_relax (virialize)
        recenter_command += ' '+str(round(evolution_time,2)) # t_evolve
        recenter_command += ' ./parms/dxv.param' # position-velocity parameters
        recenter_command += ' >> log.txt 2>&1' # Write to log.txt

        # Execute
        sim_timing0 = time.time()
        print('\nRunning GADGET-2...')
        print('Command is: '+recenter_command+'\n')
        # subprocess.call('rm -f ./log.txt', shell=True)
        subprocess.call(recenter_command, shell=True)
        print('Done running GADGET-2, took '+str(round(time.time()-sim_timing0))+' seconds\n')

        # Check the recentering.
        ophsetup.CheckRecentering('out/recentered_999.hdf5')
        print('Expected (code):')
        print( 'x: '+str( apk[0] ) )
        print( 'y: '+str( apk[1] ) )
        print( 'z: '+str( apk[2] ) )
        print( 'vx:'+str( apk[3] ) )
        print( 'vy:'+str( apk[4] ) )
        print( 'vz:'+str( apk[5] )+'\n' )

        # # Check the velocity dispersion.
        # ophsetup.CheckSigma('out/recentered_999.hdf5')
        # print('Expected (code):')
        # print( 'x: '+str( '?' ) ) # Determine what the velocity dispersion is?
        # print( 'y: '+str( '?' ) ) # Determine what the velocity dispersion is?
        # print( 'z: '+str( '?' )+'\n' ) # Determine what the velocity dispersion is?

        # Check the size.
        ophsetup.GetRhalf('out/recentered_999.hdf5', code2kpc=self.code2kpc)
        print('Expected: '+str(self.radius_3d)+' (pc)')

        # Put parameters away
        subprocess.call('mv ./parms/dxv.param '+storage_path+'/dxv_'+self.filename+'.param', shell=True)
        subprocess.call('mv ./log.txt '+storage_path+'/log_'+self.filename+'.txt', shell=True)

        # Now do something with the end result. Get home first
        os.chdir(home_path)
        subprocess.call('mv '+sim_path+'/out/recentered_999.hdf5 '+dest_path+'/'+self.filename+'.hdf5', shell=True)

        # Cleanup the results of the recentering run.
        subprocess.call('rm -rf '+sim_path+'/out/*', shell=True)

        # Now done
        print('\nDone Recentering for '+self.filename)

    #def

    def Simulate(   self,
                    nprocs=8,
                    sim_dir='simulate',
                    init_dir='./initial_conditions',
                    dest_dir='./results',
                    snap_dir='./snaps',
                    storage_dir='./storage/simulate'):
        '''
        Simulate:
            Run a simulation of the ophiuchus stream.

        Args:
            nprocs (int) - Number of processors to use when simulating [5]
            sim_dir (str) - The directory to perform the simulation in
                [./simulate]
            init_dir (str) - The directory in which the initial conditions are
                stored [./recenter]
            dest_dir (str) - The directory to store the results in [./results]
            snap_dir (str) - The directory to store the snaps in [./snaps]
            storage_dir (str) - The directory to store misc. parameter files
                [./sim_storage]

        Returns:
            None

        Outputs:
            Plots and analysis of the results.
        '''

        print('\n=================================')
        print('Simulating '+self.filename)
        print('=================================')

        # Append self.filename to directories automatically
        if init_dir[-len(self.filename):] != self.filename:
            init_dir += '/'+self.filename
        if dest_dir[-len(self.filename):] != self.filename:
            dest_dir += '/'+self.filename
        if snap_dir[-len(self.filename):] != self.filename:
            snap_dir += '/'+self.filename
        if storage_dir[-len(self.filename):] != self.filename:
            storage_dir += '/'+self.filename
        ##fi

        # Get the argument directory root paths
        init_path = os.path.abspath(init_dir)
        sim_path = os.path.abspath(sim_dir)
        dest_path = os.path.abspath(dest_dir)
        snap_path = os.path.abspath(snap_dir)
        storage_path = os.path.abspath(storage_dir)

        # Make sure directory infrastructure exists
        if os.path.isdir(snap_path) == False:
            os.makedirs(snap_path)
        ##fi

        if os.path.isdir(storage_path) == False:
            os.makedirs(storage_path)
        ##fi

        # pdb.set_trace()
        if os.path.isdir(init_path) == False:
            sys.exit('Need to Recenter first and place results in init_dir')
        ##fi

        # Set home directory, descend into the simulating directory.
        home_path = os.getcwd()
        os.chdir(sim_path)

        # Prepare the simulation. First ensure the directory structure exists.
        if os.path.isdir('./out') == False:
            os.mkdir('out')
        ##

        # Place the initial condition file in the simulation directory.
        print('\nWriting parameters...')
        subprocess.call('cp -f '+init_path+'/'+self.filename+'.hdf5 ./',
                        shell=True)

        # Get timing information:
        simulation_time_raw = ophsetup.GetTimesIC4( self.apocenter,
                                                    code2Myr=self.code2Myr)
        # Snapshot timescale
        snapshot_deltat = 0.2/self.code2Myr
        # Lengthen the simulation time by 12 snapshot intervals
        simulation_time = simulation_time_raw + (15*snapshot_deltat)
        # Begin writing snapshots 12 snapshot intervals before the raw
        # simulation time
        snapshot_t0 = simulation_time_raw - (12*snapshot_deltat)

        # Write the anchor file
        anchor_file = open('./anchors.txt', 'w')
        anchor_file.write('2\n')
        anchor_file.write('0.0000000 1.0\n')
        anchor_file.write(str(simulation_time)+' 1.0')
        anchor_file.close()

        # Prepare the simulation, first get all the parameters for the
        # potential.
        potparms = ophsetup.GenPotParmsIC4( output=True, code2kpc=self.code2kpc,
                                            code2Msol=self.code2Msol)
        # Unpack
        Md,Mb,Mhvir,Ad,Bd,Rc,Rs,conc,alpha = potparms

        # Now write the parameter file.
        params_in = open('./parms/autogadget_template.param','r')
        params_out = open('./autogadget.param','w')
        params_out.write('% James parameters\n')
        params_out.write('InitCondFile '+self.filename+'\n')
        params_out.write('TimeBetSnapshot '+str(snapshot_deltat)+'\n')
        params_out.write('TimeOfFirstSnapshot '+str(snapshot_t0)+'\n')
        params_out.write('\n% Potential parameters\n')
        params_out.write('MdExt '+str(Md)+'\n')
        params_out.write('MbExt '+str(Mb)+'\n')
        params_out.write('MhvirExt '+str(Mhvir)+'\n')
        params_out.write('aExt '+str(Ad)+'\n')
        params_out.write('bExt '+str(Bd)+'\n')
        params_out.write('cExt '+str(Rc)+'\n')
        params_out.write('dExt '+str(Rs)+'\n')
        params_out.write('concExt '+str(conc)+'\n')
        params_out.write('alphaExt '+str(alpha)+'\n')

        # Write the rest, which is the same:
        for line_in in params_in:
            params_out.write(line_in)
        ##li
        params_out.close()

        # Write the command to run GADGET-2
        gadget_command = 'mpirun -np '+str(nprocs) # Run mpi with nprocs processors
        gadget_command += ' ./source/Gadget2-MOD-'+os.environ['HOSTNAME'] # Gadget-2 binary
        gadget_command += ' ./autogadget.param' # GADGET-2 parameters
        gadget_command += ' >> log.txt 2>&1' # Write stdout/stderr to log.txt

        # Now run the simulation
        sim_timing0 = time.time()
        print('\nRunning GADGET-2...')
        print('Command is: '+gadget_command+'\n')
        # subprocess.call('rm -f log.txt')
        subprocess.call(gadget_command, shell=True)
        print('Done running GADGET-2, took '+str(round(time.time()-sim_timing0))+' seconds')

        # Put snaps in the snap_dir
        subprocess.call('mv out/snapshot*.hdf5 '+snap_path, shell=True)

        # Put parameters files away
        subprocess.call('mv ./anchors.txt '+storage_path+'/anchors_'+self.filename+'.txt', shell=True)
        subprocess.call('mv ./autogadget.param '+storage_path+'/autogadget_'+self.filename+'.param', shell=True)
        subprocess.call('mv ./log.txt '+storage_path+'/log_'+self.filename+'.txt', shell=True)

        # Cleanup the results of the recentering run.
        subprocess.call('rm -rf '+sim_path+'/out/*', shell=True)
        subprocess.call('rm -f '+sim_path+'/'+self.filename+'.hdf5', shell=True)
        subprocess.call('rm -f '+sim_path+'/autogadget.param-usedvalues', shell=True)

        # Go home
        os.chdir(home_path)

        print('\nDone Simulating for '+self.filename)
    #def

    def FindBestSnap(   self,
                        snap_dir='./snaps'):
        '''
        FindBestSnap:
            Get the snap that best matches the Ophiuchus stream data.

        Args:

        Returns:

        Outputs:

        '''

        print('\n===============================================')
        print('Finding best snap for '+self.filename)
        print('===============================================\n')

        # Handle directories. Append filename to the snap_dir directory.
        # Get the paths to the snap_dir and home directory.
        if snap_dir[-len(self.filename):] != self.filename:
            snap_dir += '/'+self.filename
        ##fi
        snap_path = os.path.abspath(snap_dir)
        home_path = os.getcwd()

        # Make sure directory infrastructure exists
        if os.path.isdir(snap_path) == False:
            sys.exit(snap_dir+' does not exist, must run the simulation first')
        ##fi

        # Get a list of the snaps and sort the list
        snap_files = glob.glob(snap_path+'/*.hdf5')
        ophmisc.sort_nicely(snap_files)

        # Declare the minimum and maximum galactic longitude that will be used
        # to represent the extent of the stream.
        stream_min_glong = 3.81 # lmax from Sesar+ 2016
        stream_max_glong = 5.85 # lmin from Sesar+ 2016

        # Declare the array that will be used to hold the number of particles
        # between the minimum and maximum galactic longitude.
        np_in_stream = np.zeros(len(snap_files))

        # Loop over all of the snapshots
        for i,snap_file in enumerate(snap_files):

            time0 = time.time()

            # Get the data
            print('Analyzing snap '+str(int(i))+' ...')
            data = Snapdata( os.path.abspath(snap_file) )

            # Get the galactocentric x,y, and z components
            x = data.x * self.code2kpc
            y = data.y * self.code2kpc
            z = data.z * self.code2kpc

            # Convert positions to l,b
            galactocen_coords = coords.Galactocentric(x=x * units.kpc,
                                                      y=y * units.kpc,
                                                      z=z * units.kpc)
            galactic_coords = galactocen_coords.transform_to(coords.Galactic)
            gl = np.array(galactic_coords.l)

            # Now find the number of particles between the minimum and
            # maximum galactic longitudes
            np_in_stream[i] = len( np.where(    ( gl > stream_min_glong ) &
                                                ( gl < stream_max_glong ) )[0] )

            del data,x,y,z,galactic_coords,galactocen_coords,gl

            # process = psutil.Process(os.getpid())
            # print( str(int(process.memory_info().vms / (1024)))+' kB' )
            # print( str( process.memory_percent() ) )

            print( str(time.time()-time0) )

        #for

        # Determine the snap with the most number of particles in the window.
        max_snap_ind = np.argmax( np_in_stream )
        max_snap_n = np_in_stream[ max_snap_ind ]

        # Find all snapshots within 5% of this maximum value.
        where_snap_close_max = np.where( np_in_stream > ( 0.95*max_snap_n ) )[0]
        avg_max_snap_ind = int(round(np.average( where_snap_close_max )))


        # Fill in the best snap
        self.matched_snap = snap_files[ avg_max_snap_ind ]
        matched_snap_n = avg_max_snap_ind
        matched_snap_out = open(snap_dir+'/'+self.filename+'_matched_snap.txt', 'w')
        matched_snap_out.write(snap_files[ matched_snap_n ])
        matched_snap_out.close()

        # Make a small plot of the number of particles in the stream.
        print('Plotting...')
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.scatter( np.arange(len(snap_files)), np_in_stream, color='Blue' )
        ax.scatter( matched_snap_n, np_in_stream[matched_snap_n], color='Red' )
        ax.axhline( 0.95*max_snap_n, linestyle='dashed', color='Black',
                    label='95% of Max' )
        ax.axvline( np.average( where_snap_close_max ), linestyle='dotted',
                    color='Black', label='Avg at 5% Max')
        leg = ax.legend(loc='lower left')
        leg.get_frame().set_alpha(0.2)
        ax.set_ylabel('N')
        ax.set_xlabel('Snapshot N')
        plt.savefig(snap_path+'/'+self.filename+'_snap_n_in_stream.pdf')
        fig.clf()
        plt.close()

        #pdb.set_trace()

    #def

    def Analyze(    self,
                    snapshot='Best',
                    output_dir='./stream_params',
                    snap_dir='./snaps',
                    n_mc=50,
                    background='uniform'):
        '''
        Analyze:
            Perform the present day analysis and comparison to the Bernard+14
            data.

        Args:
            snapshot (string) - The name of the snapshot to use. If 'Best' then
                use the self.matched_snap [Best]
            output_dir (string) - The directory in which to store the analysis
                outputs
            snap_dir (string) - The directory in which to find the snap file.
                If './snaps' then append self.filename to the path [./snaps]
            n_mc (int) - The number of times to do Monte-Carlo for errors.

        Returns:
            None

        Outputs:
            parameter file (.FIT and .txt) - A file with the output parameters.
            plots (.pdf) - Plots of the analysis.

        '''

        print('\n===============================================')
        print('Performing analysis for '+self.filename)
        print('===============================================\n')

        ############################################################################
        # Intro
        ############################################################################

        # Get the conversion factors for this progenitor.
        code2Myr = self.code2Myr
        code2kpc = self.code2kpc
        code2kms = self.code2kms
        code2Msol = self.code2Msol

        # Get the output absolute path, automatically add self.filename
        if output_dir[-len(self.filename):] != self.filename:
            output_dir += '/'+self.filename
        output_path = os.path.abspath(output_dir)

        # Make sure directory structures exist
        if os.path.isdir(output_path) == False:
            os.makedirs(output_path)
        ##fi

        ########################################################################
        # Hardcoded keywords
        ########################################################################

        cm = matplotlib.cm.get_cmap('Blues') # Colormap
        salpeter = False # Use a Salpeter IMF

        # Matplotlib parameters.
        mpl.rcParams['font.family'] = 'serif'
        mpl.rcParams['xtick.major.size'] = 6
        mpl.rcParams['ytick.major.size'] = 6
        fs = 18 # Fontsize

        # Background noise, in MSTO stars per square arcsecond.
        noiseval = 8.5E-5
        noiseval_deg = noiseval * 3600 * 3600 # Also per degree

        ############################################################################
        # Read data, perform calculations
        ############################################################################

        # Read
        if snapshot == 'Best':
            data = Snapdata(os.path.abspath(self.matched_snap))
        else:
            data = Snapdata(os.path.abspath(snapshot))
        n_part = data.npart
        p_mass = data.p_mass * code2Msol
        x = data.x * code2kpc
        y = data.y * code2kpc
        z = data.z * code2kpc
        vx = data.vx * code2kms
        vy = data.vy * code2kms
        vz = data.vz * code2kms
        pot = data.p
        time = data.time * code2Myr

        # Define the ratio of number of GADGET-2 particles to number of MSTO stars
        NtoN = 0.232 * p_mass
        if salpeter == True:
            NtoN = 0.131 * p_mass
        ##fi

        # Get line of sight velocity by dotting radial vector with velocity vector
        vrad = ophtasks.CalcGalacticVRad(x,y,z,vx,vy,vz)

        # Convert positions to l,b
        galactocen_coords = coords.Galactocentric(x=x * units.kpc,
                                                  y=y * units.kpc,
                                                  z=z * units.kpc)
        galactic_coords = galactocen_coords.transform_to(coords.Galactic)
        gl = np.array(galactic_coords.l)
        gb = np.array(galactic_coords.b)

        # Calculate distance to sun in kpc, manually
        dist = np.sqrt(np.square(x+8.3)+np.square(y)+np.square(z-0.027))

        # Get surface density for weighting the fits
        surfdens = ophtasks.CalcDensity2D(gl,gb,p_mass)
        # Inverse surface density weights low surface density (tails) higher
        # for the fit
        surfdens_weights = np.power(surfdens,-3)
        surfdens_argsort = np.argsort(surfdens)
        # Get only the 50%
        # where_low_surfdens = surfdens_argsort[:int(n_part/2)]

        # Find particles bound to the globular cluster
        # First get particle velocities with respect to the potential weighted
        # center of mass.
        vcm = np.sqrt(  np.square(data.vx-np.average(data.vx, weights=-data.p)) + \
                        np.square(data.vy-np.average(data.vy, weights=-data.p)) + \
                        np.square(data.vz-np.average(data.vz, weights=-data.p)) )
        ke = np.square(vcm)/2
        e_tot = data.p + ke
        where_unbound = np.where(e_tot > 0)[0]
        bound_core = False
        f_unbound = float(len(where_unbound)) / float(n_part)
        if f_unbound < 0.9:
            bound_core = True
        ##fi



        # Find particles that are more than two half-mass radii from the
        # potential-weighted center of the cluster
        # Only do this for streams where there is still a large bound component
        if bound_core == True:
            pwcx = np.average(x, weights=-pot)
            pwcy = np.average(y, weights=-pot)
            pwcz = np.average(z, weights=-pot)
            pwdist = np.sqrt( np.square(x-pwcx)+np.square(y-pwcy)+np.square(z-pwcz) )
            where_fit_selec = np.where( 1000*pwdist > 10*self.radius_3d )[0]
            # pdb.set_trace()
        ##fi
        else:
            # Just take everything
            where_fit_selec = np.where( x )[0]
        ##ie

        # pdb.set_trace()
        # print('Updated to low surface density fitting')

        # Get Ophiuchus stream star coordinates
        ophx,ophy,ophz = ophstream_props.GetGalCenXYZ()
        ophdist = np.sqrt(np.square(ophx+8.3)+np.square(ophy)+np.square(ophz-0.027))
        ophl,ophb = ophstream_props.Getlb()
        ophvrad = ophstream_props.GetVRad()
        n_oph = len(ophx)

        ################################################################################

        # Fit the stream with a weighted quadratic:

        # The quadratic function
        def quad_func(x,a,b,c):
            return a*np.power(x,2)+b*x+c

        # Weighted quadratic
        where_quad_fit = np.where( (gl[where_fit_selec] > 2) &
                                   (gl[where_fit_selec] < 8) &
                                   (gb[where_fit_selec] > 28) &
                                   (gb[where_fit_selec] < 33) )[0]
        # Common indicator of a simulation gone wrong.
        if len(where_quad_fit) == 0:
            sys.exit('Warning, no particles in quadratic window! Exiting...')
        ##fi
        wquad_popt,_ = curve_fit(quad_func,gl[where_fit_selec][where_quad_fit],
                                gb[where_fit_selec][where_quad_fit],
                                sigma=np.power(pot,-2)[where_fit_selec][where_quad_fit],
                                p0=[-0.15,1.0,30],
                                bounds=([-0.15,0.5,28],[-0.09,1.0,32]))
        wquada, wquadb, wquadc = wquad_popt
        # if bound_core == False:
        #     print('a, b, c = '+str(wquada)+', '+str(wquadb)+', '+str(wquadc))
        print('N unbound = '+str(len(where_unbound)))
        print('frac unbound = '+str(f_unbound))
        print('Bound Core = '+str(bound_core))
        #print('N in tail = '+str(len(where_fit_selec)))
        wquadgl = np.arange(0,10,0.1)
        wquadgb = quad_func(wquadgl,wquada,wquadb,wquadc)

        # if bound_core == True:
        #     gl = gl[where_fit_selec]
        #     gb = gb[where_fit_selec]
        #     dist = dist[where_fit_selec]
        #     vrad = vrad[where_fit_selec]
        #     n_part = len(where_fit_selec)
        # ##fi

        #pdb.set_trace()

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
        close_l_oph = np.zeros(n_oph)
        for i in range(n_part):
            close_l[i] = newton(stream_quad_dot, gl[i], tol=0.001, maxiter=500,
                                    args=(gl[i],gb[i],wquada,wquadb,wquadc))
        ###i
        for i in range(n_oph):
            close_l_oph[i] = newton(stream_quad_dot, ophl[i], tol=0.001, maxiter=500,
                                args=(ophl[i],ophb[i],wquada,wquadb,wquadc))
        ###i
        close_b = wquada*np.square(close_l) + wquadb*close_l + wquadc
        close_b_oph = wquada*np.square(close_l_oph) + wquadb*close_l_oph + wquadc

        # Latitude in the new coordinate system is the distance from the curve:
        stream_lat = np.sqrt(np.square(close_l-gl)+np.square(close_b-gb))
        stream_lat_oph = np.sqrt(np.square(close_l_oph-ophl)+np.square(close_b_oph-ophb))

        # Get the sign of the latitude by seeing if the new point is lower than the
        # curve.
        for i in range(n_part):
            if gb[i] < ( wquada*(gl[i])**2 + wquadb*gl[i] + wquadc ):
                stream_lat[i] *= -1
        ###i
        for i in range(n_oph):
            if ophb[i] < ( wquada*(ophl[i])**2 + wquadb*ophl[i] + wquadc ):
                stream_lat_oph[i] *= -1
        ###i

        # Longitude in the new coordinate system is the arc length along the curve
        adjusted_gl = np.array(gl)
        adjusted_gl[np.where(adjusted_gl>180)[0]] = adjusted_gl[np.where(adjusted_gl>180)[0]]-360
        initial_long = np.median(adjusted_gl)
        stream_long = np.zeros(n_part)
        stream_long_oph = np.zeros(n_oph)
        arclen_args = (wquada,wquadb)
        for i in range(n_part):
            stream_long[i] = quad_integrate(stream_quad_arc, initial_long, gl[i],
                                            args=arclen_args )[0]
        ###i
        for i in range(n_oph):
            stream_long_oph[i] = quad_integrate(stream_quad_arc, initial_long, ophl[i],
                                            args=arclen_args )[0]
        ###is

        ################################################################################
        # Plot images of streams, in (l,b), (L,B) and sigma.
        ################################################################################

        # Make the output files.
        outplots = PdfPages(output_path+'/plots_'+self.filename+'.pdf')
        stat_file = open(output_path+'/stream_params_'+self.filename+'.txt','w')

        # Make the figure.
        fig = plt.figure(figsize=(8,8))

        # Choose the center of the distribution.
        center_l = np.average(ophl)
        center_b = np.average(ophb)

        # Set the width of the frame.
        ax_w = 5 #Degrees.
        n_bins = 100
        gb_lim_lo = 25
        gb_lim_hi = 35
        gl_lim_lo = 0
        gl_lim_hi = 10

        # Adjust area for the different latitude of each bin.
        bin_area = ((ax_w*2*3600)/n_bins)**2 # Area in arcseconds.
        bin_lats = (np.arange(0,1,1/float(n_bins))*ax_w*2) + gb_lim_lo # Latitudes.
        bin_cos = np.cos( np.radians(bin_lats) ) # Cosine of the latitudes.
        bin_area = np.multiply(bin_area,bin_cos) # Gives a solid angle dW = cosbdldb

        # Ranges, not prettiest right now, but works with rotation.
        ax_range = [[gb_lim_lo,gb_lim_hi],
                    [gl_lim_lo,gl_lim_hi]]
        ax_ext = [gl_lim_hi,gl_lim_lo,
                    gb_lim_lo,gb_lim_hi]

        # Make the image.
        ax = fig.add_subplot(111)
        if bound_core == True:
                # pdb.set_trace()
                ax_hist,_,_ = np.histogram2d(   gb[where_fit_selec],
                                                gl[where_fit_selec],
                                                bins=n_bins, range=ax_range)
        else:
                ax_hist,_,_ = np.histogram2d(gb, gl, bins=n_bins, range=ax_range)
        ##ie
        ax_hist = np.rot90(ax_hist,k=2) # Not pretty but it works!
        ax_hist = np.divide(ax_hist,bin_area) # Convert to surface density.
        ax_hist *= NtoN # Turn to Number of MSTO stars.
        ax_im = ax.imshow(ax_hist, cmap=cm, norm=LogNorm(vmin=1E-5, vmax=1E-3),
                            extent=ax_ext, interpolation='nearest')
        divider_ax = make_axes_locatable(ax)
        cax = divider_ax.append_axes('right', size='5%', pad=0.05)
        ax_cbar = fig.colorbar(ax_im, cax=cax)

        # Mask out values below noise.
        ax_hist_belownoise = np.where((ax_hist > 0) & (ax_hist < noiseval))
        ax_mask = np.zeros((100,100,4))
        ax_mask[:,:,3] = 0.0
        ax_mask[ax_hist_belownoise[0],ax_hist_belownoise[1],:3] = 128 # Grey triplet
        ax_mask[ax_hist_belownoise[0],ax_hist_belownoise[1],3] = 1.0

        # Mask out the values if requested.
        if background == 'mask':
            ax_mask_im = ax.imshow(ax_mask, extent=ax_ext, interpolation='nearest')
        ##fi

        # Add a uniform noise if requested.
        if background == 'uniform':
            ax_bg = np.random.normal(loc=8E-5, scale=1.5E-5, size=ax_hist.shape)
            ax_hist += ax_bg
            ax_bg_im = ax.imshow(ax_hist, cmap=cm, norm=LogNorm(vmin=5E-5, vmax=5E-4),
                                    extent=ax_ext, interpolation='nearest')
        ##fi

        # Label the plot.
        ax_cbar.set_label(r'Surface density ($N_{MSTO}/\prime\prime ^{2}$)', fontsize=16)
        ax.set_xlabel(r' $l$ [deg]', fontsize=fs)
        ax.set_ylabel(r' $b$ [deg]', fontsize=fs)
        ax.scatter(ophl,ophb,edgecolor='Red',facecolor='None',s=40,label='Stream Members')
        ax.plot(wquadgl,wquadgb,'-r')
        ax.set_xlim(8,2)
        ax.set_ylim(28.5,32.5)

        # Denote if probably still bound
        if bound_core == True:
            ax.annotate('This be farkd', xy=(0.2,0.8), xycoords='axes fraction', color='red')
        ##fi

        plt.savefig(outplots, format='pdf')
        plt.clf()
        dispersion_weights = ax_hist

        ################################################################################

        # Plot stream in spline coordinates.

        # Set the width of the frame
        ax_w = 5 #Degrees
        n_bins = 100
        gb_lim_lo = -5
        gb_lim_hi = 5
        gl_lim_lo = -5
        gl_lim_hi = 5

        # Adjust area for the different latitude of each bin.
        bin_area = ((ax_w*2*3600)/n_bins)**2 # Area in arcseconds
        fixed_bin_area = bin_area # Save for later
        bin_lats = (np.arange(0,1,1/float(n_bins))*ax_w*2) + gb_lim_lo # Latitudes.
        bin_cos = np.cos( np.radians(bin_lats) ) # Cosine of the latitudes.
        bin_area = np.multiply(bin_area,bin_cos) # Gives a solid angle dW = cosbdldb

        # Ranges, not pretty right now but it works.
        ax_range = [[gb_lim_lo,gb_lim_hi],
                    [gl_lim_lo,gl_lim_hi]]
        ax_ext = [gl_lim_hi,gl_lim_lo,
                    gb_lim_lo,gb_lim_hi]

        # Make the image.
        ax = fig.add_subplot(111)
        ax_hist,_,_ = np.histogram2d(stream_lat, stream_long, bins=n_bins, range=ax_range)
        ax_hist_partcounts = ax_hist * NtoN # Useful for later.
        ax_hist = np.rot90(ax_hist,k=2)
        ax_hist = np.divide(ax_hist,bin_area)
        ax_hist *= NtoN
        ax_im = ax.imshow(ax_hist, cmap=cm, norm=LogNorm(vmin=1E-6,
                   vmax=1E-3), extent=ax_ext, interpolation='nearest')
        divider_ax = make_axes_locatable(ax)
        cax = divider_ax.append_axes('right', size='5%', pad=0.05)
        ax_cbar = fig.colorbar(ax_im, cax=cax)

        # Mask out values below noise.
        ax_im_spline_belownoise = np.where((ax_hist > 0) & (ax_hist < noiseval))
        ax_spline_mask = np.zeros((100,100,4))
        ax_spline_mask[:,:,3] = 0.0
        ax_spline_mask[ax_im_spline_belownoise[0],ax_im_spline_belownoise[1],:3] = 128
        ax_spline_mask[ax_im_spline_belownoise[0],ax_im_spline_belownoise[1],3] = 1.0
        ax_spline_mask_im = ax.imshow(ax_spline_mask, extent=ax_ext, interpolation='nearest')

        # Label the figure.
        ax_cbar.set_label(r'Surface density ($N_{MSTO}/\prime\prime ^{2}$)', fontsize=16)
        ax.set_xlabel(r' $\Lambda$ [deg]', fontsize=fs)
        ax.set_ylabel(r' $B$ [deg]', fontsize=fs)
        ax.scatter(stream_long_oph,stream_lat_oph,c='Red',s=40)
        ax.axhline(0,color='r')
        ax.set_xlim(3,-3)
        ax.set_ylim(-2,2)
        plt.savefig(outplots, format='pdf')
        plt.clf()

        # Evaluate the total number of stars in this projection. Useful later,
        # in number of MSTO stars!
        N_noise_per_pixel = fixed_bin_area * noiseval
        where_above_noise =  np.where(ax_hist_partcounts-N_noise_per_pixel > 0)
        N_above_noise = np.sum( ax_hist_partcounts[where_above_noise] - N_noise_per_pixel )
        guess_sdmax = np.amax(ax_hist)

        ################################################################################

        # Plot the velocity dispersion.

        # Set the width of the frame
        ax_w = 5 #Degrees
        n_bins = 100
        gb_lim_lo = 25
        gb_lim_hi = 35
        gl_lim_lo = 0
        gl_lim_hi = 10

        # Ranges, not pretty right now but it works.
        ax_range = [[gb_lim_lo,gb_lim_hi],
                    [gl_lim_lo,gl_lim_hi]]
        ax_ext = [gl_lim_hi,gl_lim_lo,
                    gb_lim_lo,gb_lim_hi]

        # Make the image.
        ax = fig.add_subplot(111)
        ax_hist,_,_,_ = bin2d(gb, gl, vrad, statistic=np.std, bins=n_bins, range=ax_range)
        ax_hist = np.rot90(ax_hist,k=2)
        ax_im = ax.imshow(ax_hist, cmap=matplotlib.cm.get_cmap('Reds'), vmin=0.1,
                        vmax=1.5, extent=ax_ext, interpolation='nearest')
        divider_ax = make_axes_locatable(ax)
        cax = divider_ax.append_axes('right', size='5%', pad=0.05)
        ax_cbar = fig.colorbar(ax_im, cax=cax)

        # Mask out bad values, same as before.
        ax_mask_im = ax.imshow(ax_mask, extent=ax_ext, interpolation='nearest')

        # Label the figure.
        ax_cbar.set_label('Velocity Dispersion (km/s)', fontsize=16)
        ax.set_xlabel(r' $l$ [deg]', fontsize=fs)
        ax.set_ylabel(r' $b$ [deg]', fontsize=fs)
        ax.set_xlim(8,2)
        ax.set_ylim(28.5,32.5)
        ax.plot(wquadgl,wquadgb,'-r')
        ax.scatter(ophl,ophb,c='DodgerBlue',s=40)
        plt.savefig(outplots, format='pdf')
        plt.clf()

        # Measure the velocity dispersion of the stream.
        masked_vdisp = np.ma.masked_invalid(ax_hist)
        masked_vdisp[ax_hist_belownoise] = np.ma.masked
        stream_vdisp = np.ma.average(masked_vdisp, weights=dispersion_weights)
        stream_vdisp_err = np.sqrt(np.ma.average(np.square(masked_vdisp-stream_vdisp), weights=dispersion_weights))
        #stat_file.write('Total dispersion: '+str(stream_vdisp)+'\n\n')

        ################################################################################
        # Analyze and plot the length and width of the stream.
        ################################################################################

        # We will require the nominal noise for the measured properties in order
        # to decide how long we need to run the Monte-Carlo evaluator for
        lsum_err_obs = 38
        wsum_err_obs = 43
        len_err_obs = 0.25
        fwhm_err_obs = 0.8

        # Decide the size of chunks that will be evalauted through the
        # Monte-Carlo method.
        mc_chunk = 10

        ### First length

        # Noise
        lnoise = np.sqrt(50)

        # Spatial density histograms for length of stream. Make the axis
        ax = fig.add_subplot(111)

        # Include latitudes -0.1 < B < 0.1 arcminutes, following Bernard+14.
        l_in_stream = np.where((stream_lat > -0.1) & (stream_lat < 0.1))[0]
        stream_longl = stream_long[l_in_stream]

        # Do Monte-Carlo evaluation, create the arrays to hold data.
        llen_mc = np.array([]) # Empty, to be filled
        lsum_mc = np.array([])
        llen_mc_temp = np.empty(mc_chunk) # Chunk sized
        lsum_mc_temp = np.empty(mc_chunk)

        # Create variables to track if each measurable has reached it's noise
        # requirement.
        llen_noise_good = False
        lsum_noise_good = False

        # While loop if either of the measureables have not reached their
        # noise requirements
        while (llen_noise_good == False) or (lsum_noise_good == False):

            # Reset the counter for the Monte-Carlo evaluation to make sure we
            # get mc_chunk good simulations
            mc_counter = 0

            # Loop over 10 evaluations of the Monte-Carlo loop
            while mc_counter < mc_chunk:

                # try-catch every Exception (doesn't include SystemExit or
                # KeyboardInterrupt) for ease so it doesn't crash due to weird
                # measurement errors:
                try:

                    # Make the length histogram.
                    lhist_bins = 80
                    lhist, lbin_edges = np.histogram(stream_longl,
                                                        bins=lhist_bins,
                                                        range=(-10,10) )
                    lbin_centers = lbin_edges + (lbin_edges[1]-lbin_edges[0])/2

                    # Add noise.
                    lhist = lhist.astype(float)
                    lhist *= NtoN # Convert from GADGET-2 to MSTO particles
                    lhist_n = lhist + np.random.normal(0,np.sqrt(50),lhist_bins)
                    lhist_n += np.linspace(20,80,lhist_bins) # background decreases with Long

                    # In order to get the length start at the peak and crawl left and right
                    # until you drop below the noise.
                    start_lpeak = np.argmax(lhist)
                    where_lgtn = np.array([start_lpeak])
                    stop_lr = 0
                    stop_ll = 0
                    lgtn_ind = 0
                    # Crawl rightl
                    while stop_lr == 0:
                        lgtn_ind += 1
                        # Include a breakout condition:
                        if (start_lpeak + lgtn_ind) == 79:
                            stop_lr = 1
                            break
                        ##fi
                        llim = 2*lbin_centers[start_lpeak+lgtn_ind]+50+lnoise
                        if lhist_n[start_lpeak + lgtn_ind] > llim:
                            where_lgtn = np.append(where_lgtn,
                                                    start_lpeak+lgtn_ind)
                        ##fi
                        else: stop_lr = 1
                    ##wh
                    lgtn_ind = 0
                    # Crawl_left
                    while stop_ll == 0:
                        lgtn_ind += 1
                        if (start_lpeak - lgtn_ind) == 0:
                            stop_ll = 1
                            break
                        ##fi
                        llim = 2*lbin_centers[start_lpeak-lgtn_ind]+50+lnoise
                        if lhist_n[start_lpeak-lgtn_ind] > llim:
                            where_lgtn = np.append(where_lgtn,
                                                    start_lpeak-lgtn_ind)
                        else: stop_ll = 1
                    ##wh
                    where_lgtn = np.sort(where_lgtn)

                    # Calculate the length and the sum:
                    llen_temp = max(lbin_centers[where_lgtn])-min(lbin_centers[where_lgtn])
                    lsum_excess_temp = len(where_lgtn)*(2*np.average(lbin_centers[where_lgtn])+50)
                    lsum_temp = np.sum(lhist_n[where_lgtn])-lsum_excess_temp

                    # Make sure the errors and all values are real:
                    if np.isfinite(llen_temp) == False:
                        continue
                    ##fi
                    if np.isfinite(lsum_temp) == False:
                        continue
                    ##fi

                # If we find an error continue to the next loop iteration
                except Exception:
                    continue
                ##te

                # If it worked add to the temporary chunk array
                llen_mc_temp[mc_counter] = llen_temp
                lsum_mc_temp[mc_counter] = lsum_temp

                # Tick up the Monte-Carlo counter
                mc_counter += 1
            ##wh

            # Now that we have mc_chunk new good measurements add them to the
            # existing array and see if the reduced errors meet the criterion:

            if llen_noise_good == False:
                llen_mc = np.concatenate((llen_mc, llen_mc_temp))
                if (np.std(llen_mc) / np.sqrt(len(llen_mc)-1) ) < (0.1*len_err_obs):
                    llen_noise_good = True
                ##fi
            ##fi

            if lsum_noise_good == False:
                lsum_mc = np.concatenate((lsum_mc, lsum_mc_temp))
                if (np.std(lsum_mc) / np.sqrt(len(lsum_mc)-1) ) < (0.1*lsum_err_obs):
                    lsum_noise_good = True
                ##fi
            ##fi
        ##wh

        # Average the Monte Carlo results and get some errors.
        llen = np.average(llen_mc)
        llen_err = np.std(llen_mc) / np.sqrt(len(llen_mc)-1) # Reduced Error!
        lsum = np.average(lsum_mc)
        lsum_err = np.std(lsum_mc) / np.sqrt(len(lsum_mc)-1) # Reduced Error!

        # Plot.
        ax.step(lbin_edges[:-1], lhist_n, where='post', color='k')
        ax.plot((lbin_edges[-1],lbin_edges[-2]), (lhist_n[-1],lhist_n[-1]), '-k')
        ax.step(lbin_edges[where_lgtn+1], lhist_n[where_lgtn+1], where='post', color='r',
                label='In Stream')
        ax.step(lbin_edges[where_lgtn-1], lhist_n[where_lgtn-1], where='post', color='r')
        ax.step(lbin_edges[:-1], lhist, where='post', color='b', label='Stream')
        ax.plot((lbin_edges[-1],lbin_edges[-2]), (lhist[-1],lhist[-1]), '-b')
        ax.set_xlabel(r' $\Lambda$ [deg]', fontsize=fs)
        ax.set_ylabel('$N_{MSTO}$', fontsize=fs)
        ax.set_xlim(5,-5)
        ax.legend()
        plt.savefig(outplots, format='pdf')
        plt.clf()

        ################################################################################

        # Same for width:

        # Noise
        wback = 35
        wnoise = np.sqrt(wback)

        # Define the gaussian function
        def gaussian_func(x, a, x0, sigma, d):
            return a*np.exp(-(x-x0)**2/(2*sigma**2))+d
        #def

        # Make spatial density histograms for length of stream.
        ax = fig.add_subplot(111)

        # Include longitudes -1 < L < 1 degrees, following Bernard+14.
        w_in_stream = np.where((stream_long > -1) & (stream_long < 1))[0]
        stream_latw_am = stream_lat[w_in_stream] * 60 # Latitude in arcminutes
        whist_bins = 133 # Number of bins for the histogram

        # Do Monte-Carlo, create the arrays to hold data.
        wfwhm_mc = np.array([]) # Empty to be filled
        wfwhm_err_fit_mc = np.array([])
        wsum_mc = np.array([])
        wfwhm_mc_temp = np.empty(mc_chunk) # Temporary chunk size
        wsum_mc_temp = np.empty(mc_chunk)
        wfwhm_err_fit_mc_temp = np.empty(mc_chunk)

        # Create variables to track if each measurable has reached it's noise
        # requirement.
        wfwhm_noise_good = False
        wsum_noise_good = False

        # While loop if either of the measureables have not reached their
        # noise requirements
        while (wfwhm_noise_good == False) or (wsum_noise_good == False):

            # Reset the counter for the Monte-Carlo evaluation to make sure we
            # get mc_chunk good simulations
            mc_counter = 0

            # Loop over 10 evaluations of the Monte-Carlo loop
            while mc_counter < mc_chunk:

                # try-catch every Exception (doesn't include SystemExit or
                # KeyboardInterrupt) for ease so it doesn't crash due to weird
                # measurement errors:
                try:

                    whist, wbin_edges = np.histogram(stream_latw_am,
                                                        bins=whist_bins,
                                                        range=(-60,60))
                    wbin_centers = wbin_edges + (wbin_edges[1]-wbin_edges[0])/2

                    # Add noise, values follow Bernard+14
                    whist = whist.astype(float)
                    whist *= NtoN # Convert to MSTO stars.
                    whist_n = whist + np.random.normal(0,wnoise,whist_bins)
                    whist_n += wback

                    # Fit Gaussian to the width component.
                    wpopt, wcov = curve_fit(gaussian_func, wbin_centers[:-1],
                                                whist_n, p0=(100,0,20,wback))
                    wpopt[2] = np.absolute(wpopt[2]) # Gaussian allows +/-, force +
                    werr = np.sqrt(np.diag(wcov))
                    wgaussx = np.linspace(-60,60,num=1000)
                    wgaussy = gaussian_func(wgaussx, wpopt[0], wpopt[1],
                                                wpopt[2], wpopt[3])
                    wfwhm_temp = wpopt[2]*2.35482
                    wfwhm_err_fit_temp = werr[2]*2.35482

                    # Determine the Full Width at three times the Noise, then integrate the
                    # stream members between those limits.
                    wfwn_r = wpopt[1]+wpopt[2]*np.sqrt(2*np.log(wpopt[0]/(np.sqrt(35))))
                    wfwn_l = wpopt[1]-wpopt[2]*np.sqrt(2*np.log(wpopt[0]/(np.sqrt(35))))
                    where_wfwn = np.where(  ( wbin_centers > wfwn_l ) &
                                            ( wbin_centers < wfwn_r ) )[0]

                    wsum_temp = np.sum(whist_n[where_wfwn])-(len(where_wfwn)*wback )

                    # Make sure the fitting error is smaller than twice the
                    # observational error, otherwise continue
                    if wfwhm_err_fit_temp > 2*fwhm_err_obs:
                        continue
                    ##fi

                    # Make sure the errors and all values are real:
                    if np.isfinite(wfwhm_err_fit_temp) == False:
                        continue
                    ##fi
                    if np.isfinite(wfwhm_temp) == False:
                        continue
                    ##fi
                    if np.isfinite(wsum_temp) == False:
                        continue
                    ##fi

                except Exception:
                    continue
                ##te

                # If it worked add to the array
                wfwhm_mc_temp[mc_counter] = wfwhm_temp
                wsum_mc_temp[mc_counter] = wsum_temp
                wfwhm_err_fit_mc_temp[mc_counter] = wfwhm_err_fit_temp

                # Tick up the Monte-Carlo counter
                mc_counter += 1
            ##wh

            # Now that we have mc_chunk new good measurements add them to the
            # existing array and see if the reduced errors meet the criterion:
            if wfwhm_noise_good == False:
                wfwhm_mc = np.concatenate((wfwhm_mc, wfwhm_mc_temp))
                wfwhm_err_fit_mc = np.concatenate((wfwhm_err_fit_mc, wfwhm_err_fit_mc_temp))
                if (np.std(wfwhm_mc) / np.sqrt(len(wfwhm_mc)-1) ) < (0.1*fwhm_err_obs):
                    wfwhm_noise_good = True
                ##fi
            ##fi

            if wsum_noise_good == False:
                wsum_mc = np.concatenate((wsum_mc, wsum_mc_temp))
                if (np.std(wsum_mc) / np.sqrt(len(wsum_mc)-1) ) < (0.1*wsum_err_obs):
                    wsum_noise_good = True
                ##fi
            ##fi
        ##wh

        # Compute errors
        wfwhm_err = np.std(wfwhm_mc) / np.sqrt(len(wfwhm_mc)-1)
        wfwhm_err_tot = np.sqrt(np.square(np.average(wfwhm_err_fit_mc))+np.square(wfwhm_err)) # Add fitting error
        wsum_err = np.std(wsum_mc) / np.sqrt(len(wsum_mc)-1)

        # Compute properties of the width histogram, In arcminutes!
        wfwhm = np.average(wfwhm_mc)
        wsum = np.average(wsum_mc)
        wsigma = (wfwhm / 2.35482)
        wmean = np.average(stream_latw_am)

        # Plot
        ax.step(wbin_edges[:-1], whist_n, where='post', color='k')
        ax.plot((wbin_edges[-1],wbin_edges[-2]), (whist_n[-1],whist_n[-1]), '-k')
        ax.step(wbin_edges[where_wfwn], whist_n[where_wfwn], where='post', color='r',
                label='In Stream')
        ax.step(wbin_edges[:-1], whist, where='post', color='b', label='Stream')
        ax.plot((wbin_edges[-1],wbin_edges[-2]), (whist[-1],whist[-1]), '-b')
        ax.plot(wgaussx, wgaussy, '-b')
        ax.set_xlabel(r' $B$ [arcmin]', fontsize=fs)
        ax.set_ylabel(r'$N_{MSTO}$', fontsize=fs)
        ax.set_xlim(-60,60)
        ax.legend()
        plt.savefig(outplots, format='pdf')
        plt.clf()

        # Also look at all the stars in the MSTO window, in the right coordinates.
        in_iso_window = np.where(   (stream_lat > -4.5/60) &
                                    (stream_lat < 4.5/60) &
                                    (stream_long < 1) &
                                    (stream_long > -1) )[0]
        n_msto_iso_window = len(in_iso_window) *NtoN

        ################################################################################

        ### Kinematics

        # Plot radial velocity against longitude
        ax1 = fig.add_subplot(211)
        pts = ax1.scatter(gl,vrad,s=1,c='k',alpha=0.5)
        pts.set_rasterized(True)
        ax1.scatter(ophl,ophvrad,s=40,c='DodgerBlue')
        ax1.set_ylabel(r'$V_{helio}$ [km/s]', fontsize=fs)
        ax1.set_xlabel(r' $l$ [deg]', fontsize=fs)
        ax1.set_xlim(10,0)
        ax1.set_ylim(267,302)
        ax1.tick_params(labelbottom='off')

        # Plot distance against longitude
        ax2 = fig.add_subplot(212)
        pts = ax2.scatter(gl,dist,s=1,c='k',alpha=0.5)
        pts.set_rasterized(True)
        ax2.scatter(ophl,ophdist,s=40,c='DodgerBlue')
        ax2.set_ylabel('Distance [kpc]', fontsize=fs)
        ax2.set_xlabel(r' $l$ [deg]', fontsize=fs)
        ax2.set_xlim(10,0)
        ax2.set_ylim(6.7,10.3)
        fig.subplots_adjust(hspace=0)
        plt.savefig(outplots, format='pdf', dpi=300)

        ################################################################################

        # Close matplotlib
        fig.clf()
        plt.close()

        stat_file.write('# Length + err, length sum + err, width FWHM + err, width sum + err, Velocity dispersion\n\n')
        #stat_file.write('# Central density = '+str(newt_sdmax)+', warned = '+str(newton_warning)+'\n')
        #stat_file.write(str(newt_n_total)+'\n')
        stat_file.write(str(llen)+'\n')
        stat_file.write(str(llen_err)+'\n\n')
        stat_file.write(str(lsum)+'\n')
        stat_file.write(str(lsum_err)+'\n\n')
        stat_file.write(str(wfwhm)+'\n')
        stat_file.write(str(wfwhm_err)+'\n\n')
        stat_file.write(str(wsum)+'\n')
        stat_file.write(str(wsum_err)+'\n\n')
        stat_file.write(str(stream_vdisp)+'\n')
        stat_file.write(str(stream_vdisp_err))

        outplots.close()

        ################################################################################

        # Now write output file containing all information

        out_names = (   "Name",
                        "Mass",
                        "Radius_3D",
                        "LSum",
                        "LSum_Err",
                        "WSum",
                        "WSum_Err",
                        "Length",
                        "Length_Err",
                        "Width",
                        "Width_Err",
                        "Width_Err_Rand",
                        "Sigma",
                        "Sigma_Err",
                        "Frac_Unbound",
                        "Bound_Core")

        out_data_dtypes = (  "str","float","float","float","float","float","float",
                            "float","float","float","float","float","float",
                            "float","float","str")

        out_cols = np.array([   self.filename,
                                self.mass,
                                self.radius_3d,
                                lsum,
                                lsum_err,
                                wsum,
                                wsum_err,
                                llen,
                                llen_err,
                                wfwhm,
                                wfwhm_err_tot,
                                wfwhm_err,
                                stream_vdisp,
                                stream_vdisp_err,
                                f_unbound,
                                bound_core])

        # for i,val in enumerate(out_cols):
        #     out_cols = np.array([val])
        # #for
        # pdb.set_trace()
        out_data = table.Table( out_cols,
                                names = out_names,
                                dtype = out_data_dtypes)
        out_data.write(output_path+'/stream_params_'+self.filename+'.FIT',
                            format='fits', overwrite=True)
    #def
#cls
