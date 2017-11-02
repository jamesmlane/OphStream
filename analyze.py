# ----------------------------------------------------------------------------
#
# TITLE - analyze.py
# AUTHOR - James Lane
# PROJECT - Ophiuchus Stream
# CONTENTS:
#	1. PresentDay
#
# ----------------------------------------------------------------------------
#
# Docstrings and metadata:
'''
Primary analysis code for examining simulations of the Ophiuchus stream
as seen at the present day.
'''
__author__ = "James Lane"

#Imports
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm
import matplotlib.cm
from astropy import coordinates as coords
from astropy import units
from scipy.stats import binned_statistic_2d as bin2d
from scipy.stats import linregress
from scipy.optimize import curve_fit,newton
from scipy.integrate import quad as quad_integrate
from scipy.integrate import dblquad as quad_integrate2d
import scipy.interpolate
from snapdata import Snapdata
import ophstream.units
import ophstream.analyze
import ophstream.stream_properties
import ophstream.misc
import ophstream.tasks
import pdb
import os
import sys
import glob

def PresentDay(run_num,unitsys,choose_snap=False,background='uniform'):
    '''
    PresentDay:

    Perform the present day analysis and comparison to the Bernard+14 data.

    Args:
        run_num (string): The run identification number.
        unitsys (string): The system of units used in the simulation.
        choose_snap (Boolean): Use a selected snapshot which is found in the
            present directory [False]
        background (string): The type of noise to add (background or mask)
    '''

    ############################################################################
    # Keywords
    ############################################################################

    cm = matplotlib.cm.get_cmap('Blues')
    salpeter = False

    # Matplotlib parameters.
    mpl.rcParams['font.family'] = 'serif'
    mpl.rcParams['xtick.major.size'] = 6
    mpl.rcParams['ytick.major.size'] = 6
    fs = 18 # Fontsize

    # number of MC runs.
    n_mc = 50

    # Background noise, in MSTO stars per square arcsecond.
    noiseval = 8.5E-5
    noiseval_deg = noiseval * 3600 * 3600 # Also per degree

    ############################################################################
    # Read data, perform calculations
    ############################################################################

    #Set units
    code2Myr,code2kpc,code2kms,code2Msol = ophstream.units.GetConversions(unitsys)

    # Get a list of hdf5 output files
    if choose_snap == True:
        snap_files = sorted(glob.glob('./*.hdf5'))
    else:
        snap_files = sorted(glob.glob('../*.hdf5'))
    n_snap = -1
    ophstream.misc.sort_nicely(snap_files)

    # Read in the specified file:
    abs_path = os.path.abspath(snap_files[n_snap])
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

    # Define the ratio of number of GADGET-2 particles to number of MSTO stars
    NtoN = 0.232 * p_mass
    if salpeter == True:
        NtoN = 0.131 * p_mass
    ##fi

    # Get line of sight velocity by dotting radial vector with velocity vector
    vrad = ophstream.tasks.CalcGalacticVRad(x,y,z,vx,vy,vz)

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
    surfdens = ophstream.tasks.CalcDensity2D(gl,gb,p_mass)
    surfdens_weights = np.power(surfdens,-1/3)

    # Get Ophiuchus stream star coordinates
    ophx,ophy,ophz = ophstream.stream_properties.GetGalCenXYZ()
    ophdist = np.sqrt(np.square(ophx+8.3)+np.square(ophy)+np.square(ophz-0.027))
    ophl,ophb = ophstream.stream_properties.Getlb()
    ophvrad = ophstream.stream_properties.GetVRad()
    n_oph = len(ophx)

    ################################################################################

    # Fit the stream with a weighted quadratic:

    # The quadratic function
    def quad_func(x,a,b,c):
        return a*np.power(x,2)+b*x+c

    # Weighted quadratic
    where_quad_fit = np.where( (gl > 2) &
                               (gl < 8) &
                               (gb > 28) &
                               (gb < 33) )[0]
    # Common indicator of a simulation gone wrong.
    if len(where_quad_fit) == 0:
        sys.exit('Warning, no particles in quadratic window! Exiting...')
    ##fi
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
    close_l_oph = np.zeros(n_oph)
    for i in range(n_part):
        close_l[i] = newton(stream_quad_dot, gl[i], tol=0.001, maxiter=100,
                            args=(gl[i],gb[i],wquada,wquadb,wquadc))
    ###i
    for i in range(n_oph):
        close_l_oph[i] = newton(stream_quad_dot, ophl[i], tol=0.001, maxiter=100,
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
    initial_long = np.median(gl)
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

    # Create output files.
    outplots = PdfPages('run_'+run_num+'_plots_gallb.pdf')
    fig = plt.figure(figsize=(8,8))
    stat_file = open('run_'+run_num+'_stream_params.txt','w')

    # Plot stream in galactic coordinates.

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
    ax_hist,_,_ = np.histogram2d(gb, gl, bins=n_bins, range=ax_range)
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
        ax_bg_im = ax.imshow(ax_hist, cmap=cm, norm=LogNorm(vmin=1E-5, vmax=1E-3),
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

    # Evaluate the total number of stars in this projection. Useful later, in
    # number of MSTO stars!
    N_noise_per_pixel = fixed_bin_area * noiseval
    where_above_noise =  np.where(ax_hist_partcounts-N_noise_per_pixel > 0)
    N_above_noise = np.sum( ax_hist_partcounts[where_above_noise] - N_noise_per_pixel )
    guess_sdmax = np.amax(ax_hist)

    #pdb.set_trace()

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
    #stat_file.write('Total dispersion: '+str(stream_vdisp)+'\n\n')

    ################################################################################
    # Analyze and plot the length and width of the stream.
    ################################################################################

    # Noise
    lnoise = np.sqrt(50)

    # Make spatial density histograms for length of stream.
    ax = fig.add_subplot(111)

    # Include latitudes -0.1 < B < 0.1 arcminutes, following Bernard+14.
    l_in_stream = np.where((stream_lat > -0.1) & (stream_lat < 0.1))[0]
    stream_longl = stream_long[l_in_stream]

    # Properties of the length histogram, in degrees
    lsigma = np.std(stream_longl) # Proxy for the length.
    lfwhm = lsigma * 2.35482
    lmean = np.average(stream_longl) # Center of the stream.

    # Do Monte-Carlo, create the arrays to hold data.
    llen_mc = np.zeros(n_mc)
    lsum_mc = np.zeros(n_mc)

    for i in range(n_mc):

        # Make the length histogram.
        lhist_bins = 80
        lhist, lbin_edges = np.histogram(stream_longl,bins=lhist_bins,range=(-10,10))
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
            llim = 2*lbin_centers[start_lpeak + lgtn_ind] + 50 + lnoise
            if lhist_n[start_lpeak + lgtn_ind] > llim:
                where_lgtn = np.append(where_lgtn,start_lpeak+lgtn_ind)
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
            llim = 2*lbin_centers[start_lpeak - lgtn_ind] + 50 + lnoise
            if lhist_n[start_lpeak - lgtn_ind] > llim:
                where_lgtn = np.append(where_lgtn,start_lpeak-lgtn_ind)
            else: stop_ll = 1
        ##wh
        where_lgtn = np.sort(where_lgtn)

        # Calculate the length and the sum:
        llen_mc[i] = max(lbin_centers[where_lgtn]) - min(lbin_centers[where_lgtn])
        lsum_excess = len(where_lgtn)*(2*np.average(lbin_centers[where_lgtn]) + 50)
        lsum_mc[i] = np.sum(lhist_n[where_lgtn]) - lsum_excess
    ###i

    # Average the Monte Carlo results and get some errors.
    llen = np.median(llen_mc)
    llen_err = np.std(llen_mc)
    lsum = np.median(lsum_mc)
    lsum_err = np.std(lsum_mc)

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
    stream_latw_am = stream_lat[w_in_stream] * 60
    whist_bins = 133

    # Do Monte-Carlo, create the arrays to hold data.
    wfwhm_mc = np.zeros(n_mc)
    wfwhm_err_mc = np.zeros(n_mc)
    wsum_mc = np.zeros(n_mc)

    for i in range(n_mc):
        whist, wbin_edges = np.histogram(stream_latw_am,bins=whist_bins,range=(-60,60))
        wbin_centers = wbin_edges + (wbin_edges[1]-wbin_edges[0])/2

        # Add noise, values follow Bernard+14
        whist = whist.astype(float)
        whist *= NtoN # Convert to MSTO stars.
        whist_n = whist + np.random.normal(0,wnoise,whist_bins)
        whist_n += wback

        # Fit Gaussian to the width component.
        wpopt, wcov = curve_fit(gaussian_func, wbin_centers[:-1], whist_n,
                                p0=(100,0,20,wback))
        wpopt[2] = np.absolute(wpopt[2]) # Gaussian allows +/-, force +
        werr = np.sqrt(np.diag(wcov))
        wgaussx = np.linspace(-60,60,num=1000)
        wgaussy = gaussian_func(wgaussx, wpopt[0], wpopt[1], wpopt[2], wpopt[3])
        wfwhm = wpopt[2]*2.35482
        wfwhm_err = werr[2]*2.35482
        wfwhm_mc[i] = wfwhm
        wfwhm_err_mc[i] = wfwhm_err

        # Determine the Full Width at three times the Noise, then integrate the
        # stream members between those limits.
        wfwn_r = wpopt[1] + wpopt[2] * np.sqrt(2*np.log(wpopt[0]/(np.sqrt(35))))
        wfwn_l = wpopt[1] - wpopt[2] * np.sqrt(2*np.log(wpopt[0]/(np.sqrt(35))))
        where_wfwn = np.where(  ( wbin_centers > wfwn_l ) &
                                ( wbin_centers < wfwn_r ) )[0]
        wsum_mc[i] = np.sum(whist_n[where_wfwn]) - ( len(where_wfwn) * wback )
    ###i

    # Compute errors
    wfwhm_err = np.sqrt(np.square(np.average(wfwhm_err_mc))+np.square(np.std(wfwhm_mc)))
    wsum_err = np.std(wsum_mc)

    # Compute properties of the width histogram, In arcminutes!
    wfwhm = np.median(wfwhm_mc) # Median to avoid bad-fit biases.
    wsum = np.median(wsum_mc)
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
    plt.clf()

    ################################################################################

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
    stat_file.write(str(stream_vdisp))

    outplots.close()

    ################################################################################

    # # Do evaluation of central surface density.
    #
    # # Provide a range of central surface densities to test: in N MSTO / square
    # # arcseconds. Also convert gaussian parameters to arcseconds. do EVERYTHING in
    # # arcseconds.
    # n_test = np.array([1,10,25,50,75,100,200,500]) * 10**-4
    # lmean_as = lmean * 3600
    # lsigma_as = lsigma * 3600
    # wmean_as = wmean * 60
    # wsigma_as = wsigma * 60
    #
    # # Define a function that solves for x, given a gaussian and a target value n
    # # in order to determine the integration boundaries given a noise value.
    # # Standard gaussian parameters.
    # def solve_gaussian_x(n, a, x0, sigma):
    #     return np.sqrt(2) * sigma * np.sqrt( np.log( a/n ) )
    # #def
    #
    # # Define a function that will calculate the total number of particles above
    # # the noise and use Newtons method to calculate the ideal central density.
    # # n0 -- central density, nv -- noise in N MSTO / square arcsecond
    # # l0 -- length mean, w0 -- width mean, lsig -- length sigma,
    # # ngtn -- Number of MSTO stars above the noise in the image.
    # def calc_n_total(n0, nv, l0, w0, lsig, wsig, ngtn):
    #     #Calculate bounds.
    #     lbound = solve_gaussian_x( nv, n0, l0, lsig )
    #     wbound = solve_gaussian_x( nv, n0, w0, wsig )
    #
    #     #Calculate the total, only integrate the exponential part, no amplitude
    #     # (multiplied in at the end) or vertical offset.
    #     int_n_total =  quad_integrate(gaussian_func,l0-lbound,l0+lbound,
    #                                     args=(1,l0,lsig,0))[0]
    #     int_n_total *= quad_integrate(gaussian_func,w0-wbound,w0+wbound,
    #                                     args=(1,w0,wsig,0))[0]
    #     int_n_total *= n0
    #
    #     return int_n_total - ngtn
    # #def
    #
    # newton_warning = False
    #
    # try: newt_sdmax = newton( calc_n_total, guess_sdmax, tol=0.0001, maxiter=1000,
    #                         args=(noiseval, lmean_as, wmean_as, lsigma_as, wsigma_as, N_above_noise) )
    # except RuntimeError:
    #     newt_sdmax = 0
    #     newton_warning = True
    #     newt_n_total = N_above_noise
    #
    # if newton_warning == False:
    #     newt_n_total = calc_n_total(newt_sdmax, noiseval, lmean_as, wmean_as,
    #                                 lsigma_as, wsigma_as, 0)
    # ##fi

    ################################################################################
