#!/usr/bin/env python3

import numpy as np
import matplotlib
from scipy import interpolate
import line_plotting, h5py, argparse, os
from mpl_toolkits import mplot3d
from scipy.signal import detrend
from scipy import stats
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn import preprocessing
from scipy.signal import savgol_filter as sg

def smad(freq_time, sigma, clip=True):
    freq_time=freq_time.copy()
    #mads = stats.median_absolute_deviation(freq_time, axis=0)
    #threshold=1.4826*sigma
    #for j,k in enumerate(mads):
    #    cut = threshold*k
    #    if clip:
    #        freq_time[freq_time[:,j]>=cut,j]=cut
    #        freq_time[freq_time[:, j]<=-cut, j]=cut
    #return freq_time
    medians = np.median(freq_time, axis=0)
    sigs = 1.4826*sigma*stats.median_absolute_deviation(freq_time, axis=0)
    if clip:
        return np.clip(freq_time, a_min=medians-sigs, a_max=medians+sigs)
    else:
        for j, sig in enumerate(sigs):
            freq_time[np.absolute(freq_time[:, j] - medians[j]) >= sig, j] = 0.0
        return freq_time

def spec_sad(gulp, window=65):
    gulp = gulp.copy()
    """
    Calculates Savgol Absolute Deviations along the spectral axis

    Args:
       frame: number of time samples to calculate the SAD

       sigma: cutoff sigma

    Returns:
     
       Dynamic Spectrum with values clipped
    """
    data_type = gulp.dtype
    return sg(gulp, window, 2, axis=1).astype(data_type)

 
def main(**options):
   #import h5 to numpy array
    with h5py.File(options['file'], 'r') as f:
        dm_time = np.array(f['data_dm_time'])
        freq_time = detrend(np.array(f['data_freq_time'])[:, ::-1].T)
        dm_time[dm_time != dm_time] = 0
        freq_time[freq_time != freq_time] = 0
        freq_time -= np.median(freq_time)
        freq_time /= np.std(freq_time)
        fch1, foff, nchan, dm, cand_id, tsamp, dm_opt, snr, snr_opt, width = f.attrs['fch1'], \
                                                                             f.attrs['foff'], f.attrs['nchans'], \
                                                                             f.attrs['dm'], f.attrs['cand_id'], \
                                                                             f.attrs['tsamp'], f.attrs['dm_opt'], \
                                                                             f.attrs['snr'], f.attrs['snr_opt'], \
                                                                             f.attrs['width']

    if width > 1:
            ts = np.linspace(-128, 128, 256) * tsamp * width * 1000 / 2
    else:
            ts = np.linspace(-128, 128, 256) * tsamp * 1000

    sigma=3.0
    scut = smad(freq_time, sigma)
    sg_smooth = True
    if sg_smooth:
        scut = spec_sad(scut, window=7)

    #map = preprocessing.scale(np.array(freq_time, dtype=float), axis=1)*5000
    map = np.array(scut, dtype=float)
    print('map.shape=', map.shape)
    
    if True: #plot map
        plt.imshow(map)
        plt.show(block=True)
        plt.savefig('input.png')
        print(f"max:{np.max(map)}, min:{np.min(map)}: std:{np.std(map)}")
    #map = map[::-1, :] # reverse y axis

    # create grid
    n_lines = 256
    n_points = 1700 # number of data points
    n_int_points = 1700
    x = np.linspace(0, map.shape[1]-1, n_points, dtype=int)
    y_values = np.linspace(0, map.shape[0]-1, n_lines, dtype=int)
    #print(f"x: {x}")
    #print(f"y_values: {y_values}")
    highest_value = np.max(map)
    print('highest_value=', highest_value)
    # create figure
    cmap = plt.get_cmap('binary')
    fig, ax = plt.subplots()
    fig.set_facecolor((1.0, 1.0, 1.0))
    ax.set_aspect('equal')
    print(map.shape[1])
    print(map.shape[0]*1.3)
    ax.set_xlim(0, map.shape[1])
    ax.set_ylim(0, map.shape[0]*1.3) #add a bit to allow for mountains to flow over the figure
    ax.set_axis_off()

    # draw each individual line
    y_values_test = y_values[len(y_values)//3:len(y_values)//3+10]

    for y_value in y_values:
        print('drawing line at y_value=', y_value, end='\r')

        # get z data from map
        line = np.array(np.ones(n_points)*y_value, dtype=int)
        z = map[line, x]

        # set all values smaller than -30000 to 0
        #remove_values = z < -30000
        #z[remove_values] = 0

        smoothen = False
        # interpolate to n_int_points
        if smoothen:
            tck = interpolate.splrep(x, z, s=options['smooth'], k=1)#s=30000000, k=1)
            x_int = np.linspace(min(x), max(x), n_int_points)
            z_int = interpolate.splev(x_int, tck)
            y_int = np.array(np.ones(n_int_points)*y_value, dtype=int)
            # set all nans again (resized)
            #z_int[[remove_values[int(n_points * i / n_int_points)] for i in range(n_int_points)]] = np.nan
            
        else:
            x_int, y_int, z_int = x, line, z
            # set all nans again
            #z_int[remove_values] = np.nan

        # set dynamics colors and widths
        colors = [cmap(0.15+0.85*value/highest_value) if not np.isnan(value) else cmap(0) for value in z_int]
        widths = [0.1 + 0.2*value / highest_value for value in z_int]

        line_plotting.plot_line_2d(ax, x_int, y_int, z_int, z_fraction=options['z_frac'], linewidths=widths, linecolors=colors)

    #print('saving svg')
    #plt.savefig("graph.svg")
    print('\nsaving png')
    if options["name"]:
        plt.savefig(options["name"])#, dpi=500)
    else:
        plt.savefig(os.path.splitext(os.path.basename(options["file"]))[0]+'.png', dpi=400)
    #plt.show()


if __name__ == "__main__":
    parser=argparse.ArgumentParser(
        description="Plot pulse dynamic spectra as a line.",
        prog='line-plotter.py',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter 
    )

    required=parser.add_argument_group('required arguments:')
    required.add_argument('-f','--file',type=str,
                          help=".h5 file to analyze",
                          required=True)

    semi_opt=parser.add_argument_group('arguments set to defaults:')
    semi_opt.add_argument('-o','--output_path',type=str,
                          help='path where output is saved',
                          default=".")
    semi_opt.add_argument('-s','--smooth',type=int,
                          help='How much to smooth the spectra',
                          default=30000)
    semi_opt.add_argument('-z','--z_frac',type=float,
                          help='z fraction',
                          default=10)
        
    optional=parser.add_argument_group('other optional arguments:')
    optional.add_argument('-n','--name',type=int,nargs='+',
                          help="name of outputfile, default: input name")

    args = vars(parser.parse_args())
    main(**args) 
