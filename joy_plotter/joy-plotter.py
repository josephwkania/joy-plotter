#!/usr/bin/env python3

import argparse
import os

import h5py
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate, stats
from scipy.signal import detrend
from scipy.signal import savgol_filter as sg

import line_plotting

# from mpl_toolkits import mplot3d

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
matplotlib.use("Agg")


def smad(freq_time, sigma=3, clip=True):
    """
    Spectral Median Absoulte Deviation filter to clip rfi

    Args:
        freq_time: freq_time/dynamic spectra to be filtered

        sigma: sigma to clip at

        clip: clip the values to the given sigma

    Returns:

        Dynamic Spectra with the values clipped
    """
    # mads = stats.median_absolute_deviation(freq_time, axis=0)
    # threshold=1.4826*sigma
    # for j,k in enumerate(mads):
    #    cut = threshold*k
    #    if clip:
    #        freq_time[freq_time[:,j]>=cut,j]=cut
    #        freq_time[freq_time[:, j]<=-cut, j]=cut
    # return freq_time
    medians = np.median(freq_time, axis=0)
    sigs = sigma * stats.median_abs_deviation(freq_time, axis=0, scale='normal')
    if clip:
        return np.clip(freq_time, a_min=medians - sigs, a_max=medians + sigs)
    else:
        for j, sig in enumerate(sigs):
            freq_time[np.absolute(freq_time[:, j] - medians[j]) >= sig, j] = 0.0
        return freq_time


def spec_sad(gulp, window=7):
    """
    Uses Savgol filter to smooth along the time axis

    Args:
       gulp: dynamic spectra to be smoothed

       window: number of point to smooth

    Returns:

       Dynamic Spectrum smoothed along the time axis
    """
    data_type = gulp.dtype
    return sg(gulp, window, 2, axis=1).astype(data_type)


def main():
    """
    Produces Joy Division Unknown Pleasures like plots of dynamic spectra

    Args:
        None

    Returns:

        Nothing

    """
    parser = argparse.ArgumentParser(
        description="Plot pulse dynamic spectra as a line.",
        prog="line-plotter.py",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    required = parser.add_argument_group("required arguments:")
    required.add_argument(
        "-f", "--file", type=str, help=".h5 file to analyze", required=True
    )

    semi_opt = parser.add_argument_group("arguments set to defaults:")
    semi_opt.add_argument(
        "-o", "--output_path", type=str, help="path where output is saved", default="."
    )
    semi_opt.add_argument(
        "-s", "--smooth", type=int, help="How much to smooth the spectra", default=7
    )
    semi_opt.add_argument("-z", "--z_frac", type=float, help="z fraction", default=10)
    semi_opt.add_argument(
        "-t",
        "--taper",
        type=float,
        help="illumination taper, 0=no taper, 1 full taper",
        default=0.80,
    )
    semi_opt.add_argument(
        "--flip", action="store_false", help="flip the frequencies order", default=True
    )
    semi_opt.add_argument(
        "--digital", action="store_true", help="Use discrete levels", default=False
    )

    optional = parser.add_argument_group("other optional arguments:")
    optional.add_argument(
        "-n",
        "--name",
        type=str,
        nargs="?",
        help="name of outputfile, default: input name",
    )
    options = vars(parser.parse_args())

    with h5py.File(options["file"], "r") as f:
        dm_time = np.array(f["data_dm_time"])
        freq_time = detrend(np.array(f["data_freq_time"])[:, ::-1].T)
        dm_time[dm_time != dm_time] = 0
        freq_time[freq_time != freq_time] = 0
        freq_time -= np.median(freq_time)
        freq_time /= np.std(freq_time)
        """
        fch1, foff, nchan, dm, cand_id, tsamp, dm_opt, snr, snr_opt, width = (
            f.attrs["fch1"],
            f.attrs["foff"],
            f.attrs["nchans"],
            f.attrs["dm"],
            f.attrs["cand_id"],
            f.attrs["tsamp"],
            f.attrs["dm_opt"],
            f.attrs["snr"],
            f.attrs["snr_opt"],
            f.attrs["width"],
        )
        """

    # if width > 1:
    #     ts = np.linspace(-128, 128, 256) * tsamp * width * 1000 / 2
    # else:
    #     ts = np.linspace(-128, 128, 256) * tsamp * 1000

    scut = smad(freq_time, sigma=3)  # clip using spectral mad filter
    scut = spec_sad(
        scut, window=options["smooth"]
    )  # use a Savitzky-Golay filter to smooth along the frequency axis

    map = np.array(scut, dtype=float)

    # plt.imshow(map) #plot input data
    # plt.show(block=True)
    # plt.savefig('input.png')

    if options["flip"]:
        map = map[::-1, :]  # reverse y axis

    # create grid
    n_lines = 256
    n_points = 1700  # number of data points
    n_int_points = 1700
    x = np.linspace(0, map.shape[1] - 1, n_points, dtype=int)
    y_values = np.linspace(0, map.shape[0] - 1, n_lines, dtype=int)
    highest_value = np.max(map)

    # create figure
    plt.style.use("dark_background")
    # cmap = plt.get_cmap('binary')
    cmap = plt.get_cmap("binary_r")  # flip the colors for dark background
    # fig, ax = plt.subplots()
    fig, ax = plt.subplots(figsize=(5, 5))

    fig.set_facecolor((1.0, 1.0, 1.0))
    ax.set_aspect("equal")

    # ax.set_xlim(0, map.shape[0]*1.3)
    # ax.set_ylim(-map.shape[0]*0.4, map.shape[0]*1.4)
    # ax.set_xlim(0, map.shape[1])
    # ax.set_ylim(-40, map.shape[0]*1.3)
    # add a bit to allow for mountains to flow over the figure
    ax.set_axis_off()

    # draw each individual line
    # y_values_test = y_values[len(y_values) // 3: len(y_values) // 3 + 10]

    for y_value in y_values[::-1]:
        print(f"Drawing line at y_value={y_value}  ", end="\r")

        # get z data from map
        line = np.array(np.ones(n_points) * y_value, dtype=int)
        z = map[line, x]

        # set all values smaller than -30000 to 0
        # remove_values = z < -30000
        # z[remove_values] = 0

        # interpolate to n_int_points
        if options["digital"]:
            x_new = np.linspace(min(x), max(x), num=len(x))
            # tck = interpolate.CubicSpline(x_new, z)
            tck = interpolate.splrep(x_new, z, k=1)  # s=30000000, k=1)
            x_int = np.linspace(min(x), max(x), n_int_points)
            z_int = interpolate.splev(x_int, tck)
            y_int = np.array(np.ones(n_int_points) * y_value, dtype=int)
            # set all nans again (resized)
            # z_int[[remove_values[int(n_points * i / n_int_points)] for i in range(n_int_points)]] = np.nan

        else:
            x_int, y_int, z_int = x, line, z
            # set all nans again
            # z_int[remove_values] = np.nan

        # set dynamics colors and widths
        # colors = [cmap(0.15+0.85*value/highest_value) if not np.isnan(value) else cmap(0) for value in z_int]
        # original
        colors = [
            cmap(1.0 - options["taper"] + options["taper"] * value / highest_value)
            if not np.isnan(value)
            else cmap(0)
            for value in z_int
        ]
        widths = [0.2 + 0.2 * value / highest_value for value in z_int]

        line_plotting.plot_line_2d(
            ax,
            x_int,
            y_int,
            z_int,
            z_fraction=options["z_frac"],
            linewidths=widths,
            linecolors=colors,
        )

    # print('saving svg')
    # plt.savefig("graph.svg")

    if options["name"]:
        img_name = options["name"]
    else:
        img_name = os.path.splitext(os.path.basename(options["file"]))[0] + ".png"

    print(f"\nSaving {img_name}")
    # plt.tight_layout()
    plt.savefig(img_name, dpi=400)
    print("Done!")


if __name__ == "__main__":
    main()
