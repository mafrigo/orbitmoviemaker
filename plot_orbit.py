import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib.animation as mani
from parameters import *


def static_orbit_plot(orbit, frame=-1, rmax=None, color='r', orbit_label=None, inclination=0, alphaorbit=None, savefile=None):
    # Note: only uses x and y
    plt.style.use('dark_background')
    x, y, z, t = orbit.x, orbit.y, orbit.z, orbit.t
    if alphaorbit is None:
        alphaorbit = default_alphaorbit
    if inclination != 0:
        ang = np.pi * float(inclination) / 180.
        xnew = x * np.cos(ang) - z * np.sin(ang)
        znew = z * np.cos(ang) + x * np.sin(ang)
        x = xnew
        z = znew

    nframes = len(x)
    if rmax is None:
        rmax = np.sqrt(np.max(np.array(x) ** 2 + np.array(y) ** 2 + np.array(z) ** 2))
    if frame == -1:
        frame = nframes - 1

    # Pseudo-3D correction
    robs = obs_fac * rmax  # radius from which we are observing
    x = x * (robs / (z + robs))
    y = y * (robs / (z + robs))
    zcond = (z > -0.99 * obs_fac * rmax)  # whether orbit should be visible from POV
    frame_zcond = (z[frame] > -0.99 * obs_fac * rmax)

    # Plotting overall orbit
    if orbit_part_to_be_plotted == 'whole':
        plt.plot(x[zcond], y[zcond], color='grey', zorder=1, linewidth=1, alpha=alphaorbit)
    elif orbit_part_to_be_plotted == 'past':
        x = x[:frame + 1]
        y = y[:frame + 1]
        z = z[:frame + 1]
        t = t[:frame + 1]
        zcond = zcond[:frame + 1]
        nframes = len(x)
        plt.plot(x[zcond], y[zcond], color='grey', zorder=1, linewidth=1, alpha=alphaorbit)
    else:
        raise IOError("orbit_part_to_be_plotted can either be 'whole' or 'past'")

    # Plotting moving particle with trail
    alpha = np.zeros(nframes)
    alpha[frame - trail_length:frame] = np.arange(trail_length).astype(float) / float(trail_length)
    alpha = alpha[zcond]
    rgbcolor = colors.to_rgba(color)
    color_array = np.asarray([(rgbcolor[0], rgbcolor[1], rgbcolor[2], a ** 2) for a in alpha])
    ax0 = plt.gca()
    plt.scatter(x[zcond], y[zcond], c=color_array, zorder=2, s=1)  # plot trail
    if frame_zcond:  # if the particle is in the visible area
        ax0.scatter(x[frame], y[frame], color=color, edgecolor='k', s=20, zorder=3)  # plot particle

    # Markers and labels
    if addcenter:
        ax0.scatter(0., 0., marker='x', color='grey', s=5, zorder=1)
    if circlerad:
        circle = plt.Circle((0, 0), circlerad, color='w', zorder=1)
        ax0.add_artist(circle)
    if time_label:
        time = t[frame]
        plt.text(0.5 * rmax, -0.88 * rmax, str(int(time)) + " Myr", color='w', fontsize=12)
    if orbit_label is not None:
        plt.text(0.5 * rmax, 0.88 * rmax, orbit_label, color=color, alpha=alphaorbit / 0.6, fontsize=16)

    # Final settings
    ax0.set_xlim([-rmax, rmax])
    ax0.set_ylim([-rmax, rmax])
    ax0.axes.set_aspect('equal')
    plt.axis('off')
    plt.tight_layout()
    if savefile:
        plt.savefig("output/"+savefile)


def make_orbit_movie(orbit, output_label='testmovie',
                     startframe=0, endframe=None,
                     rmax_start=None, rmax_end=None,
                     inclination_start=0., inclination_end=None,
                     skipevery=1, color='r',  vanish=None, orbit_label=None,
                     axes_mini_plot=True, file_frame_delay=0, saveframes=False):

    # Initialize parameters
    x, y, z, t = orbit.x, orbit.y, orbit.z, orbit.t
    if saveframes and not os.path.exists(output_label):
        os.mkdir(output_label)
    vanish_length = vanish_length_base * skipevery
    if endframe is None:
        endframe = len(x)
    if rmax_start is None:
        rmax_start = np.sqrt(np.max(np.array(x) ** 2 + np.array(y) ** 2 + np.array(z) ** 2))
    if rmax_end is None:
        rmax_end = rmax_start
    if inclination_end is None:
        inclination_end = inclination_start

    # Initialize movie
    plt.rcParams['animation.ffmpeg_path'] = ffmpeg_path
    FFMpegWriter = mani.writers['ffmpeg']
    metadata = dict(title=output_label, artist='Matplotlib')
    writer = FFMpegWriter(fps=frames_per_second, metadata=metadata)
    fig = plt.figure()

    # Main loop over movie frames
    with writer.saving(fig, "output/"+output_label+".mp4", 100):
        for frame in np.arange(startframe, endframe):
            if frame % skipevery != 0:
                continue
            if vanish is not None:  # Make orbit vanish at start or end of movie
                if vanish == 'start':
                    alphaorbit = 0.6 * np.min([float(frame - startframe) / vanish_length, 1.])
                elif vanish == 'end':
                    alphaorbit = 0.6 * np.min([float(endframe - frame) / vanish_length, 1.])
                elif isinstance(vanish, float):
                    alphaorbit = vanish
                else:
                    raise IOError("Vanish type not recognized")
            else:
                alphaorbit = 0.6

            rmax = rmax_start + (rmax_end - rmax_start) * (frame - startframe) / (endframe - startframe)
            inclination = inclination_start + (inclination_end - inclination_start) * (frame - startframe) / (endframe - startframe)
            static_orbit_plot(orbit, frame=frame, inclination=inclination, rmax=rmax, color=color,
                              orbit_label=orbit_label, alphaorbit=alphaorbit)
            if axes_mini_plot:
                plot_axes(rmax, inclination)
            if saveframes:
                if file_frame_delay is None:
                    plt.savefig(output_label+'/orbit%.4i' % frame)
                else:
                    plt.savefig(output_label+'/orbit%.4i' % (frame + file_frame_delay))

            writer.grab_frame()
            plt.clf()


def plot_axes(rmax, inclination, arrowsize=0.08, arrow_offset=0.05):
    ang = np.pi * inclination / 180.
    axes_center_loc = [-0.93 * rmax, 0.93 * rmax]
    arrow_length = arrowsize * rmax
    displ = arrow_offset * arrow_length
    x_real = [np.sqrt(arrow_length ** 2 - displ ** 2), -displ, 0.]
    y_real = [displ, np.sqrt(arrow_length ** 2 - displ ** 2), 0.]
    z_real = [0., 0., arrow_length]

    def real_to_proj(vec):
        xor = vec[0]
        yor = vec[1]
        zor = vec[2]
        xnew = xor * np.cos(ang) - zor * np.sin(ang)
        znew = zor * np.cos(ang) + xor * np.sin(ang)
        robs = obs_fac * rmax  # radius from which we are observing
        x = xnew * (robs / (znew + robs))
        y = yor * (robs / (znew + robs))
        z = znew
        vecnew = [x, y, z]
        return vecnew

    x_proj = real_to_proj(x_real)
    y_proj = real_to_proj(y_real)
    z_proj = real_to_proj(z_real)
    plt.arrow(axes_center_loc[0], axes_center_loc[1], x_proj[0], x_proj[1], color='springgreen',
              width=0.000000000001)  # head_width=0.000000000001)
    plt.arrow(axes_center_loc[0], axes_center_loc[1], y_proj[0], y_proj[1], color='r',
              width=0.000000000001)  # , head_width=0.000000000001)
    plt.arrow(axes_center_loc[0], axes_center_loc[1], z_proj[0], z_proj[1], color='dodgerblue',
              width=0.000000000001)  # , head_width=0.000000000001)