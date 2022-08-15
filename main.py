from src.read_orbit import read_orbit_file
from src.plot_orbit import *

# Note: Adapt path to ffmpeg in parameters.py to your machine

if __name__ == '__main__':
    # Load orbit
    orbit_filename = "sample/sample_orbit_data.ascii"
    orbit = read_orbit_file(orbit_filename)
    # Interpolate orbit to double the number of frames
    orbit.interpolate_orbit(new_frames_per_frame=2)
    # Static orbit image
    static_orbit_plot(orbit, frame=300, savefile="testimage")
    # Movie with rotation and zooming out
    make_orbit_movie(orbit, startframe=4400, endframe=5400, skipevery=2, color="orange", output_label="testmovie",
                     inclination_start=0, inclination_end=180, rmax_start=0.3, rmax_end=0.9)
