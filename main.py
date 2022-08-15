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
    static_orbit_plot(orbit, frame=1500, savefile="testimage", color="orange")
    # Movie with rotation and zooming out
    make_orbit_movie(orbit, startframe=4400, endframe=5400, skipevery=2, color="orange", output_label="testmovie",
                     angle_start=0, angle_end=180, rad_start=0.3, rad_end=0.9)
