# Parameters
ffmpeg_path = "C:/Users/mattf/PycharmProjects/orbitmoviemaker/ffmpeg.exe"  # Change to your ffmpeg path
frames_per_second = 30
observer_relative_distance = 1.  # relative to rmax (half the size of the image)
default_alphaorbit = 0.6  # alpha of the main orbit (not particle)
vanish_length_base = 100  # length in frames of vanishing orbit (if activated)
trail_length = 50  # length of the particle trail, in number of frames
orbit_part_to_be_plotted = 'past'  # Can be "whole" or "past"
addcenter = False  # Mark center with an x
circlerad = None  # if a number, plots a circle with that radius (e.g. for BH sphere of influence)
scale_label = False  # Add spatial scale in the bottom of the image
time_label = True  # if True, adds time in the bottom of the screen
