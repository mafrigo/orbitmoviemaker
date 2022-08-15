import unittest
from src.plot_orbit import *
from src.read_orbit import *

orbit_filename = "../sample/sample_orbit_data.ascii"


class TestReadOrbit(unittest.TestCase):
    def test_read_orbit(self):
        orbit = read_orbit_file(orbit_filename)
        self.assertEqual(len(orbit.x), 4096)

    def test_interpolation(self):
        orbit = read_orbit_file(orbit_filename)
        orbit.interpolate_orbit(new_frames_per_frame=2)
        self.assertEqual(len(orbit.x), 8192)
        plt.close()


class TestPlotOrbit(unittest.TestCase):
    def test_plot(self):
        self.orbit = read_orbit_file(orbit_filename)
        static_orbit_plot(self.orbit, frame=1500)


if __name__ == '__main__':
    unittest.main()
