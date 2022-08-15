import numpy as np
from scipy.interpolate import interp1d


def read_orbit_file(filename):
    data = np.loadtxt(filename)
    x = data[:, 0]
    y = data[:, 1]
    z = data[:, 2]
    t = data[:, 3]
    orbit = Orbit(x, y, z, t, filename)
    return orbit


class Orbit:
    def __init__(self, x, y, z, t, filename):
        self.x = x
        self.y = y
        self.z = z
        self.t = t
        self.filename = filename

    def interpolate_orbit(self, new_frames_per_frame=2, savefile=None):
        x, y, z, t = self.x, self.y, self.z, self.t
        tnew = np.linspace(np.min(t), np.max(t), len(t) * (new_frames_per_frame))
        fx = interp1d(t, x, kind='cubic')
        fy = interp1d(t, y, kind='cubic')
        fz = interp1d(t, z, kind='cubic')
        xnew = fx(tnew)
        ynew = fy(tnew)
        znew = fz(tnew)
        self.x, self.y, self.z, self.t = xnew, ynew, znew, tnew
        if savefile:
            new_data = np.transpose(np.array([self.x, self.y, self.z, self.t]))
            np.savetxt("../sample/"+savefile, new_data, fmt="%f", header="x(kpc) y(kpc) z(kpc) t(Myr)")


def file_format_transform(old_orbit_filename="../sample/old_sample_orbit_data.ascii"):
    data = np.loadtxt(old_orbit_filename, skiprows=2)
    x = 3.5 * data[:, 0]  # in kpc
    y = 3.5 * data[:, 1]
    z = 3.5 * data[:, 2]
    with open(old_orbit_filename, "r") as file:
        first_line = file.readline()
    deltat = float(first_line[10:]) * 13.0552 / 0.7  # in Myr
    t = np.arange(len(x)) * deltat
    print(t)
    new_data = np.transpose(np.array([x, y, z, t]))
    np.savetxt("../sample/sample_orbit_data.ascii", new_data, fmt="%f", header="x(kpc) y(kpc) z(kpc) t(Myr)")