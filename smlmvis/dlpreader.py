import numpy as np

class DlpReader(object):
    ''' Read a .3dlp file into memory'''
    def __init__(self, filename):
        self._filename = filename
        self._values = None
        self._points = None
        self._columns = ['std_x', 'std_y', 'std_z', 'amplitude', 'framenumber']
        self._read()

    def _read(self):
        A = np.loadtxt(self._filename)
        self._points = A[:,0:3] # X Y Z
        self._values = A[:,3:] # remainder

    @property
    def points(self):
        return self._points

    @property
    def values(self):
        return self._values

    @property
    def value_names(self):
        return self._columns

    def points_generator(self):
        cnt = 0
        l = len(self.points)
        while cnt < l:
            yield self._points[cnt]

    def values_generator(self):
        cnt = 0
        l = len(self.values)
        while cnt < l:
            yield self._values[cnt]

if __name__=='__main__':
    r = DlpReader('../data/20170725_1_G.3dlp')
    p = r.points[0]
    print(type(p))
    print(p.shape)
