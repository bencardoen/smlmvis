import numpy as np
import pandas as pd

class AbbelightReader(object):
    def __init__(self, filename):
        self._filename = filename
        self._values = None
        self._points = None
        self._columns = ['id', 'frame', 'sigma1 [nm]',
       'sigma2 [nm]', 'intensity [photon]', 'offset [photon]',
       'bkgstd [photon]', 'uncertainty_xy [nm]', 'uncertainty_z [nm]',
       'ratio_SAF [photon]', 'z_SAF [nm]', 'z_Astigm [nm]']
        self._read()

    def _read(self):
        data = pd.read_csv(self._filename)
        self._points = data[['x [nm]', 'y [nm]', 'z [nm]']].values.copy()
        self._values = data[['id', 'frame','sigma1 [nm]',
       'sigma2 [nm]', 'intensity [photon]', 'offset [photon]',
       'bkgstd [photon]', 'uncertainty_xy [nm]', 'uncertainty_z [nm]',
       'ratio_SAF [photon]', 'z_SAF [nm]', 'z_Astigm [nm]']].values.copy()

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
