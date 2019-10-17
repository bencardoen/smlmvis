import numpy as np
import pandas as pd

class RainStormReader(object):
    def __init__(self, filename):
        self._filename = filename
        self._values = None
        self._points = None
        self._columns = ['idx', 'frame_idx', 'x_coord', 'y_coord', 'I', 'sig_x',
       'sig_y', 'avg_brightness', 'res', 'res_Row', 'res_Col', 'roi_min',
       'Sum_signal', 'Sum_signal_ph', 'x_std', 'y_std',
       'ellip_xy']
        self._read()

    def _read(self):
        data = pd.read_csv(self._filename)
        pxtonm = max(data['x_nm']) / max(data['x_coord'])
        self._points = data[['x_nm', 'y_nm', 'z_coord']].values.copy()
        self._points[:,2] *= pxtonm
        self._values = data[['idx', 'frame_idx', 'x_coord', 'y_coord', 'I', 'sig_x',
       'sig_y', 'avg_brightness', 'res', 'res_Row', 'res_Col', 'roi_min',
       'Sum_signal', 'Sum_signal_ph', 'x_std', 'y_std',
       'ellip_xy']].values.copy()

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
