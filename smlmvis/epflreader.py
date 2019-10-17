import numpy as np
import logging
logger = logging.getLogger('global')
import pandas as pd




class EPFLReader(object):
    def __init__(self, filename):
        '''
        Parse EPFL Challenge dataset
        :param filename: Path to file
        '''
        self._filename = filename
        self._values = None
        self._points = None
        logger.debug("Starting decode for {}".format(filename))
        self._read_ascii()
        logger.debug("Complete")


    def _read_ascii(self):
        data = pd.read_csv(self._filename)
        pts = data[['xnano', 'ynano', 'znano']]
        assert(pts.shape[1] == 3)
        vls = data[['frame', 'Ground-truth', 'intensity ']]
        self._points=pts.values.copy()
        self._values=vls.values.copy()
        self._columns = ['frame', 'Ground-truth', 'intensity ']


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
