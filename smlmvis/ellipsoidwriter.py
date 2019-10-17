from smlmvis.vtuwriter import VtuWriter
from collections import namedtuple
from numpy import sin as nsin
from numpy import cos as ncos
from numpy import arccos as narccos
from numpy import random as nrp
from numpy import sqrt
from numpy import pi
import numpy as np


class AbstractWriter:
    def __init__(self):
        pass

    @property
    def written(self):
        return self._written

    def write(self, name, points, values):
        VtuWriter(name, points, values)
        self._written = (points, values)


class Ellipsoid(AbstractWriter):
    def __init__(self, a=1, b=1, c=1, center=(0,0,0), perturbfactor = 1):
        self.a = a
        self.b = b
        self.c = c
        self.center = center
        self.perturbfactor = perturbfactor


    def _point(self, u, v):
        c = len(u)
        points = np.empty((c, 3))
        nv = nsin(v)
        points[:,0] = self.a * ncos(u) *nv + self.center[0]
        points[:,1] = self.b * nsin(u) *nv + self.center[1]
        points[:,2] = self.c * ncos(v) + self.center[2]#
        if self.perturbfactor != 1:
            points = points * (1 + np.random.random((c,3))) * self.perturbfactor
        return points

    def sample(self, rng=nrp, count=1, seed=0):
        rng.seed(seed)
        return self._point(rng.uniform(0, 1, count)*pi, rng.uniform(0, 1, count)*(2*pi))

    def write(self, count, name, values = None):
        points = self.sample(rng=nrp, count=count)
        if values is None:
            values = points
        super().write(name, points, values)

class HalfEllipsoid(AbstractWriter):
    def __init__(self, a, b, c, center, perturbfactor, zlimits, minima, maxima):
        '''
        Create an ellipsoid(a,b,c, center, perturbfactor) sliced by zlimits embedded in a plane specified by m, M.
        '''
        self.plane = Planar(minima, maxima)
        self.half_ellipse = OpenEllipsoid(a, b, c, center, perturbfactor, zlimits)

    def _in_ellipsoid(self, point):
        '''
        Test if point is contained within the member ellipsoid.
        #Fix non center intercept (i.e. also use Z intercept)
        '''
        x, y, _ = [(p - c)**2 for p,c in zip(point, self.half_ellipse.center)]
        rx, ry = self.half_ellipse.a, self.half_ellipse.b
        return (x / (rx**2)) + (y / (ry**2)) <= 1

    def write(self, count, name, values=None):
        ef_points = self.plane.sample(count=count)
        # Remove all points where x,y are within the area of the ellipsoid specified.
        # All x where d(x, center) < center[1]+-a)
        # All y where d(y, center) < b)
        mask = [~self._in_ellipsoid(ef_points[i,:]) for i in range(count)]
        e_points = ef_points[mask,:]
        masked = len(mask)
        p_points = self.half_ellipse.sample(rng=nrp, count=masked) # Not quite fair since sphere will be more depth

        points = np.concatenate((e_points, p_points))
        if values is None:
            values = points
        values = points
        super().write(name, points, values)


class OpenEllipsoid(Ellipsoid):
    def __init__(self, a=1, b=1, c=1, center=(0,0,0), perturbfactor =1, zrange=(-float('inf'), float('inf'))):
        assert(c>=1)
        super().__init__(a, b, c, center, perturbfactor)
        self.zrange = zrange

    @staticmethod
    def ztov(z, c, ctr):
        '''
        Given z coordinate, scale of Z radius and center point, return the angular threshold value correspoding with z (0-pi).
        # TODO fix c < 1 flipping range
        '''
        arg = (z-ctr)/c
        if arg >= 1:
            return 0
        if arg <= -1:
            return 0
        return narccos(arg)


    def sample(self, rng=nrp, count=1, seed=0):
        '''
        Sample count points (shape=count x 3) in 3D space using rng and seed.
        '''
        assert(count>0)
        rng.seed(seed)
        vm, vM = [OpenEllipsoid.ztov(z, self.c, self.center[2]) for z in self.zrange]
        h = int(count//2) if count > 1 else 1 # count pos int, int() floors
        fcount, scount = (h, h+1) if (count % 2) else (h, h)
        lowv = vM + rng.uniform(0, 1, fcount) * (vm - vM)      # explicit since rng api may not have range (random())
        highv = 2*pi - vm + rng.uniform(0, 1, scount) * (vm - vM)
        vpoints = np.concatenate((lowv, highv))
        upoints = rng.uniform(0,1,count)*pi
        return self._point(upoints, vpoints)


class GaussianEllipsoid(AbstractWriter):
    '''
    Create 3D Gaussian.
    '''
    def __init__(self, means, sigmas, seed = 0):
        self._rng = nrp
        nrp.seed(seed)
        self._means = means
        self._sigmas = sigmas

    def sample(self, count=1):
        points = np.empty((count, 3))
        for index, (mean, sigma) in enumerate(zip(self._means, self._sigmas)):
            points[:,index] = self._rng.normal(mean, sigma, count)
        return points

    def write(self, count, name):
        points = self.sample(count)
        values = np.arange(count).reshape(count, 1)
        super().write(name, points, values)

class Planar(AbstractWriter):
    def __init__(self, minima, maxima, seed=0):
        self._rng = nrp
        self._rng.seed(seed)
        self.minima = minima
        self.maxima = maxima

    def sample(self, count):
        points = np.empty((count, 3))
        for index, (min_ax, max_ax) in enumerate(zip(self.minima, self.maxima)):
            points[:,index] = self._rng.uniform(min_ax, max_ax, count)
        return points

    def write(self, count, name):
        points = self.sample(count)
        values = points
        super().write(name, points, values)


class Ellipsoids:
    def __init__(self, params):
        self.es = [Ellipsoid(*param) for param in params]

    def write(self, counts, names):
        k = len(self.es)
        for index, (ellipsoid, count, name) in enumerate(zip(self.es, counts, names)):
            values = np.full((count, 1), index)
            ellipsoid.write(count, name, values)

    def points(self):
        return [ell.written for ell in self.es]
