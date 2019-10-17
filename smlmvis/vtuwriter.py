import vtk
import numpy as np

class TemporalVtuWriter(object):
    def __init__(self, filename, points, values, incremental=False, limit=0, collate_frames=1):
        '''
        Writes filename<number>.vtu of points, values sliced along frames.
        The slicenumber corresponds with the last fields of the k-tuple in values.
        '''
        self._filename = filename
        self._points = points
        self._values = values
        self._incremental = incremental
        self._limit = limit
        self._skipframes = 1 if collate_frames <= 1 else collate_frames
        self._writeAll()

    def _writeAll(self):
        fname = self._filename
        # find out if doing binary split 40k times is more efficient
        lastframe = 1
        lastindex = 0
        for index, (point, value) in enumerate(zip(self._points, self._values)):
            fnumber = int(value[-1])
            if fnumber > lastframe + self._skipframes - 1:
                # print("New framenumber from {} to {} with index going from {} to {}".format(lastframe, fnumber, lastindex, index))
                l = 0 if self._incremental else lastindex
                p, v = self._points[l:index, :], self._values[l:index, :]
                sfname = "{}{}".format(fname, fnumber-1)
                VtuWriter(sfname, p, v)
                lastframe, lastindex = fnumber, index
                if self._limit != 0 and fnumber > self._limit:
                    return


class VtuWriter(object):
    def __init__(self, filename, points, values):
        '''
        Appends .vtu to filename
        '''
        self._points= vtk.vtkPoints()
        self._grid = vtk.vtkUnstructuredGrid()
        self._values = vtk.vtkDoubleArray()
        self._values.SetName('point_values_array')
        self._grid.SetPoints(self._points)
        self._grid.GetPointData().SetScalars(self._values)
        self._loadPoints(points, values)
        self._write("{}.vtu".format(filename))

    def _loadPoints(self, points, values):
        poly = vtk.vtkPolyVertex()
        poly.GetPointIds().SetNumberOfIds(points.shape[0])
        for i, (point, value) in enumerate(zip(points, values)):
            poly.GetPointIds().SetId(i, self._points.InsertNextPoint(point))
            self._values.InsertNextValue(values[i][-1])
        self._grid.InsertNextCell(poly.GetCellType(), poly.GetPointIds())

    def _write(self, filename):
        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetFileName(filename)
        writer.SetInputData(self._grid)
        writer.Write()


class GSDWriter(VtuWriter):
    def __init__(self, filename, points, values):
        '''
        Appends .vtu to filename
        '''
        self._points= vtk.vtkPoints()
        self._grid = vtk.vtkUnstructuredGrid()
        self._values = vtk.vtkDoubleArray()
        self._values.SetName('point_values_array')
        self._grid.SetPoints(self._points)
        self._grid.GetPointData().SetScalars(self._values)
        self._loadPoints(points, values)
        self._write("{}.vtu".format(filename))

    def _loadPoints(self, points, values):
        poly = vtk.vtkPolyVertex()
        poly.GetPointIds().SetNumberOfIds(points.shape[0])
        for i, (point, value) in enumerate(zip(points, values)):
            poly.GetPointIds().SetId(i, self._points.InsertNextPoint(point))
            self._values.InsertNextValue(values[i][-3]) # fix arbitrary data
        self._grid.InsertNextCell(poly.GetCellType(), poly.GetPointIds())


class KChannelVtuWriter(VtuWriter):
    def __init__(self, filename, points, values, channels):
        self._channels = channels
        self._points = vtk.vtkPoints()
        self._grid = vtk.vtkUnstructuredGrid()
        self._values = vtk.vtkDoubleArray()
        self._values.SetName('point_values_array')
        self._grid.SetPoints(self._points)
        self._grid.GetPointData().SetScalars(self._values)
        self._loadPoints(points, values)
        self._write("{}.vtu".format(filename))

    def _loadPoints(self, pointset, valueset): # values is set
        for index, (points, _, values) in enumerate(zip(pointset, self._channels, valueset)):
            poly = vtk.vtkPolyVertex()
            poly.GetPointIds().SetNumberOfIds(points.shape[0])
            for i, (point, value) in enumerate(zip(points, values)):
                poly.GetPointIds().SetId(i, self._points.InsertNextPoint(point))
                self._values.InsertNextValue(index)
            self._grid.InsertNextCell(poly.GetCellType(), poly.GetPointIds())


class NNWriter(VtuWriter):
    def __init__(self, filename, points, values, channels, edges):
        '''
        Writes a k-way 1-nearest neighbour graph between points using edges.
        '''
        self._channels = channels
        self._points = vtk.vtkPoints()
        self._grid = vtk.vtkUnstructuredGrid()
        self._values = vtk.vtkDoubleArray()
        self._edges = edges # Dict-> {offset_src, offset_trg} -> [(src_index, trg_index)]
        self._values.SetName('point_values_array')
        self._grid.SetPoints(self._points)
        self._grid.GetPointData().SetScalars(self._values)
        self._loadPoints(points, values)
        self._write("{}.vtu".format(filename))

    def _loadPoints(self, pointset, valueset): # values is set
        # Points
        for index, (points, _, values) in enumerate(zip(pointset, self._channels, valueset)):
            poly = vtk.vtkPolyVertex()
            poly.GetPointIds().SetNumberOfIds(points.shape[0])
            for i, (point, value) in enumerate(zip(points, values)):
                poly.GetPointIds().SetId(i, self._points.InsertNextPoint(point))
                self._values.InsertNextValue(index)
            self._grid.InsertNextCell(poly.GetCellType(), poly.GetPointIds())
        # Edges
        for offset, edges in self._edges.items():
            sourceoffset, targetoffset = offset
            for source, target in edges:
                self._grid.InsertNextCell(vtk.VTK_LINE, 2, [sourceoffset+source, targetoffset+target])


class GraphWriter(VtuWriter):
    def __init__(self, filename, midpoints, values, tetras):
        '''
        Writes a k-way 1-nearest neighbour graph between points using edges.
        '''
        # self._channels=channels
        self._points= vtk.vtkPoints()
        self._grid = vtk.vtkUnstructuredGrid()
        self._values = vtk.vtkDoubleArray()
        # self._edges = edges # Dict-> {offset_src, offset_trg} -> [(src_index, trg_index)]
        self._values.SetName('point_values_array')
        self._grid.SetPoints(self._points)
        self._grid.GetPointData().SetScalars(self._values)
        self._loadPoints(midpoints, values, tetras)
        self._write("{}.vtu".format(filename))

    def _loadPoints(self, points, values, tetras):
        # Points
        poly = vtk.vtkPolyVertex()
        poly.GetPointIds().SetNumberOfIds(points.shape[0])
        for i, (point, value) in enumerate(zip(points, values)):
            poly.GetPointIds().SetId(i, self._points.InsertNextPoint(point))
            self._values.InsertNextValue(value)
        self._grid.InsertNextCell(poly.GetCellType(), poly.GetPointIds())
        # Tetras
        # print('Adding {} tetras'.format(len(tetras)))
        for tetra in tetras:
            self._grid.InsertNextCell(vtk.VTK_TETRA, 4, tetra)


class AlphaGraphWriter(VtuWriter):
    def __init__(self, filename, midpoints, values, tetras):
        '''
        Writes a k-way 1-nearest neighbour graph between points using edges.
        '''
        # self._channels=channels
        self._points= vtk.vtkPoints()
        self._grid = vtk.vtkUnstructuredGrid()
        self._values = vtk.vtkDoubleArray()
        # self._edges = edges # Dict-> {offset_src, offset_trg} -> [(src_index, trg_index)]
        self._values.SetName('point_values_array')
        self._grid.SetPoints(self._points)
        self._grid.GetPointData().SetScalars(self._values)
        self._loadPoints(midpoints, values, tetras)
        self._write("{}.vtu".format(filename))

    def _loadPoints(self, points, values, tetras):
        # Points
        poly = vtk.vtkPolyVertex()
        poly.GetPointIds().SetNumberOfIds(points.shape[0])
        for i, (point, value) in enumerate(zip(points, values)):
            poly.GetPointIds().SetId(i, self._points.InsertNextPoint(point))
            self._values.InsertNextValue(value)
        self._grid.InsertNextCell(poly.GetCellType(), poly.GetPointIds())
        # Tetras
        # print('Adding {} tetras'.format(len(tetras)))
        for tetra in tetras:
            self._grid.InsertNextCell(vtk.VTK_TETRA, 4, tetra)


class MatWriter(VtuWriter):
    def __init__(self, filename, mat, flip = True):
        '''
        Write a square matrix to VTK.
        Flip results in [0,0] to [n,n] being mapped as top left, lower right.
        '''
        assert(mat is not None)
        self._points = vtk.vtkPoints()
        self._grid = vtk.vtkUnstructuredGrid()
        self._values = vtk.vtkDoubleArray()
        self._values.SetName('point_values_array')
        self._grid.SetPoints(self._points)
        self._grid.GetPointData().SetScalars(self._values)
        self._loadPoints(mat)
        self._write("{}.vtu".format(filename))

    def _loadPoints(self, mat):
        n = mat.shape[0]
        nsq = n**2
        points = np.empty((nsq, 3))
        values = np.empty((nsq, 1))
        for k in range(nsq):
            i, j = divmod(k, n)
            points[k] = [n-i,j,0]
            values[k] = mat[i,j] if (i,j) in mat else 0

        poly = vtk.vtkPolyVertex()
        poly.GetPointIds().SetNumberOfIds(points.shape[0])
        for i, (point, value) in enumerate(zip(points, values)):
            poly.GetPointIds().SetId(i, self._points.InsertNextPoint(point))
            self._values.InsertNextValue(value)
        self._grid.InsertNextCell(poly.GetCellType(), poly.GetPointIds())


class MidPointWriter(VtuWriter):
    def __init__(self, filename, midpoints):
        '''
        Writes a point cloud of k midpoints
        '''
        self._k = len(midpoints)
        self._points = vtk.vtkPoints()
        self._grid = vtk.vtkUnstructuredGrid()
        self._values = vtk.vtkDoubleArray()
        self._midpoints = midpoints
        self._values.SetName('point_values_array')
        self._grid.SetPoints(self._points)
        self._grid.GetPointData().SetScalars(self._values)
        self._loadPoints()
        self._write("{}.vtu".format(filename))

    def _loadPoints(self):
        for (src_offset, trg_offset),midpoints in self._midpoints.items():
            p, v = zip(*midpoints)
            points = np.empty((len(p), len(p[0])))
            for index, coord in enumerate(p):
                points[index, :] = coord
            values = np.empty((len(v), len(v[0])))
            for index, value in enumerate(v):
                values[index, :] = value
            poly = vtk.vtkPolyVertex()
            poly.GetPointIds().SetNumberOfIds(points.shape[0])
            for i, (point, value) in enumerate(zip(points, values)):
                poly.GetPointIds().SetId(i, self._points.InsertNextPoint(point))
                self._values.InsertNextValue(value[-1]) # Use distance as value
            self._grid.InsertNextCell(poly.GetCellType(), poly.GetPointIds())
