import numpy as np
import struct
import logging
logger = logging.getLogger('global')


def type_to_unpack(name):
    d = {'UINT32':'I', 'float':'f'}
    return d[name]

def bytelength(name):
    d = {'I':4, 'f': 4}
    return d[name]

def desctotypes(filename):
    with open(filename, 'r') as f:
        return [(word[:word.index('(')].strip(), type_to_unpack(word[word.index('(')+1:word.index(')')])) for word in f.readline().split(',')]


def mask_z(points, values, threshold):
    '''
    :param points: np array 3xN
    :param values: np array 6xN
    :param threshold: Z values, all rows where abs(points[:,2]) exceeds treshold are omitted
    :return: truncated points, values
    '''
    n = len(points)
    mask = abs(points[:,2]) < threshold
    trimmedp, trimmedv = points[mask].copy(), values[mask].copy()
    diff = n - len(trimmedp)
    pct = (diff / n) * 100
    # logger.info('Removed {} points, {} % of data'.format(diff , pct))
    return trimmedp, trimmedv, pct


def gsd_preprocess(points, values):
    p, v, pct = mask_z(points, values, 1e15) # GSD has corrupt Z values (~-500 to 500 with outliers 1e36
    return p, v, pct

def filter_photoncount_outliers(points, values, k_sd=3):
    '''
    Photoncount filter for GSD
    :param points: N x 3 numpy array
    :param values: N x 6 numpy array where col [3] is photoncount
    :param k_sd: all points < mean + k_sd * sd are filtered out
    :return: p', v'
    '''
    pc = values[:,3]
    pm, ps = np.mean(pc), np.std(pc)
    threshold = pm + k_sd*ps
    logger.info('Treshold {}'.format(threshold))
    m = values[:,3] <= threshold
    pk, vk = points[m].copy(), values[m].copy()
    logger.info('Removed {} % points'.format(100 - len(pk)/len(points) * 100))
    return pk, vk

def filter_z_plane(points, values, z, Z):
    '''
    :param points: N x 3 np array
    :param values: N x 6 np array
    :param z: min z value tolerated
    :param Z: max z value tolerated
    :return: A filtered copy of points, values
    '''
    assert(points.shape[1] == 3)
    logger.info('Min {} Max {}'.format(z, Zb))
    m = (points[:,2] > z) & (points[:,2] <= Z)
    return points[m].copy(), values[m].copy()


def slice_roi(points, values, roimin, roimax):
    mask = None
    # True if value in accepted range
    for i, (m, M) in enumerate(zip(roimin, roimax)):
        mm = np.logical_and(points[:,i] > m, points[:,i] < M)
        if mask is not None:
            mask = np.logical_and(mask, mm)
        else:
            mask = mm
    fp = points[mask].copy()
    fv = values[mask].copy()
    N = len(points)
    d = len(fp)
    logger.info('Removed {} out of {} points'.format(N-d, N))
    return fp, fv


def split_gsd_data(fpoints, fvalues):
    x, y, z = np.split(fpoints, 3, 1)
    stack_id, frame_id, eventid, pcount, sigma_x, sigma_y = np.split(fvalues, 6, 1)
    return (x, y, z), (stack_id, frame_id, eventid, pcount, sigma_x, sigma_y)


def dilutions(tree, imnr):
    return list(tree[imnr].keys())


def antibodies(tree, imnr, dilution):
    return list(tree[imnr][dilution].keys())


def files(tree, imnr, dilution, ab):
    fs = tree[imnr][dilution][ab]
    for filename, fullname in fs.items():
        if filename.endswith('.bin'):
            return fullname


def loadcell(n, tree):
    '''
    :param n: String integer cell number
    :param tree:  Nested dict of cell number, dilution, antibody , filename
    :return: dict of data[dilution][antiobody] = (x,y,z), (stack_id, frame_id, eventid, pcount, sigma_x, sigma_y)
    '''
    data = {}
    dis = dilutions(tree, n)
    for di in dis:
        data[di] = {}
        ab_s = antibodies(tree, n, di)
        for ab in ab_s:
            bin_filename = files(tree, n, di, ab)
            r = GSDReader(bin_filename)
            raw_points, raw_values = r.points, r.values
            fpoints, fvalues = mask_z(raw_points, raw_values, 1e30)
            (x,y,z), (stack_id, frame_id, eventid, pcount, sigma_x, sigma_y) = split_gsd_data(fpoints, fvalues)
            data[di][ab]=(x,y,z), (stack_id, frame_id, eventid, pcount, sigma_x, sigma_y)
    return data


def loadcellfile(fname):
    r = GSDReader(fname)
    raw_points, raw_values = r.points, r.values
    fpoints, fvalues = mask_z(raw_points, raw_values, 1e30)
    (x, y, z), (stack_id, frame_id, eventid, pcount, sigma_x, sigma_y) = split_gsd_data(fpoints, fvalues)
    return (x, y, z), (stack_id, frame_id, eventid, pcount, sigma_x, sigma_y)


class GSDReader(object):
    '''
        Read a GSD file into memory
        Filename is the name of the binary data file, with a corresponding filename.desc file in the same location.
        This header file is parsed to get the alignment (int, float, etc).
    '''
    def __init__(self, filename, preprocess=True, binary=True):
        '''
        Parse GSD files.
        :param filename: Path to file (if binary, expects filename.desc with headers for encoding
        :param preprocess: If true, remove invalid values
        :param binary: If true, reads binary files. If False, ascii.
        '''
        self._filename = filename
        self._values = None
        self._points = None
        self._preprocess = preprocess
        logger.debug("Starting decode for {}".format(filename))
        self._binary = binary
        if self._binary:
            self._GSD_Encoding=desctotypes(filename+'.desc')
            logger.debug('Decoded types: {}'.format(self._GSD_Encoding))
            self._unpack_string = ''.join(t for _,t in self._GSD_Encoding)
            # logger.info('Unpack String: {}'.format(self._unpack_string))
            self._columns = [n for n, _ in self._GSD_Encoding]
            self._read_binary()
        else:
            self._read_ascii()
        self._post_read()

    def _post_read(self):
        if self._points is None:
            logger.error("Decoding failed, can't preprocess")
            return
        N = self._points.shape[0]
        if self._preprocess:
            # logger.debug("Pre Processing data for invalid Z coordinates")
            fp, fv, pct = gsd_preprocess(self._points, self._values)
            # logger.info('PCt {}'.format(pct))
            self._pct = pct
            self._points = fp
            self._values = fv
            logger.debug("Retained {} out {} points".format(self._points.shape[0], N))
        else:
            pass
        logger.debug('Mean X {} Mean Y {} Mean Z{}'.format(np.mean(self._points[:,0]), np.mean(self._points[:,1]), np.mean(self._points[:,2])))


    def _read_ascii(self):
        data = None
        with open(self._filename) as f:
            line = f.readline()
            if '#' in line:
                line = line[1:]
            self._GSD_Encoding = [(word[:word.index('(')].strip(), type_to_unpack(word[word.index('(')+1:word.index(')')])) for word in line.split(',')]
            self._columns = [n for n, _ in self._GSD_Encoding]
        try:
            decoded = np.loadtxt(self._filename, delimiter=',', dtype='float64')
            assert(decoded.shape[1] == 9)
            self._values = np.empty((decoded.shape[0], decoded.shape[1]-3))
            self._points = decoded[:,3:6]
            self._values[:,0:3] = decoded[:,0:3]
            self._values[:,3:] = decoded[:,6:]

        except Exception as e:
            logger.error('Failed parsing file with exception {}'.format(e))


    def _read_binary(self):
        data = None
        linelength = sum(bytelength(t) for t in self._unpack_string)

        with open(self._filename, "rb") as f:
            data = f.read()

        bytecount = len(data)
        assert(bytecount % linelength == 0)
        linecount = int(bytecount / linelength)
        # logger.debug('Have {} bytes, {} bytes per line, {} lines.'.format(bytecount, linelength, linecount))

        points = np.empty((linecount, 3))
        values = np.empty((linecount, 6))
        if data:
            for linenumber in range(linecount):
                line = data[ linenumber * linelength : (linenumber+1)*linelength]
                assert(len(line)== linelength)
                decoded = struct.unpack(self._unpack_string, line)
                x, y, z = decoded[3:6]
                points[linenumber] = x, y, z
                values[linenumber] = decoded[0:3] + decoded[6:]
            logger.info('Decoding complete of file {}'.format(self._filename))
            self._points = points
            self._values = values
        else:
            logger.error('Failed reading data.')
        # logger.debug('Mean X {} Mean Y {} Mean Z{}'.format(np.mean(self._points[:,0]), np.mean(self._points[:,1]), np.mean(self._points[:,2])))



    @property
    def points(self):
        return self._points

    @property
    def invalid_data_pct(self):
        return self._pct

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
#
# if __name__=='__main__':
#     # r = DlpReader('../data/20170725_1_G.3dlp')
#     # p = r.points[0]
#     # print(type(p))
#     # print(p.shape)
#     # word_types = desctotypes("../../data/1000_Image1_cy3b_eventlist_July252018.bin.desc")
#     r = GSDReader("../../data/600_Image2_CY3B_eventlist_July242018.bin")
#     logger.info('Have {} points'.format(len(r.points)))