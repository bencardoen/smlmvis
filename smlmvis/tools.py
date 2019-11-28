import os
import numpy as np
from scipy.spatial.ckdtree import cKDTree
from PIL import Image
import logging
FORMAT = "[@ %(asctime)s %(filename)s : %(lineno)s - %(funcName)20s() ] %(message)s"
logging.basicConfig(format=FORMAT, datefmt='%H:%M:%S')
lgr = logging.getLogger('global')
lgr.setLevel(logging.INFO)


def leaves(node, res):
    """
    Traverse a tree rooted at in postorder
    :param node: current node
    :param res: accumulator where 3d points stored in leaf node are appended to
    :return: None
    """
    leaf = True
    if node.lesser:
        leaf = False
        leaves(node.lesser, res)
    if node.greater:
        leaf = False
        leaves(node.greater, res)
    if leaf:
        res.append(node.indices)


def treedir(pathname):
    """
    Create a recursive dictionary mapping of a directory structure starting at pathname
    A leaf is filename->full path of file
    A internal node is dirname : {dict of contents}
    :param pathname:
    :return: {dir : {contents} or filename : filepath}
    """
    (t, d, f) = next(os.walk(pathname))
    return {**{_dir : treedir(os.path.join(pathname, _dir)) for _dir in d}, **{file:os.path.join(t, file) for file in f}}



def topix(x, y, pxsz):
    """
    Map x,y coordinates to pixel size
    :param x:
    :param y:
    :param pxsz:
    :return:
    """
    assert (x >= 0)
    assert (y >= 0)
    return np.array([int(round(x / pxsz)), int(round(y / pxsz))], dtype=np.uint)


def computerecondensity(d3d, label, leafs=16, PIX=10, IMGMAX=40000):
    """
    Compute Local Effective Resolution
    :param d3d:  3D points
    :param label:
    :param leafs: Number of leafs (SNR = np.sqrt(leafs/2))
    :param PIX: nm to pixel (e.g. 10nm per pixel = 100nm^2
    :param IMGMAX: Upper limit of the image ROI
    :return: The CKDTree of the points (with leaf size), the image array, the points in leafs, and the image array where im[pix_x, pix_y] = LER
    """
    points = d3d[:, :2].copy()
    mx, MX = np.min(points[:, 0]), np.max(points[:, 0])
    my, MY = np.min(points[:, 1]), np.max(points[:, 1])
    if mx < 0:
        points[:, 0] += np.abs(mx)
        print("Neg pos for {}".format(label))
    if my < 0:
        points[:, 1] += np.abs(my)
        print("Neg pos for {}".format(label))

    tr = cKDTree(points, leafsize=leafs)
    lfs = []
    root = tr.tree
    leaves(root, lfs)

    assert (MX < IMGMAX)
    assert (MY < IMGMAX)
    pixels = np.empty((len(lfs), 2), dtype=np.float64)
    imgsize = int(np.round(IMGMAX / PIX))
    imarray = np.zeros((imgsize, imgsize), dtype=np.float32)
    for ip, p in enumerate(lfs):
        ps = points[p]
        #         assert(ps.shape[0] <= leafs)
        assert (ps.shape[1] == 2)
        _mx = min(ps[:, 0])
        _Mx = max(ps[:, 0])
        _my = min(ps[:, 1])
        _My = max(ps[:, 1])
        _x, _y = topix(_mx, _my, PIX)  # Lower left
        _X, _Y = topix(_Mx, _My, PIX)  # Upper right
        # ELR = nr of points over area xy - XY
        #         imarray[_x:_X, _y:_Y] += ps.shape[0]
        pixels[ip, 0] = ((_X - _x) + 1) * ((_Y - _y) + 1)  # Superpixel area
        pixels[ip, 1] = ps.shape[0] / pixels[ip, 0]  # Locs / px^2
        for pri in range(ps.shape[0]):
            pti = ps[pri, :2]
            xp, yp = topix(pti[0], pti[1], PIX)
            assert (xp < IMGMAX)
            assert (
                        yp < IMGMAX)  # You'd assume this check is not nec. as an index error would follow, but you'd be wrong in interesting cases
            try:
                imarray[xp, yp] += pixels[ip, 1]
            except OverflowError as e:
                print("OE {} , {} <- {} {}".format(xp, yp, pixels[ip, 1], pixels[ip, 0]))
                print("OE {} , {} <- {:.2f} {:.2f}".format(xp, yp, pti[0], pti[1]))
                raise
    return tr, imarray, pixels


def computeSNRLE(leafs=2, pix=10, MAX=60000, data=None, outpath="."):
    """
    Compute the local effective resolution for data
    :param leafs: nr of leafs (SNR = sqrt(leafs/2)
    :param pix: Pixel size in nm (e.g. 10nm/pixel = 100nm^2 in 2D pixels
    :param MAX: Maximum image ROI range in pixels
    :param data: Structured dict {cell : {channel : reader.points}}
    :param outpath: writeable directory to create tiff files in
    :return: {"cell_channel" : (pixels, imagearray)} where pixels are the nonnegative pixels with LRE value
    """
    lgr.info("Computing LRE for leaf size {}, {} nm/pixel".format(leafs, pix))
    if data is None:
        return
    pixelmap ={}
    SNR = np.sqrt(leafs/2)
    for cell in data:
        lgr.info("Cell {}".format(cell))
        for channel in data[cell]:
            lgr.info("Channel {}".format(channel))
            d3d = data[cell][channel].points
            label = "{}_{}".format(cell, channel)
            tr, imarray, pixels = computerecondensity(d3d, label, leafs, pix, MAX)
            img = Image.fromarray((imarray/np.max(imarray)*255).astype(np.uint8), mode='L')
            img.save(os.path.join(outpath, '{}_oct_{:.2f}.tiff'.format(label, SNR)))
#                 sns.distplot(pixels[:,1])
            pixelmap[label] = pixels, imarray
    return pixelmap