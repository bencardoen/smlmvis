import numpy as np
from scipy.spatial.ckdtree import cKDTree


def leaves(node, res):
    leaf = True
    if node.lesser:
        leaf = False
        leaves(node.lesser, res)
    if node.greater:
        leaf = False
        leaves(node.greater, res)
    if leaf:
        res.append(node.indices)


def topix(x, y, pxsz):
    assert (x >= 0)
    assert (y >= 0)
    return np.array([int(round(x / pxsz)), int(round(y / pxsz))], dtype=np.uint)


def computerecondensity(d3d, label, leafs=16, PIX=10, IMGMAX=40000):
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