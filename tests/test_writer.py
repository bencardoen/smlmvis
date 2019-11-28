import smlmvis.vtuwriter as v
import smlmvis.epflreader as e
import pandas as pd
import os
import filecmp
import urllib
import tempfile

def test_wrt():
    with tempfile.TemporaryDirectory() as indir:

        url = "http://bigwww.epfl.ch/smlm/challenge2016/datasets/MT0.N1.HD/sample/activations.csv"
        urllib.request.urlretrieve(url, os.path.join(indir, 'test.csv'))
        challenge_data = e.EPFLReader(os.path.join(indir, 'test.csv'))

        v.VtuWriter(os.path.join(indir, 'test'), challenge_data.points, challenge_data.values)

        url2 = 'http://vault.sfu.ca/index.php/s/NMhjJmaAWfw7sAa/download'
        urllib.request.urlretrieve(url2, os.path.join(indir, 'testref.vtu'))
        assert (filecmp.cmp(os.path.join(indir, 'test.vtu'), os.path.join(indir, 'testref.vtu'), shallow=False))
