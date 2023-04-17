import sys

import astropandas as apd
import pandas as pd


infiles = sys.argv[1:-1]
outfile = sys.argv[-1]

indata = [apd.read_auto(path) for path in infiles]
outdata = pd.concat(indata, axis=1)
apd.to_auto(outdata, outfile)
