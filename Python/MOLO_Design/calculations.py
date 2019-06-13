__author__ = "Eivind Sonju"
__copyright__ = "Copyright (C) 2017-2019 Verbun AS. All rights reserved."
__version__ = "2.0"

import numpy as np

from sea_loads import HydroCoefficients


class TransferFunctions(HydroCoefficients,object):
    def __init__(self,settings):
        super().__init__(settings)

    def rao(self, fe, m, ma, c, k, w):
        return np.absolute(fe / (-w ** 2 * (m + ma) + 1j * w * (c) + k))


    def get_rao(self, idof, idir):
        fe = self._fe[:, idir, idof]
        m = self._m[idof][idof]
        ma = self._ma[:, idof, idof]
        c = self._c_hyd[:, idof, idof]
        k = self._k[idof][idof]
        return self.rao(fe, m, ma, c, k, self._w)
