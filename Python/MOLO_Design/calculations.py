__author__ = "Eivind Sonju"
__copyright__ = "Copyright (C) 2017-2019 Verbun AS. All rights reserved."
__version__ = "2.0"

import numpy as np
import scipy

from loads import Sea_and_Inertia_Loads


class TransferFunctions(Sea_and_Inertia_Loads, object):
    def __init__(self, settings):
        super().__init__(settings)

    # def rao(self, fe, m, ma, c, k, w):
    #     return np.absolute(fe / (-w ** 2 * (m + ma) + 1j * w * (c) + k))
    #
    #
    # def get_rao(self, idof, idir):
    #     fe = self._fe[:, idir, idof]
    #     m = self._m[idof][idof]
    #     ma = self._ma[:, idof, idof]
    #     c = self._c_hyd[:, idof, idof]
    #     k = self._k[idof][idof]
    #     return self.rao(fe, m, ma, c, k, self._w)

    def get_h(self, idir):
        fe = self._fe[:, idir, :]
        m = self._m
        ma = self._ma
        c = self._c_hyd
        k = self._k
        return self.H(fe, m, ma, c, k, self._w)

    def H(self, f, m, ma, c, k, vw):
        container = np.zeros([len(vw), 6], dtype=complex)
        for i, w in enumerate(vw):
            this_ma=ma[i, :, :]
            denom = np.asarray(-w ** 2 * (m + this_ma) + 1j * w * c[i, :, :] + k, dtype=complex)
            daf = scipy.linalg.inv(denom)
            # daf =np.asarray([[1/x for x in col]for col in denom])
            x = daf @ f[i, :]
            container[i, :] = x

    # def H(self, f, m, ma, c, k, vw):
    #     container = np.zeros([len(vw), 6], dtype=complex)
    #     for iw, w in enumerate(vw):
    #         for i in [2,4]:
    #             this_ma=ma[iw, i, i]
    #             denom = np.asarray(-w ** 2 * (m[i,i] + ma[iw, i, i]) + 1j * w * c[iw, i, i] + k[i,i], dtype=complex)
    #             daf = 1/denom
    #             # daf =np.asarray([[1/x for x in col]for col in denom])
    #             x = daf * f[iw, i]
    #             container[iw, i] = x

        return container

