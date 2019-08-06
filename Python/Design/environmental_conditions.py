import numpy as np


class Short_Term_Wave_Conditions():
    def __init__(self, hs, tp, gamma=None):
        self._hs = hs
        self._tp = tp

        if gamma == None:
            self._gamma = self.gamma()
        else:
            self._gamma = gamma

        self._tz = self.tp2tz()

    def gamma(self):
        x = self._tp / np.sqrt(self._hs)
        if x <= 3.6:
            return 5
        elif x < 5:
            return np.exp(5.75 - 1.15 * x)
        else:
            return 1

    def s_jonswap(self, w):
        sig_a = 0.07
        sig_b = 0.09
        delta_sig = sig_b - sig_a
        wp = 2 * np.pi / self._tp

        a_gamma = 1 - 0.287 * np.log(self._gamma)

        def spec_pm(w):
            return (5 / 16) * (self._hs ** 2) * (wp ** 4) * (w ** (-5)) * np.exp(-(5 / 4) * ((w / wp) ** (-4)))

        def spec_j(w):
            def sig(w):
                return sig_a if w <= wp else sig_b

            sig_ab = np.array(list(map(sig, w)))

            return a_gamma * spec_pm(w) * self._gamma ** np.exp(-0.5 * ((w - wp) / sig_ab * wp))

        if self._gamma == 1:
            return spec_pm(w)
        else:
            return spec_j(w)

    def tp2tz(self):
        return (0.6673 + 0.05037 * self._gamma - 0.006230 * self._gamma ** 2 + 0.0003341 * self._gamma ** 3) * self._tp

    def expected_largest_maximum(self, h, w):
        r = np.abs(h ** 2) * self.s_jonswap(w)[:, np.newaxis]
        dw = w[1] - w[0]
        sig_r_m0 = sum(r) * dw
        nz = 3 * 3600 / self._tz
        return np.sqrt(sig_r_m0) * (np.sqrt(2 * np.log(nz)) + 0.5772 / np.sqrt(2 * np.log(nz)))
