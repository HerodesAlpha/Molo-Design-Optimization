import numpy as np
from scipy.stats import norm, lognorm,weibull_min
from data_tables import world_wide_distribution_parameters as wwdp


class Short_Term_Wave_Conditions():
    def __init__(self, hs, tp=None, tz=None, gamma=None):
        self._hs = hs
        self._tp = tp
        self._tz = tz

        if gamma != None:
            self._gamma = gamma
        else:
            self._gamma = None

        if (tp == None and tz == None) or (tp != None and tz != None):
            print('Either Tp or Tz must be given')
            exit()
        elif tp != None:
            self._tz = self.tp2tz(self._tp, self._hs)
        else:
            self._tp = self.tz2tp(self._tz, self._hs)

        if self._gamma == None:
            self._gamma = self.gamma(self._tp, self._hs)

    def gamma(self, tp, hs):
        if self._gamma != None:
            return self._gamma
        else:
            x = tp / np.sqrt(hs)
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

    def tp2tz(self, tp, hs):
        return (0.6673 + 0.05037 * self.gamma(tp, hs) - 0.006230 * self.gamma(tp, hs) ** 2 + 0.0003341 * self.gamma(tp,
                                                                                                                    hs) ** 3) * tp

    def tz2tp(self, tz, hs):
        reltol = 0.01
        tp = tz
        while not np.isclose(self.tp2tz(tp, hs), tz, reltol):
            tp *= 1 + reltol
        return tp

    def expected_largest_maximum(self, h, w):
        r = np.abs(h ** 2) * self.s_jonswap(w)
        dw = w[1] - w[0]
        sig_r_m0 = sum(r) * dw
        nz = 3 * 3600 / self._tz
        return np.sqrt(sig_r_m0) * (np.sqrt(2 * np.log(nz)) + 0.5772 / np.sqrt(2 * np.log(nz)))


class Long_Term_Wave_Conditions():
    def __init__(self, area):
        self._alpha_s = wwdp[area - 1, 1]
        self._beta_s = wwdp[area - 1, 2]
        self._a1 = wwdp[area - 1, 3]
        self._a2 = wwdp[area - 1, 4]
        self._b1 = wwdp[area - 1, 5]
        self._b2 = wwdp[area - 1, 6]

    def logninv(self, x, mu, sigma):
        return lognorm(s=sigma, scale=np.exp(mu)).ppf(x)

    def contour_line(self, return_period):  # Return period in years, statistics conditioned for 3hr storms
        pf = 1 / (return_period * 365 * 8)
        beta = norm.ppf((1 - pf), 0, 1)
        phi = np.linspace(0, 2 * np.pi, 100, endpoint=True)
        u1 = np.cos(phi) * beta
        u2 = np.sin(phi) * beta

        # hs = self._beta_s * (-np.log(1 - norm.cdf(u1, 0, 1))) ** (1 / self._alpha_s)
        x_hs = norm.cdf(u1, 0, 1)
        hs = weibull_min(self._beta_s).ppf(x_hs)*self._alpha_s
        x = norm.cdf(u2, 0, 1)
        mu = 0.70 + self._a1 * hs ** self._a2
        sigma = 0.07 + self._b1 * np.exp(self._b2 * hs)
        tz = self.logninv(x, mu, sigma)

        return hs, tz
