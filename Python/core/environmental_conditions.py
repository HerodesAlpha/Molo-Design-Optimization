"""
Environmental conditions module for wave statistics and spectra.

This module provides classes for short-term and long-term wave conditions,
including JONSWAP spectrum calculations and contour line generation.
"""

import numpy as np
from scipy.stats import lognorm, norm, weibull_min
from typing import Optional

from core.data_tables import world_wide_distribution_parameters as wwdp


class Short_Term_Wave_Conditions:
    """
    Short-term wave conditions with JONSWAP spectrum.
    
    Handles wave parameters (Hs, Tp, Tz) and spectrum calculations.
    """
    
    def __init__(
        self,
        hs: float,
        tp: Optional[float] = None,
        tz: Optional[float] = None,
        gamma: Optional[float] = None
    ) -> None:
        self._hs = hs
        self._tp = tp
        self._tz = tz
        self._gamma = gamma

        if (tp is None and tz is None) or (tp is not None and tz is not None):
            raise ValueError('Either Tp or Tz must be given (not both, not neither)')
        elif tp is not None:
            self._tz = self.tp2tz(self._tp, self._hs)
        else:
            self._tp = self.tz2tp(self._tz, self._hs)

        if self._gamma is None:
            self._gamma = self.gamma_func(self._tp, self._hs)

    def gamma_func(self, tp: float, hs: float) -> float:
        """
        Calculate JONSWAP peak enhancement factor.
        
        Args:
            tp: Peak period
            hs: Significant wave height
            
        Returns:
            Peak enhancement factor gamma
        """
        if self._gamma is not None:
            return self._gamma
        
        x = tp / np.sqrt(hs)
        if x <= 3.6:
            return 5.0
        elif x < 5:
            return np.exp(5.75 - 1.15 * x)
        else:
            return 1.0

    def s_jonswap(self, w: np.ndarray) -> np.ndarray:
        """
        Calculate JONSWAP wave spectrum.
        
        Args:
            w: Angular frequency array [rad/s]
            
        Returns:
            Spectral density array
        """
        sig_a = 0.07
        sig_b = 0.09
        delta_sig = sig_b - sig_a
        wp = 2 * np.pi / self._tp

        a_gamma = 1 - 0.287 * np.log(self._gamma)

        def spec_pm(w):
            v1 = (5 / 16) * (self._hs ** 2) * (wp ** 4) * (w ** (-5))
            v2 = np.exp(-(5 / 4) * ((w / wp) ** (-4)))
            return v1 * v2

        def spec_j(w):
            def sig(w):
                return sig_a if w <= wp else sig_b

            sig_ab = np.array(list(map(sig, w)))
            pm = spec_pm(w)

            v3 = np.exp(-0.5 * np.square((w - wp) / (sig_ab * wp)))
            return a_gamma * pm * np.power(self._gamma, v3)

        if self._gamma == 1:
            return spec_pm(w)
        else:
            return spec_j(w)

    def tp2tz(self, tp, hs):
        return (0.6673
                + 0.05037 * self.gamma_func(tp, hs)
                - 0.006230 * self.gamma_func(tp, hs) ** 2
                + 0.0003341 * self.gamma_func(tp, hs) ** 3) * tp

    def tz2tp(self, tz, hs):
        reltol = 0.001
        tp = tz
        while not np.isclose(self.tp2tz(tp, hs), tz, reltol):
            tp *= 1 + reltol
        return tp

    def expected_largest_maximum(self, h, w):
        s = self.s_jonswap(w)
        r = np.abs(h ** 2) * s
        dw = w[1] - w[0]
        sig_r_m0 = sum(r) * dw
        nz = 3 * 3600 / self._tz
        return np.sqrt(sig_r_m0) * (np.sqrt(2 * np.log(nz)) + 0.5772 / np.sqrt(2 * np.log(nz)))

    @property
    def hs(self):
        return self._hs

    @property
    def tp(self):
        return self._tp

    @property
    def tz(self):
        return self._tz

    @property
    def wp(self):
        return 2*np.pi/self._tp



    @property
    def gamma(self):
        return self._gamma

class Long_Term_Wave_Conditions:
    """
    Long-term wave conditions for contour line generation.
    
    Uses worldwide distribution parameters to generate environmental contours.
    """
    
    def __init__(self, area: int) -> None:
        """
        Initialize long-term wave conditions for a specific area.
        
        Args:
            area: Area number (1-104) corresponding to worldwide distribution parameters
        """
        self._alpha_s = wwdp[area - 1, 1]
        self._beta_s = wwdp[area - 1, 2]
        self._a1 = wwdp[area - 1, 3]
        self._a2 = wwdp[area - 1, 4]
        self._b1 = wwdp[area - 1, 5]
        self._b2 = wwdp[area - 1, 6]

    def logninv(self, x: float, mu: float, sigma: float) -> float:
        """
        Inverse lognormal distribution.
        
        Args:
            x: Probability value
            mu: Location parameter
            sigma: Scale parameter
            
        Returns:
            Inverse lognormal value
        """
        return lognorm(s=sigma, scale=np.exp(mu)).ppf(x)

    def contour_line(self, return_period: float) -> np.ndarray:
        """
        Generate environmental contour line for given return period.
        
        Args:
            return_period: Return period in years (statistics conditioned for 3hr storms)
            
        Returns:
            Array of [Hs, Tz] pairs defining the contour
        """
        pf = 1 / (return_period * 365 * 8)
        beta = norm.ppf((1 - pf), 0, 1)
        phi = np.linspace(-np.pi / 2, np.pi / 2, 10, endpoint=True)
        u1 = np.cos(phi) * beta
        u2 = np.sin(phi) * beta

        # hs = self._beta_s * (-np.log(1 - norm.cdf(u1, 0, 1))) ** (1 / self._alpha_s)
        x_hs = norm.cdf(u1, 0, 1)
        hs = weibull_min(self._beta_s).ppf(x_hs) * self._alpha_s
        x = norm.cdf(u2, 0, 1)
        mu = 0.70 + self._a1 * hs ** self._a2
        sigma = 0.07 + self._b1 * np.exp(self._b2 * hs)
        tz = self.logninv(x, mu, sigma)

        return np.transpose(np.vstack((hs, tz)))
