import numpy as np


class rv_histogramdd:
    """
    Takes the output of numpy.histogramdd (bins, density) and returns an object that
    provides a pdf(x) function.
    Assumes that numpy.histogramdd was called with density=True 
    """

    def __init__(self, histogram):
        self._hpdf = histogram[0]
        self._hbins = histogram[1]
        self.nvars = len(self._hbins)
        self.nbins = [len(self._hbins[i]) - 1 for i in range(self.nvars)]

    def pdf(self, x):
        if len(x) != self.nvars:
            raise ValueError(
                f"Sample has wrong dimensions {len(x)}, should be {self.nvars}"
            )
        # Find index in histogram array by searching bin array for x
        idx = np.zeros(self.nvars, dtype=int)
        for i in range(self.nvars):
            idx[i] = np.searchsorted(self._hbins[i], x[i]) - 1  # get bin
        if (min(idx) < 0) or (
            np.any(idx >= self.nbins)
        ):  # at least one idx is outside of bin bounds
            return 0  # zero pdf
        return self._hpdf[tuple(idx)]

    def maxp(self):
        return self._hpdf.max()

    def minp(self):
        return self._hpdf.min()

