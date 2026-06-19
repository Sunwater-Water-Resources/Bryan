from scipy import optimize
import numpy as np
from scipy.special import ndtri, ndtr


class Curve:
    def __init__(self):
        self.aep_2000 = 0
        self.aep_1000 = 0
        self.pmp_depth = 0
        self.pmp_aep = 0
        self.quantiles = None
        self.aeps = None
        self.validity = True

    def setup_rainfall_boundaries(self, aep_2000, aep_1000, pmp_depth, pmp_aep):
        self.aep_2000 = aep_2000
        self.aep_1000 = aep_1000
        self.pmp_depth = pmp_depth
        self.pmp_aep = pmp_aep
        self.quantiles = np.array([self.aep_1000, self.aep_2000, self.pmp_depth])
        self.aeps = np.array([1000.0, 2000.0, self.pmp_aep])

    def test_validity(self):
        return True


class CoercedQuadratic(Curve):
    def __init__(self, method='SiriwardenaWeinmann1998'):
        super(CoercedQuadratic, self).__init__()
        self.method = method
        self.method_exception = f'The {self.method} method is not recognised, use either: SiriwardenaWeinmann1998 | HillAndOthers2000'
        self.Sgc = 0
        self.xd = 0
        self.Sgap = 0
        self.a1 = 0
        self.a2 = 0
        self.inflection_aep = 0

    def fit_curve(self):
        # Get the transformed AEPs
        if self.method == 'SiriwardenaWeinmann1998':
            Xpmp = np.log10(self.pmp_aep)
            Xy1 = np.log10(1000)
            Xy2 = np.log10(2000)
        elif self.method == 'HillAndOthers2000':
            Xpmp = ndtri(1-1/self.pmp_aep)
            Xy1 = ndtri(1-1/1000)
            Xy2 = ndtri(1-1/2000)
        else:
            raise Exception(self.method_exception)

        # Get the standardised quantiles
        Rpmp = np.log10(self.pmp_depth) / np.log10(self.aep_2000)
        Ry1 = np.log10(self.aep_1000) / np.log10(self.aep_2000)
        Ry2 = 1.0

        # Get the slopes
        self.Sgc = (Ry2 - Ry1) / (Xy2 - Xy1)
        self.xd = Xpmp - Xy2
        self.Sgap = (Rpmp - Ry2) / self.xd

        # Get quadratic parameters
        self.a1 = self.Sgc * self.xd
        self.a2 = (self.Sgap - self.Sgc) * self.xd

    def get_quantile(self, aep):
        Ty2 = np.log10(self.aep_2000)  # Transformed depth for 1 in 2000
        if self.method == 'SiriwardenaWeinmann1998':
            Xy2 = np.log10(2000)
            Xy = np.log10(aep)
        elif self.method == 'HillAndOthers2000':
            Xy2 = ndtri(1-1/2000)
            Xy = ndtri(1-1/aep)
        else:
            raise Exception(self.method_exception)
        x = Xy - Xy2
        Ry = 1 + self.a1 * (x / self.xd) + self.a2 * (x / self.xd) ** 2
        Ty = Ry * Ty2
        return 10 ** Ty

    def test_validity(self):
        if not self.Sgc <= 2 * self.Sgap:
            self.validity = False
            print(f'Failed test... Sgc of {self.Sgc} > 2Sgap of {2 * self.Sgap}')
        else:
            self.validity = True
        return self.validity

    def test_validity_by_aep(self, aep=-99):
        if aep < 0:
            aep = self.pmp_aep
        # Compute the inflection point
        inflection_point = -self.a1 / (2 * self.a2) * self.xd
        if self.method == 'SiriwardenaWeinmann1998':
            self.inflection_aep = 10 ** (inflection_point + np.log10(2000))
        elif self.method == 'HillAndOthers2000':
            inflection_p = inflection_point + ndtri(1-1/2000)
            self.inflection_aep = 1 / (1 - ndtr(inflection_p))
        print(f'The fitted curve has an inflection at an AEP of 1 in {np.around(self.inflection_aep, 0)}')
        # Compare the aep
        self.validity = (aep < self.inflection_aep)
        print(f'Comparing this against an AEP of 1 in {aep}: {self.validity}')
        return self.validity


class GEV(Curve):
    def __init__(self):
        super(GEV, self).__init__()
        self.location = 0
        self.scale = 0
        self.shape = 0

    def fit_curve(self):
        # print(f'Fitting GEV to aeps:', self.aeps)
        solution = optimize.root(self.root_shape, x0=0.1)
        # print(solution)
        self.shape = solution.x[0]

    def root_shape(self, shape):
        q = 1 - (-np.log(1 - 1 / self.aeps)) ** shape
        q_1_2 = self.quantiles[0] - self.quantiles[1]
        q_1_3 = self.quantiles[0] - self.quantiles[2]
        q_div = q_1_2 / q_1_3
        numerator = q[0] - q[1]
        denominator = q[0] - q[2]
        root = numerator / denominator - q_div
        self.scale = q_1_2 * shape[0] / numerator
        self.location = self.quantiles[0] - (self.scale / shape[0]) * q[0]
        return root

    def get_quantile(self, aep):
        factor = 1 - (-np.log(1 - 1 / aep)) ** self.shape
        quantile = self.location + self.scale / self.shape * factor
        return quantile

    def print_parameters(self):
        location = np.around(self.location, 1)
        scale = np.around(self.scale, 2)
        shape = np.around(self.shape, 3)
        print('\nGEV curve parameters:')
        print(f'Location: {location} | Scale: {scale} | Shape: {shape}')
