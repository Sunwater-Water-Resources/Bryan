def set_aep_of_pmp(self, area):
        if area < 100:  # 100 km²
            aep_of_pmp = 1e-7
        elif area > 1e5:  # 100,000 km²
            aep_of_pmp = 1e-4
        else:
            aep_of_pmp = 10 ** (np.log10(area) - 9)

        aep_of_pmp = 1 / aep_of_pmp  # Express as 1 in Y
        self.aep_of_pmp = aep_of_pmp