import os.path
import random
from pathlib import Path
import scipy.stats
from scipy.special import ndtri, ndtr
from scipy import stats
from scipy import interpolate
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


class SampleScheme:
    """
    This class handles the Monte Carlo scheme. The main dataframe is a log of the samples generated for
    each simulation.
    """
    def __init__(self, lower_aep, upper_aep, number_of_main_divisions, number_of_sub_divisions,
                 number_of_temporal_patterns=10, sample_method='random', output_folder=''):

        self.aep_of_pmp = None
        self.m = number_of_main_divisions
        self.n = number_of_sub_divisions
        self.number_of_temporal_patterns = number_of_temporal_patterns
        self.output_folder = output_folder
        self.df = pd.DataFrame(index=range(self.m * self.n),
                               columns=['m', 'n', 'rain_z', 'rain_aep', 'mean_rain_mm', 'tp', 'storm_method',
                                        'tp_frequency', 'il_p', 'il_scaling', 'preburst_p', 'preburst_proportion',
                                        'preburst_mm', 'initial_loss', 'cl_p', 'cl_scaling', 'continuing_loss',
                                        'residual_depth', 'lake_z', 'embedded_bursts', 'ADV',
                                        'inflow', 'level', 'outflow'])

        # container for the tpt results
        self.quantiles = {'inflow': pd.DataFrame(columns=['probability', 'inflow']),
                          'level': pd.DataFrame(columns=['probability', 'level']),
                          'outflow': pd.DataFrame(columns=['probability', 'outflow'])}

        self.upper_aep = upper_aep
        self.lower_aep = lower_aep
        self.sample_method = sample_method

        self.z_upper = ndtri(1 - 1 / self.upper_aep)
        self.z_lower = ndtri(1 - 1 / self.lower_aep)

        print(f'\n Rainfall sampling from 1 in {self.lower_aep} AEP to 1 in {self.upper_aep} AEP')

        print(f'\n Setting up {self.m} main divisions from {self.z_lower} to {self.z_upper}')
        self.main_divisions = np.linspace(self.z_lower, self.z_upper, self.m + 1)
        self.p_divisions = ndtr(self.main_divisions)
        print('\nExceedance probability divisions:')
        print(self.p_divisions)

        # Create empty sample spaces
        self.rain_samples = pd.DataFrame()
        self.temporal_pattern_samples = pd.DataFrame()

    def random_std_norm_variates(self):
        n = self.m * self.n
        # numbers = np.random.rand(n)
        # sample = ndtri(numbers)
        rng = np.random.default_rng()
        sample = rng.standard_normal(n)
        # sample = np.random.randn(n)
        return sample

    def random_std_norm_variates_truncated(self, lower_p=0.0001, upper_p=0.9999):
        # This is a thought - not used for anything
        low_z = ndtri(lower_p)
        high_z = ndtri(upper_p)
        n = self.m * self.n
        rv = stats.truncnorm(low_z, high_z)
        sample = rv.rvs(n)
        return sample

    def sample_percentiles(self, lower=0.0, upper=1.0, loc=0.5, scale=0.34135, method='uniform'):
        n = self.m * self.n
        if method == 'normal':
            a, b = (lower - loc) / scale, (upper - loc) / scale
            rv = stats.truncnorm(a, b, loc=loc, scale=scale)
            sample = rv.rvs(n)
        else:
            sample = np.random.uniform(low=0.0, high=1.0, size=n)
            if upper < 1.0:
                sample = np.where(sample < upper, sample, upper)
            if lower > 0.0:
                sample = np.where(sample > lower, sample, lower)
        return sample

    def setup_rain_sample_space(self):
        # Set up the sample space matrix for rainfall depths
        print('\nSetting up the sample space matrix for rainfall depths')
        print(f'Using the {self.sample_method} sampling method within each main division.')
        sample_space = np.zeros(shape=(self.m, self.n))
        for i in range(self.m):
            div_0 = self.main_divisions[i]
            div_1 = self.main_divisions[i + 1]
            if self.sample_method == 'stratified':
                sample_space[i] = np.linspace(div_0, div_1, self.n, endpoint=False)
            elif self.sample_method == 'normally distributed':
                rv = stats.truncnorm(div_0, div_1)
                sample_space[i] = rv.rvs(self.n)
            elif self.sample_method == 'uniformly distributed':
                sample_space[i] = np.random.rand(self.n) * (div_1 - div_0) + div_0
        self.rain_samples = pd.DataFrame(sample_space)
        # print(self.rain_samples)

    def get_rain_sample(self, main_division, sub_division):
        z = self.rain_samples.loc[main_division, sub_division]
        print(f'Sampling rainfall from {main_division} main division and {sub_division} sub division with z of {z}')
        return z

    def get_temporal_pattern_sample(self, front_w, front_patterns, m=3):
        # the 'm' parameter is the number of more forward loaded events
        # used in the weightings: top third (3 out of 10) was used in the Sharpe analysis
        # returns integer between 0 and 9
        if front_w < 0.1001:
            sample = np.random.randint(self.number_of_temporal_patterns)
        else:
            # back_w = (self.number_of_temporal_patterns - 3 * front_w) / (self.number_of_temporal_patterns - m)
            back_w = (1 - 3 * front_w) / (self.number_of_temporal_patterns - m)
            # get the set to sample from
            patterns = []
            probabilities = []
            check = 0.0
            for i in range(self.number_of_temporal_patterns):
                patterns.append(i)
                pattern_number = i + 1
                if pattern_number in front_patterns:
                    probabilities.append(front_w)
                    check += front_w
                else:
                    check += back_w
                    probabilities.append(back_w)
            # print('Front pattern numbers', front_patterns)
            # print(f'Front weighting: {front_w} | Back weighting: {back_w} | Sum check: {check}')
            # print(patterns)
            # print(probabilities)
            sample = np.random.choice(patterns, p=probabilities)
        return sample

    def store_simulations(self, filename='simulation_out.csv'):
        # filename = Path(filename).stem + '__mcdf.csv'                # Add common suffix
        # output_path = os.path.join(self.output_folder, filename)
        filename = f'{filename}__mcdf.csv'
        try:
            os.makedirs(os.path.dirname(filename), exist_ok=True)
            print('Storing the monte carlo method analysis file:', filename)
            self.df.to_csv(filename)
        except IOError:
            input("Could not save the simulation file. The file may be open in Excel. Please close the file and press enter.")
            self.df.to_csv(filename)

    def compute_quantiles(self, start_q, end_q, step_q, result_type, output_filename):
        number_of_q = int((end_q - start_q) / step_q) + 1
        quantiles = np.linspace(start=start_q, stop=end_q, num=number_of_q)
        compute_cols = ['z_min', 'z_max', 'p_min', 'p_max', 'pMi', 'num', 'pH', 'pH x pMi']
        compute_df = pd.DataFrame(index=range(self.m), columns=compute_cols)
        compute_df['z_max'] = self.main_divisions[1:]
        compute_df['z_min'] = self.main_divisions[0:self.m]
        compute_df['p_min'] = ndtr(compute_df['z_min'])
        compute_df['p_max'] = ndtr(compute_df['z_max'])
        compute_df['pMi'] = compute_df['p_max'] - compute_df['p_min']
        print(compute_df)
        probability = []
        for quantile in quantiles:
            for div in range(self.m):
                div_df = self.df[self.df.m == div]
                div_df.set_index('n', inplace=True)
                # print(f'\nProcessing division {div}')
                # print(div_df)
                count = 0
                for sim in range(self.n):
                    result = div_df[result_type].loc[sim]
                    # print(f'\nDiv: {div} | Sim: {sim} | Quantile: {quantile} | Threshold: {result}')
                    if quantile < result:
                        count += 1
                compute_df.loc[div, 'num'] = count

            compute_df['pH'] = compute_df['num'] / self.n
            compute_df['pH x pMi'] = compute_df['pH'] * compute_df['pMi']
            # print(f'\n Quantile: {quantile} | Threshold: {result}')
            # print(compute_df)
            quant_p = compute_df['pH x pMi'].sum()
            probability.append(quant_p)

        # self.quantile_df = pd.DataFrame(columns=['probability', result_type])
        quant_df = self.quantiles[result_type]
        quant_df['probability'] = probability
        quant_df[result_type] = quantiles
        quant_df['aep (1 in x)'] = 1/quant_df['probability']
        # output_path = os.path.join(self.output_folder, output_filename)
        print('\nWriting analysis results:', output_filename)
        quant_df.to_csv(output_filename)
        print(quant_df)

    def smooth_percentiles(self, percentiles, output_filename):
        print('\nDeveloping smooth models of the percentiles')
        # Set up the dataframe container
        std_aeps = np.array(self.get_standard_aeps())
        std_z = ndtri(1 - 1 / std_aeps)
        column_titles = percentiles.columns.values.tolist()
        lower_title = column_titles[2]
        upper_title = column_titles[3]
        df = pd.DataFrame(index=std_aeps, columns=['z', lower_title, upper_title])
        df.index.name = 'AEP'
        df['z'] = std_z

        # Build the lower bound model
        lower_df = percentiles[['z', 'AEP', lower_title]]
        lower_df = lower_df.loc[lower_df[lower_title] > 0]
        z = lower_df['z'].to_numpy()
        y = np.log10(lower_df[lower_title].to_numpy())
        lower_model = np.poly1d(np.polyfit(z, y, 6))

        # Build the upper bound model
        upper_df = percentiles[['z', 'AEP', upper_title]]
        upper_df = upper_df.loc[upper_df[upper_title] > 0]
        z = upper_df['z'].to_numpy()
        y = np.log10(upper_df[upper_title].to_numpy())
        upper_model = np.poly1d(np.polyfit(z, y, 6))

        # Get the estimates
        percentiles.set_index('AEP', inplace=True)
        for ind, z in enumerate(std_z):
            aep = std_aeps[ind]

            # Get the lower bound estimate
            try:
                result = np.around(10 ** lower_model(z), 3)
            except:
                result = current_lower
            if np.isinf(result) or result > 1000000:
                df.loc[aep, lower_title] = 0
            else:
                df.loc[aep, lower_title] = result

            # Get the upper bound estimate
            try:
                result = np.around(10 ** upper_model(z), 3)
            except:
                result = current_upper
            if np.isinf(result) or result > 1000000:
                df.loc[aep, upper_title] = 0
            else:
                df.loc[aep, upper_title] = result

        df.to_csv(output_filename)
        df.reset_index(inplace=True)
        print(df)
        return df
        
    def main_division_percentiles(self, result_type, output_filename, percentiles=(5, 95)):
        all_intervals = []
        std_nrm_vars = []
        header = [f'{i}th percentile' for i in percentiles]

        # extract the percentiles
        for division in range(self.m):
            estimates = self.df.loc[self.df['m'] == division, result_type]
            estimates = estimates.to_numpy()
            intervals = np.percentile(estimates, percentiles)
            all_intervals.append(intervals)
            std_nrm_vars.append((self.main_divisions[division] + self.main_divisions[division + 1]) / 2)

        # store the results in a dataframe
        df = pd.DataFrame(std_nrm_vars, columns=['z'])
        aeps = 1 / (1 - ndtr(std_nrm_vars))
        df['AEP'] = aeps
        df[header] = all_intervals
        print(f'\nConfidence intervals for {result_type}')
        print(df)
        df.to_csv(output_filename, index=False)
        return df

    def compute_std_quantiles(self, result_type, output_filename=None, ):
        # m = self.m
        mcdf = self.df
        
        # compute_cols = ['z_min', 'z_max', 'p_min', 'p_max', 'pMi', 'num', 'pH', 'pH x pMi']
        # compute_df = pd.DataFrame(index=range(m), columns=compute_cols, dtype=float)
        # compute_df['z_max'] = self.main_divisions[1:]
        # compute_df['z_min'] = self.main_divisions[0:m]
        # compute_df['p_min'] = ndtr(compute_df['z_min'])
        # compute_df['p_max'] = ndtr(compute_df['z_max'])
        # compute_df['pMi'] = compute_df['p_max'] - compute_df['p_min']
        
        # # Special treatment for first and last intervals as per Nathan and Weinmann 2013
        # compute_df.loc[0, 'pMi'] = compute_df.loc[0, 'p_max']           # Non-exceedance probability of upper bound
        # compute_df.loc[m-1, 'pMi'] = 1 - compute_df.loc[m-1, 'p_min']   # Exceedance probability of lower bound
        
        # print(compute_df)
        
        # mcdf[f'{result_type}_aep'] = 0.0
        # for i in mcdf.index:
        #     value = mcdf.loc[i, result_type]
            
        #     compute_df['num'] = mcdf.groupby('m').apply(lambda div: sum(div[result_type] > value))
        #     compute_df['pH'] = compute_df['num'] / self.n
        #     compute_df['pH x pMi'] = compute_df['pH'] * compute_df['pMi']
            
        #     mcdf.loc[i, f'{result_type}_aep'] = compute_df['pH x pMi'].sum()

        # Get aep of quantiles for each simulation
        print(f'\nRunning TPT analysis for {result_type}...')
        tpt = TotalProbTheorem(self.m, self.n, self.main_divisions, mcdf[['m', result_type]])
        # tpt.tpt_limits = {'lower': mcdf.loc[mcdf['m'] == 0, result_type].mean(),
        #                   'upper': mcdf.loc[mcdf['m'] == self.m - 1, result_type].mean()}
        tpt.tpt_limits = {'lower': self.lower_aep * 0.98, 'upper': self.upper_aep * 1.02}
        print('Setting limits for the TPT analysis:', tpt.tpt_limits)
        mcdf[f'{result_type}_aep'] = 0.0
        mcdf[f'{result_type}_aep'] = mcdf[result_type].apply(tpt.assign_aep, result_type=result_type)
        # mcdf.dropna(subset=[f'{result_type}_aep'], inplace=True)
        print(mcdf)
        
        # Get the standard set of quantiles
        std_aeps = self.get_standard_aeps()
        std_aeps_frac = 1 / np.array(std_aeps)
        ascending = mcdf[[f'{result_type}_aep', result_type]].sort_values(f'{result_type}_aep')
        # quantiles = np.interp(std_aeps_frac,
        #                       ascending[f'{result_type}_aep'],
        #                       ascending[result_type],
        #                       left=np.nan,
        #                       right=np.nan)

        # change the interpolation to being in log-normal space
        std_norm_var = ndtri(1 - std_aeps_frac)
        ascending['z'] = ndtri(1 - ascending[f'{result_type}_aep'])
        ascending.sort_values('z', inplace=True)
        quantiles = np.interp(std_norm_var,
                              ascending['z'],
                              np.log10(ascending[result_type]),
                              left=None,
                              right=np.nan)
        quantiles = 10 ** quantiles

        quant_df = pd.DataFrame(std_aeps, columns=['aep (1 in x)'])
        # quant_df['aep (1 in x)'] = std_aeps
        quant_df['probability'] = std_aeps_frac
        quant_df[result_type] = quantiles
        # quant_df = quant_df[['aep (1 in x)', 'probability', result_type]]       # Re-arrange columns
        
        if output_filename:
            print('\nWriting analysis results:', output_filename)
            quant_df.to_csv(output_filename, index=False)
        print(quant_df)
        self.quantiles[result_type] = quant_df
        
        return

    def get_standard_aeps(self):
        # Create a list of standard AEPS, e.g.: 2, 5, 10, 20, ...
        lower_aep = self.lower_aep
        if lower_aep < 2:  # limit the lower AEP to no less than 1 in 2
            lower_aep = 2
        aep = lower_aep
        std_aeps = []
        # while aep < self.upper_aep/4:
        while aep <= self.upper_aep:
            # if aep >= 4 * self.lower_aep:
            if aep >= lower_aep:
                std_aeps.append(aep)

            if np.log10(aep / 2) % 1 == 0.0:  # is multiple of 20
                aep = aep * 5 // 2
            else:
                aep = aep * 2

        if self.aep_of_pmp is not None:
            if self.aep_of_pmp < self.upper_aep:
                std_aeps.append(self.aep_of_pmp)
                std_aeps.sort()

        return std_aeps

    def plot_tpt_results_2(self, result_type, output_filename, show_fig=False, percentiles=None):
        fig, ax = plt.subplots()

        # inserting the main data points
        aeps = self.df['rain_aep'].to_numpy()
        z = ndtri(1 - 1 / aeps)
        results = self.df[result_type]
        ax.plot(z, results, 'o', alpha=0.25, markeredgewidth=0, markersize=2.5)

        # Inserting the percentiles if they were provided
        if percentiles is not None:
            column_titles = percentiles.columns.values.tolist()
            upper_title = column_titles[3]
            lower_title = column_titles[2]
            percentiles = percentiles.astype(float)
            # z = percentiles['z'].to_numpy()
            # lower = percentiles[lower_title].to_numpy()
            # upper = percentiles[upper_title].to_numpy()
            # print(percentiles[['z', lower_title, upper_title]])
            # input()
            # ax.fill_between(z, lower, upper, alpha=0.2, linewidth=1, color='k')
            ax.fill_between(percentiles['z'],
                            percentiles[lower_title],
                            percentiles[upper_title],
                            alpha=0.35, linewidth=1, color='k')
            ax.plot(percentiles['z'], percentiles[lower_title], '--', color='k', alpha=0.9)
            ax.plot(percentiles['z'], percentiles[upper_title], '--', color='k', alpha=0.9)

        # Inserting the TPT result
        df = self.quantiles[result_type].replace(0, 1)
        x = df['aep (1 in x)']
        z = ndtri(1 - 1 / x)
        y = df[result_type]
        ax.plot(z, y, '-', color='k')

        # Formatting the plot
        min_flow = y.min()
        if min_flow > 1000:
            bottom = 1000
        elif min_flow > 100:
            bottom = 100
        elif min_flow > 10:
            bottom = 10
        else:
            bottom = 1
        if not result_type == 'level':
            ax.set_yscale('log')
            ax.set_ylim(bottom=bottom)
        plot_labels = {'inflow': 'Inflow (m³/s)',
                       'outflow': 'Outflow (m³/s)',
                       'level': 'Lake level (m AHD)'}
        if result_type[0: 3] == 'Vol':
            duration = result_type[3: 6]
            duration = duration.replace('h', '')
            plot_labels = {result_type: f'{duration} hr Volume (m³)'}
        std_aeps = np.array(self.get_standard_aeps())
        std_z = ndtri(1 - 1 / std_aeps)
        ax.set_xticks(std_z)
        ax.set_xticklabels(std_aeps, rotation=90)
        ax.set_xlabel("AEP (1 in X)")
        ax.set_ylabel(plot_labels[result_type])
        plt.tight_layout()
        if show_fig:
            plt.show()
        print('\nWriting analysis plots:', output_filename)
        plt.savefig(output_filename)
        plt.clf()

    def plot_tpt_results(self, result_type, output_filename, show_fig=False):
        df = self.quantiles[result_type]
        ax1 = df.plot(x='aep (1 in x)', y=result_type, logx=True, color='black')
        mc_result = self.df[['rain_aep', result_type]]
        # mc_result['aep'] = 1 / (1 - mc_result['rain_p'])
        mc_result.plot(ax=ax1, x='rain_aep', y=result_type, logx=True, kind='scatter', color='b', alpha=0.2)
        if show_fig:
            plt.show()
        # plt.savefig(os.path.join(self.output_folder, output_filename))
        print('\nWriting analysis plots:', output_filename)
        plt.savefig(output_filename)


class TotalProbTheorem:
    # See Table 4.4.3 in Chapter 4 of Book 4 in ARR2019.
    # Though the ARR example does not get the pMi for the upper and lower bound correct.
    def __init__(self, m, n, main_divisions, mcdf):
        self.m = m 
        self.n = n
        self.mcdf = mcdf
        self.tpt_limits = None

        # computation dataframe for sampled domain
        compute_cols = ['z_min', 'z_max', 'p_min', 'p_max', 'pMi', 'num', 'pH', 'pH x pMi']
        compute_df = pd.DataFrame(index=range(m), columns=compute_cols, dtype=float)
        compute_df['z_max'] = main_divisions[1:]
        compute_df['z_min'] = main_divisions[0:m]
        compute_df['p_min'] = ndtr(compute_df['z_min'])
        compute_df['p_max'] = ndtr(compute_df['z_max'])
        compute_df['pMi'] = compute_df['p_max'] - compute_df['p_min']

        # Special treatment for first and last intervals (unsampled domain) as per Nathan and Weinmann 2013
        edge_df = pd.DataFrame(index=[-1, m], columns=compute_cols, dtype=float)
        edge_df.loc[-1, 'pMi'] = compute_df.loc[0, 'p_min']  # Non-exceedance probability of upper bound
        edge_df.loc[m, 'pMi'] = 1 - compute_df.loc[m-1, 'p_max']   # Exceedance probability of lower bound

        check_space = compute_df['pMi'].sum() + edge_df['pMi'].sum()
        print('\nSample Space:')
        print(compute_df)
        print('\nEdge space:')
        print(edge_df)
        print(f'\nCheck on the probability space (should be 1.0):', check_space)

        # Removed by RGS: treating the first and last intervals as being on either side of the sampled domain.
        # compute_df.loc[0, 'pMi'] = compute_df.loc[0, 'p_max']           # Non-exceedance probability of upper bound
        # compute_df.loc[m-1, 'pMi'] = 1 - compute_df.loc[m-1, 'p_min']   # Exceedance probability of lower bound
        
        # print(compute_df)
        self.compute_df = compute_df
        self.edge_df = edge_df

        self.upper_factor_assumption = edge_df.loc[-1, 'pMi']
        print('\nAn assumed factor is needed to account for the upper side of the unsampled probability space...')
        print(f'Assuming a factor of {np.around(self.upper_factor_assumption, 3)} of the bounding interval\n')
    
    def assign_aep(self, peak_value, result_type):
        # m = self.m
        n = self.n
        compute_df = self.compute_df
        edge_df = self.edge_df
        compute_df['num'] = self.mcdf.groupby('m').apply(lambda div: sum(div[result_type] > peak_value))
        compute_df['pH'] = compute_df['num'] / self.n

        # # Special treatment for first and last intervals as per Nathan and Weinmann 2013
        # pQr0 = compute_df.loc[0, 'pH']
        # compute_df.loc[0, 'pH'] = np.sqrt(0.1 * pQr0 * pQr0)        # Geometric mean of (1 + 0.1 ) * p[Q>q|Ri*]
        # compute_df.loc[m-1, 'pH'] = np.sqrt(compute_df.loc[m-1, 'pH'])    # Geometric mean of p[Q>q|Ri] and 1.0
        # edge_df.loc[0, 'pH'] = np.sqrt(0.1 * compute_df.loc[0, 'pH'] ** 2)
        edge_df.loc[-1, 'pH'] = np.sqrt(self.upper_factor_assumption * compute_df.loc[0, 'pH'] ** 2)
        edge_df.loc[self.m, 'pH'] = np.sqrt(compute_df.loc[self.m-1, 'pH'] * 1)

        compute_df['pH x pMi'] = compute_df['pH'] * compute_df['pMi']
        edge_df['pH x pMi'] = edge_df['pH'] * edge_df['pMi']
        aep = compute_df['pH x pMi'].sum() + edge_df['pH x pMi'].sum()

        # print('Level:', peak_value)
        # print('Sample space probability:', compute_df['pH x pMi'].sum())
        # print('Outer space probability (upper):')
        # print(edge_df)

        # # Test boundary interval as fraction of probability
        # if aep != 0.0:
        #     p0 = compute_df.loc[0,'pH x pMi'] / aep
        #     pM = compute_df.loc[m-1,'pH x pMi'] / aep
        #     if p0 > 0.2 or pM > 0.2:
        #         aep = np.nan

        # Conservative condition: test all sims in first two intervals < peak_value; and all sims in last interval exceed peak_value
        # if np.all(compute_df.iloc[np.r_[0:2, -1:0]]['num'].to_numpy() == np.array([0, 0, n])):
        #     return aep
        # else:
        #     return np.nan

        # Flow based limits (mean in first and last interval). Abandoned to stretch the TPT curve used aep limits
        # if self.tpt_limits is not None:
        #     if peak_value < self.tpt_limits['lower'] or peak_value > self.tpt_limits['upper']:
        #        aep = np.nan

        # limit extent of TPT curve by AEP bounds
        # if self.tpt_limits is not None:
        #     inverse_aep = 1 / aep
        #     if inverse_aep < self.tpt_limits['lower'] or inverse_aep > self.tpt_limits['upper']:
                # print(f'AEP of 1 in {inverse_aep} is outside the limits')
        #        aep = np.nan
        return aep

