"""
lib/ReservoirRouting.py

Fast routing-only rerunner for the TFD Monte Carlo design flood study.

ReservoirRoutingSimulator re-routes a fixed inflow hydrograph from an existing
URBS run under an alternative dam rating curve (ELS / SQ), without re-running
hydrology. It uses the Storage-Indication form of the Modified Puls method
vectorised across all simulations on each timestep -- no scipy.optimize, no Numba.

FastTPT vectorises the Total Probability Theorem (counts of exceedances per
main-division) via one-time per-group sorts plus np.searchsorted. Math mirrors
lib/MCScheme.py::TotalProbTheorem.

Invoked from Main.py when an Excel sim-list row has Method == 'reservoir routing'.
"""
import os
import re
import sys
import time
import json
import platform
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.special import ndtr, ndtri

from lib.Lake import LakeConditions
from lib.EnbAnalysis import analyse_ensemble
from lib.Simulator import Logger


# ---------------------------------------------------------------------------
# Rating curve I/O (mirrors lib/Routing.py:set_dam_properties_from_urbs)
# ---------------------------------------------------------------------------

def _read_els(els_path):
    """Read URBS ELS file (CSV with EL, V columns). Returns sorted (EL, V) arrays."""
    df = pd.read_csv(els_path)
    el = df['EL'].astype(float).to_numpy()
    v = df['V'].astype(float).to_numpy()
    order = np.argsort(el)
    return el[order], v[order]


def _read_sq(sq_path, fsv):
    """Read URBS SQ file: 5 header lines + 'storage_above_FSV(ML) flow(m^3/s)' pairs.
    Returns (S_abs, O) where S_abs is absolute storage in ML and O in m^3/s."""
    storages, flows = [], []
    with open(sq_path) as f:
        for _ in range(6):
            line = f.readline()
        while line:
            line = line.strip()
            if line:
                s, q = re.split(r'\s+', line)
                storages.append(float(s))
                flows.append(float(q))
            line = f.readline()
    s_abs = np.array(storages, dtype=float) + fsv
    o = np.array(flows, dtype=float)
    return s_abs, o


def _read_indexed(path, description, sims_list_column):
    """Read a .csv or .parquet table whose first column is its index.

    Both the inflow hydrographs and the MCDF are written in either format by a
    Monte Carlo run, so both are read through here. Parquet needs pyarrow.

    ``description`` and ``sims_list_column`` only shape the error messages.
    """
    extension = os.path.splitext(path)[1].lower()
    if extension not in ('.csv', '.parquet'):
        raise ValueError(
            f'{description} files must be .csv or .parquet, got "{extension}": {path}'
        )

    if not os.path.isfile(path):
        # The same file is often kept in both formats, so point at the
        # counterpart rather than just reporting the missing name.
        counterpart = os.path.splitext(path)[0] + (
            '.parquet' if extension == '.csv' else '.csv')
        hint = f'\n  A {os.path.splitext(counterpart)[1]} version does exist: {counterpart}' \
               f'\n  Update the "{sims_list_column}" column in the sims list to match.' \
               if os.path.isfile(counterpart) else ''
        raise FileNotFoundError(f'{description} file not found: {path}{hint}')

    if extension == '.csv':
        return pd.read_csv(path, index_col=0)

    df = pd.read_parquet(path)
    # Parquet written with index=False carries the index as an ordinary first
    # column, so it must be promoted explicitly - there is no positional
    # index_col here.
    df.set_index(df.columns[0], inplace=True)
    return df


def _read_inflows(inflow_path):
    """Read an inflow hydrograph file, indexed by time (hours), one column per sim.

    Accepts .csv or .parquet. Parquet needs pyarrow installed.
    """
    return _read_indexed(inflow_path, 'Inflow', 'Inflow')


def _read_mcdf(mcdf_path):
    """Read a Monte Carlo dataframe, one row per realisation.

    Accepts .csv or .parquet. Parquet needs pyarrow installed.
    """
    return _read_indexed(mcdf_path, 'MCDF', 'Input MCDF')


def _build_routing_grid(el_els, v_els, fsv, sq_storage, sq_flow, n_grid=20000):
    """Build the storage/outflow/level lookup grids.

    psi(S) = 2000*S/dt_s + O(S) is built later (depends on dt). O(S) is 0 below
    FSV and interpolated from the SQ curve above. Level is from the ELS curve.
    """
    s_max = max(float(v_els.max()), float(sq_storage.max())) * 1.05
    s_grid = np.linspace(0.0, s_max, n_grid)
    o_grid = np.where(
        s_grid <= fsv,
        0.0,
        np.interp(s_grid, sq_storage, sq_flow, left=0.0, right=sq_flow[-1])
    )
    l_grid = np.interp(s_grid, v_els, el_els,
                       left=el_els[0], right=el_els[-1])
    return {'S': s_grid, 'O': o_grid, 'L': l_grid}


# ---------------------------------------------------------------------------
# Vectorised Storage-Indication Modified Puls
# ---------------------------------------------------------------------------

def _route_storage_indication(inflows_arr, adv, grids, dt_hours):
    """Vectorised Modified Puls in Storage-Indication form.

    Continuity (one timestep):  (S2-S1)*1000/dt_s = 0.5*(I1+I2) - 0.5*(O1+O2)
    Rearranged:                  psi(S2) = 2000*S1/dt_s - O1 + (I1+I2)
                              =  psi(S1) - 2*O1 + (I1+I2)

    where psi(S) = 2000*S/dt_s + O(S) is strictly monotone in S, so a single
    np.interp per timestep solves for S2 across all sims simultaneously.

    Parameters
    ----------
    inflows_arr : (n_t, n_s) m^3/s
    adv         : (n_s,) initial storage in ML
    grids       : dict from _build_routing_grid (keys S, O, L)
    dt_hours    : timestep in hours

    Returns S, O, L each shape (n_t, n_s) -- ML, m^3/s, m AHD.
    """
    n_t, n_s = inflows_arr.shape
    dt_s = 3600.0 * dt_hours
    S_grid = grids['S']
    O_grid = grids['O']
    L_grid = grids['L']
    psi_grid = 2000.0 * S_grid / dt_s + O_grid

    S = np.empty((n_t, n_s))
    O = np.empty((n_t, n_s))
    L = np.empty((n_t, n_s))

    S[0] = adv
    O[0] = np.interp(adv, S_grid, O_grid)
    L[0] = np.interp(adv, S_grid, L_grid)

    psi_prev = 2000.0 * S[0] / dt_s + O[0]
    for t in range(1, n_t):
        rhs = psi_prev - 2.0 * O[t - 1] + inflows_arr[t - 1] + inflows_arr[t]
        S[t] = np.interp(rhs, psi_grid, S_grid)
        O[t] = np.interp(S[t], S_grid, O_grid)
        L[t] = np.interp(S[t], S_grid, L_grid)
        psi_prev = 2000.0 * S[t] / dt_s + O[t]
    return S, O, L


# ---------------------------------------------------------------------------
# Vectorised Total Probability Theorem
# ---------------------------------------------------------------------------

class FastTPT:
    """Vectorised TPT -- AEP for all K peaks in O((m+1)*K log n) ops.

    Matches the math of lib/MCScheme.py::TotalProbTheorem (lines 459-547):
      pH_i  = (count in group i strictly greater than peak) / n
      pH_-1 = sqrt(upper_factor_assumption * pH_0**2)
      pH_m  = sqrt(pH_{m-1} * 1)
      aep   = sum_i (pH_i * pMi_i) + pH_-1 * pMi_-1 + pH_m * pMi_m
    """

    def __init__(self, m, n, main_divisions):
        self.m = int(m)
        self.n = int(n)
        main_divisions = np.asarray(main_divisions, dtype=float)
        z_min = main_divisions[:-1]
        z_max = main_divisions[1:]
        p_min = ndtr(z_min)
        p_max = ndtr(z_max)
        self.pMi = p_max - p_min                # (m,)
        self.pMi_low = float(p_min[0])          # below lowest z
        self.pMi_high = float(1.0 - p_max[-1])  # above highest z
        self.upper_factor_assumption = self.pMi_low

    def assign_aep(self, peaks, group_ids):
        """Vectorised assign_aep equivalent to TotalProbTheorem.assign_aep applied
        to every row in mcdf.

        peaks     : (K,) array
        group_ids : (K,) array of m-indices (0..m-1)
        Returns   : (K,) array of exceedance probabilities.
        """
        peaks = np.asarray(peaks, dtype=float)
        group_ids = np.asarray(group_ids, dtype=int)

        grouped = [np.sort(peaks[group_ids == i]) for i in range(self.m)]
        sizes = np.array([len(g) for g in grouped])

        K = peaks.size
        num = np.empty((self.m, K))
        for i in range(self.m):
            num[i] = sizes[i] - np.searchsorted(grouped[i], peaks, side='right')
        pH = num / self.n                        # (m, K)
        pH_low = np.sqrt(self.upper_factor_assumption * pH[0] ** 2)
        pH_high = np.sqrt(pH[-1])
        aep = pH.T @ self.pMi + pH_low * self.pMi_low + pH_high * self.pMi_high
        return aep


def standard_aeps(lower_aep, upper_aep, aep_of_pmp=None):
    """Mirror of SampleScheme.get_standard_aeps (lib/MCScheme.py:350-373)."""
    if lower_aep < 2:
        lower_aep = 2
    aep = lower_aep
    out = []
    while aep <= upper_aep:
        if aep >= lower_aep:
            out.append(aep)
        if np.log10(aep / 2) % 1 == 0.0:
            aep = aep * 5 // 2
        else:
            aep = aep * 2
    if aep_of_pmp is not None and aep_of_pmp < upper_aep and aep_of_pmp not in out:
        out.append(int(aep_of_pmp))
        out.sort()
    return out


_PLOT_LABEL = {
    'inflow':  'Peak inflow (m³/s)',
    'outflow': 'Peak outflow (m³/s)',
    'level':   'Peak lake level (m AHD)',
}


def _plot_tpt(result_type, rain_aeps, peaks, std_q, out_path,
              lower_aep, upper_aep, aep_of_pmp=None):
    """Scatter of Monte Carlo realisations (peak vs rainfall AEP) overlaid with
    the TPT quantile curve. Mirrors lib/MCScheme.py:plot_tpt_results_2 style."""
    fig, ax = plt.subplots(figsize=(7, 5))

    # Scatter: x = z-score of the *rainfall* AEP (one point per realisation)
    valid_pts = (rain_aeps > 1) & np.isfinite(peaks)
    if valid_pts.any():
        z_pts = ndtri(1.0 - 1.0 / rain_aeps[valid_pts])
        ax.plot(z_pts, peaks[valid_pts], 'o',
                markersize=2.5, markeredgewidth=0, alpha=0.25,
                color='C0', label='Monte Carlo realisations')

    # TPT line: x = z-score of the *result* AEP at each standard quantile
    line_df = std_q.dropna(subset=[result_type])
    if not line_df.empty:
        z_line = ndtri(1.0 - line_df['probability'].to_numpy())
        ax.plot(z_line, line_df[result_type].to_numpy(),
                '-', color='k', linewidth=1.5, label='TPT')

    # AEP-labelled ticks
    ticks_aeps = np.array(standard_aeps(lower_aep, upper_aep, aep_of_pmp), dtype=float)
    tick_z = ndtri(1.0 - 1.0 / ticks_aeps)
    ax.set_xticks(tick_z)
    ax.set_xticklabels([f'{int(a):,}' if a >= 1 else f'{a:g}' for a in ticks_aeps],
                       rotation=90)
    ax.set_xlabel('AEP (1 in X)')
    ax.set_ylabel(_PLOT_LABEL.get(result_type, result_type))
    if result_type != 'level':
        ax.set_yscale('log')
        if valid_pts.any():
            ymin = max(np.nanmin(peaks[valid_pts & (peaks > 0)]) if (peaks[valid_pts] > 0).any() else 0.1, 0.1)
            ax.set_ylim(bottom=ymin * 0.5)
    ax.grid(True, which='both', alpha=0.3)
    ax.legend(loc='best', framealpha=0.9)
    plt.tight_layout()
    fig.savefig(out_path, dpi=120)
    plt.close(fig)


def compute_std_quantiles(peaks, aeps, lower_aep, upper_aep, aep_of_pmp=None):
    """Log-normal-space interpolation to evaluate the peak at each standard AEP.
    Mirrors lib/MCScheme.py:326-334."""
    std_aeps = np.array(standard_aeps(lower_aep, upper_aep, aep_of_pmp), dtype=float)
    std_p = 1.0 / std_aeps
    valid = (aeps > 0) & (peaks > 0)
    if not valid.any():
        return pd.DataFrame({'aep (1 in x)': std_aeps,
                             'probability': std_p,
                             'value': np.full_like(std_aeps, np.nan)})
    z_data = ndtri(1.0 - aeps[valid])
    log_y = np.log10(peaks[valid])
    order = np.argsort(z_data)
    z_data, log_y = z_data[order], log_y[order]
    z_std = ndtri(1.0 - std_p)
    interp_log_y = np.interp(z_std, z_data, log_y, left=np.nan, right=np.nan)
    return pd.DataFrame({'aep (1 in x)': std_aeps,
                         'probability': std_p,
                         'value': 10.0 ** interp_log_y})


# ---------------------------------------------------------------------------
# Main simulator (called by Main.py via dispatch)
# ---------------------------------------------------------------------------

class ReservoirRoutingSimulator:
    """Re-route fixed inflow hydrographs from a prior URBS run under an alternative
    rating curve, then optionally re-analyse the result.

    Takes the output of either scheme. Which one is in front of it is worked out
    from the columns of the input database (_detect_scheme), and it decides both
    how the antecedent dam volume is resolved and how the routed peaks are
    analysed:

      monte carlo - one row per realisation, with the m / n sample position.
                    Re-analysed with the Total Probability Theorem (FastTPT).
      ensemble    - one row per AEP x duration x temporal pattern. Re-analysed
                    with lib/EnbAnalysis.py: the median pattern of each duration
                    and the critical duration for each AEP.

    Required Excel sim-list columns:
      Include, Method (= 'reservoir routing'), Output file,
      Run models, Analyse results, Input MCDF, Inflow,
      ELS file, SQ file, FSL, Hydrographs folder, Results folder

    Optional columns:
      Store hydrographs (default 'yes'), Output suffix, Config file,
      ADV, ADV source, Lake config, Log file, Comment.

    'Input MCDF' names the input database for both schemes ('Input database' is
    accepted as an alias). 'Config file' should be the Monte Carlo config the
    input mcdf was generated with: the TPT parameters come from its
    scheme_config, and the mcdf sample is checked against them before the results
    are analysed (_validate_sample). It is not used by the ensemble scheme.

    -- Monte Carlo input: 'ADV source' selects where the initial storage comes from:
      mcdf              - the ADV column of the input mcdf (default, and what
                          happens if the column is absent)
      lake_z            - the ADV is recomputed from the lake_z column of the
                          input mcdf using the distribution in the 'Lake config'
                          file, so the same sample can be re-routed under a
                          different antecedent storage distribution
      lake_z correlated - as above, but the correlation layers in the lake config
                          are (re)applied to lake_z against the mcdf rain_z first.
                          Only for mcdf files sampled without a correlation - the
                          lake_z written by a correlated Monte Carlo run already
                          has the correlation in it.

    -- Ensemble input: the ADV is a single value for the whole run, so it comes
    from the 'ADV' column of the sims list, read the same way the ensemble method
    itself reads it (lib/Lake.py):
      <a number>        - that many ML
      fsv               - the full supply volume of the *new* rating curve. This
                          is the one to use when re-routing under a different
                          dam: taking the ADV from the input database instead
                          would silently start every event at the old dam's FSV.
      mav               - the volume the 'Lake config' exceedance curve gives at
                          z = 0, i.e. the median of the distribution a Monte Carlo
                          run samples. Capped against the new curve's FSV like fsv.
      database          - the ADV column of the input database, i.e. whatever the
                          source run started from. For reproducing that run.
    'ADV source' is a Monte Carlo concept and is ignored for ensemble input.
    """

    ADV_SOURCES = ('mcdf', 'lake_z', 'lake_z correlated')

    def __init__(self, sim_row, filepaths, test_runs=0):
        self.start = time.time()
        self.sim_row = sim_row
        self.filepaths = filepaths
        self.test_runs = test_runs
        self.basename = str(sim_row['Output file'])

        self.hydrographs_folder = str(sim_row['Hydrographs folder'])
        self.results_folder = str(sim_row['Results folder'])
        os.makedirs(self.hydrographs_folder, exist_ok=True)
        os.makedirs(self.results_folder, exist_ok=True)

        self._start_log()
        print(f'\n=== ReservoirRoutingSimulator: {self.basename} ===')
        self._log_inputs()

        self.run_models = str(sim_row['Run models']).strip().lower() == 'yes'
        self.do_analysis = str(sim_row['Analyse results']).strip().lower() == 'yes'
        self.adv_source = self._get_adv_source()
        self.scheme = None  # set by _detect_scheme once the database is read

        if self.run_models:
            self._load_curves()
            self._load_inflows()
            self._load_mcdf()
            self._route()
            self._write_hydrographs()
            self._write_mcdf()

        if self.do_analysis:
            self._ensure_mcdf_loaded()
            self._analyse()

        elapsed = time.time() - self.start
        print(f'\nReservoirRoutingSimulator finished in {elapsed:.2f} s '
              f'({elapsed / 60:.2f} min)')

    # ----------------------------------------------------------------- log

    def _log_path(self):
        """Where this simulation's log goes.

        The sims-list 'Log file' if there is one, else a file beside the results
        so a reservoir routing run always leaves a record - QA needs to be able to
        tie a set of results back to the inputs that produced it, and this method
        used to write nothing at all.
        """
        if 'Log file' in self.sim_row.index and pd.notna(self.sim_row['Log file']):
            path = str(self.sim_row['Log file']).strip()
            if path:
                return path
        return os.path.join(self.results_folder, f'{self.basename}{self._suffix()}_log.txt')

    def _start_log(self):
        """Tee the console output to the log file, as the other two methods do."""
        self.logger = None
        log_path = self._log_path()
        folder = os.path.dirname(log_path)
        if folder:
            os.makedirs(folder, exist_ok=True)
        try:
            self.logger = Logger(log_path)
            sys.stdout = self.logger
        except OSError as error:
            # A log that cannot be opened is not a reason to lose the simulation.
            print(f'WARNING: could not open the log file "{log_path}": {error}')

    @staticmethod
    def _file_details(path):
        """Size and modification time of an input file, for the log.

        Which file was read matters as much as its name when results are being
        checked months later - two runs can name the same relative path and get
        different data.
        """
        try:
            stat = os.stat(path)
            when = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(stat.st_mtime))
            size = stat.st_size
            for unit, scale in (('MB', 1e6), ('kB', 1e3)):
                if size >= scale:
                    return f'{size / scale:,.1f} {unit}, modified {when}'
            return f'{size:,d} bytes, modified {when}'
        except OSError:
            return 'NOT FOUND'

    def _log_inputs(self):
        """Write everything that determines this simulation's results.

        Printed rather than returned so it lands in both the log file and the
        console, and written before anything is loaded so that a run which fails
        on a missing input still records what it was looking for.
        """
        row = self.sim_row
        print(f'\n{"-" * 78}\nSIMULATION INPUTS\n{"-" * 78}')
        print(f'  run at            {time.strftime("%Y-%m-%d %H:%M:%S")} '
              f'on {platform.node()} ({platform.system()} {platform.release()})')
        print(f'  python            {platform.python_version()}')
        print(f'  working folder    {os.getcwd()}')
        print(f'  log file          {self._log_path()}')

        print('\n  Sims list row:')
        # Everything the sims list carries, so the log is a complete statement of
        # what was asked for - including any columns added since this was written.
        for key in row.index:
            value = row[key]
            if pd.isna(value):
                value = ''
            print(f'    {str(key):<22} {value}')

        print('\n  Input files:')
        for label, column in (('input database', 'Input MCDF'),
                              ('input database', 'Input database'),
                              ('inflows', 'Inflow'),
                              ('ELS (storage)', 'ELS file'),
                              ('SQ (rating)', 'SQ file'),
                              ('lake config', 'Lake config'),
                              ('method config', 'Config file')):
            if column not in row.index or pd.isna(row[column]) or not str(row[column]).strip():
                continue
            path = str(row[column])
            print(f'    {label:<16} {os.path.abspath(path)}')
            print(f'    {"":<16}   {self._file_details(path)}')

        print('\n  Outputs:')
        print(f'    {"results folder":<16} {os.path.abspath(self.results_folder)}')
        print(f'    {"hydrographs":<16} {os.path.abspath(self.hydrographs_folder)}')
        print(f'{"-" * 78}')

    # ------------------------------------------------------------------ IO

    def _load_curves(self):
        els_path = str(self.sim_row['ELS file'])
        sq_path = str(self.sim_row['SQ file'])
        self.fsl = float(self.sim_row['FSL'])
        print(f'Loading ELS: {els_path}')
        self.el_els, self.v_els = _read_els(els_path)
        self.fsv = float(np.interp(self.fsl, self.el_els, self.v_els))
        print(f'  FSL = {self.fsl} m AHD  ->  FSV = {self.fsv:.1f} ML')
        print(f'Loading SQ:  {sq_path}')
        self.sq_storage, self.sq_flow = _read_sq(sq_path, self.fsv)
        self.grids = _build_routing_grid(self.el_els, self.v_els, self.fsv,
                                         self.sq_storage, self.sq_flow)

    def _load_inflows(self):
        inflow_path = str(self.sim_row['Inflow'])
        print(f'Loading inflows: {inflow_path}')
        t0 = time.time()
        df = _read_inflows(inflow_path)
        df.interpolate(method='slinear', inplace=True)
        self.time_index = df.index.to_numpy(dtype=float)
        self.sim_columns = df.columns.tolist()
        self.inflows_arr = df.to_numpy(dtype=float)
        self.dt_hours = float(self.time_index[1] - self.time_index[0])
        print(f'  {self.inflows_arr.shape[0]} timesteps x '
              f'{self.inflows_arr.shape[1]} sims, dt = {self.dt_hours} h '
              f'(read in {time.time() - t0:.2f} s)')

    def _get_adv_source(self):
        """Where the antecedent dam volume comes from - see the class docstring."""
        source = 'mcdf'
        for column in ('ADV source', 'ADV sample source'):
            if column in self.sim_row.index and pd.notna(self.sim_row[column]):
                source = ' '.join(str(self.sim_row[column]).split()).lower()
                break
        if source in ('', 'adv'):
            source = 'mcdf'
        if source not in self.ADV_SOURCES:
            raise ValueError(
                f'"ADV source" of "{source}" is not recognised. '
                f'Use one of: {" | ".join(self.ADV_SOURCES)}'
            )
        # Not reported here: the scheme of the input database is not known yet, and
        # "ADV source" applies to Monte Carlo input only. _load_mcdf says what was
        # actually used once it does know.
        return source

    def _database_path(self):
        """The input database, under either of the two accepted column names."""
        for column in ('Input database', 'Input MCDF'):
            if column in self.sim_row.index and pd.notna(self.sim_row[column]):
                return str(self.sim_row[column]), column
        raise ValueError('The sims list needs an "Input MCDF" (or "Input database") '
                         'filepath for the reservoir routing method.')

    def _detect_scheme(self):
        """Work out whether the input database came from a Monte Carlo or an
        ensemble run, from the columns each scheme writes.

        The Monte Carlo scheme writes the sample position (m, n) of every
        realisation; the ensemble scheme writes the duration and temporal pattern
        of every event. Neither writes the other's, so the two are distinguishable
        without another sims-list column to get out of step with the data.
        """
        columns = set(self.mcdf.columns)
        is_mc = {'m', 'n'} <= columns
        is_enb = {'duration', 'tp'} <= columns and not is_mc
        if is_mc:
            scheme = 'monte carlo'
        elif is_enb:
            scheme = 'ensemble'
        else:
            raise ValueError(
                'The input database is neither a Monte Carlo nor an ensemble result:\n'
                '  a Monte Carlo database has "m" and "n" columns (the sample position),\n'
                '  an ensemble database has "duration" and "tp" columns.\n'
                f'  This one has: {", ".join(sorted(columns))}'
            )
        print(f'Input database is from the {scheme} scheme')
        return scheme

    def _load_mcdf(self):
        mcdf_path, column = self._database_path()
        print(f'Loading input database ({column}): {mcdf_path}')
        self.mcdf = _read_mcdf(mcdf_path)
        if len(self.mcdf) != self.inflows_arr.shape[1]:
            raise ValueError(
                f'Input database row count ({len(self.mcdf)}) does not match number of '
                f'inflow columns ({self.inflows_arr.shape[1]}).'
            )
        self.scheme = self._detect_scheme()
        if self.scheme == 'ensemble':
            self.adv = self._ensemble_adv()
        elif self.adv_source == 'mcdf':
            print('ADV source: the ADV column of the input MCDF')
            self.adv = self.mcdf['ADV'].to_numpy(dtype=float)
        else:
            print(f'ADV source: {self.adv_source} - the ADV will be resampled '
                  f'from the lake config file')
            self.adv = self._sample_adv_from_lake_z()

    def _ensemble_adv(self):
        """The antecedent dam volume for an ensemble re-route.

        An ensemble run holds the lake at one starting volume for every event, so
        the ADV comes from the sims list rather than per row. The sims-list values
        are read by LakeConditions exactly as the ensemble method itself reads
        them, so a number and the "fsv" and "mav" keywords behave identically here
        - except that both keywords resolve against the rating curve being routed,
        not the one the input database was generated with. That is the point:
        taking the input database's own ADV under a new curve would start every
        event at the old dam's full supply volume.
        """
        if 'ADV source' in self.sim_row.index and pd.notna(self.sim_row['ADV source']):
            print('NOTE: "ADV source" applies to Monte Carlo input only and is ignored '
                  'for an ensemble - the ADV comes from the "ADV" column of the sims list.')
        if 'ADV' not in self.sim_row.index or pd.isna(self.sim_row['ADV']):
            raise ValueError('The sims list needs an "ADV" for an ensemble re-route: '
                             'a volume in ML, "fsv", "mav", or "database".')

        adv_value = self.sim_row['ADV']
        n_sims = self.inflows_arr.shape[1]

        if isinstance(adv_value, str) and adv_value.strip().lower() == 'database':
            if 'ADV' not in self.mcdf.columns:
                raise ValueError('The input database has no "ADV" column to take the '
                                 'antecedent dam volume from.')
            adv = self.mcdf['ADV'].to_numpy(dtype=float)
            print(f'ADV taken from the input database: {np.unique(adv)} ML')
            return adv

        if isinstance(adv_value, str) and adv_value.strip().lower() == 'varying':
            # LakeConditions only blocks this for Method == 'ensemble', and this row
            # says 'reservoir routing', so it would otherwise fall through to the
            # lake config and sample 130 different starting volumes.
            raise ValueError('An "ADV" of "varying" cannot be used with ensemble input - '
                             'an ensemble holds one starting volume for the whole run. '
                             'Use a volume in ML, "fsv", "mav", or "database".')

        lake = LakeConditions(self.sim_row)
        lake.set_full_supply_volume(self.fsv)
        if lake.antecedent_volume is None:
            raise ValueError(f'The "ADV" of "{adv_value}" did not resolve to a volume.')
        volume = float(lake.antecedent_volume)
        print(f'ADV for all {n_sims} events: {volume:.1f} ML '
              f'(FSV of the routed curve = {self.fsv:.1f} ML)')
        self.mcdf['ADV_input'] = self.mcdf['ADV'] if 'ADV' in self.mcdf.columns else np.nan
        self.mcdf['ADV'] = volume
        return np.full(n_sims, volume)

    def _sample_adv_from_lake_z(self):
        """Recompute the ADV for every realisation from the sampled lake z values
        in the MCDF, using the distribution in the lake config file.

        This lets the one set of inflow hydrographs be re-routed under a different
        antecedent storage distribution without re-running the hydrology.
        """
        lake_config = self.sim_row.get('Lake config')
        if lake_config is None or pd.isna(lake_config):
            raise ValueError(
                'A "Lake config" filepath is needed in the sims list when the '
                f'ADV source is "{self.adv_source}" - check the sims list!'
            )
        if 'lake_z' not in self.mcdf.columns:
            raise ValueError(
                f'The input MCDF has no "lake_z" column, so the ADV cannot be '
                f'resampled. Use an ADV source of "mcdf" or an MCDF from a Monte '
                f'Carlo run with a varying ADV.'
            )

        lake = LakeConditions.from_lake_config(str(lake_config))
        # The FSV comes from the *new* rating curve, so any capping is against the
        # storage being routed rather than the one the MCDF was generated with.
        lake.set_full_supply_volume(self.fsv)

        lake_z = self.mcdf['lake_z'].to_numpy(dtype=float)
        if self.adv_source == 'lake_z correlated':
            if not lake.is_correlated:
                raise ValueError(
                    'The ADV source is "lake_z correlated" but the lake config '
                    'file has no "correlation_layer_info" key.'
                )
            if 'rain_z' not in self.mcdf.columns:
                raise ValueError('The input MCDF has no "rain_z" column to correlate against.')
            print('\nApplying correlation for lake initial conditions...')
            rain_z = self.mcdf['rain_z'].to_numpy(dtype=float)
            lake_z = lake.correlation.apply_correlations_to_array(rain_z, lake_z)
            self.mcdf['lake_z_correlated'] = lake_z
        elif lake.is_correlated:
            print('\nNOTE: the lake config file has correlation layers, but the lake z '
                  'values in the MCDF are used as sampled.\n'
                  '      Use an ADV source of "lake_z correlated" to apply the '
                  'correlation (only valid if the MCDF was sampled without one).')

        adv = lake.get_lake_volumes(lake_z)
        print(f'\nResampled ADV for {len(adv)} realisations: '
              f'min {adv.min():.1f}, mean {adv.mean():.1f}, max {adv.max():.1f} ML '
              f'(FSV = {self.fsv:.1f} ML, volume cap = {lake.volume_cap})')
        # Record what was actually routed - the ADV column from the source run no
        # longer applies. The original is kept alongside it for checking.
        if 'ADV' in self.mcdf.columns:
            self.mcdf['ADV_input'] = self.mcdf['ADV']
        self.mcdf['ADV'] = adv
        return adv

    # ----------------------------------------------------------------- core

    def _route(self):
        print('Routing (Storage-Indication, vectorised)...')
        t0 = time.time()
        self.S, self.O, self.L = _route_storage_indication(
            self.inflows_arr, self.adv, self.grids, self.dt_hours
        )
        print(f'  Routed in {time.time() - t0:.2f} s')

    # ----------------------------------------------------------------- IO

    def _suffix(self):
        if 'Output suffix' in self.sim_row.index and pd.notna(self.sim_row['Output suffix']):
            return '_' + str(self.sim_row['Output suffix'])
        return ''

    def _write_hydrographs(self):
        store = 'yes'
        if 'Store hydrographs' in self.sim_row.index and pd.notna(self.sim_row['Store hydrographs']):
            store = str(self.sim_row['Store hydrographs']).strip().lower()
        if store != 'yes':
            print('Store hydrographs = no -> skipping hydrograph CSVs')
            return
        suffix = self._suffix()
        out_specs = {'outflows': self.O, 'levels': self.L, 'volumes': self.S}
        for name, arr in out_specs.items():
            path = os.path.join(self.hydrographs_folder,
                                f'{self.basename}_{name}{suffix}.csv')
            print(f'Writing {path}')
            t0 = time.time()
            df = pd.DataFrame(arr, index=self.time_index, columns=self.sim_columns)
            df.index.name = 'time'
            df.to_csv(path)
            print(f'  ({time.time() - t0:.2f} s)')

    def _output_base(self):
        """The output basename the results are written under.

        The two schemes name their databases differently, because each has to be
        found again by the code that analyses it: the Monte Carlo mcdf keeps the
        '__mcdf' tag the rest of Bryan expects, while the ensemble database is
        '<base>.csv' beside the plots/ and csv/ folders that lib/EnbAnalysis.py
        writes into - the same layout EnsembleSimulator produces.
        """
        if self.scheme == 'ensemble':
            return os.path.join(self.results_folder, f'{self.basename}{self._suffix()}')
        return os.path.join(self.results_folder, f'{self.basename}__mcdf')

    def _write_mcdf(self):
        # nanmax is critical: inflow series interpolated via 'slinear' can have
        # trailing NaNs (sims with shorter URBS runs, and every duration but the
        # longest in an ensemble hydrograph file), and the routing arrays
        # propagate those NaNs to the tail of the level/outflow series.
        self.mcdf['inflow'] = np.nanmax(self.inflows_arr, axis=0)
        self.mcdf['outflow'] = np.nanmax(self.O, axis=0)
        self.mcdf['level'] = np.nanmax(self.L, axis=0)
        out_path = f'{self._output_base()}.csv'
        print(f'Writing new {self.scheme} database: {out_path}')
        self.mcdf.to_csv(out_path)

    # --------------------------------------------------------------- analyse

    def _ensure_mcdf_loaded(self):
        if hasattr(self, 'mcdf'):
            return
        # The scheme is normally set when the database is read, but an
        # analysis-only row (Run models = no) has not read one yet, so fall back
        # to the input and detect from that.
        mcdf_path, column = self._database_path()
        for scheme, tag in (('monte carlo', '__mcdf'), ('ensemble', self._suffix())):
            out_path = os.path.join(self.results_folder, f'{self.basename}{tag}.csv')
            if os.path.isfile(out_path):
                print(f'Loading routed {scheme} database for analysis: {out_path}')
                self.mcdf = pd.read_csv(out_path, index_col=0)
                self.scheme = self._detect_scheme()
                return
        print(f'Loading database for analysis (fallback to input {column}): {mcdf_path}')
        self.mcdf = _read_mcdf(mcdf_path)
        self.scheme = self._detect_scheme()

    def _tpt_config(self):
        """Pull TPT params from the row's Config file JSON, with sensible defaults
        derived from the MCDF if the JSON is missing or partial.

        The Monte Carlo config file holds these keys inside "scheme_config", so
        look there first and fall back to the top level for a config written just
        for the routing. Where a key comes from matters: a wrong lower/upper AEP
        silently shifts the main division boundaries and so the TPT weights, which
        is why the source of each value is reported and the sample is then checked
        against them in _validate_sample.
        """
        cfg = {}
        if 'Config file' in self.sim_row.index and pd.notna(self.sim_row['Config file']):
            cfg_path = str(self.sim_row['Config file'])
            if os.path.isfile(cfg_path):
                with open(cfg_path) as f:
                    cfg = json.load(f)
                print(f'\nReading the TPT parameters from: {cfg_path}')
            else:
                print(f'\nWARNING: the "Config file" in the sims list was not found: {cfg_path}')
        else:
            print('\nWARNING: no "Config file" given in the sims list for this simulation.')

        scheme = cfg.get('scheme_config', {})
        self.config_sources = {}

        def from_config(key, default):
            if key in scheme:
                self.config_sources[key] = 'scheme_config'
                return scheme[key]
            if key in cfg:
                self.config_sources[key] = 'config file'
                return cfg[key]
            self.config_sources[key] = 'DEFAULT (not in the config file)'
            return default

        m_count = int(from_config('number_of_main_divisions', self.mcdf['m'].nunique()))
        n_count = int(from_config('number_of_sub_divisions', len(self.mcdf) // max(m_count, 1)))
        lower_aep = float(from_config('lower_aep', 2))
        upper_aep = float(from_config('upper_aep', 1e7))
        aep_of_pmp = from_config('aep_of_pmp', None)
        if aep_of_pmp is not None:
            aep_of_pmp = float(aep_of_pmp)

        print('\nTPT parameters:')
        for key, value in (('number_of_main_divisions', m_count),
                           ('number_of_sub_divisions', n_count),
                           ('lower_aep', lower_aep),
                           ('upper_aep', upper_aep),
                           ('aep_of_pmp', aep_of_pmp)):
            print(f'  {key:<26} {str(value):<12} from the {self.config_sources[key]}')
        scheme_keys = ('number_of_main_divisions', 'number_of_sub_divisions',
                       'lower_aep', 'upper_aep')
        if any(self.config_sources[key].startswith('DEFAULT') for key in scheme_keys):
            print('  NOTE: defaulted values are not necessarily those the MCDF was sampled with.\n'
                  '        The AEP bounds set the main division boundaries used by the TPT.')

        return m_count, n_count, lower_aep, upper_aep, aep_of_pmp

    def _validate_sample(self, m_count, n_count, main_divisions):
        """Check that the MCDF sample matches the scheme it is about to be analysed
        with - the sampled rainfall must sit inside the AEP range, in the right main
        division, with the expected number of samples per division.

        The TPT weights each main division by its own probability slice, so an MCDF
        analysed against a scheme it was not sampled with gives quantiles that are
        wrong without looking wrong. Nothing downstream would catch it.
        """
        print('\nChecking the MCDF sample against the TPT parameters...')
        problems = []
        n_rows = len(self.mcdf)

        # --- structure of the sample space --------------------------------
        if 'm' not in self.mcdf.columns:
            raise Exception('The MCDF has no "m" column - the main division of each realisation.')

        expected_rows = m_count * n_count
        if n_rows != expected_rows:
            problems.append(
                f'The MCDF has {n_rows} realisations but the scheme describes '
                f'{m_count} x {n_count} = {expected_rows}.')

        group_ids = self.mcdf['m'].to_numpy(dtype=int)
        found_ids = np.unique(group_ids)
        expected_ids = np.arange(m_count)
        if not np.array_equal(found_ids, expected_ids):
            problems.append(
                f'The main divisions in the MCDF are {found_ids.min()} to {found_ids.max()} '
                f'({found_ids.size} of them), not 0 to {m_count - 1} as the scheme expects.')

        counts = np.array([(group_ids == i).sum() for i in expected_ids])
        odd = expected_ids[counts != n_count]
        if odd.size:
            problems.append(
                f'{odd.size} main division(s) do not hold {n_count} realisations, e.g. '
                + ', '.join(f'm={i} has {counts[i]}' for i in odd[:5])
                + ('...' if odd.size > 5 else '')
                + '. The TPT divides the exceedance counts by the number of sub divisions, '
                  'so the assigned AEPs would be biased.')

        if 'n' in self.mcdf.columns and not problems:
            # Only worth reporting once the counts are right - it locates duplicates
            # (the same sub division twice) that the counts alone would not show.
            pairs = set(zip(group_ids.tolist(),
                            self.mcdf['n'].to_numpy(dtype=int).tolist()))
            if len(pairs) != expected_rows:
                problems.append(
                    f'The MCDF holds {len(pairs)} distinct (m, n) sample positions rather '
                    f'than {expected_rows} - there are repeated realisations.')

        # --- the sampled rainfall itself ----------------------------------
        if 'rain_z' not in self.mcdf.columns:
            problems.append('The MCDF has no "rain_z" column, so the sampled rainfall '
                            'cannot be checked against the AEP range of the scheme.')
        else:
            rain_z = self.mcdf['rain_z'].to_numpy(dtype=float)
            z_low, z_up = main_divisions[0], main_divisions[-1]
            tol = 1e-6
            blank = ~np.isfinite(rain_z)
            if blank.any():
                # e.g. an MCDF from a run stopped early by the test_runs key - the
                # empty realisations would count as non-exceedances in the TPT.
                problems.append(
                    f'{blank.sum()} of {n_rows} realisations have no rainfall z. The MCDF looks '
                    f'incomplete (a run stopped part way through?), and the TPT needs the whole '
                    f'sample.')
            outside = (rain_z < z_low - tol) | (rain_z > z_up + tol)
            if outside.any():
                problems.append(
                    f'{outside.sum()} of {n_rows} sampled rainfall z values fall outside the '
                    f'{z_low:.4f} to {z_up:.4f} range of the scheme (1 in {1 / (1 - ndtr(z_low)):.0f} '
                    f'to 1 in {1 / (1 - ndtr(z_up)):.0f} AEP): sampled z runs from '
                    f'{rain_z.min():.4f} to {rain_z.max():.4f}. The MCDF was almost certainly '
                    f'sampled over a different AEP range to the one being analysed.')
            elif not problems:
                # Right range and right structure, so the divisions can be checked
                misplaced = 0
                for i in expected_ids:
                    in_group = group_ids == i
                    z_group = rain_z[in_group]
                    misplaced += ((z_group < main_divisions[i] - tol)
                                  | (z_group > main_divisions[i + 1] + tol)).sum()
                if misplaced:
                    problems.append(
                        f'{misplaced} of {n_rows} realisations hold a rainfall z from outside '
                        f'their own main division. The {m_count} divisions being analysed span '
                        f'z {z_low:.4f} to {z_up:.4f}, which is not how the MCDF was divided up: '
                        f'either the number of main divisions or the AEP bounds differ from the '
                        f'ones it was sampled with.')

            if 'rain_aep' in self.mcdf.columns:
                rain_aep = self.mcdf['rain_aep'].to_numpy(dtype=float)
                with np.errstate(divide='ignore', over='ignore'):
                    expected_aep = 1.0 / (1.0 - ndtr(rain_z))
                comparable = np.isfinite(rain_aep) & np.isfinite(expected_aep) & (expected_aep > 0)
                if comparable.any():
                    error = np.abs(rain_aep[comparable] / expected_aep[comparable] - 1.0).max()
                    if error > 1e-3:
                        problems.append(
                            f'The "rain_aep" column does not match "rain_z" (up to {error:.1%} '
                            f'out) - the two columns did not come from the one sample.')

        if 'tp_w' in self.mcdf.columns:
            print('  WARNING: the MCDF has a "tp_w" column (temporal pattern probability '
                  'weights).\n           The reservoir routing TPT does not apply these weights, '
                  'so the AEPs\n           will differ from the Monte Carlo analysis of the same '
                  'MCDF.')

        if problems:
            raise Exception(
                '\nThe MCDF sample does not match the TPT parameters being used:\n  - '
                + '\n  - '.join(problems)
                + '\n\nCheck that the "Config file" column in the sims list points at the '
                  'config\nfile the MCDF was generated with (the scheme keys are read from its '
                  '"scheme_config").'
            )

        print(f'  {n_rows} realisations in {m_count} main divisions of {n_count}, '
              f'rainfall z within each division: consistent.')

    def _analyse(self):
        if self.scheme == 'ensemble':
            return self._analyse_ensemble()
        return self._analyse_tpt()

    def _analyse_ensemble(self):
        """Median pattern per duration and critical duration per AEP - the same
        analysis EnsembleSimulator runs, over the re-routed peaks."""
        out_base = self._output_base()
        print(f'\nEnsemble analysis of the routed results: {out_base}')
        analyse_ensemble(self.mcdf, out_base)

    def _analyse_tpt(self):
        m_count, n_count, lower_aep, upper_aep, aep_of_pmp = self._tpt_config()

        z_low = ndtri(1.0 - 1.0 / lower_aep)
        z_up = ndtri(1.0 - 1.0 / upper_aep)
        main_divisions = np.linspace(z_low, z_up, m_count + 1)
        self._validate_sample(m_count, n_count, main_divisions)

        tpt = FastTPT(m_count, n_count, main_divisions)
        group_ids = self.mcdf['m'].to_numpy(dtype=int)
        suffix = self._suffix()

        rain_aeps = self.mcdf['rain_aep'].to_numpy(dtype=float)

        for result_type in ('inflow', 'outflow', 'level'):
            print(f'\nVectorised TPT: {result_type}')
            t0 = time.time()
            peaks = self.mcdf[result_type].to_numpy(dtype=float)
            aeps = tpt.assign_aep(peaks, group_ids)
            self.mcdf[f'{result_type}_aep'] = aeps
            std_q = compute_std_quantiles(peaks, aeps, lower_aep, upper_aep, aep_of_pmp)
            std_q = std_q.rename(columns={'value': result_type})
            out_q = os.path.join(self.results_folder,
                                 f'{self.basename}__{result_type}_quantiles{suffix}.csv')
            std_q.to_csv(out_q, index=False)
            print(f'  TPT + quantiles in {time.time() - t0:.2f} s -> {out_q}')

            png_path = os.path.join(self.results_folder,
                                    f'{self.basename}__{result_type}_tpt{suffix}.png')
            _plot_tpt(result_type, rain_aeps, peaks, std_q, png_path,
                      lower_aep, upper_aep, aep_of_pmp)
            print(f'  Plot -> {png_path}')

        out_path = f'{self._output_base()}.csv'
        print(f'\nUpdating MCDF with new *_aep columns: {out_path}')
        self.mcdf.to_csv(out_path)
