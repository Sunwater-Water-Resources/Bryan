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
import time
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.special import ndtr, ndtri


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
    rating curve, then optionally re-apply the Total Probability Theorem.

    Required Excel sim-list columns:
      Include, Method (= 'reservoir routing'), Output file,
      Run models, Analyse results, Input MCDF, Inflow,
      ELS file, SQ file, FSL, Hydrographs folder, Results folder

    Optional columns:
      Store hydrographs (default 'yes'), Output suffix, Config file,
      Log file, Comment.
    """

    def __init__(self, sim_row, filepaths, test_runs=0):
        self.start = time.time()
        self.sim_row = sim_row
        self.filepaths = filepaths
        self.test_runs = test_runs
        self.basename = str(sim_row['Output file'])

        print(f'\n=== ReservoirRoutingSimulator: {self.basename} ===')

        self.hydrographs_folder = str(sim_row['Hydrographs folder'])
        self.results_folder = str(sim_row['Results folder'])
        os.makedirs(self.hydrographs_folder, exist_ok=True)
        os.makedirs(self.results_folder, exist_ok=True)

        self.run_models = str(sim_row['Run models']).strip().lower() == 'yes'
        self.do_analysis = str(sim_row['Analyse results']).strip().lower() == 'yes'

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
        df = pd.read_csv(inflow_path, index_col=0)
        df.interpolate(method='slinear', inplace=True)
        self.time_index = df.index.to_numpy(dtype=float)
        self.sim_columns = df.columns.tolist()
        self.inflows_arr = df.to_numpy(dtype=float)
        self.dt_hours = float(self.time_index[1] - self.time_index[0])
        print(f'  {self.inflows_arr.shape[0]} timesteps x '
              f'{self.inflows_arr.shape[1]} sims, dt = {self.dt_hours} h '
              f'(read in {time.time() - t0:.2f} s)')

    def _load_mcdf(self):
        mcdf_path = str(self.sim_row['Input MCDF'])
        print(f'Loading MCDF: {mcdf_path}')
        self.mcdf = pd.read_csv(mcdf_path, index_col=0)
        if len(self.mcdf) != self.inflows_arr.shape[1]:
            raise ValueError(
                f'MCDF row count ({len(self.mcdf)}) does not match number of '
                f'inflow columns ({self.inflows_arr.shape[1]}).'
            )
        self.adv = self.mcdf['ADV'].to_numpy(dtype=float)

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

    def _write_mcdf(self):
        # nanmax is critical: inflow series interpolated via 'slinear' can have
        # trailing NaNs (sims with shorter URBS runs), and the routing arrays
        # propagate those NaNs to the tail of the level/outflow series.
        self.mcdf['inflow'] = np.nanmax(self.inflows_arr, axis=0)
        self.mcdf['outflow'] = np.nanmax(self.O, axis=0)
        self.mcdf['level'] = np.nanmax(self.L, axis=0)
        out_path = os.path.join(self.results_folder, f'{self.basename}__mcdf.csv')
        print(f'Writing new MCDF: {out_path}')
        self.mcdf.to_csv(out_path)

    # --------------------------------------------------------------- analyse

    def _ensure_mcdf_loaded(self):
        if hasattr(self, 'mcdf'):
            return
        out_path = os.path.join(self.results_folder, f'{self.basename}__mcdf.csv')
        if os.path.isfile(out_path):
            print(f'Loading MCDF for analysis: {out_path}')
            self.mcdf = pd.read_csv(out_path, index_col=0)
        else:
            mcdf_path = str(self.sim_row['Input MCDF'])
            print(f'Loading MCDF for analysis (fallback to input): {mcdf_path}')
            self.mcdf = pd.read_csv(mcdf_path, index_col=0)

    def _tpt_config(self):
        """Pull TPT params from the row's Config file JSON, with sensible defaults
        derived from the MCDF if the JSON is missing or partial."""
        cfg = {}
        if 'Config file' in self.sim_row.index and pd.notna(self.sim_row['Config file']):
            cfg_path = str(self.sim_row['Config file'])
            if os.path.isfile(cfg_path):
                with open(cfg_path) as f:
                    cfg = json.load(f)
        m_count = int(cfg.get('number_of_main_divisions', self.mcdf['m'].nunique()))
        n_count = int(cfg.get('number_of_sub_divisions', len(self.mcdf) // max(m_count, 1)))
        lower_aep = float(cfg.get('lower_aep', 2))
        upper_aep = float(cfg.get('upper_aep', 1e7))
        aep_of_pmp = cfg.get('aep_of_pmp', None)
        if aep_of_pmp is not None:
            aep_of_pmp = float(aep_of_pmp)
        return m_count, n_count, lower_aep, upper_aep, aep_of_pmp

    def _analyse(self):
        m_count, n_count, lower_aep, upper_aep, aep_of_pmp = self._tpt_config()
        print(f'\nTPT config: m={m_count}, n={n_count}, '
              f'lower=1in{lower_aep}, upper=1in{upper_aep}, aep_of_pmp={aep_of_pmp}')

        z_low = ndtri(1.0 - 1.0 / lower_aep)
        z_up = ndtri(1.0 - 1.0 / upper_aep)
        main_divisions = np.linspace(z_low, z_up, m_count + 1)

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

        out_path = os.path.join(self.results_folder, f'{self.basename}__mcdf.csv')
        print(f'\nUpdating MCDF with new *_aep columns: {out_path}')
        self.mcdf.to_csv(out_path)
