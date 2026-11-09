"""The world's maximum observed floods, and the PMF and PMPF laid against them.

A reasonableness check for the report: the dam's peak PMF and PMPF inflows
plotted on the catchment area - peak discharge plane beside the largest floods
ever recorded and the envelope drawn over them.

**The source** is Li, Wang and Li (2013), *Maximum observed floods in China*,
Hydrological Sciences Journal 58(3), 728-735 (doi:10.1080/02626667.2013.772299):

- ``WORLD_FLOODS`` is its Table 1, "World's greatest known floods": 54 floods
  from 0.9 to 4,640,300 km². Eleven are *surveyed* rather than gauged (its
  asterisk) and the Amazon at Obidos is a *modified* value, 303,000 m³/s in
  place of the original 370,000 (its double asterisk). Transcribed from the
  table and checked against its Fig. 1.
- ``wmo_envelope`` is its equation (2), the envelope over Table 1 published in
  the WMO *Manual on Estimation of Probable Maximum Precipitation* (WMO-No. 1045,
  2009): Qm = 1830 A^0.316 for 300 < A < 3,000,000 km², 154 A^0.738 below
  300 km². The two pieces step at 300 km² (10,367 to 11,097 m³/s) because the
  equation does, and it is drawn as written. It is a fitted envelope, not a
  strict bound: six of the 54 floods sit up to 12% above it, as in the paper's
  Fig. 1, so a PMF just over the line is not by itself out of the ordinary.
- ``herschy_envelope`` is its equation (1), Herschy's 2003 envelope (IAHS
  2003), kept for comparison with older reports.

**What is plotted for each PMF group** (core/ensemble.py's PMF groups):

- the **PMF inflow** is the largest inflow in the group's ensemble database -
  the peak of the largest flood, which is what a maximum observed flood is.
  That is not always the inflow of the event that gives the PMF *level*
  (``ensemble.pick``), which the PMF table quotes; the two are shown side by
  side on the page so the difference is visible.
- the **PMPF inflow** is the Monte Carlo group's inflow envelope - the largest
  inflow over its durations - at the AEP of the PMP, as the "Peaks at one AEP"
  table reads it with ``flows_at = own``.

pandas, numpy and the standard library only; nothing here imports nicegui.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

SOURCE = ("Li, C., Wang, G. and Li, R. (2013). Maximum observed floods in China. "
          "Hydrological Sciences Journal 58(3), 728-735.")
WMO_LABEL = "WMO (2009) envelope"
HERSCHY_LABEL = "Herschy (2003) envelope"


@dataclass(frozen=True)
class Flood:
    number: int
    country: str
    station: str
    area_km2: float
    peak_m3s: float
    date: str
    surveyed: bool = False
    modified: bool = False


_S, _M = "surveyed", "modified"

# Li et al. (2013) Table 1: (no., country, station, area km², peak m³/s, date, note)
_TABLE_1 = [
    (1, "China", "Sunjiazhai, Tamigou, Taohe", 0.9, 159, "10/08/1979", _S),
    (2, "Puerto Rico (USA)", "Q. de Los Cedros, Isabela", 1.79, 212, "07/05/1970", ""),
    (3, "USA", "Nevada, Humboldt River trib", 2.2, 251, "31/05/1973", ""),
    (4, "China", "Nangou, Shandong", 3.27, 366, "22/07/1989", _S),
    (5, "USA", "Big Creek, nr Waynesville, North Carolina", 3.4, 368, "30/08/1940", ""),
    (6, "China", "Chengjiaqu, Chengjiaqugou, Fenhe", 5.6, 457, "31/07/1971", _S),
    (7, "USA", "Halawa, Halawa Stream, Halawa (Molokai)", 12, 762, "04/02/1965", ""),
    (8, "USA", "Lane Canyon, Oregon", 13.1, 807, "26/07/1965", ""),
    (9, "USA", "El Rancho Arroyo nr Pojoazue, New Mexico", 17.4, 1_250, "22/08/1952", ""),
    (10, "USA", "Meyers Canyon, Oregon", 32.9, 1_540, "13/07/1956", ""),
    (11, "USA", "Bronco CK, near Wikieup, Arizona", 49.2, 2_080, "18/08/1971", ""),
    (12, "USA", "Colorado, El Dorado Canyon, Nelson Landing", 56.5, 2_150, "14/09/1974", ""),
    (13, "USA", "Wailua, South Fork Wailua, Lihue", 58, 2_470, "15/04/1963", ""),
    (14, "China", "Baisha, Nandujiang, Hainan", 75.3, 3_420, "1894", _S),
    (15, "Mexico", "Arroyo San Bartolo, Arroyo San Bartolo, San Bartolo", 81, 3_000,
     "30/09/1976", ""),
    (16, "China", "Shaohuyingzi, Tangtouhe, Dalinghe", 97.2, 4_000, "01/08/1930", ""),
    (17, "China", "Xiyuangong, Xiyuanxi, Minjiang", 142, 4_600, "19/09/1909", _S),
    (18, "China", "Waziyu, Waziyuhe, Dalinghe", 154, 5_320, "01/08/1930", ""),
    (19, "China", "Xiaoshankou, Shihe, Bohaiwan", 171, 7_000, "07/1894", _S),
    (20, "Taiwan, China", "Tungton, Choshu, Choshui", 259, 7_780, "24/08/1979", ""),
    (21, "New Caledonia", "Quaième, Quaième, Derniers Rapides", 330, 10_400, "24/12/1981", ""),
    (22, "China", "Banqiao, Ruhe, Huaihe", 762, 13_100, "07/08/1975", ""),
    (23, "USA", "West Nueces River, Texas", 1_041, 16_400, "14/06/1935", ""),
    (24, "China", "Foziling, Pihe, Huaihe", 1_840, 17_600, "14/07/1969", ""),
    (25, "Taiwan, China", "Daduqiao, Wuxi, Taiwan", 1_980, 18_300, "08/08/1959", ""),
    (26, "Taiwan, China", "Jiji, Choshui, Taiwan", 2_310, 20_000, "01/08/1996", ""),
    (27, "Japan", "Shingu Oga", 2_350, 19_025, "26/09/1959", ""),
    (28, "Philippines", "Cagayan, Cagayan, Echaque Isabella", 4_244, 17_550, "1959", ""),
    (29, "China", "Baoqiao, Changhuajiang, Hainan", 4_634, 28_300, "1887", _S),
    (30, "India", "Madhopur, Ravi, Indus", 6_087, 26_050, "01/09/1988", ""),
    (31, "USA", "Pecos Comstock, Texas", 9_100, 26_800, "28/06/1954", ""),
    (32, "China", "Linyi, Yihe", 10_315, 30_000, "09/08/1730", _S),
    (33, "Democratic People's Republic of Korea", "Taedong-gang, Taedong-gang, Mirom",
     12_175, 29_000, "29/08/1967", ""),
    (34, "China", "Sanjiangkou, Lishui, Changjiang", 15_242, 31_100, "05/07/1935", ""),
    (35, "China", "Fengtan, Qujiang, Changjiang", 16_595, 32_300, "22/09/1847", _S),
    (36, "India", "Jhalawaw, Chambal, Ganga", 22_584, 37_000, "08/1969", ""),
    (37, "Republic of Korea", "Han, Han, Goan", 23_880, 37_000, "18/07/1925", ""),
    (38, "India", "Betwa Sahijna", 44_870, 43_800, "26/07/1971", ""),
    (39, "Madagascar", "Mangoky, Mangoky, Banian", 50_000, 38_000, "05/02/1933", ""),
    (40, "China", "Huanggou, Yalujiang", 55_420, 44_800, "11/08/1888", ""),
    (41, "India", "Gujarat, Wuyika", 62_225, 42_475, "08/1968", ""),
    (42, "India", "Narmada, Gurdeshwar", 88_000, 69_400, "06/09/1970", ""),
    (43, "India", "Namada", 98_420, 70_790, "1954", ""),
    (44, "India", "Godavari", 299_320, 78_690, "17/09/1959", ""),
    (45, "India", "Godavari, Polaeshwaram", 307_800, 87_250, "1907", ""),
    (46, "India", "Godavari Dowlaishwaram", 315_000, 88_400, "1959", ""),
    (47, "Bangladesh", "Brahmaputra Jiamuna", 530_000, 93_500, "01/08/1955", ""),
    (48, "Bangladesh", "Brahmaputra, Bahadurabad", 580_000, 98_500, "30/08/1988", ""),
    (49, "China", "Wanxian, Changjiang", 974_900, 108_000, "18/07/1870", _S),
    (50, "China", "Yichang, Changjiang", 1_005_500, 105_000, "20/07/1870", _S),
    (51, "Bangladesh", "Ganges, Baluliya", 1_650_000, 132_000, "01/09/1988", ""),
    (52, "Russia", "Laptev Sea, Lena, Kusur", 2_430_000, 189_000, "08/06/1967", ""),
    (53, "Russia", "Lena, Kusur", 2_430_000, 194_000, "11/06/1944", ""),
    (54, "Brazil", "Obidos, Amazon", 4_640_300, 303_000, "06/1953", _M),
]

WORLD_FLOODS = [Flood(number, country, station, float(area), float(peak), date,
                      surveyed=note == _S, modified=note == _M)
                for number, country, station, area, peak, date, note in _TABLE_1]


def wmo_envelope(area_km2: float) -> float:
    """WMO (2009) envelope, Li et al. (2013) eq. (2); NaN outside 0-3,000,000 km²."""
    area = float(area_km2)
    if not 0 < area < 3_000_000:
        return math.nan
    return 154.0 * area ** 0.738 if area < 300 else 1830.0 * area ** 0.316


def herschy_envelope(area_km2: float) -> float:
    """Herschy (2003) envelope, Li et al. (2013) eq. (1)."""
    area = float(area_km2)
    if not area > 0:
        return math.nan
    return 100.0 * area ** 0.8 if area <= 90 else 500.0 * area ** 0.43


def envelope_curve(function, low=0.5, high=3_000_000, per_decade=24) -> list:
    """(area, peak) pairs along an envelope, evenly spaced in log area, with the
    WMO curve's step at 300 km² kept as a vertical rather than smoothed over."""
    points = []
    count = int(round(math.log10(high / low) * per_decade))
    for step in range(count + 1):
        area = low * (high / low) ** (step / count)
        points.append(area)
    if function is wmo_envelope:
        points += [299.999, 300.0]
    out = []
    for area in sorted(set(points)):
        value = function(min(area, high * 0.999999))
        if math.isfinite(value):
            out.append((area, value))
    return out


# -- the dam's own floods --------------------------------------------------------

@dataclass
class DamFlood:
    """One PMF group's two floods on the world plane."""

    label: str
    pmf: float = math.nan            # m³/s, the largest ensemble inflow
    pmf_level_event: float = math.nan  # m³/s, the inflow of the PMF-level event
    pmpf: float = math.nan           # m³/s, the inflow envelope at the AEP of the PMP
    problems: list = field(default_factory=list)


def ratio_to_envelope(area_km2, peak) -> float:
    """A peak as a fraction of the WMO envelope at the catchment's area."""
    limit = wmo_envelope(area_km2) if area_km2 else math.nan
    if not (math.isfinite(limit) and peak and math.isfinite(peak)):
        return math.nan
    return float(peak) / limit


def dam_floods(study, section, durations=None) -> list:
    """A DamFlood per PMF group in the study's PMF section. Reads results, so the
    page calls it off the event loop."""
    from . import ensemble, reporttables
    from .study import StudyError

    pmp_aep = float(section.get("pmp_aep") or ensemble.DEFAULT_PMP_AEP)
    out = []
    for entry in section.get("groups") or []:
        flood = DamFlood(label=entry.get("label") or entry["ensemble"].get("group", ""))
        try:
            frame = ensemble.load(study, entry["ensemble"]["run"], entry["ensemble"]["group"])
            inflow = frame["inflow"].dropna()
            if inflow.empty:
                flood.problems.append(f"{flood.label}: no inflows in the PMF database")
            else:
                flood.pmf = float(inflow.max())
            chosen = ensemble.pick(frame, ensemble.HIGHEST)
            if chosen.found:
                flood.pmf_level_event = float(chosen.inflow)
        except StudyError as exc:
            flood.problems.append(f"{flood.label} PMF: {exc}")
        mc = entry.get("mc") or {}
        if mc.get("group"):
            curves = reporttables.group_curves(study, mc.get("run", ""), mc["group"],
                                               durations)
            value = curves.envelope("inflow").get(pmp_aep, math.nan)
            if math.isfinite(float(value)):
                flood.pmpf = float(value)
            else:
                flood.problems.append(f"{flood.label} PMPF: no inflow quantile at "
                                      f"1 in {pmp_aep:,.0f} in {mc['group']}")
        else:
            flood.problems.append(f"{flood.label} PMPF: no Monte Carlo group")
        out.append(flood)
    return out
