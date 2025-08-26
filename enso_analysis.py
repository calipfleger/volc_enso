"""ENSO analysis utilities for Samalas eruption studies."""

from __future__ import annotations

from typing import Iterable, Dict, Tuple

from itertools import combinations
import numpy as np
import xarray as xr
from scipy import stats

from nino_indices import nino34
from samalas_setup import FILEPATH, VAR_NAME


def load_ts_ensemble(onset: str, members: Iterable[int], variable: str = "TS") -> xr.DataArray:
    """Load TS ensemble members for a given eruption onset month.

    Parameters
    ----------
    onset : str
        Name of the eruption onset directory, e.g. ``'January_1x'``.
    members : iterable of int
        Member numbers to load.
    variable : str, optional
        Variable name in the NetCDF files, by default ``"TS"``.

    Returns
    -------
    xr.DataArray
        Concatenated DataArray with a new ``member`` dimension.
    """
    paths = [f"{FILEPATH}{onset}/{variable}_m{m:02d}.nc" for m in members]
    datasets = [xr.open_dataarray(p) for p in paths]
    return xr.concat(datasets, dim="member")


def compute_ts_anomaly(exp: xr.DataArray, control: xr.DataArray) -> xr.DataArray:
    """Return experiment minus control monthly anomalies."""
    control_clim = control.groupby("time.month").mean("time")
    return exp.groupby("time.month") - control_clim


def classify_pre_eruption_phase(
    sst: xr.DataArray,
    months: int = 12,
    threshold: float = 1.0,
    eruption_index: int | None = None,
) -> xr.DataArray:
    """Label ENSO phase based on the pre-eruption Nino3.4 index.

    The input data must contain at least ``months`` time steps prior to the
    eruption.  By default the eruption is assumed to occur ``months`` steps
    into the record; override ``eruption_index`` if the eruption occurs at a
    different position.
    """
    if eruption_index is None:
        eruption_index = months
    n34 = nino34(sst)
    pre = n34.isel(time=slice(eruption_index - months, eruption_index))
    mean = pre.mean("time")
    std = pre.std("time")
    phase = xr.full_like(mean, "Neutral", dtype=object)
    phase = xr.where(mean > threshold * std, "El Nino", phase)
    phase = xr.where(mean < -threshold * std, "La Nina", phase)
    return phase


def composite_post_eruption(
    ts_anom: xr.DataArray,
    phases: xr.DataArray,
    months: int = 24,
    eruption_index: int | None = None,
) -> Dict[str, xr.DataArray]:
    """Composite post-eruption TS anomalies by ENSO phase."""
    if eruption_index is None:
        eruption_index = 0
    post = ts_anom.isel(time=slice(eruption_index, eruption_index + months))
    comps: Dict[str, xr.DataArray] = {}
    for phase in ["El Nino", "Neutral", "La Nina"]:
        mask = phases == phase
        if mask.any():
            comps[phase] = post.sel(member=mask).mean("member")
    return comps


def analyze_eruption_seasonality(
    months: Iterable[str] = VAR_NAME,
    members: Iterable[int] = range(1, 11),
    variable: str = "TS",
    pre_months: int = 12,
    post_months: int = 24,
    control: xr.DataArray | None = None,
    eruption_index: int | None = None,
) -> Dict[str, Dict[str, xr.DataArray]]:
    """Analyze how pre-eruption ENSO phase modulates post-eruption TS response.

    Parameters
    ----------
    months : iterable of str, optional
        Eruption onset directories to examine.
    members : iterable of int, optional
        Ensemble member numbers to load.
    variable : str, optional
        Variable name in the NetCDF files.
    pre_months : int, optional
        Number of months before eruption for ENSO classification.
    post_months : int, optional
        Number of months after eruption for compositing.
    control : xr.DataArray, optional
        Control run for anomaly calculation. If ``None``, anomalies are
        computed relative to the ensemble mean of the pre-eruption period.
    eruption_index : int, optional
        Index of the eruption month within the time dimension. Defaults to
        ``pre_months`` assuming the record starts ``pre_months`` months before
        the eruption.

    Returns
    -------
    dict
        Nested dictionary mapping onset month -> phase -> composite DataArray.
    """
    if eruption_index is None:
        eruption_index = pre_months
    results: Dict[str, Dict[str, xr.DataArray]] = {}
    for onset in months:
        ts = load_ts_ensemble(onset, members, variable=variable)
        if control is not None:
            ts_anom = compute_ts_anomaly(ts, control)
        else:
            base = ts.isel(time=slice(eruption_index - pre_months, eruption_index)).mean("time")
            ts_anom = ts - base
        phase_labels = classify_pre_eruption_phase(ts, months=pre_months, eruption_index=eruption_index)
        composites = composite_post_eruption(ts_anom, phase_labels, months=post_months, eruption_index=eruption_index)
        results[onset] = composites
    return results


def _global_mean(data: xr.DataArray, lat_name: str = "lat", lon_name: str = "lon") -> xr.DataArray:
    """Cosine-latitude weighted global mean."""
    weights = np.cos(np.deg2rad(data[lat_name]))
    return data.weighted(weights).mean(dim=[lat_name, lon_name])


def _nino34_anomaly(
    ts: xr.DataArray,
    control: xr.DataArray | None = None,
    pre_months: int = 12,
    eruption_index: int | None = None,
) -> xr.DataArray:
    """Return Niño 3.4 index anomalies from a TS ensemble.

    Parameters
    ----------
    ts : xr.DataArray
        TS ensemble with a ``time`` dimension.
    control : xr.DataArray, optional
        Control run for anomaly calculation.
    pre_months : int, optional
        Number of months prior to the eruption used for the baseline.
    eruption_index : int, optional
        Index of the eruption month within the time dimension. Defaults to
        ``pre_months``.
    """
    if eruption_index is None:
        eruption_index = pre_months
    if control is not None:
        ts_anom = compute_ts_anomaly(ts, control)
    else:
        base = ts.isel(time=slice(eruption_index - pre_months, eruption_index)).mean("time")
        ts_anom = ts - base
    return nino34(ts_anom)


def ttest_onset_differences(
    onsets: Iterable[str] = VAR_NAME,
    members: Iterable[int] = range(1, 11),
    variable: str = "TS",
    pre_months: int = 12,
    post_months: int = 24,
    control: xr.DataArray | None = None,
    eruption_index: int | None = None,
) -> Dict[Tuple[str, str], xr.Dataset]:
    """Pairwise t-tests comparing eruption months.

    Computes global-mean TS and Niño 3.4 anomalies for each onset month and
    returns Student's t statistics and p-values for all month pairs.

    Parameters
    ----------
    onsets : iterable of str, optional
        Eruption onset directories to compare.
    members : iterable of int, optional
        Ensemble member numbers to load.
    variable : str, optional
        Variable name in the NetCDF files.
    pre_months : int, optional
        Months before the eruption used for baseline.
    post_months : int, optional
        Months after the eruption included in the test.
    control : xr.DataArray, optional
        Control run for anomaly calculation.
    eruption_index : int, optional
        Index of the eruption month within the time dimension. Defaults to
        ``pre_months``.
    """
    if eruption_index is None:
        eruption_index = pre_months
    ts_means: Dict[str, xr.DataArray] = {}
    n34_anoms: Dict[str, xr.DataArray] = {}
    for onset in onsets:
        ts = load_ts_ensemble(onset, members, variable=variable)
        if control is not None:
            ts_anom = compute_ts_anomaly(ts, control)
        else:
            base = ts.isel(time=slice(eruption_index - pre_months, eruption_index)).mean("time")
            ts_anom = ts - base
        ts_mean = _global_mean(ts_anom).isel(time=slice(eruption_index, eruption_index + post_months))
        n34 = _nino34_anomaly(ts, control, pre_months, eruption_index).isel(time=slice(eruption_index, eruption_index + post_months))
        ts_means[onset] = ts_mean
        n34_anoms[onset] = n34

    results: Dict[Tuple[str, str], xr.Dataset] = {}
    for a, b in combinations(onsets, 2):
        ts_t, ts_p = stats.ttest_ind(ts_means[a], ts_means[b], axis=0, equal_var=False)
        n34_t, n34_p = stats.ttest_ind(n34_anoms[a], n34_anoms[b], axis=0, equal_var=False)
        results[(a, b)] = xr.Dataset(
            {
                "ts_t": ("time", ts_t),
                "ts_p": ("time", ts_p),
                "n34_t": ("time", n34_t),
                "n34_p": ("time", n34_p),
            },
            coords={"time": ts_means[a].time},
        )
    return results

