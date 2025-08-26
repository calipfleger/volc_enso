"""Utilities for computing ENSO-related sea surface temperature indices.

These functions mirror the approach used in ``V1_1x_Sam.ipynb`` for
selecting latitude/longitude regions and applying cosine-latitude
weights prior to averaging.
"""

from __future__ import annotations

import numpy as np
import xarray as xr


def _area_weighted_mean(
    data: xr.DataArray,
    lat_name: str = "lat",
    lon_name: str = "lon",
) -> xr.DataArray:
    """Return cosine(latitude) weighted mean over ``lat`` and ``lon``."""
    weights = np.cos(np.deg2rad(data[lat_name]))
    return data.weighted(weights).mean(dim=[lat_name, lon_name])


def nino_region(
    data: xr.DataArray,
    lat_bounds: tuple[float, float],
    lon_bounds: tuple[float, float],
    lat_name: str = "lat",
    lon_name: str = "lon",
) -> xr.DataArray:
    """Compute area-mean SST index over a latitude/longitude box.

    Parameters
    ----------
    data : xr.DataArray
        Sea surface temperature field with ``lat`` and ``lon`` dimensions.
    lat_bounds, lon_bounds : tuple
        (min, max) bounds in degrees.
    """
    subset = data.sel({lat_name: slice(*lat_bounds), lon_name: slice(*lon_bounds)})
    return _area_weighted_mean(subset, lat_name=lat_name, lon_name=lon_name)


def nino34(
    data: xr.DataArray,
    lat_name: str = "lat",
    lon_name: str = "lon",
) -> xr.DataArray:
    """Niño 3.4 index: 5°S-5°N, 190°E-240°E."""
    return nino_region(data, (-5, 5), (190, 240), lat_name, lon_name)


def nino3(
    data: xr.DataArray,
    lat_name: str = "lat",
    lon_name: str = "lon",
) -> xr.DataArray:
    """Niño 3 index: 5°S-5°N, 210°E-270°E."""
    return nino_region(data, (-5, 5), (210, 270), lat_name, lon_name)


def nino4(
    data: xr.DataArray,
    lat_name: str = "lat",
    lon_name: str = "lon",
) -> xr.DataArray:
    """Niño 4 index: 5°S-5°N, 160°E-210°E."""
    return nino_region(data, (-5, 5), (160, 210), lat_name, lon_name)


# -----------------------------------------------------------------------------
# Notebook-style helpers copied from V1_1x_Sam.ipynb
# -----------------------------------------------------------------------------

def calculate_nino34(data: xr.DataArray) -> xr.DataArray:
    """Calculate NINO3.4 index from temperature data."""
    nino34 = data.sel(lat=slice(-5, 5), lon=slice(190, 240))
    weights = np.cos(np.deg2rad(nino34.lat))
    nino34_avg = nino34.weighted(weights).mean(dim=["lat", "lon"])
    return nino34_avg


def season_nino34(m: xr.DataArray) -> tuple[xr.DataArray, xr.DataArray]:
    """Split a NINO3.4 index into DJF and JJA components."""
    DA_DJF = m.sel(time=m.time.dt.season == "DJF")
    DA_JJA = m.sel(time=m.time.dt.season == "JJA")
    return DA_DJF, DA_JJA

