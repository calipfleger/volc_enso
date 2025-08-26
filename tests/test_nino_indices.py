import pytest

np = pytest.importorskip("numpy")
xr = pytest.importorskip("xarray")
cftime = pytest.importorskip("cftime")

from nino_indices import calculate_nino34, season_nino34


def _synthetic_sst():
    time = [cftime.DatetimeNoLeap(1255 + (i // 12), (i % 12) + 1, 1) for i in range(24)]
    lat = np.linspace(-10, 10, 3)
    lon = np.linspace(190, 240, 4)
    data = np.random.rand(len(time), len(lat), len(lon))
    return xr.DataArray(data, coords={"time": time, "lat": lat, "lon": lon},
                        dims=("time", "lat", "lon"))


def test_nino34_and_seasons():
    sst = _synthetic_sst()
    n34 = calculate_nino34(sst)
    djf, jja = season_nino34(n34)
    assert n34.dims == ("time",)
    # Ensure seasons are filtered correctly
    assert set(djf['time'].dt.season.values) == {"DJF"}
    assert set(jja['time'].dt.season.values) == {"JJA"}
