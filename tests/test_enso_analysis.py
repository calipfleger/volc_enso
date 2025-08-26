import pytest

xr = pytest.importorskip("xarray")
np = pytest.importorskip("numpy")
cftime = pytest.importorskip("cftime")

import enso_analysis as ea


def _synthetic_ts():
    """Return synthetic TS data shaped like the ensembles."""
    # 10 members, 99 months, 96 lat, 144 lon
    time = [cftime.DatetimeNoLeap(1255 + (i + 11)//12, ((i + 11) % 12) + 1, 1) for i in range(99)]
    lat = np.linspace(-90, 90, 96)
    lon = np.linspace(0, 357.5, 144)
    data = np.random.rand(10, 99, 96, 144)
    return xr.DataArray(data, coords={"member": range(10), "time": time, "lat": lat, "lon": lon},
                        dims=("member", "time", "lat", "lon"), name="TS")


def test_analyze_eruption_seasonality_runs():
    fake = _synthetic_ts()
    ea.load_ts_ensemble = lambda onset, members, variable="TS": fake
    res = ea.analyze_eruption_seasonality(
        months=["January_1x", "April_1x", "July_1x", "October_1x"],
        members=range(1, 11),
        eruption_index=12,
    )
    assert set(res.keys()) == {"January_1x", "April_1x", "July_1x", "October_1x"}
    for comps in res.values():
        for phase in comps:
            assert comps[phase].shape[0] == 24  # time dimension


def test_ttest_onset_differences_runs():
    fake = _synthetic_ts()
    ea.load_ts_ensemble = lambda onset, members, variable="TS": fake
    res = ea.ttest_onset_differences(
        onsets=["January_1x", "April_1x", "July_1x", "October_1x"],
        members=range(1, 11),
        eruption_index=12,
    )
    assert len(res) == 6  # number of pairwise combinations
