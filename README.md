# volc_enso

Utilities and analyses for exploring how major volcanic eruptions influence
ENSO behavior.

## Requirements

Install the scientific Python stack before running the modules or tests:

```
pip install -r requirements.txt
```

## Nino index utilities

The repository includes a ``nino_indices.py`` module that provides functions
for computing Niño 3, Niño 3.4, and Niño 4 indices using cosine‑latitude
weighting.  It also exposes ``calculate_nino34`` and ``season_nino34`` helper
functions that mirror the original ``V1_1x_Sam.ipynb`` notebook API so the
notebook can import these routines directly.

## ENSO response analysis

``enso_analysis.py`` offers a workflow to relate pre-eruption ENSO phase to
post-eruption surface-temperature anomalies.  It loads ensemble members for
each eruption month, classifies El Niño/La Niña/Neutral conditions using the
Niño 3.4 index, and composites the subsequent temperature response to evaluate
seasonal sensitivity of the 1258 Samalas eruption.  The module also provides
pairwise t-tests comparing January, April, July and October onset experiments,
quantifying whether global-mean TS and Niño 3.4 anomalies differ across
eruption seasons.  Analysis functions assume the record begins a fixed number
of months before the eruption but expose an ``eruption_index`` parameter to
override this for arbitrary time axes.

## Repository structure

- ``Samalas_ENSO_analysis.ipynb`` – consolidated notebook that demonstrates the
  full analysis workflow for one eruption month.
- ``V1_1x_Sam.ipynb`` – original notebook from which helper routines and setup
  code were ported.
- ``nino_indices.py`` – reusable cosine‑weighted ENSO index calculators.
- ``enso_analysis.py`` – functions for classifying pre‑eruption ENSO phase and
  compositing post‑eruption anomalies.
- ``samalas_setup.py`` – paths, plotting utilities, and constants shared across
  scripts.
- ``tests/`` – lightweight unit tests verifying index calculations and analysis
  helpers.

## Data

The utilities expect NetCDF files containing ensemble simulations of the 1258
Samalas eruption and a corresponding control run.  Each ensemble should supply
the ``TS`` variable with dimensions ``member × time × lat × lon``.  The
analysis has been exercised with 10 members and 99 months of data, but other
shapes are supported via the ``eruption_index`` argument.

## Usage

1. Install dependencies with ``pip install -r requirements.txt``.
2. Edit ``samalas_setup.py`` to point ``DATA_DIR`` at your NetCDF directory.
3. Open ``Samalas_ENSO_analysis.ipynb`` and run the cells, or call functions in
   ``enso_analysis.py`` from your own scripts to automate experiments across
   January, April, July, and October eruption timings.

## Running tests

Basic unit tests exercise the ENSO index calculators and the high‑level
analysis functions.  Run them with ``pytest``:

```bash
pytest tests
```

Some tests require optional packages (e.g., ``xarray``); missing dependencies
will cause those tests to be skipped.

## License

This project is released under the MIT License.  See ``LICENSE`` for details.
