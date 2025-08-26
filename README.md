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
