# CSC3107 Week 12 - Raster Maps of SST Anomalies

This repository contains the completed Quarto analysis for the Week 12 lab on NOAA sea-surface temperature anomalies during the December 2009 El NiÃ±o and December 2010 La NiÃ±a phases.

## Submission artifacts

- `lab12_sst_anomalies.qmd` - reproducible R/Quarto source with all required verifier checks and seven rendered plots.
- `lab12_sst_anomalies.html` - self-contained rendered submission.
- `README.txt` - provenance supplied with the NOAA dataset.

The analysis uses the Mollweide equal-area projection for the final maps. The pooled 99th percentile of absolute open-ocean anomaly is 2.76 Â°C, giving shared symmetric limits of -3 Â°C to +3 Â°C.

## Data

The 62 daily NetCDF files are intentionally excluded from Git because they total about 102 MB and are published unchanged by NOAA/NCEI. Download the December 2009 and December 2010 AVHRR files from:

<https://www.ncei.noaa.gov/data/sea-surface-temperature-optimum-interpolation/v2.1/access/avhrr/>

Place them in this layout before rendering:

```text
oisst_daily/
  el_nino_2009-12/
    oisst-avhrr-v02r01.20091201.nc
    ...
    oisst-avhrr-v02r01.20091231.nc
  la_nina_2010-12/
    oisst-avhrr-v02r01.20101201.nc
    ...
    oisst-avhrr-v02r01.20101231.nc
```

Dataset DOI: <https://doi.org/10.25921/RE9P-PT57>

## Render

Tested with R 4.6.0, Quarto 1.9.38, and these R packages: `tidyverse` 2.0.0, `stars` 0.7-3, `sf` 1.1.1, `scales` 1.4.0, `ggnewscale` 0.5.2, `units` 1.0.1, and `lwgeom` 0.2-16.

```powershell
quarto render lab12_sst_anomalies.qmd --cache-refresh
```
