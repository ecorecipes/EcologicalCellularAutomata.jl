# High-water events at North Core Banks

`synthetic_storms.csv` and `storm_series.csv` are taken from the input files of MEEB, the Mesoscale Explicit
Ecogeomorphic Barrier model, version 1.0.1 (I. R. B. Reeves, U.S. Geological Survey,
https://code.usgs.gov/spcmsc/meeb/-/tree/1.0.1; Reeves et al. 2025, Geoscientific Model Development 18,
9319-9348, https://doi.org/10.5194/gmd-18-9319-2025). MEEB is in the public domain in the United States as
a work of the U.S. Geological Survey, and its authors waive copyright and related rights worldwide through
the CC0 1.0 Universal dedication. No warranty is made by the USGS or the U.S. Government as to the data. These
files are not covered by this package's licence.

## What was done

The two NumPy arrays of MEEB 1.0.1's `Input/` folder were read and written as CSV, rounded to six decimals:

- `SyntheticStorms_NCB-CE_10k_1979-2020_Beta0pt039_BermEl1pt78.npy` (9999 events from the copula sea-storm
  model of Wahl et al. 2016, fitted to the 1979-2020 hindcast offshore of North Core Banks, NC) became
  `synthetic_storms.csv`, keeping the four columns MEEB uses (its first, second, fourth and fifth); the
  total water level, runup and lowest runup it also holds are recomputed by MEEB on the local beach slope.
- `StormTimeSeries_1979-2020_NCB-CE_Beta0pt039_BermEl1pt78.npy` (222 events of the hindcast) became
  `storm_series.csv` unchanged.

MEEB 1.1.0 ships the same synthetic events.

## Format

`synthetic_storms.csv`: a header, then 9999 rows of `Hs` (significant wave height, m), `duration` (h), `SL`
(still water level, m NAVD88) and `Tp` (peak wave period, s).

`storm_series.csv`: a header, then 222 rows of `iteration` (the week of the event, counted in fiftieths of a
year from 1 January 1979, always even), `Rhigh` (total water level, m NAVD88), `Rlow` (lowest runup, m
NAVD88) and `duration` (h).

# The Florence test site

`florence_site.bin` is cropped from three more arrays of MEEB 1.0.1's `Input/` folder (same provenance and
dedication): `Init_NCB-NewDrum-Ocracoke_2017_PreFlorence_18400-23400.npy` (elevation and plant cover from 2017
lidar and imagery), `Init_NCB-NewDrum-Ocracoke_2018_PostFlorence_18400-23400.npy` (elevation after Hurricane
Florence, October 2018) and `Mask_NCB-NewDrum-Ocracoke_2018_Florence_18400-23400.npy` (the overwash extent
digitised from imagery and lidar). It holds the test site of the paper's Fig. 5a and Table 1, the one
MEEB's `Tools/test_overwash.py` runs: rows 850 to 1349 (0-based) of those arrays, 19250 to 19750 m alongshore
on North Core Banks, and columns 149 to 849 (0-based), which hold the whole barrier; the 149 columns seaward
of them are all below MHW. Elevations were rounded to the millimetre and plant covers to 1/250.

## Format

Little-endian binary: the 8 bytes `MEEBSITE`; the number of rows and of columns (`Int32` each, 500 and 701);
then five arrays of rows × columns, column-major: the elevation before and after the storm (`Int16`,
millimetres NAVD88), the cover of the burial-tolerant and of the burial-intolerant type (`UInt8`, in 1/250),
and the overwash extent (`UInt8`, 0 or 1).
