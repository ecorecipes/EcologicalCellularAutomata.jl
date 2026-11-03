# Asian monthly climate (CRU CL 2.0)

`cru_asia.bin` contains information from CRU CL v2.0, the 10-minute climatology of mean monthly surface
climate over global land areas (1961–1990) of New, Lister, Hulme & Makin (2002), *Climate Research* 21,
1–25, published by the Climatic Research Unit, University of East Anglia
(https://crudata.uea.ac.uk/cru/data/hrg/tmc/). It is a derived database made available under the Open
Database License (ODbL) 1.0, https://opendatacommons.org/licenses/odbl/1-0/; any rights in its individual
contents are licensed under the Database Contents License, https://opendatacommons.org/licenses/dbcl/1-0/.
Attribution: Climatic Research Unit, University of East Anglia.

This file is distributed under the ODbL, not under this package's licence.

## What was done

`make_cru_asia.jl`, beside this file, made it from `grid_10min_tmp.dat.gz` (mean temperature) and
`grid_10min_reh.dat.gz` (relative humidity). The 10' cells were grouped in blocks of 2 × 2 into cells of
1/3 degree over 26° E–180° E and 11° S–78° N (267 rows from the north, 462 columns from the west); a cell
was kept if at least two of its four 10' cells are land and it lies in Asia, within a rough boundary (no
Europe west of the Turkish straits, the Black Sea, the Caucasus and the Ural river and mountains; no Africa
west of Suez and of a line along the Red Sea to Bab-el-Mandeb, nor the Horn of Africa; no New Guinea), and its values are the means over its land cells, rounded to 0.1.

## Format

49,628 records of 26 little-endian 16-bit integers: the row and the column of the cell, its 12 monthly
mean temperatures in tenths of a degree Celsius and its 12 monthly mean relative humidities in tenths of a
per cent, January first.
