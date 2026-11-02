# Simiatug landscape

`simiatug.csv` holds four layers of the landscape published with the model "Pest dispersion (simplified
representation)", version 1.1.0, by François Rebaudo (IRD), 2013, https://www.comses.net/codebases/2278/releases/1.1.0/,
the code of Rebaudo, Crespo-Pérez, Silvain & Dangles (2011), Journal of Artificial Societies and Social
Simulation 14(3), 7 (doi 10.18564/jasss.1802).

The original files (`elevation_simiatug.txt`, `temperature_simiatug.txt`, `village_simiatug.txt`,
`population_simiatug.txt`) are licensed under the Academic Free License 3.0 (AFL-3.0); they were joined
into one table, one row per cell, with the columns `x, y` (the source's patch coordinates, 0 to 79),
`elevation` (m), `temperature` (°C), `village` (an identifier, 0 for none) and `population`. The values
are unchanged. The source's code describes the village populations as invented: real population sizes
were not available to its authors.

This file is distributed under the AFL-3.0, not under this package's licence.
