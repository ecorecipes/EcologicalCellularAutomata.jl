# Water levels at Scheveningen

`scheveningen_levels.csv` is derived from the water levels measured at Scheveningen by Rijkswaterstaat, the
Dutch directorate for public works and water management, as published through its WaterWebservices
(https://rwsdata.nl/waterdata/; the data can also be viewed at https://waterinfo.rws.nl). Taken on
6 November 2026. Source: Rijkswaterstaat, https://rwsdata.nl.

Rijkswaterstaat states that anyone may view its open data and use it without attribution ("Iedereen is vrij
deze open data in te zien en zonder bronvermelding te gebruiken"), and the Dutch national data catalogue
lists Waterinfo Rijkswateren as public domain ("Publiek domein",
https://data.overheid.nl/dataset/9442-waterinfo-rijkswateren). Rijkswaterstaat accepts no liability for
errors in the data. Its website asks those who reuse its information to cite it, give the date, and let it
know through its contact form. This file is not covered by this package's licence.

## What was done

`make_scheveningen_levels.jl`, beside this file, made it from the measured high and low waters (`WATHTE` in
cm against NAP, grouping `GETETM2`, process type `meting`, location `scheveningen`). Each row is the highest
level in a block of 14.6 days from 1 January 1970, one spring-neap cycle, for the first 1250 blocks: 50
years of 25 levels, the forcing Teixeira et al. (2023) describe. The values are otherwise unchanged. No
measurement was excluded. Some recent values are marked by Rijkswaterstaat as not yet checked.

## Format

A header, then 1250 rows: `cycle_start`, the date a block starts, and `level_m`, its highest level in metres
against NAP, the Dutch datum. Teixeira et al. write their levels against mean sea level; no correction was
made for the difference.
